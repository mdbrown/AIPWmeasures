#' Calculate accuracy measures using augmented inverse probability weighted estimates
#'
#' Estimates AIPW
#'
#' @param
#'
#'
#' @return data.frame with estimates and standard errors.
#'
#' @examples
#'
#'
#' @export
#'
AIPWmeasures <- function(N, cohortdata,
                          vp,
                          typex,
                          typey,
                          NBp,
                          predict.time,
                          method, design=NULL, AugPos=NULL, CalVar)
{
  if (method == 'TrueW'){
  	if (design == 'CCH')  {
  		cohortdata$wi = 1/P0HAT.CCH.FUN(cohortdata, type=2, ncch0 = sum(cohortdata$vi==1 & cohortdata$di == 0), ncch1 = sum(cohortdata$vi==1 & cohortdata$di == 1))

       }

  	if (design == 'NCC')  {
  		cohortdata$wi = ifelse(cohortdata$di==1,1 ,1/P0HAT.NCC.FUN(cohortdata, NULL, nmatch)) }

    subdata <- cohortdata[cohortdata$vi == 1,]
   #print(sum(subdata$wi))
   }

   if (method == 'TrueWZ'){
  		if (design == 'CCH')  {

  		cohortdata$wi = 1/P0HAT.CCH.Z.FUN(cohortdata, type=2);
      subdata <- cohortdata[cohortdata$vi == 1,]
  		}
  		if (design=='NCC') {
  			subdata <- cohortdata[cohortdata$vi == 1,]
  		}
  	}
   #print(sum(subdata$wi))


   if (method == 'AugW') {

   	    if (is.null(AugPos)) {print("Need the positions of the variables used in estimating weights")}
    	tmp.data0 <- cohortdata[cohortdata$di==0,]
		tmp.data1 <- cohortdata[cohortdata$di==1,]
#browser()
#    input.sm <- tmp.data0[,AugPos]

#		tmp.data0$wi = 1/sm.regression(input.sm,tmp.data0$vi,eval.points=input.sm,eval.grid=FALSE, display =   			"none")$estimate
    if(design == "NCC"){
    if(myAugPos[2]==4){
    tmp.mod <- locfit(vi~lp(xi, y1, nn=nn.par), data = tmp.data0, family = "binomial", link = "logit")
		}else{
		  tmp.mod <- locfit(vi~lp(xi, zi, nn = nn.par), data = tmp.data0, family = "binomial", link = "logit")

		}

		tmp.data1$wi = 1
		tmp.data0$wi <- 1/fitted(tmp.mod, data = tmp.data0)

    }else if(design == "CCH"){
      if(myAugPos[1]==4){

        tmp.mod <- locfit(vi~lp(y1, nn=nn.par), data = tmp.data0, family = "binomial", link = "logit")
        tmp.data0$wi <- 1/fitted(tmp.mod, data = tmp.data0)

        tmp.mod <- locfit(vi~lp(y1, nn=nn.par), data = tmp.data1, family = "binomial", link = "logit")
        tmp.data1$wi <- 1/fitted(tmp.mod, data = tmp.data1)

      }else{
        tmp.mod <- locfit(vi~lp(zi, nn = nn.par), data = tmp.data0, family = "binomial", link = "logit")
        tmp.data0$wi <- 1/fitted(tmp.mod, data = tmp.data0)

        tmp.mod <- locfit(vi~lp(zi, nn = nn.par), data = tmp.data1, family = "binomial", link = "logit")
        tmp.data1$wi <- 1/fitted(tmp.mod, data = tmp.data1)
      }
    }


		cohortdata <- rbind(tmp.data0, tmp.data1)
		subdata = cohortdata[cohortdata$vi==1,]

    #print( sum(subdata$wi))
   }


   junk = GetRTdata.SEMIPW.sorted.FUN(subdata, predict.time)  ## data.RT and data are sorted by linear predictor
   data.RT   <- junk$data.RT

   cutoff = data.RT[,1]
   ncutoff = length(cutoff)
   cutoff = sort(cutoff)[-c(1:4,(ncutoff-4):ncutoff)]
   #cutoff = sort(data.RT[,1])
   cutpos = sum.I(cutoff,">=", data.RT[,1])

   if (sum(cutpos<=0)>0) {
   cutoff = cutoff[-c(1:sum(cutpos<=0)) ]}
   RT.out    <- EstRTall.fixcutoff.FUN(cutoff,data.RT)

  betahat      <- junk$beta

#  TG        <- RT.out[[2]]
  rho       <- RT.out[[3]]
  AUC       <- RT.out$AUC
  IDI       <- RT.out$IDI
  ITPR      <- RT.out$ITPR
  IFPR      <- RT.out$IFPR

  RT.out    <- RT.out$RT.out.c
  RTvp.out  <- NULL


  if (!is.null(vp)) {
  		nvp = length(vp)
  		RTvp.out = rep(0,nvp)
    	for (pp in 1:nvp) {
			RTvp.out[pp] = EstRTvp.FUN(RT.out,vp[pp],typex[pp],typey[pp])
		}
  }

  if (!is.null(vv)) {
	   pcf = EstRTvp.FUN(RT.out,vv,'v','TPR')
	   pnf = EstRTvp.FUN(RT.out,mytpr,'TPR','v')
  }


NBp.out = NULL

if (!is.null(NBp)) {

  for (jj in 1:length(NBp)) {

    NBp.out=rbind(NBp.out,EstNB.FUN(RT.out,NBp[jj],rho))
  }
}


  est = c(betahat, AUC, IDI, ITPR,IFPR,RTvp.out,pcf,pnf, NBp.out )

  if (CalVar)  {
  	    subdata = cbind(junk$data,data.RT[,c(2,1)])
  		names(subdata)=c("times","status","zi","y1","y2","vi","weights","Sy","linearY")
        jjunk = Est.Wexp(cutoff,subdata,N,RT.out,predict.time,vp,typex,typey,vv,mytpr)
        Wexp = cbind(jjunk$Wexp.beta,jjunk$Wexp.AUC,jjunk$Wexp.IDI,jjunk$Wexp.ITPR,jjunk$Wexp.IFPR,jjunk$Wexp.vp,jjunk$Wexp.pcf,jjunk$Wexp.pnf)
        if (method == 'TrueW') {
        	if (design == 'CCH') {

        		sd = sqrt(Est.Var.CCH.trueweights(N,Wexp,subdata,subdata$status)[[3]])
        	}
        	if (design == 'NCC') {
        		 tk = sort(subdata$times[cohortdata$status == 1])
            	 pi.tk = sum.I(tk, "<=", cohortdata$xi)/N
        		sd = sqrt(Est.Var.NCC.trueweights(N,Wexp,subdata,subdata$status,tk,pi.tk)[[3]])
        	}
        }

        if (method == 'TrueWZ') {
        	    stratum = subdata$zi
        	    stratum = ifelse(subdata$status==1,stratum+3,stratum)
              tmpestvar <- Est.Var.CCH.trueweights(N,Wexp,subdata,stratum)
        		sd = sqrt(tmpestvar[[3]])


        }

		if (method == 'AugW') {
			if (design =='CCH')
        	{


        sd = sqrt(Est.Var.CCH.augweights(N,Wexp,subdata,AugPos,subdata$status)[[3]]) }
        	if (design =='NCC')
        	{sd = sqrt(Est.Var.NCC.augweights(N,Wexp,subdata,AugPos)[[3]]) }
        	     	    }

   c(est,sd,jjunk$err)
  } else {c(est)}
}

## function for calculating expensions

###########################################################################
#  Estimates the base expansion from which the variance
#     estimates will be made; This is for Sty, St, and Fc.
#
#    N is the cohort size
#    data is the dataset, and must have columns:
#        "Y", "status", "times", "linearY" and "Sy" from Est.Sy.SEMIPW
#    sometimes with 'zi' the matching variable and 'vi' used
#      SHOULD ALSO HAVE
#        "weights" if relevant...filled in with 1's below
#    predict.time is a scalar time used to divide cases from controls

#
Est.Wexp<-
function(cutoff,data,N,RT.out,predict.time,uu0Vec,typexVec,typeyVec,vv=NULL,mytpr=NULL) {

  if(missing(data))      { stop("Est.Wexp0: data not specified") }

  if( !("status" %in% names(data)) )  { stop(sprintf(errFormat,"status")) }
  if( !("times" %in% names(data)) )  { stop(sprintf(errFormat,"times")) }


  if(!"weights" %in% names(data)) {
    data$weights=1
  }

  # First, fit the survival model
    data = data[order(data$linearY),]   ## it is sorted it before this function;
    #data = data[order(data$Sy),]
    Y  <- as.matrix(data[,!is.element(names(data), c("times", "zi", "status", "weights", "vi","Sy","linearY"))])

    cutpos = sum.I(cutoff,">=", data$linearY)
    ncut = length(cutpos)
    nr = nrow(data)
    np = dim(Y)[2]

    fit  = coxph(Surv(data$times,data$status)~Y,
                 method="breslow", weight=data$weights)

    # Doing riskmat, haz0 and time by hand since coxph.detail appears
    #  to be a newer R feature & some users may have not updated their R.
    #    Note: this hazard is frequently normalized,
    #    by multiplying by exp(mean(data$Y)*fit$coef), but that is
    #    not necessary here, as our haz0 below doesn't want it.
    rrk      =  exp(data$linearY)
    dataD    =  subset(data[order(data$times),],status==1)
    riskmat  =  t(sapply(data$times,function(x) x >= dataD$times))

    s0   = t(riskmat) %*% (rrk*data$weights) ## length of nt
    s1   = t(riskmat) %*% t(VTM(rrk*data$weights,np)*t(Y))  ## nt *np
    haz0      = dataD$weights / apply(riskmat*rrk*data$weights,2,sum)
    cumhaz0   = cumsum(haz0)
    cumhaz.t0 = sum.I(predict.time, ">=", dataD$times, haz0)
   ## CondSyk   = exp(-cumhaz.t0*rrk) ## check it is the same as Sy
  # CondSmyY = exp(-cumhaz.t0*exp(myY%*%fit$coef))
    tmpind    = (data$times<=predict.time)&(data$status==1)
    tmpind.t  = sum.I(data$times[tmpind], ">=", dataD$times)

    resid.sco         = resid(fit, type="score")
    Wexp.beta         = resid.sco %*% fit$var * N
    Wexp.Lam1         = rep(0, nr)

    Wexp.Lam1[tmpind] = N/s0[tmpind.t]
    Wexp.Lam1 = Wexp.Lam1 -
               sum.I(pmin(predict.time,data$times), ">=", dataD$times,
               haz0/s0)*rrk*N
    Wexp.Lam2 = Wexp.beta %*% sum.I(predict.time, ">=", dataD$times,
                  haz0*s1/t(VTM(s0,np)))
    Wexp.Lam  = Wexp.Lam1 - Wexp.Lam2

    # end of most basic expansions...
    #    next calcs are derived expansions of performance measures
   # Fyk  = sum.I(data$Sy, ">=", data$Sy, data$weights)/sum(data$weights)  ##Fyk = P(linearY <c)
    Fyk  = sum.I(data$linearY, ">=", data$linearY, data$weights)/sum(data$weights)  # Fyk is the distribution of linear predictor under the cox model, same if use -Sy but not Sy  ## P(linearY<= cutoff)
    #Fyk   = rank(data$Y,ties="max")/nrow
    dFyk = Fyk - c(0,Fyk[-nr])

    St0.Fyk   = cumsum(data$Sy*dFyk)  ## St0.Fyk = P(T> t0,linearY<=c) for c at each linearY
    St0       = max(St0.Fyk)          ## St0 = P(T>t0)
    St0.Syk   = St0-St0.Fyk           ## St0.Syk = P(T>t0,Sy>c) for c at each linearY

   ## Wexp.Cond.Stc = -VTM(data$Sy*rrk, nr)*(t(VTM(Wexp.Lam,nr))+cumhaz.t0*Wexp.beta%*%t(Y)) ## iid expansion of Sy at all available Y;
    Wexp.Cond.Stc = -VTM(data$Sy*rrk, nr)*t(VTM(Wexp.Lam,nr))   ## for beta0
    Wexp.Cond.Stc.withB = -VTM(data$Sy*rrk, nr)*(t(VTM(Wexp.Lam,nr))+cumhaz.t0*Wexp.beta%*%t(Y)) ## for betahat
    Wexp.Stc  = t(sum.I(cutoff, "<", data$linearY, t(Wexp.Cond.Stc)*dFyk)) +
                data$Sy*(data$linearY > VTM(cutoff,nr)) -
                VTM(St0.Syk[cutpos], nr)     ## iid expansion of St0.Syk

    Wexp.St   = colSums(t(Wexp.Cond.Stc)*dFyk) + data$Sy - St0  ## iid expansion of St0, which is the same as the first column of Wexp.Stc: Wexp.Stc[,1] = Wexp.St

    Wexp.Fc   = 1*((data$linearY <= VTM(cutoff,nr)) - VTM(Fyk[cutpos],nr))  ## iid expansion of Fyc = P(Sty<c) for c known; nr*ncut
    ## below I calculate the derivative of Fc with respect to beta


  #  list(Wexp.Cond.Stc  = Wexp.Cond.Stc, Wexp.Stc  = Wexp.Stc, Wexp.St = Wexp.St, Wexp.Fc = Wexp.Fc) }

  ## Assemble for classic performance measures: given linear predictor;
   Wexp.all = as.list(1:8);
   names(Wexp.all)=c("RiskT","v","FPR","TPR","rho","NPV","PPV","RiskTWithB")
   Wexp.all[[1]] = -Wexp.Cond.Stc[,cutpos]; ## iid expansion of conditional risk: Fty = P(T<t|Y) for each person n*n ### note in this version this does not include variation in beta;
   Wexp.all[[2]] = Wexp.Fc; ## iid expansion of v=Fyc = P(linearY<c)
   Wexp.all[[3]] = ( - Wexp.St * VTM(RT.out[,4], nr) + Wexp.Stc)/St0
   Wexp.all[[4]] = (Wexp.St*VTM(RT.out[,5], nr) -Wexp.Fc - Wexp.Stc)/(1-St0)
   Wexp.all[[5]] = -Wexp.St;  ## iid expansion of P(T>t)
   Wexp.all[[6]] = (Wexp.St-Wexp.Stc-VTM(RT.out[,6], nr)*Wexp.Fc)/VTM(Fyk[cutpos],nr)
   Wexp.all[[7]]= (VTM(RT.out[,7]-1, nr)*Wexp.Fc-Wexp.Stc)/VTM(1-Fyk[cutpos],nr)
   Wexp.all[[8]] = Wexp.Cond.Stc.withB[,cutpos];

   # Wexp.all[[3]] = ( - Wexp.St * VTM(RT.out[,4], nr) + Wexp.Stc)/St0 +Wexp.beta%*%DFPR.beta  ## iid expansion of P(Lphat>c|T>t)
 # Wexp.all[[4]] = (Wexp.St*VTM(RT.out[,5], nr) -Wexp.Fc - Wexp.Stc)/(1-St0) + Wexp.beta%*%DTPR.beta  ## iid expansion of P(Lphat>c|T<t)
 #  Wexp.all[[6]] = (Wexp.St-Wexp.Stc-VTM(RT.out[,6], nr)*Wexp.Fc)/VTM(Fyk[cutpos],nr)+ Wexp.beta%*%DNPV.beta
 #  Wexp.all[[7]]= (VTM(RT.out[,7]-1, nr)*Wexp.Fc-Wexp.Stc)/VTM(1-Fyk[cutpos],nr)+ Wexp.beta%*%DPPV.beta

   np = length(fit$coef)
 	ncut = length(cutoff)
    newbeta1=newbeta2=rep(0,np)
    DFc.beta = 	DSt.beta =DTPR.beta = DFPR.beta= DPPV.beta=DNPV.beta=DF.beta=matrix(0,np,ncut)
    DAUC.beta = DITPR.beta = DIFPR.beta = DIDI.beta = DPCF.beta = DPNF.beta = DSt.beta = rep(0,np)
    DRTvp.beta = matrix(0,np,length(vp))
    err = 0;
    for (j in c(1:np)) {
    	newbeta = newbeta1 = newbeta2 = fit$coef
    	newbeta1[j] = newbeta[j] + delta.beta[j]
    	newbeta2[j] = newbeta[j] - delta.beta[j]

    	junk1 = Cal.All.Beta(data,cumhaz.t0,cutoff,newbeta1,Y,vv,mytpr,vp,typex,typey)
    	junk2 = Cal.All.Beta(data,cumhaz.t0,cutoff,newbeta2,Y,vv,mytpr,vp,typex,typey)
    	DTPR.beta[j,] = (junk1$TPR-junk2$TPR)/(2*delta.beta[j])
    	DTPR.beta[j,] =predict(loess(DTPR.beta[j,]~cutoff),cutoff)

       	DFPR.beta[j,] = (junk1$FPR-junk2$FPR)/(2*delta.beta[j])
       	DFPR.beta[j,] =predict(loess(DFPR.beta[j,]~cutoff),cutoff)

       	DPPV.beta[j,] = (junk1$PPV-junk2$PPV)/(2*delta.beta[j])
    	DPPV.beta[j,] =predict(loess(DPPV.beta[j,]~cutoff),cutoff)

    	DNPV.beta[j,] = (junk1$NPV-junk2$NPV)/(2*delta.beta[j])
    	DNPV.beta[j,] =predict(loess(DNPV.beta[j,]~cutoff),cutoff)

        DFc.beta[j,] = (junk1$Fyk-junk2$Fyk)/(2*delta.beta[j])
    	DFc.beta[j,] =predict(loess(DF.beta[j,]~cutoff),cutoff)
        err = err + ifelse(length(junk1$TPR)==length(junk2$TPR),0,1)

        DSt.beta[j] = (junk1$St-junk2$St)/(2*delta.beta[j])
        DAUC.beta[j] = (junk1$AUC-junk2$AUC)/(2*delta.beta[j])
        DIDI.beta[j] = (junk1$IDI-junk2$IDI)/(2*delta.beta[j])
        DITPR.beta[j] = (junk1$ITPR-junk2$ITPR)/(2*delta.beta[j])
        DIFPR.beta[j] = (junk1$IFPR-junk2$IFPR)/(2*delta.beta[j])
        DPCF.beta[j] = (junk1$PCF-junk2$PCF)/(2*delta.beta[j])
        DPNF.beta[j] = (junk1$PNF-junk2$PNF)/(2*delta.beta[j])
        DRTvp.beta[j,] = (junk1$RTvp-junk2$RTvp)/(2*delta.beta[j])
    }

   Wexp.ITPR = cbind(0,Wexp.all[[4]])%*%(c(RT.out$RiskT,1)-c(0,RT.out$RiskT))+
	              (cbind(Wexp.all[[1]],0)-cbind(0,Wexp.all[[1]]))%*%c(1,RT.out$TPR)+Wexp.beta%*%DITPR.beta
   Wexp.IFPR = cbind(0,Wexp.all[[3]])%*%(c(RT.out$RiskT,1)-c(0,RT.out$RiskT))+
	              (cbind(Wexp.all[[1]],0)-cbind(0,Wexp.all[[1]]))%*%c(1,RT.out$FPR)+Wexp.beta%*%DIFPR.beta
   Wexp.IDI= Wexp.ITPR - Wexp.IFPR
   Wexp.AUC = cbind(0,Wexp.all[[4]])%*%(c(1,RT.out$FPR)-c(RT.out$FPR,0))+
	             (cbind(0,Wexp.all[[3]])-cbind(Wexp.all[[3]],0))%*%c(1,RT.out$TPR)+Wexp.beta%*%DAUC.beta
   nvp = length(uu0Vec)
   Wexp.vp  = matrix(0,nr,nvp)
   for (pp in 1:nvp) {
   		uu0 = uu0Vec[pp]
   		if (typexVec[pp] == "cutoff") {
   			uuk = sort(RT.out[,1]);
   			tmpind = sum.I(uu0,">=",uuk)
   			ind0.y = match(typeyVec[pp],c("RiskT","v","FPR","TPR","rho","NPV","PPV","RiskTWithB"))
   			if (ind0.y<8) {
   			myD = cbind(DFc[,tmpind],DFPR[,tmpind],DTPR[,tmpind],DSt,DNPV[,tmpind],DPPV[,tmpind])
   			Wexp.vp[,pp] = Wexp.all[[ind0.y]][,tmpind] +Wexp.beta%*%myD[,(ind0.y-1)]} else { Wexp.vp[,pp] = Wexp.all[[ind0.y]][,tmpind] }
   		} else {
   			ind0.x = match(typexVec[pp],c("cutoff","RiskT","v","FPR","TPR","NPV","PPV"))
   			ind0.y = match(typeyVec[pp],c("cutoff","RiskT","v","FPR","TPR","NPV","PPV"))
   		 ### need to add check for permissible values
    		dYY.hat = dAcc.FUN(uu0, uu=RT.out[,ind0.x], A.uu=RT.out[,ind0.y],bw=NULL)
    		uuk = RT.out[,ind0.x]
            tmpind <- which.min(abs(uuk - uu0))[1]

           #  uuk = sort(RT.out[,ind0.x])
   		   #	tmpind = sum.I(uu0,">=",uuk)
      		ind0.x = match(typexVec[pp],c("RiskT","v","FPR","TPR","rho","NPV","PPV","RiskTWithB"))
      		ind0.y = match(typeyVec[pp],c("RiskT","v","FPR","TPR","rho","NPV","PPV","RiskTWithB"))
	  		Wexp.vp[,pp] = Wexp.all[[ind0.y]][,tmpind]-dYY.hat*Wexp.all[[ind0.x]][,tmpind] + Wexp.beta%*%DRTvp.beta[,pp]
   		}
   }
   Wexp.pcf = Wexp.pnf = NULL
   if (!is.null(vv)) {

   	    dYY.hat = dAcc.FUN(vv, uu=RT.out[,3], A.uu=RT.out[,5],bw=NULL)
	    uuk = RT.out[,3]
        tmpind <- which.min(abs(uuk - vv))[1]
        Wexp.pcf = Wexp.all[[4]][,tmpind]-dYY.hat*Wexp.all[[2]][,tmpind]+Wexp.beta%*%DPCF.beta
        uuk = RT.out[,5]
        tmpind <- which.min(abs(uuk - mytpr))[1]
 		dYY.hat = dAcc.FUN(mytpr, uu=RT.out[,5], A.uu=RT.out[,3],bw=NULL)
		Wexp.pnf = Wexp.all[[2]][,tmpind]-dYY.hat*Wexp.all[[4]][,tmpind]+Wexp.beta%*%DPNF.beta
  }

   list(err=err, Wexp.beta = Wexp.beta, Wexp.AUC = Wexp.AUC,Wexp.IDI=Wexp.IDI,Wexp.ITPR=Wexp.ITPR,Wexp.IFPR=Wexp.IFPR,Wexp.vp=Wexp.vp,Wexp.pcf=Wexp.pcf,Wexp.pnf=Wexp.pnf)
 }

 Est.Var.CCH.trueweights = function(N,Wexp,data,stratum) {
 	Wexp = as.matrix(Wexp)
 	cohort.variance = colSums(data$weights*(Wexp/N)^2)
 	robust.variance = colSums((data$weights*Wexp/N)^2)
 	## strata by cases and conrols
 	stra = sort(unique(stratum))
	nstra=length(stra);
	np = dim(Wexp)[2]
	strvar = rep(0,np);
	for (i in 1:nstra) {
		straWt = data$weights[stratum==stra[i]]
		straWexp = as.matrix(Wexp[stratum==stra[i],])
		ns = length(straWt)
		tempstratavar = (ns-1)*(straWt[1]-1)*straWt[1]*apply(straWexp/N,2,var)
		strvar = strvar + tempstratavar
	}
	list(cohort.variance=cohort.variance,robust.variance=robust.variance,variance = cohort.variance+strvar)
}

Est.Var.NCC.trueweights = function(N, Wexp,data,stratum,tk,pi.tk) {
	Wexp = as.matrix(Wexp)
 	cohort.variance = colSums(data$weights*(Wexp/N)^2)
 	robust.variance = colSums((data$weights*Wexp/N)^2)
    adj.tk = as.matrix(sum.I(tk,"<=",data$times,data$weights*Wexp/N*(data$weights-1))/N)
    var.adj = apply(adj.tk^2/pi.tk^2,2,sum)
    list(cohort.variance=cohort.variance,robust.variance=robust.variance, variance = robust.variance-var.adj*nmatch)
}

EstNB.FUN<-function(RT.out,NBp,rho) {

  FPR.p = EstRTvp.FUN(RT.out, NBp, "RiskT", "FPR")
  TPR.p = EstRTvp.FUN(RT.out,NBp,"RiskT","TPR")

  rho*TPR.p-NBp/(1-NBp)*(1-rho)*FPR.p
}
#####################3

#################
##  for one dimension smoothing: try my own smooth function similar to the above;
Est.Var.CCH.augweights = function(N,Wexp,data,posXZ,stratum) {
 	Wexp = as.matrix(Wexp)
 	cohort.variance = colSums(data$weights*(Wexp/N)^2)
 	robust.variance = colSums((data$weights*Wexp/N)^2)
 	## strata by cases and conrols
 	stra = sort(unique(stratum))
	nstra=length(stra);
	strvar =0;
	np = dim(Wexp)[2]
	strvar = rep(0,np); straWexp=NULL;
	for (i in 1:nstra) {
		straWt = data$weights[stratum==stra[i]];		straWexp = as.matrix(Wexp[stratum==stra[i],])
		ns = length(straWt)
		straX = data[stratum==stra[i],posXZ]

        bw=1.06*min(sd(straX),IQR(straX)/1.34)*ns^(-bw.power)
        junk = Kern.FUN(straX,straX,bw)
        junk = junk/VTM(colSums(junk),dim(junk)[2])
        straWexp = straWexp/N - junk%*%straWexp/N
		tempstratavar = colSums(straWt*(straWt-1)*(straWexp-colSums(straWexp)/ns)^2)
		#tempstratavar = colSums(straWt*(straWt-1)*(straWexp)^2)
		strvar = strvar + tempstratavar
	}
	list(cohort.variance=cohort.variance,robust.variance=robust.variance,variance = cohort.variance+strvar)
}

Kern.FUN <- function(zz,zi,bw,kern0="gauss") ## returns an (n x nz) matrix ##
{
out = (VTM(zz,length(zi))- zi)/bw
switch(kern0,
"epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
"gauss"= dnorm(out)/bw)
}

Est.Var.CCH.augweights = function(N,Wexp,data,posXZ,stratum){
 	Wexp = as.matrix(Wexp)
 	cohort.variance = colSums(data$weights*(Wexp/N)^2)
 	robust.variance = colSums((data$weights*Wexp/N)^2)
 	## strata by cases and conrols
 	stra = sort(unique(stratum))
	nstra=length(stra);
	strvar =0;
	np = dim(Wexp)[2]
	strvar = rep(0,np); straWexp=NULL;
	for (i in 1:nstra) {
		straWt = data$weights[stratum==stra[i]];		straWexp = as.matrix(Wexp[stratum==stra[i],])
		ns = length(straWt)
		straX = data[stratum==stra[i],posXZ]
		for (j in 1:np) {


			#musR = sm.regression(straX,straWexp[,j]/N,eval.points=straX,eval.grid=FALSE,display = "none")$estimate

      musR = fitted(locfit(straWexp[,j]/N~lp(straX, nn = nn.par)))
      straWexp[,j] = straWexp[,j]/N-musR

	    		    }
		tempstratavar = colSums(straWt*(straWt-1)*(straWexp-colSums(straWexp)/ns)^2)
		#tempstratavar = colSums(straWt*(straWt-1)*(straWexp)^2)
		strvar = strvar + tempstratavar
	}
	list(cohort.variance=cohort.variance,robust.variance=robust.variance,variance = cohort.variance+strvar)
}

################
##################
Est.Var.NCC.augweights = function(N,Wexp,data,posXZ) {
 	Wexp = as.matrix(Wexp)
 	cohort.variance = colSums(data$weights*(Wexp/N)^2)
 	robust.variance = colSums((data$weights*Wexp/N)^2)
 	## strata by cases and conrols

	np = dim(Wexp)[2]
	strvar = rep(0,np); straWexp=NULL;
		straWt = data$weights[data$status==0];
		straWexp = as.matrix(Wexp[data$status==0,])
		ns = length(straWt)
		straX = data[data$status==0,posXZ]
		for (j in 1:np) {

      tmp <- data.frame("Y" = straWexp[,j]/N, straX)
		  musR = fitted(locfit(straWexp[,j]/N~lp(times, y1, nn = nn.par), data = tmp), data = tmp)

		#	musR = sm.regression(straX,straWexp[,j]/N,eval.points=straX,eval.grid=FALSE,display = "none")$estimate
	    	straWexp[,j] = straWexp[,j]/N-musR
	    }
		strvar = colSums(straWt*(straWt-1)*(straWexp-colSums(straWexp)/ns)^2)
		#tempstratavar = colSums(straWt*(straWt-1)*(straWexp)^2)


	list(cohort.variance=cohort.variance,robust.variance=robust.variance,variance = cohort.variance+strvar)
}





########################################
#### function for point estimation


GetRTdata.SEMIPW.sorted.FUN = function(data,predict.time) {


  Sy = EST.Sy.SEMIPW(data,predict.time)

  beta    <- Sy[[1]]
  linearY <- Sy[[3]]
  ooo = order(linearY)
  data.RT <- cbind(linearY, Sy[[2]], data$wi)[ooo,]

  Fck <- sum.I(data.RT[,1], ">=", data.RT[,1], data.RT[,3])/sum(data.RT[,3])

  data.RT <- cbind(data.RT[,-c(3)],Fck)

  list(beta = beta, data.RT = data.RT, data=data[ooo,])
}


EST.Sy.SEMIPW<-function(data, predict.time) {
  Y.old  <- as.matrix(data[,!is.element(names(data), c("di", "zi", "xi", "wi", "vi"))])

  data <- data[order(-data$xi),]
  Y  <- as.matrix(data[,!is.element(names(data), c("di", "zi", "xi", "wi", "vi"))])

  fit<-coxph( Surv(data$xi,data$di) ~ Y, weight = data$wi)

  beta<-fit$coef
  linearY <- Y%*%beta

  r.riskset <-data$wi/cumsum(exp(linearY)*data$wi)

  Lambda0t  <- sum(r.riskset[(data$xi <= predict.time)&(data$di == 1)])

  linearY <- Y.old%*%beta
  Sy <- exp(-Lambda0t*exp(linearY))

  list(beta,Sy, linearY)
}

Cal.All.Beta=function(data,cumhaz.t0,cutoff,beta,Y,v,tpr,vp=NULL,typex=NULL,typey=NULL) {
 	    linearY = c(Y%*%beta)
 		rrk     =  exp(linearY)
    	Fyk  = sum.I(linearY, ">=", linearY, data$weights)/sum(data$weights)
    	nY = length(linearY)
    	dFyk = Fyk - c(0,Fyk[-nY])
    	CondSyk   = exp(-cumhaz.t0*rrk)
    	St0.Fyk   = cumsum(CondSyk*dFyk)  ## St0.Fyk = P(T> t0,Sy<=c)
    	St0       = max(St0.Fyk)          ## St0 = P(T>t0)
    	Ft0     = 1-St0
    	Ft0.Fyk = Fyk-St0.Fyk


    	RiskT = 1-CondSyk

    	cutpos = sum.I(cutoff,">=", linearY)
    	FPR   = (St0-St0.Fyk)/St0     ## P(Y> ck|T> t0)  values at empirical ck
        TPR   = (Ft0-Ft0.Fyk)/Ft0     ## P(Y> ck|T<=t0)
        NPV = St0.Fyk/Fyk
    	PPV = (Ft0-Ft0.Fyk)/(1-Fyk)

    	tmpind <- which.min(abs(Fyk - v))[1]
    	pcf = TPR[tmpind]
    	tmpind <- which.min(abs(TPR - tpr))[1]
    	pnf = Fyk[tmpind]

		FPR = FPR[cutpos]
        TPR = TPR[cutpos]
        PPV=PPV[cutpos]
        NPV=NPV[cutpos]
    	ITPR = sum(c(1,TPR)*(c(RiskT[cutpos],1)-c(0,RiskT[cutpos])))
  		IFPR = sum(c(1,FPR)*(c(RiskT[cutpos],1)-c(0,RiskT[cutpos])))
  		IDI  = ITPR - IFPR
  		AUC = sum(c(1,TPR)*(c(1,FPR)-c(FPR,0)))

  		RT.out= data.frame("cutoff" = cutoff,"RiskT" = 1-CondSyk[cutpos],"v"= Fyk[cutpos], "FPR" = FPR,"TPR" = TPR,"NPV" = NPV, "PPV" = PPV)

  		if (!is.null(vp)) {
  			nvp = length(vp)
  			RTvp.out = rep(0,nvp)
    		for (pp in 1:nvp) {
				RTvp.out[pp] = EstRTvp.FUN(RT.out,vp[pp],typex[pp],typey[pp])
			}
  		}

  		list(TPR=TPR,FPR=FPR,PPV=PPV,NPV=NPV,Fyk=Fyk[cutpos],St = St0, IDI = IDI,ITPR = ITPR, IFPR=IFPR,AUC=AUC,PCF=pcf,PNF = pnf,RTvp = RTvp.out)
 }

 Cal.All.Dbeta = function(data,cumhaz.t0,cutoff,betahat,delta.beta,Y,v,tpr,vp=NULL,typex=NULL,typey=NULL) {

 	np = length(betahat)
 	ncut = length(cutoff)
    newbeta1=newbeta2=rep(0,np)
    DFc.beta = 	DSt.beta =DTPR.beta = DFPR.beta= DPPV.beta=DNPV.beta=DF.beta=matrix(0,np,ncut)
    DAUC.beta = DITPR.beta = DIFPR.beta = DIDI.beta = DPCF.beta = DPNF.beta = DSt.beta = rep(0,np)
    DRTvp.beta = matrix(0,np,length(vp))
    err = 0;
    for (j in c(1:np)) {
    	newbeta = newbeta1 = newbeta2 = betahat
    	newbeta1[j] = newbeta[j] + delta.beta[j]
    	newbeta2[j] = newbeta[j] - delta.beta[j]

    	junk1 = Cal.All.Beta(data,cumhaz.t0,cutoff,newbeta1,Y,vv,mytpr,vp,typex,typey)
    	junk2 = Cal.All.Beta(data,cumhaz.t0,cutoff,newbeta2,Y,vv,mytpr,vp,typex,typey)
    	DTPR.beta[j,] = (junk1$TPR-junk2$TPR)/(2*delta.beta[j])
    	DTPR.beta[j,] =predict(loess(DTPR.beta[j,]~cutoff),cutoff)

       	DFPR.beta[j,] = (junk1$FPR-junk2$FPR)/(2*delta.beta[j])
       	DFPR.beta[j,] =predict(loess(DFPR.beta[j,]~cutoff),cutoff)

       	DPPV.beta[j,] = (junk1$PPV-junk2$PPV)/(2*delta.beta[j])
    	DPPV.beta[j,] =predict(loess(DPPV.beta[j,]~cutoff),cutoff)

    	DNPV.beta[j,] = (junk1$NPV-junk2$NPV)/(2*delta.beta[j])
    	DNPV.beta[j,] =predict(loess(DNPV.beta[j,]~cutoff),cutoff)

        DFc.beta[j,] = (junk1$Fyk-junk2$Fyk)/(2*delta.beta[j])
    	DFc.beta[j,] =predict(loess(DF.beta[j,]~cutoff),cutoff)
        err = err + ifelse(length(junk1$TPR)==length(junk2$TPR),0,1)

        DSt.beta[j] = (junk1$St-junk2$St)/(2*delta.beta[j])
        DAUC.beta[j] = (junk1$AUC-junk2$AUC)/(2*delta.beta[j])
        DIDI.beta[j] = (junk1$IDI-junk2$IDI)/(2*delta.beta[j])
        DITPR.beta[j] = (junk1$ITPR-junk2$ITPR)/(2*delta.beta[j])
        DIFPR.beta[j] = (junk1$IFPR-junk2$IFPR)/(2*delta.beta[j])
        DPCF.beta[j] = (junk1$PCF-junk2$PCF)/(2*delta.beta[j])
        DPNF.beta[j] = (junk1$PNF-junk2$PNF)/(2*delta.beta[j])
        DRTvp.beta[j,] = (junk1$RTvp-junk2$RTvp)/(2*delta.beta[j])
    }
    list(DTPR.beta=DTPR.beta,DFPR.beta=DFPR.beta,DPPV.beta=DPPV.beta,DNPV.beta=DNPV.beta,DF.beta=DF.beta,DSt.beta = DSt.beta,DIDI.beta = DIDI.beta,DITPR.beta = DITPR.beta, DIFPR.beta=DIFPR.beta,DAUC.beta=AUC.beta,DPCF.beta=DPCF.beta,DPNF.beta = DPNF.betea,DRTvp.beta = DRTvp.beta)
}


################ individual components
 Cal.pcf.pnf=function(data,cumhaz.t0,cutoff,beta,Y,v,tpr) {
 	    linearY = c(Y%*%beta)
 	    cutpos = sum.I(cutoff,">=", linearY)
 		rrk     =  exp(linearY)

    	CondSyk   = exp(-cumhaz.t0*rrk)

    	Fyk  = sum.I(linearY, ">=", linearY, data$weights)/sum(data$weights)

    	nY = length(linearY)
    	dFyk = Fyk - c(0,Fyk[-nY])
    	St0.Fyk   = cumsum(CondSyk*dFyk)  ## St0.Fyk = P(T> t0,Sy<=c)
    	St0       = max(St0.Fyk)          ## St0 = P(T>t0)

    	TPR = (1-St0-Fyk+St0.Fyk)/(1-St0)
    	tmpind <- which.min(abs(Fyk - v))[1]
    	pcf = TPR[tmpind]
    	tmpind <- which.min(abs(TPR - tpr))[1]
    	pnf = Fyk[tmpind]
    	list(PCF=pcf,PNF = pnf )
 }





EstRTall.fixcutoff.FUN<-function(cutoff,data)
{
  ck      = data[,1]   ## aka linearY
  nc      = length(ck)
  CondSck = data[,2]   ## aka Sy
  Fck     = data[,3]
  ncut = length(cutoff)

  dFck    = Fck - c(0,Fck[-nc])
  St0.Fck = cumsum(CondSck*dFck)## St0.Fck = P(T> t0,Y<=ck)
  Ft0.Fck = Fck-St0.Fck         ## Ft0.Fck = P(T<=t0,Y<=ck)
  St0     = max(St0.Fck)        ## St0     = P(T> t0      )
  Ft0     = 1-St0               ## Ft0     = P(T<=t0      )
  FPR.e   = (St0-St0.Fck)/St0     ## P(Y> ck|T> t0)  values at empirical ck
  TPR.e   = (Ft0-Ft0.Fck)/Ft0     ## P(Y> ck|T<=t0)
  NPV.e   = St0.Fck/Fck           ## P(T> t0|Y<=ck)
  PPV.e   = (Ft0-Ft0.Fck)/(1-Fck) ## P(T<=t0|Y> ck)

  cutpos = sum.I(cutoff,">=", ck)
  FPR.c   = FPR.e[cutpos]     ## P(Y> ck|T> t0)  values at empirical ck
  TPR.c   = TPR.e[cutpos]     ## P(Y> ck|T<=t0)
  NPV.c   = NPV.e[cutpos]          ## P(T> t0|Y<=ck)
  PPV.c   = PPV.e[cutpos] ## P(T<=t0|Y> ck)

#  RT.out  = data.frame("cutoff" = ck,"RiskT" = 1-CondSck,"v"= Fck, "FPR" = FPR.e,"TPR" = TPR.e, "NPV" = NPV.e,"PPV" = PPV.e,)
  RT.out.c  = data.frame("cutoff" = cutoff,"RiskT" = 1-CondSck[cutpos],"v"= Fck[cutpos], "FPR" = FPR.c,"TPR" = TPR.c,"NPV" = NPV.c, "PPV" = PPV.c)

  TG    = sum(abs(1-CondSck-Ft0)*(-Fck+c(Fck[-1],1)))
 # AUC   = sum(TPR.c*(FPR.c-c(FPR.c[-1],0)))  # original
 # AUC = sum(TPR.c*(FPR.c-c(FPR.c[-1],0)))+(1-FPR.c[1]) # for c
  AUC = sum(c(1,TPR.c)*(c(1,FPR.c)-c(FPR.c,0)))
  RiskT = 1-CondSck
  #mmm  = length(RiskT)
  #ITPR = sum(TPR.e*(RiskT-c(0,RiskT[-mmm])))
  #IFPR = sum(FPR.e*(RiskT-c(0,RiskT[-mmm])))
  #ITPR = sum(TPR.c*(RiskT[cutpos]-c(0,RiskT[cutpos][-ncut])))
  #IFPR = sum(FPR.c*(RiskT[cutpos]-c(0,RiskT[cutpos][-ncut])))
  ITPR = sum(c(1,TPR.c)*(c(RiskT[cutpos],1)-c(0,RiskT[cutpos])))
  IFPR = sum(c(1,FPR.c)*(c(RiskT[cutpos],1)-c(0,RiskT[cutpos])))
  IDI  = ITPR - IFPR

 # list(RT.out = RT.out, RT.out.c = RT.out.c, TG = TG, rho = Ft0, AUC = AUC, IDI = IDI, ITPR = ITPR , IFPR = IFPR)
 list(RT.out.c = RT.out.c, TG = TG, rho = Ft0, AUC = AUC, IDI = IDI, ITPR = ITPR , IFPR = IFPR)
}

EstRTvp.FUN<-function(RT.out, uu0, typex, typey) {

  ind0.x = match(typex,c("cutoff","RiskT","v","FPR","TPR","NPV","PPV"))
  ind0.y = match(typey,c("cutoff","RiskT","v","FPR","TPR","NPV","PPV"))

  uuk = RT.out[,ind0.x]
  tmpind <- which.min(abs(uuk - uu0))[1]
  RT.out[tmpind, ind0.y]
}


GetRTdata.SEMIPW.FUN<-function(data,predict.time) {

  Sy = EST.Sy.SEMIPW(data,predict.time)

  beta    <- Sy[[1]]
  linearY <- Sy[[3]]

  data.RT <- cbind(linearY, Sy[[2]], data$wi)[order(linearY),]

  Fck <- sum.I(data.RT[,1], ">=", data.RT[,1], data.RT[,3])/sum(data.RT[,3])

  data.RT <- cbind(data.RT[,-c(3)],Fck)

  list(beta = beta, data.RT = data.RT)
}

sum.I<-function(yy,FUN,Yi,Vi=NULL){
	if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
	pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
	if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
	if (!is.null(Vi)) {
      	if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
		Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
		return(rbind(0,Vi)[pos+1,])
	} else return(pos)
}

VTM <- function(vc, dm)
{
    matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}

bw.power = 0.30
Kern.FUN <- function(zz,zi,bw,kern0="gauss") ## returns an (n x nz) matrix ##
{
out = (VTM(zz,length(zi))- zi)/bw
switch(kern0,
"epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
"gauss"= dnorm(out)/bw)
}



## Sampling Weights
P0HAT.CCH.Z.FUN <- function(cohortdata, type = 2)
{

  di <- cohortdata$di
  zi <- cohortdata$zi
  vi <- cohortdata$vi
  z.levels <- sort(unique(zi))


  p0hat <- rep(0, length(di))

  if(type == 1){


   for( myz in z.levels){

     p0hat[zi == myz] <- sum((zi==myz)*(vi==1))/sum(zi==myz)

   }
   #cases all have weight 1
   p0hat[di == 1] <- 1

  }else {

   for( myz in z.levels){
     #case
     p0hat[zi == myz & di == 1] <-    sum((zi==myz)*(di==1)*(vi==1))/sum((zi==myz)*(di==1))
     #p0hat[zi == myz & di == 1] <-    ncch.case.z[myz+1]/sum((zi==myz)*(di==1))
     #control
     p0hat[zi == myz & di == 0] <-    sum((zi==myz)*(di==0)*(vi==1))/sum((zi==myz)*(di==0))
     #p0hat[zi == myz & di == 0] <-    ncch.control.z[myz+1]/sum((zi==myz)*(di==0))
     }

  }

  p0hat
}
## Sampling Weights
#  include both finite and quota sample
P0HAT.NCC.FUN <- function(data, Rstar.size = NULL, nmatch)
{
  xi  <- data$xi
  di  <- data$di
  tj  <- xi[di==1];
  N   <- length(xi);

  ooo <- order(tj)
  tj  <- tj[ooo]

  if (is.null(Rstar.size)) {

    nj <- sum.I(tj, "<=", xi);
    wj <- c(1,cumprod(1-nmatch/(nj-1)))

  }else {

    Rstar.size <- Rstar.size[ooo];
    #wj = c(1,cumprod(1-1/(N/Rstar.size-1/nmatch)  ))
    wj <- c(1,cumprod(1-Rstar.size/(N-1)))  # results are very similar to above
  }

  1-wj[sum.I(xi,">",tj)+1]
}

SIM.NCC.subcase.FUN<- function(data0, quota=F, nmatch)
{

  data0=cohortdata
  nn <- nrow(data0)
  xi <- data0$xi
  di <- data0$di
  vi <- rep(0, nn)
  #zi <- data0$zi
  zi <- rep(1,nn)

  vi <- di #choose all cases

  ncase=sum(di)
  ind.case <- (1:nn)[di==1]
  vi[ind.case] = rbinom(ncase,1,0.1)
  err <- 0

  matchid      <- NULL
  matchgrp     <- NULL
  case.ind     <- NULL
  riskset.size.nz <- NULL
  riskset.size.z  <- NULL
  Rstar.size   <- NULL  #used for Borgan method

    ind.case = (1:nn)[vi==1&di==1]
    #for each case
    for(j in 1:sum(vi[ind.case])){

      tmpind <- ind.case[j]
      risksetind   <- (1:nn)[xi>xi[tmpind]]
      risksetind.z <- (1:nn)[(xi>xi[tmpind])&(zi==zi[tmpind])]

         # if riskset is empty, no control would be selected
         # if riskset size < nmatch, only select # of available

      err = err + 1*(length(risksetind)<nmatch)

      if(length(risksetind)>0){

        controlind <- sample(risksetind,min(nmatch,length(risksetind)))
	      vi[controlind] <- 1

        matchid      <- c(matchid, c(tmpind, controlind))
      	matchgrp     <- c(matchgrp, rep(j, length(controlind)+1))

        case.ind     <- c(case.ind, c(1,rep(0,length(controlind))))

        riskset.size.nz <- c(riskset.size.nz, length(risksetind)+1)
	      riskset.size.z  <- c(riskset.size.z, length(risksetind.z)+1)
      } #end if
    } #end for

  cohortdata    <- data0
  cohortdata$vi <- vi

  list(cohortdata   = cohortdata,
       matchinfo    = cbind(matchgrp,case.ind,matchid),
       riskset.size = cbind(riskset.size.nz,riskset.size.z),
       Rstar.size   = Rstar.size)
}

## Sampling Weights
#  include both finite and quota sample
P0HAT.NCC.subcase.FUN <- function(data, Rstar.size = NULL, nmatch,pi)
{
  xi  <- data$xi
  di  <- data$di
  tj  <- xi[di==1];
  N   <- length(xi);

  ooo <- order(tj)
  tj  <- tj[ooo]

  if (is.null(Rstar.size)) {

    nj <- sum.I(tj, "<=", xi);
    wj <- c(1,cumprod(1-nmatch/(nj-1)*pi))

  }else {

    Rstar.size <- Rstar.size[ooo];
    #wj = c(1,cumprod(1-1/(N/Rstar.size-1/nmatch)  ))
    wj <- c(1,cumprod(1-Rstar.size/(N-1)))  # results are very similar to above
  }

  1-wj[sum.I(xi,">",tj)+1]
}
##Calculate weights for NCC samples with or without matching on Z
P0HAT.NCC.Z.FUN <- function(data, Riskset.size, nmatch)
{
  xi <- data$xi
  di <- data$di
  zi <- data$zi
  tj <- xi[di==1]
  zj <- zi[di==1]

  ooo <- order(tj)
  tj  <- tj[ooo]
  zj  <- zj[ooo];

  p0hati <- rep(1, length(xi))

  Riskset.size <- Riskset.size[ooo];

  for (i in 1:length(xi)) {

    is.zi <- zj == zi[i]

    wj <- cumprod(1-(nmatch*is.zi/(Riskset.size-1)))

    if (sum(tj<xi[i])>0 & sum(tj<xi[i])<=length(wj) ) p0hati[i] = 1-wj[sum(tj<xi[i])]



  }

  p0hati
}

## Sampling Weights
#  include both finite and quota sample
P0HAT.NCC.FUN <- function(data, Rstar.size = NULL, nmatch)
{
  xi  <- data$xi
  di  <- data$di
  tj  <- xi[di==1];
  N   <- length(xi);

  ooo <- order(tj)
  tj  <- tj[ooo]

  if (is.null(Rstar.size)) {

    nj <- sum.I(tj, "<=", xi);
    wj <- c(1,cumprod(1-nmatch/(nj-1)))

  }else {

    Rstar.size <- Rstar.size[ooo];
    #wj = c(1,cumprod(1-1/(N/Rstar.size-1/nmatch)  ))
    wj <- c(1,cumprod(1-Rstar.size/(N-1)))  # results are very similar to above
  }

  1-wj[sum.I(xi,">",tj)+1]
}



dAcc.FUN <- function(uu0, uu=SE.yy, A.uu = Sp.yy, bw=NULL)
  {
    data = cbind(uu,A.uu); data=na.omit(data); uu=data[,1]; A.uu=data[,2]; n.uu = length(uu)
    A.uu = A.uu[order(uu)]; uu = sort(uu)
    if(is.null(bw)){bw=1.06*min(sd(uu),IQR(uu)/1.34)*n.uu^(-bw.power)}
    Ki.u0 = Kern.FUN(uu0, uu[-1], bw) ## n.uu x n.u0
    c(t(A.uu[-1]-A.uu[-n.uu])%*%Ki.u0)
  }
