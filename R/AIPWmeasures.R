#' Estimate standard measures of predictive accuracy for two-phase designs using AIPW
#'
#' Estimate measures of predictive accuracy using augmented inverse probability weights (ipw) or true ipw for two-phase biomarker validation studies
#'
#' @param time vector of time to event
#' @param event status indicator (1 for observed failure, 0 for censoring)
#' @param X matrix or data frame of covariates
#' @param subcohort indicator for selection into the subcohort (1 for selected, 0 for not selected)
#' @param aug.weights.x vector of variable to be used to derive augmented inverse probability weights.
#' @param risk.threshold vector of risk thresholds on the absolute risk scale used to calcuate cutoff-based summary measures.
#' @param landmark.time numeric value specifying landmark time at which to estimate the performance measures.
#' @param weight.method either "Aug" for augmented ipw or "True" for true ipw.
#' @param design either "CCH" for case-cohort or "NCC" for nested case control specifying the subcohort design used.
#' @param calculate.sd  should analytic standard errors be calculated.
#' @param smoothing.par nearest neighbor smoothing parameter used in locfit (default = 0.7)
#' @param pnf.threshold  numeric threshold value for which to calculate PNF(bb)
#' @param pcf.threshold  numeric threshold value for which to calculate PCF(vv)
#' @param ncc.nmatch For design = "NCC", specify the number of controls matched per case.
#'
#' @note variance calculations are unavailable for the NB measure.
#'
#' @return data.frame with estimates of AUC, IDI, ITPR, IFPR, TPR, FPR, PPV, NPV and net benefit (NB) and standard errors.
#'
#' @examples
#'
#'data(CCHsimdata)
#'
#'predict.time <- 0.75
#'
#'## augmented ipw
#'AIPWmeasures( time = CCHsimdata$xi, event = CCHsimdata$di,
#'              X = cbind(CCHsimdata$y1, CCHsimdata$y2),
#'             subcohort = CCHsimdata$vi,
#'              aug.weights.x =   CCHsimdata$y1,
#'              risk.threshold  =  c(.05, .3),
#'              landmark.time = predict.time,
#'              weight.method = 'Aug',
#'              design = "CCH",
#'              smoothing.par = 0.9,
#'              calculate.sd  = TRUE,
#'              pnf.threshold = 0.85,
#'              pcf.threshold = 0.8)
#'
#'
#'#simulated data from a ncc design with nmatch = 2
#'data("NCCsimdata")
#'AIPWmeasures( time = NCCsimdata$xi, event = NCCsimdata$di,
#'              X = cbind(NCCsimdata$y1, NCCsimdata$y2),
#'              subcohort = NCCsimdata$vi,
#'              aug.weights.x =   NCCsimdata$y1,
#'              risk.threshold  =  c(.01, .03),
#'              landmark.time = predict.time,
#'              weight.method = 'Aug',
#'              design = "NCC",
#'              smoothing.par = 0.9,
#'              calculate.sd  = TRUE,
#'              pnf.threshold = 0.85,
#'              pcf.threshold = 0.8,
#'              ncc.nmatch = 2)
#'
#' @import locfit
#' @import survival
#' @export
#'
AIPWmeasures <- function(time,event, X, subcohort,
                         aug.weights.x = NULL,
                         risk.threshold,
                         landmark.time,
                         weight.method = c("Aug", "True"),
                         design= c("CCH", "NCC"),
                         calculate.sd = TRUE,
                         smoothing.par = 0.7,
                         pnf.threshold = NULL,
                         pcf.threshold = NULL,
                         ncc.nmatch  = NULL)
{

  ###### checks

  weight.method <- match.arg(weight.method)
  design  <- match.arg(design)



  ###### end checks

  z = rep(0, length(time))

  cohortdata = data.frame(xi = time, di = event, zi = z, X, vi = subcohort )
  vp = risk.threshold

  typey = rep(c("TPR","FPR","PPV","NPV"), length(vp))
  typex = c(rep("RiskT",length(typey)))
  vp = rep(vp, rep(4, length(vp))) # four outcome measures

  NBp = risk.threshold
  predict.time = landmark.time
  method = weight.method
  N = nrow(cohortdata)
  nn.par = smoothing.par
  CalVar = calculate.sd
  mytpr = pnf.threshold # for PNF V(TPR = 0.85)
  vv = pcf.threshold # for PCF TPR(v = 0.8)
  AugWeightsX = aug.weights.x
  nmatch = ncc.nmatch

  if (method == 'True'){
    if (design == 'CCH')  {

      cohortdata$wi = 1/P0HAT.CCH.Z.FUN(cohortdata, type=2);
      subdata <- cohortdata[cohortdata$vi == 1,]
    }
    if (design=='NCC') {

      cohortdata$wi = ifelse(cohortdata$di==1,1 ,1/P0HAT.NCC.FUN(cohortdata, NULL, nmatch))
      subdata <- cohortdata[cohortdata$vi == 1,]
    }

    AugWeightsX.sub <- NULL
  }
  #print(sum(subdata$wi))


  if (method == "Aug") {

    #if (is.null(AugPos)) {print("Need the positions of the variables used in estimating weights")}
    tmp.data0 <- cohortdata[cohortdata$di==0,]
    tmp.data1 <- cohortdata[cohortdata$di==1,]
    #browser()
    #    input.sm <- tmp.data0[,AugPos]

    #		tmp.data0$wi = 1/sm.regression(input.sm,tmp.data0$vi,eval.points=input.sm,eval.grid=FALSE, display =   			"none")$estimate
    if(design == "NCC"){

     tmp.mod <- locfit(vi~lp(xi, AugWeightsX[event==0], nn=nn.par), data = tmp.data0, family = "binomial", link = "logit")


      tmp.data1$wi = 1
      tmp.data0$wi <- 1/fitted(tmp.mod, data = tmp.data0)

    }else{
     # design == "CCH"
        tmp.mod <- locfit(vi~lp(AugWeightsX[event==0], nn=nn.par), data = tmp.data0, family = "binomial", link = "logit")
        tmp.data0$wi <- 1/fitted(tmp.mod, data = tmp.data0)

        tmp.mod <- locfit(vi~lp(AugWeightsX[event==1], nn=nn.par), data = tmp.data1, family = "binomial", link = "logit")
        tmp.data1$wi <- 1/fitted(tmp.mod, data = tmp.data1)

      }


    augw0 <- aug.weights.x[cohortdata$di==0]
    augw1 <- aug.weights.x[cohortdata$di==1]
    AugWeightsX <- c(augw0, augw1)


    cohortdata <- rbind(tmp.data0, tmp.data1)

    subdata = cohortdata[cohortdata$vi==1,]
    AugWeightsX.sub <- AugWeightsX[cohortdata$vi==1]
}
    #print( sum(subdata$wi))



  junk = GetRTdata.SEMIPW.sorted.FUN(subdata, predict.time, AugWeightsX.sub)  ## data.RT and data are sorted by linear predictor
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


  #est = data.frame(c(betahat), AUC, IDI, ITPR,IFPR,RTvp.out,pcf,pnf, NB = NBp.out )
  out = data.frame(measure = c(paste0("coef.x", 1:length(betahat)),
                               "AUC", "IDI", "ITPR", "IFPR", typey, "PCF", "PNF", rep("NB", length(risk.threshold))),
                   threshold = c(rep(NA, length(betahat)), NA, NA, NA, NA, vp, pcf.threshold, pnf.threshold, risk.threshold),
                   estimate = c(c(betahat), AUC, IDI, ITPR,IFPR,RTvp.out,pcf,pnf, NBp.out ))



  if (CalVar)  {

    subdata = cbind(junk$data,data.RT[,c(2,1)])

    names(subdata)=c("times","status", "zi",paste0("y", 1:length(betahat)), "vi", "weights","Sy","linearY")

    #cutoff,data,N,RT.out,predict.time,uu0Vec,typexVec,typeyVec,
    #vv=NULL,mytpr=NULL, vp

    jjunk = Est.Wexp(cutoff = cutoff, data = subdata,N = N, RT.out = RT.out,
                     predict.time =predict.time,
                     uu0Vec = vp,typexVec = typex, typeyVec = typey,
                     vv = vv, mytpr = mytpr, vp = vp)
    Wexp = cbind(jjunk$Wexp.beta,jjunk$Wexp.AUC,jjunk$Wexp.IDI,jjunk$Wexp.ITPR,jjunk$Wexp.IFPR,jjunk$Wexp.vp,jjunk$Wexp.pcf,jjunk$Wexp.pnf)

    if (method == 'True') {
      if (design == 'CCH') {

        sd = sqrt(Est.Var.CCH.trueweights(N,Wexp,subdata,subdata$status)[[3]])
      }
      if (design == 'NCC') {
        tk = sort(subdata$times[cohortdata$status == 1])
        pi.tk = sum.I(tk, "<=", cohortdata$xi)/N
        sd = sqrt(Est.Var.NCC.trueweights(N,Wexp,subdata,subdata$status,tk,pi.tk, nmatch)[[3]])
      }
    }


    if (method == 'Aug') {
      if (design =='CCH')
      {
        sd = sqrt(Est.Var.CCH.augweights(N,Wexp,subdata,AugWeightX = junk$AugWeightX,subdata$status, nn.par = nn.par)[[3]])
        }else if(design =='NCC')
      {
        sd = sqrt(Est.Var.NCC.augweights(N,Wexp,subdata,AugWeightX = junk$AugWeightX, nn.par = nn.par)[[3]])
      }
    }

    #Wexp = cbind(jjunk$Wexp.beta,jjunk$Wexp.AUC,jjunk$Wexp.IDI,jjunk$Wexp.ITPR,jjunk$Wexp.IFPR,
    #jjunk$Wexp.vp,jjunk$Wexp.pcf,jjunk$Wexp.pnf)
    #cbind(c(betahat, AUC, IDI, ITPR,IFPR,RTvp.out,pcf,pnf, NBp.out), sd)

    out$sd = c(sd, rep(NA, length(risk.threshold)))

  }

  out
}
