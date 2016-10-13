
#simulated data from a case-cohort design

data(CCHsimdata)

predict.time <- 0.75

## augmented ipw
AIPWmeasures( time = CCHsimdata$xi, event = CCHsimdata$di,
              X = cbind(CCHsimdata$y1, CCHsimdata$y2),
               subcohort = CCHsimdata$vi,
                   aug.weights.x =   CCHsimdata$y1,
                   risk.threshold  =  c(.05, .3),
                   landmark.time = predict.time,
                   weight.method = 'Aug',
                   design = "CCH",
                   smoothing.par = 0.9,
                   calculate.sd  = TRUE,
                   pnf.threshold = 0.85,
                   pcf.threshold = 0.2)

# true ipw

AIPWmeasures( time = CCHsimdata$xi, event = CCHsimdata$di,
              X = cbind(CCHsimdata$y1, CCHsimdata$y2),
              subcohort = CCHsimdata$vi,
              #aug.weights.x =   CCHsimdata$y1,
              risk.threshold  =  c(.05, .3),
              landmark.time = predict.time,
              weight.method = 'True',
              design = "CCH",
              smoothing.par = 0.9,
              calculate.sd  = TRUE,
              pnf.threshold = 0.85,
              pcf.threshold = 0.8,
              ncc.nmatch = 2)

#simulated data from a ncc design with nmatch = 2
data("NCCsimdata")
AIPWmeasures( time = NCCsimdata$xi, event = NCCsimdata$di,
                    X = cbind(NCCsimdata$y1, NCCsimdata$y2),
                    subcohort = NCCsimdata$vi,
                    aug.weights.x =   NCCsimdata$y1,
                    risk.threshold  =  c(.01, .03),
                    landmark.time = predict.time,
                    weight.method = 'Aug',
                    design = "NCC",
                    smoothing.par = 0.9,
                    calculate.sd  = TRUE,
                    pnf.threshold = 0.85,
                    pcf.threshold = 0.8,
                    ncc.nmatch = 2)

