October 13, 2016  

## AIPWmeasures

This package includes a function `AIPWmeasures` that calculates standard measures of predictive accuracy commonly used in biomarker validation studies. More specifically, this function calculates estimates using data from two-phase sampling designs ('case-cohort' and 'nested case-control') using two different methods, a standard ipw estimator (true ipw) and a novel method (augmented ipw) that has been shown to be more efficient in some contexts. See the manuscript "Improving Efficiency in Biomarker Incremental Value Evaluation under Two-phase Study Designs" by Zheng et. al. for more details.

## Install the package 

You can install the development version of the package directly from github using the R package `devtools`. The package will be added to CRAN after further testing. 

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("AIPWmeasures", "mdbrown")
```

## Quick Example 


```r
library(AIPWmeasures)
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
                   pcf.threshold = 0.8)
```

```
##    measure threshold   estimate          sd
## 1  coef.x1        NA 1.05978364 0.145055914
## 2  coef.x2        NA 0.71931257 0.142795037
## 3      AUC        NA 0.87401754 0.009552598
## 4      IDI        NA 0.35754277 0.024312293
## 5     ITPR        NA 0.46734730 0.022720905
## 6     IFPR        NA 0.10980453 0.007838267
## 7      TPR      0.05 0.95288855 0.006929574
## 8      FPR      0.05 0.53276478 0.041303221
## 9      PPV      0.05 0.26937735 0.013421904
## 10     NPV      0.05 0.97963812 0.002090129
## 11     TPR      0.30 0.62682586 0.030088946
## 12     FPR      0.30 0.09445823 0.010898073
## 13     PPV      0.30 0.57769278 0.018163560
## 14     NPV      0.30 0.92170129 0.006369866
## 15     PCF      0.80 0.65228402 0.022046510
## 16     PNF      0.85 0.61286419 0.022557534
## 17      NB      0.05 0.13960886          NA
## 18      NB      0.30 0.07356653          NA
```

```r
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
```

```
##    measure threshold   estimate          sd
## 1  coef.x1        NA 1.18479929 0.163974680
## 2  coef.x2        NA 0.66773798 0.140826774
## 3      AUC        NA 0.88616251 0.011032333
## 4      IDI        NA 0.38700949 0.029590607
## 5     ITPR        NA 0.49032068 0.026694264
## 6     IFPR        NA 0.10331119 0.008204509
## 7      TPR      0.05 0.95150194 0.005591640
## 8      FPR      0.05 0.49785154 0.043204502
## 9      PPV      0.05 0.27922295 0.014815301
## 10     NPV      0.05 0.98079949 0.001857965
## 11     TPR      0.30 0.64722244 0.031923770
## 12     FPR      0.30 0.08794821 0.011257735
## 13     PPV      0.30 0.59865952 0.023809816
## 14     NPV      0.30 0.92729881 0.006711179
## 15     PCF      0.80 0.67823896 0.025454076
## 16     PNF      0.85 0.63985556 0.026438899
## 17      NB      0.05 0.13857292          NA
## 18      NB      0.30 0.07773870          NA
```

```r
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
```

```
##    measure threshold   estimate          sd
## 1  coef.x1        NA 0.71388346 0.208118358
## 2  coef.x2        NA 1.15691575 0.271633156
## 3      AUC        NA 0.88616115 0.020995366
## 4      IDI        NA 0.18682668 0.113202687
## 5     ITPR        NA 0.20556047 0.114756791
## 6     IFPR        NA 0.01873379 0.009894893
## 7      TPR      0.01 0.89471342 0.022890509
## 8      FPR      0.01 0.32808139 0.065253683
## 9      PPV      0.01 0.05820036 0.010594222
## 10     NPV      0.01 0.99646180 0.008799836
## 11     TPR      0.03 0.76411311 0.076181378
## 12     FPR      0.03 0.16240075 0.037477153
## 13     PPV      0.03 0.09634659 0.017867684
## 14     NPV      0.03 0.99365883 0.008587344
## 15     PCF      0.80 0.79915680 0.038859288
## 16     PNF      0.85 0.73365163 0.044845399
## 17      NB      0.01 0.01658467          NA
## 18      NB      0.03 0.01201992          NA
```

```r
#true ipw
AIPWmeasures( time = NCCsimdata$xi, event = NCCsimdata$di,
                    X = cbind(NCCsimdata$y1, NCCsimdata$y2),
                    subcohort = NCCsimdata$vi,
                    risk.threshold  =  c(.01, .03),
                    landmark.time = predict.time,
                    weight.method = 'True',
                    design = "NCC",
                    smoothing.par = 0.9,
                    calculate.sd  = TRUE,
                    pnf.threshold = 0.85,
                    pcf.threshold = 0.8,
                    ncc.nmatch = 2)
```

```
##    measure threshold   estimate           sd
## 1  coef.x1        NA 0.47387905 0.3277693989
## 2  coef.x2        NA 1.27659158 0.3174745046
## 3      AUC        NA 0.87132718 0.0261367899
## 4      IDI        NA 0.15487564 0.0640699239
## 5     ITPR        NA 0.17521745 0.0634192962
## 6     IFPR        NA 0.02034181 0.0036090417
## 7      TPR      0.01 0.89313631 0.0271053204
## 8      FPR      0.01 0.36187459 0.0745445204
## 9      PPV      0.01 0.05484130 0.0093359570
## 10     NPV      0.01 0.99607842 0.0008082396
## 11     TPR      0.03 0.73906260 0.0736501276
## 12     FPR      0.03 0.16633111 0.0592791358
## 13     PPV      0.03 0.09458037 0.0269468452
## 14     NPV      0.03 0.99269531 0.0019549727
## 15     PCF      0.80 0.76516131 0.0581184633
## 16     PNF      0.85 0.70053556 0.0633486323
## 17      NB      0.01 0.01694355           NA
## 18      NB      0.03 0.01194980           NA
```
