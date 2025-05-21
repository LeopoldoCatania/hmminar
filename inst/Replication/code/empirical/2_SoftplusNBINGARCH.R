library(parallel)
library(hmminar)

## parallelization
iW = 6
cluster = makeCluster(iW)

sPath = "XXX"

#Load Data
load(file = paste(sPath, "data/vNT_1min.rdata", sep = ""))
#Load seasonal matrix
load(file = paste(sPath, "data/mD_1min.rdata", sep = ""))
# load seasonal indicator
load(file = paste(sPath, "data/vSeason_1min.rdata", sep = ""))
#load extra functions
source(paste(sPath, "code/functions.R", sep = ""))

iT = length(vN)
## in sample dataset
vN_is = vN[1:(iT/2)]
mD_is = mD[1:(iT/2), ]
## out of sample dataset
vN_oos = vN[(iT/2+1):length(vN)]
mD_oos = mD[(iT/2+1):length(vN), ]

set.seed(123)
Fit_is = Estimate_NBsoftplusINGARCH_DummySeason_DEoptim(vN_is, mD_is, dC = 1,
                                                        control = DEoptim.control(itermax = 5000, NP = 5000, steptol = 200,
                                                                                  cluster = cluster))


#Estimated parameters
vPar = Fit_is$vPar
dP = vPar["p"]

#Filtered mean and variance
vMean_is = Fit_is$Filter$vDeltaPlusM
vVar_is = Fit_is$Filter$vN*(1 - dP)/(dP^2)

## residuals
Res = (vN_is - vMean_is)/sqrt(vVar_is)

## Randomized PIT
#this is parameterized in terms of beta = 1-p
PIT_is = RandomizedPIT_NB_Softplus(vN_is, Fit_is$Filter$vN, Fit_is$vPar["p"], seed = 123)$vZ

## Collect results
NBsoftplusINGARCH_IS = list(Fit_is = Fit_is, Res = Res, PIT_is = PIT_is)

save(NBsoftplusINGARCH_IS, file = paste(sPath, "output/NBsoftplusINGARCH_IS.rdata", sep = ""))

vU = PIT_is
vU[vU > 0.9999] = 0.9999
vU[vU < 1.0 - 0.9999] = 1.0 - 0.9999

## Plot pit histogram
pdf(file = paste(sPath, "output/PIT_NBsoftplusINGARCH_is_1m.pdf", sep = ""), width = 5, height = 5)

foo =  BinTest(vU, plot = TRUE, main = "", g = 10, mar = c(1.5, 1.5, 0, 1.5), las = 1)

dev.off()

### Out of sample
iH = Fit_is$iH
vPar = Fit_is$vPar
dM0 = vPar["m0"]
dAlpha0 = vPar["alpha0"]
dAlpha1 = vPar["alpha1"]
dBeta1 = vPar["beta1"]
dP = vPar["p"]

vBeta_s = Fit_is$Filter$vBeta_s

## Filter using full sample
Filter_full = NBSoftplusSeasonalINGARCH_DummySeason_filter(vN, dM0,
                                                           dAlpha0,
                                                           dAlpha1,  dBeta1,
                                                           vSeason, vBeta_s, dC = 1.0,
                                                           dP)

#compute oos PIT
PIT_oos = RandomizedPIT_NB_Softplus(vN_oos, Filter_full$vN[(iT/2+1):length(vN)], Fit_is$vPar["p"], seed = 123)$vZ

vMean_oos = Filter_full$vDeltaPlusM[(iT/2+1):length(vN)]
vVar_oos = Filter_full$vN[(iT/2+1):length(vN)]*(1 - dP)/(dP^2)

# Forecast error
FE_oos = (vN_oos - vMean_oos)

## Simulate NB
SimNB <- function(n, dN, dP) {
  dTheta = (1-dP)/dP
  vLambda = rgamma(n, shape = dN, scale = dTheta)
  rpois(n, vLambda)
}

## Compute median predictions
set.seed(123)
Median_oos = sapply(Filter_full$vN[(iT/2+1):length(vN)], function(dN) {

  vFoo = SimNB(1e4, dN, Fit_is$vPar["p"])

  return(round(median(vFoo)))

})

## Collect results
NBsoftplusINGARCH_OOS = list(Filter_full=Filter_full, PIT_oos=PIT_oos, FE_oos=FE_oos, Mean_oos = vMean_oos,
                 Median_oos = Median_oos)

save(NBsoftplusINGARCH_OOS, file = paste(sPath, "output/NBsoftplusINGARCH_OOS_1m.rdata", sep = ""))

stopCluster(cluster)
rm(list = ls())
