library(hmminar)
library(parallel)

sPath = "XXX"

iL = 6 #Number of cores
cluster = makeCluster(iL)
clusterEvalQ(cluster, library(hmminar))

### INGARCH parametric seasonal

#Load Data
load(file = paste(sPath, "data/vNT_1min.rdata", sep = ""))
#Load seasonal matrix
load(file = paste(sPath, "data/mD_1min.rdata", sep = ""))
source(paste(sPath, "code/functions.R", sep = ""))

iT = length(vN)
## in sample dataset
vN_is = vN[1:(iT/2)]
## out of sample dataset
vN_oos = vN[(iT/2+1):length(vN)]
mD_oos = mD[(iT/2+1):length(vN),]
# Estimate the model for all H = 1, ..., H_max
lFit = list()
iH_max = 7

for (iH in 0:iH_max) {
  Fit_foo = Estimate_NBsoftplusINGARCH(vN_is, iD = 390, dC = 1, iH = iH)
  lFit[[paste(iH)]] = Fit_foo
}

# Selects H (the number of nonparametric seasons)
mIC = sapply(lFit, function(x) x$iC)

iH = paste(0:iH_max)[which.min(mIC["BIC",])] #H = 3

Fit_is =  lFit[[paste(iH)]]
#Estimated parameters
vPar = Fit_is$vPar
dP = vPar["p"]

#Filtered mean and variance
vMean_is = Fit_is$Filter$vDeltaPlusM
vVar_is = Fit_is$Filter$vN*(1 - dP)/(dP^2)

## residuals
Res = (vN_is - vMean_is[-(length(vN_is) - 1)])/sqrt(vVar_is[-(length(vN_is) - 1)])

## Randomized PIT
#this is parameterized in terms of beta = 1-p
PIT_is = RandomizedPIT_NB_Softplus(vN_is, Fit_is$Filter$vN[-(length(vN_is) - 1)], Fit_is$vPar["p"], seed = 123)$vZ

## Collect results
NBsoftplusINGARCH_IS = list(Fit_is = Fit_is, Res = Res, PIT_is = PIT_is)

save(NBsoftplusINGARCH_IS, file = paste(sPath, "output/NBsoftplusINGARCH_ps_IS.rdata", sep = ""))

### HMM-INAR

iH = 3 #number of non parametric seasons at the beginning

lComb = list(c(1, 1, 1), #iJ = 1, iL = 1, iK = 1
             c(4, 8, 3), #iJ = 4, iL = 3, iK = 8
             c(1, 7, 5), #iJ = 1, iL = 5, iK = 7,
             c(5, 6, 1), #iJ = 5, iL = 1, iK = 6,
             c(1, 9, 1), #iJ = 1, iL = 1, iK = 9
             c(5, 1, 1)) #iJ = 5, iL = 1, iK = 1

clusterExport(cluster, c("sPath", "iH", "FromNPStoPS"))

parLapplyLB(cluster, lComb, function(vComb) {

  iJ = vComb[1]
  iK = vComb[2]
  iL = vComb[3]

  sFile_NPS = paste(sPath, "output/Fit",
                    iJ, "_", iL, "_", iK, ".rdata", sep = "")

  load(sFile_NPS)

  Fit_PS = FromNPStoPS(Fit, iH = iH)

  save(Fit_PS, file = paste(sPath, "output/Fit",
                            iJ, "_", iL, "_", iK, "_ps.rdata", sep = ""))

})

stopCluster(cluster)

### Parametric seasonal results

# IN SAMPLE PIT PLOT
iT = length(vN)/2
mPIT_ps = matrix(0, iT-1, 7, dimnames = list(NULL, c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")))
mRes_ps = matrix(0, iT-1, 7, dimnames = list(NULL, c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")))

########################################### HMM INAR ##########################################################
load(paste(sPath, "output/Fit4_3_8_ps.rdata", sep = ""))
Fit_PS = AdjustFitForFiltering(Fit_PS)
mPIT_ps[, "HMMINAR"] = RandomizedPIT(Fit_PS, seed = 123)$vZ
mRes_ps[, "HMMINAR"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit_PS)$vRes

########################################### INAR ##########################################################
load(paste(sPath, "output/Fit1_1_1_ps.rdata", sep = ""))
Fit_PS = AdjustFitForFiltering(Fit_PS)
mPIT_ps[, "INAR"] = RandomizedPIT(Fit_PS, seed = 123)$vZ
mRes_ps[, "INAR"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit_PS)$vRes

########################################### Mix INAR ##########################################################
load(paste(sPath, "output/Fit1_1_9_ps.rdata", sep = ""))
Fit_PS = AdjustFitForFiltering(Fit_PS)
mPIT_ps[, "MixINAR"] = RandomizedPIT(Fit_PS, seed = 123)$vZ
mRes_ps[, "MixINAR"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit_PS)$vRes

########################################### HMM INAR L=1 ##########################################################
load(paste(sPath, "output/Fit5_1_6_ps.rdata", sep = ""))
Fit_PS = AdjustFitForFiltering(Fit_PS)
mPIT_ps[, "HMMINAR_l1"] = RandomizedPIT(Fit_PS, seed = 123)$vZ
mRes_ps[, "HMMINAR_l1"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit_PS)$vRes

########################################### HMM INAR J=1 ##########################################################
load(paste(sPath, "output/Fit1_5_7_ps.rdata", sep = ""))
Fit_PS = AdjustFitForFiltering(Fit_PS)
mPIT_ps[, "HMMINAR_j1"] = RandomizedPIT(Fit_PS, seed = 123)$vZ
mRes_ps[, "HMMINAR_j1"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit_PS)$vRes

########################################### HMM INAR K=1 ##########################################################
load(paste(sPath, "output/Fit5_1_1_ps.rdata", sep = ""))
Fit_PS = AdjustFitForFiltering(Fit_PS)
mPIT_ps[, "HMMINAR_k1"] = RandomizedPIT(Fit_PS, seed = 123)$vZ
mRes_ps[, "HMMINAR_k1"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit_PS)$vRes

########################################### NBsoftplusINGARCH ##########################################################
load(file = paste(sPath, "output/NBsoftplusINGARCH_ps_IS.rdata", sep = ""))
mPIT_ps[, "INGARCH"] = NBsoftplusINGARCH_IS$PIT_is[-1]
mRes_ps[, "INGARCH"] = NBsoftplusINGARCH_IS$Res[-1]

save(mPIT_ps, file = paste(sPath, "output/mPIT_ps.rdata", sep = ""))
save(mRes_ps, file = paste(sPath, "output/mRes_ps.rdata", sep = ""))

vModelNames = colnames(mPIT_ps)
sPath_Save = paste(sPath, "output/", sep = "")

for (sModel in vModelNames) {

  vU = mPIT_ps[,sModel]

  vU[vU > 0.9999] = 0.9999
  vU[vU < 1.0 - 0.9999] = 1.0 - 0.9999

  pdf(file = paste(sPath_Save, "PIT_", sModel, "_ps_is_1m.pdf", sep = ""), width = 5, height = 5)
  par(mar = c(2,4,1,1))
  foo =  BinTest(vU, plot = TRUE, main = "", g = 10, mar = c(1.5, 1.5, 0, 1.5), las = 1, cex.axis = 1.5)

  dev.off()
}

## Residuals ACF

iLag = 50
vModels = c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")
mACF_ps = matrix(0, iLag, length(vModels), dimnames = list(NULL, vModels))

for (sModel in vModels) {
  vRes = mRes_ps[, sModel]
  mACF_ps[,sModel] = acf(vRes, lag.max = iLag, plot = FALSE)$acf[-1,1,1]
}

iT = length(vN)/2

for (sModel in vModels) {

  vACF = mACF_ps[,sModel]
  pdf(file = paste(sPath_Save, "ACF_residuals_", sModel ,"_ps_1m.pdf", sep = ""), width = 5, height = 5)
  par(mar = c(2,4,1,1))
  plot(1:length(vACF), vACF, ylim = range(mACF_ps),
       las = 1, ylab = "", xlab = "", type = "n", cex.axis = 1.5)
  grid(10, 10, "gray80")

  lines(vACF, lwd = 2, lty = 1)

  abline(h = 1.96/sqrt(iT), col = "red")
  abline(h = -1.96/sqrt(iT), col = "red")

  dev.off()
}
# SEASONALITY PLOT

mBeta_PS = matrix(0, 390, 7, dimnames = list(NULL, c("HMMINAR", "INAR", "MixINAR", "HMMINARl1", "HMMINARj1", "HMMINARk1", "INGARCH")))

########################################### HMM INAR ##########################################################
load(paste(sPath, "output/Fit4_3_8_ps.rdata", sep = ""))
mBeta_PS[,"HMMINAR"] = Fit_PS$vDelta_S

########################################### INAR ##########################################################
load(paste(sPath, "output/Fit1_1_1_ps.rdata", sep = ""))
mBeta_PS[,"INAR"] = Fit_PS$vDelta_S

########################################### Mix INAR ##########################################################
load(paste(sPath, "output/Fit1_1_9_ps.rdata", sep = ""))
mBeta_PS[,"MixINAR"] = Fit_PS$vDelta_S

########################################### HMM INAR L=1 ##########################################################
load(paste(sPath, "output/Fit5_1_6_ps.rdata", sep = ""))
mBeta_PS[,"HMMINARl1"] = Fit_PS$vDelta_S

########################################### HMM INAR J=1 ##########################################################
load(paste(sPath, "output/Fit1_5_7_ps.rdata", sep = ""))
mBeta_PS[,"HMMINARj1"] = Fit_PS$vDelta_S

########################################### HMM INAR K=1 ##########################################################
load(paste(sPath, "output/Fit5_1_1_ps.rdata", sep = ""))
mBeta_PS[,"HMMINARk1"] = Fit_PS$vDelta_S

########################################### INGARCH ##########################################################
load(paste(sPath, "output/NBsoftplusINGARCH_ps_IS.rdata", sep = ""))
mBeta_PS[,"INGARCH"] = NBsoftplusINGARCH_IS$Fit_is$Filter$vDelta_S

for (i in 1:ncol(mBeta_PS)) {
  mBeta_PS[, i] = mBeta_PS[, i]/mBeta_PS[1,i]
}

iD = 390

vSeq = seq(0, 390, 30)

vLabels = c("09:30", "10:00", "10:30", "11:00", "11:30",
            "12:00", "12:30", "13:00", "13:30", "14:00",
            "14:30", "15:00", "15:30", "16:00")

vModels = colnames(mBeta_PS)

for(sModel in vModels){

  pdf(paste(sPath_Save, sModel, "_beta_coeff_ps_1m.pdf", sep = ""), width = 5, height = 5)

  par(mar = c(2,4,1,1))
  plot(1:390, log(mBeta_PS[, sModel]), type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 1.5)
  grid(10, 10, col = "gray80")
  lines(1:390, log(mBeta_PS[, sModel]))

  axis(1, at = vSeq, labels = vLabels, cex.axis = 1.5)
  axis(1, at = seq(0, 390, 5), labels = FALSE, tcl = -0.25)
  axis(1, at = seq(0, 390, 1), labels = FALSE, tcl = -0.25/2)

  dev.off()

}

#### Out of sample results

###############################################################################################################
########################################### HMM INAR ##########################################################

load(paste(sPath, "output/Fit4_3_8_ps.rdata", sep = ""))
Fit_PS$vBeta = Fit_PS$vDelta_S
Fit_PS$mGamma_Eta = Fit_PS$lPn$mGamma_Eta
Fit_PS$mGamma_Alpha = Fit_PS$lPn$mGamma_Alpha
Fit_PS$mOmega = Fit_PS$lPn$mOmega
Fit_PS$vLambda = Fit_PS$lPn$vLambda
Fit_PS$vAlpha = Fit_PS$lPn$vAlpha
Fit_PS$dVarPhi = Fit_PS$lPn$dVarPhi

iJ = 4
iL = 3
iK = 8

iMax = max(vN)

vLogFactorial_K = lfactorial(0:(iMax))

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit_PS$lPn$mGamma_Eta,
                                             Fit_PS$lPn$mGamma_Alpha,
                                             Fit_PS$lPn$mOmega,
                                             Fit_PS$lPn$vLambda,
                                             Fit_PS$lPn$vAlpha,
                                             Fit_PS$vDelta_S,
                                             Fit_PS$lPn$dVarPhi,
                                             vSeason,
                                             vLogFactorial_K,
                                             vO,
                                             maxIter = 1, tol = 1e-5, bSeasonal, bFilter = TRUE)


dLLK = Run[["dLLK"]]
np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD + 1

IC = ICfun(dLLK, np, iT)

Run[["iT"]]   = iT
Run[["IC"]]   = IC
Run[["iJ"]]   = iJ
Run[["iK"]]   = iK
Run[["iL"]]   = iL
Run[["iD"]]   = iD
Run[["mD"]]   = mD

Filter = Filtering_MSMixInar_TwoChains_AlfaCor(Run, ComputeMedian = TRUE)

lPred = list(Filter = Filter, Run = Run, Fit_is = Fit_PS)

save(lPred, file = paste(sPath, "output/lPred_4_3_8_ps.rdata", sep = ""))

###############################################################################################################
########################################### INAR ##########################################################


load(paste(sPath, "output/Fit1_1_1_ps.rdata", sep = ""))
Fit_PS$vBeta = Fit_PS$vDelta_S
Fit_PS$mGamma_Eta = Fit_PS$lPn$mGamma_Eta
Fit_PS$mGamma_Alpha = Fit_PS$lPn$mGamma_Alpha
Fit_PS$mOmega = Fit_PS$lPn$mOmega
Fit_PS$vLambda = Fit_PS$lPn$vLambda
Fit_PS$vAlpha = Fit_PS$lPn$vAlpha
Fit_PS$dVarPhi = Fit_PS$lPn$dVarPhi

iJ = 1
iL = 1
iK = 1

###

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit_PS$lPn$mGamma_Eta,
                                             Fit_PS$lPn$mGamma_Alpha,
                                             Fit_PS$lPn$mOmega,
                                             Fit_PS$lPn$vLambda,
                                             Fit_PS$lPn$vAlpha,
                                             Fit_PS$vDelta_S,
                                             Fit_PS$lPn$dVarPhi,
                                             vSeason,
                                             vLogFactorial_K,
                                             vO,
                                             maxIter = 1, tol = 1e-5, bSeasonal, bFilter = TRUE)


dLLK = Run[["dLLK"]]
np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD + 1

IC = ICfun(dLLK, np, iT)

Run[["iT"]]   = iT
Run[["IC"]]   = IC
Run[["iJ"]]   = iJ
Run[["iK"]]   = iK
Run[["iL"]]   = iL
Run[["iD"]]   = iD
Run[["mD"]]   = mD

Filter = Filtering_MSMixInar_TwoChains_AlfaCor(Run, ComputeMedian = TRUE)

lPred_INAR = list(Filter = Filter, Run = Run, Fit_is = Fit_PS)

save(lPred_INAR, file = paste(sPath, "output/lPred_INAR_ps.rdata", sep = ""))

gc()
###############################################################################################################
########################################### Mix INAR ##########################################################

load(paste(sPath, "output/Fit1_1_9_ps.rdata", sep = ""))
Fit_PS$vBeta = Fit_PS$vDelta_S
Fit_PS$mGamma_Eta = Fit_PS$lPn$mGamma_Eta
Fit_PS$mGamma_Alpha = Fit_PS$lPn$mGamma_Alpha
Fit_PS$mOmega = Fit_PS$lPn$mOmega
Fit_PS$vLambda = Fit_PS$lPn$vLambda
Fit_PS$vAlpha = Fit_PS$lPn$vAlpha
Fit_PS$dVarPhi = Fit_PS$lPn$dVarPhi

iJ = 1
iL = 1
iK = 9
iJ*(iJ - 1) + iL*(iL-1) + (iK-1)*iL + iJ + iK + 80 + 1

###
Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit_PS$lPn$mGamma_Eta,
                                             Fit_PS$lPn$mGamma_Alpha,
                                             Fit_PS$lPn$mOmega,
                                             Fit_PS$lPn$vLambda,
                                             Fit_PS$lPn$vAlpha,
                                             Fit_PS$vDelta_S,
                                             Fit_PS$lPn$dVarPhi,
                                             vSeason,
                                             vLogFactorial_K,
                                             vO,
                                             maxIter = 1, tol = 1e-5, bSeasonal, bFilter = TRUE)


dLLK = Run[["dLLK"]]
np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD + 1

IC = ICfun(dLLK, np, iT)

Run[["iT"]]   = iT
Run[["IC"]]   = IC
Run[["iJ"]]   = iJ
Run[["iK"]]   = iK
Run[["iL"]]   = iL
Run[["iD"]]   = iD
Run[["mD"]]   = mD

Filter = Filtering_MSMixInar_TwoChains_AlfaCor(Run, ComputeMedian = TRUE)

lPred_MixINAR = list(Filter = Filter, Run = Run, Fit_is = Fit_PS)

save(lPred_MixINAR, file = paste(sPath, "output/lPred_MixINAR_ps.rdata", sep = ""))

###############################################################################################################
########################################### HMM INAR L=1 ##########################################################

load(paste(sPath, "output/Fit5_1_6_ps.rdata", sep = ""))
Fit_PS$vBeta = Fit_PS$vDelta_S
Fit_PS$mGamma_Eta = Fit_PS$lPn$mGamma_Eta
Fit_PS$mGamma_Alpha = Fit_PS$lPn$mGamma_Alpha
Fit_PS$mOmega = Fit_PS$lPn$mOmega
Fit_PS$vLambda = Fit_PS$lPn$vLambda
Fit_PS$vAlpha = Fit_PS$lPn$vAlpha
Fit_PS$dVarPhi = Fit_PS$lPn$dVarPhi

iJ = 5
iL = 1
iK = 6

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit_PS$lPn$mGamma_Eta,
                                             Fit_PS$lPn$mGamma_Alpha,
                                             Fit_PS$lPn$mOmega,
                                             Fit_PS$lPn$vLambda,
                                             Fit_PS$lPn$vAlpha,
                                             Fit_PS$vDelta_S,
                                             Fit_PS$lPn$dVarPhi,
                                             vSeason,
                                             vLogFactorial_K,
                                             vO,
                                             maxIter = 1, tol = 1e-5, bSeasonal, bFilter = TRUE)


dLLK = Run[["dLLK"]]
np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD + 1

IC = ICfun(dLLK, np, iT)

Run[["iT"]]   = iT
Run[["IC"]]   = IC
Run[["iJ"]]   = iJ
Run[["iK"]]   = iK
Run[["iL"]]   = iL
Run[["iD"]]   = iD
Run[["mD"]]   = mD

Filter = Filtering_MSMixInar_TwoChains_AlfaCor(Run, ComputeMedian = TRUE)

lPred_HMMINAR_l1 = list(Filter = Filter, Run = Run, Fit_is = Fit_PS)

save(lPred_HMMINAR_l1, file = paste(sPath, "output/lPred_HMMINAR_l1_ps.rdata", sep = ""))

gc()
###############################################################################################################
########################################### HMM INAR J=1 ##########################################################

load(paste(sPath, "output/Fit1_5_7_ps.rdata", sep = ""))
Fit_PS$vBeta = Fit_PS$vDelta_S
Fit_PS$mGamma_Eta = Fit_PS$lPn$mGamma_Eta
Fit_PS$mGamma_Alpha = Fit_PS$lPn$mGamma_Alpha
Fit_PS$mOmega = Fit_PS$lPn$mOmega
Fit_PS$vLambda = Fit_PS$lPn$vLambda
Fit_PS$vAlpha = Fit_PS$lPn$vAlpha
Fit_PS$dVarPhi = Fit_PS$lPn$dVarPhi

iJ = 1
iL = 5
iK = 7

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit_PS$lPn$mGamma_Eta,
                                             Fit_PS$lPn$mGamma_Alpha,
                                             Fit_PS$lPn$mOmega,
                                             Fit_PS$lPn$vLambda,
                                             Fit_PS$lPn$vAlpha,
                                             Fit_PS$vDelta_S,
                                             Fit_PS$lPn$dVarPhi,
                                             vSeason,
                                             vLogFactorial_K,
                                             vO,
                                             maxIter = 1, tol = 1e-5, bSeasonal, bFilter = TRUE)


dLLK = Run[["dLLK"]]
np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD + 1

IC = ICfun(dLLK, np, iT)

Run[["iT"]]   = iT
Run[["IC"]]   = IC
Run[["iJ"]]   = iJ
Run[["iK"]]   = iK
Run[["iL"]]   = iL
Run[["iD"]]   = iD
Run[["mD"]]   = mD

Filter = Filtering_MSMixInar_TwoChains_AlfaCor(Run, ComputeMedian = TRUE)

lPred_HMMINAR_j1 = list(Filter = Filter, Run = Run, Fit_is = Fit_PS)

save(lPred_HMMINAR_j1, file = paste(sPath, "output/lPred_HMMINAR_j1_ps.rdata", sep = ""))

gc()

###############################################################################################################
########################################### HMM INAR K=1 ##########################################################

load(paste(sPath, "output/Fit5_1_1_ps.rdata", sep = ""))
Fit_PS$vBeta = Fit_PS$vDelta_S
Fit_PS$mGamma_Eta = Fit_PS$lPn$mGamma_Eta
Fit_PS$mGamma_Alpha = Fit_PS$lPn$mGamma_Alpha
Fit_PS$mOmega = Fit_PS$lPn$mOmega
Fit_PS$vLambda = Fit_PS$lPn$vLambda
Fit_PS$vAlpha = Fit_PS$lPn$vAlpha
Fit_PS$dVarPhi = Fit_PS$lPn$dVarPhi

iJ = 5
iL = 1
iK = 1

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit_PS$lPn$mGamma_Eta,
                                             Fit_PS$lPn$mGamma_Alpha,
                                             Fit_PS$lPn$mOmega,
                                             Fit_PS$lPn$vLambda,
                                             Fit_PS$lPn$vAlpha,
                                             Fit_PS$vDelta_S,
                                             Fit_PS$lPn$dVarPhi,
                                             vSeason,
                                             vLogFactorial_K,
                                             vO,
                                             maxIter = 1, tol = 1e-5, bSeasonal, bFilter = TRUE)


dLLK = Run[["dLLK"]]
np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD + 1

IC = ICfun(dLLK, np, iT)

Run[["iT"]]   = iT
Run[["IC"]]   = IC
Run[["iJ"]]   = iJ
Run[["iK"]]   = iK
Run[["iL"]]   = iL
Run[["iD"]]   = iD
Run[["mD"]]   = mD

Filter = Filtering_MSMixInar_TwoChains_AlfaCor(Run, ComputeMedian = TRUE)

lPred_HMMINAR_k1 = list(Filter = Filter, Run = Run, Fit_is = Fit_PS)

save(lPred_HMMINAR_k1, file = paste(sPath, "output/lPred_HMMINAR_k1_ps.rdata", sep = ""))

gc()

###############################################################################################################
########################################### NBsoftplusINGARCH ##########################################################

load(file = paste(sPath, "output/NBsoftplusINGARCH_ps_IS.rdata", sep = ""))
Fit_is = NBsoftplusINGARCH_IS$Fit_is

iH = Fit_is$iH
vPar = Fit_is$vPar
dM0 = vPar["m0"]
dAlpha0 = vPar["alpha0"]
dAlpha1 = vPar["alpha1"]
dBeta1 = vPar["beta1"]
dP = vPar["p"]

vDelta_S = Fit_is$Filter$vDelta_S

## Filter using full sample
Filter_full = NBSoftplusSeasonalINGARCH_filter(vN, dM0,
                                               dAlpha0,
                                               dAlpha1,  dBeta1,
                                               vDelta_S, iS = iD, dC = 1,
                                               dP)

#compute oos PIT
PIT_oos = RandomizedPIT_NB_Softplus(vN_oos, Filter_full$vN[(length(vN)/2+1):length(vN)], Fit_is$vPar["p"], seed = 123)$vZ

vMean_oos = Filter_full$vDeltaPlusM[(length(vN)/2+1):length(vN)]
vVar_oos = Filter_full$vN[(length(vN)/2+1):length(vN)]*(1 - dP)/(dP^2)

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
Median_oos = sapply(Filter_full$vN[(length(vN)/2+1):length(vN)], function(dN) {

  vFoo = SimNB(1e4, dN, Fit_is$vPar["p"])

  return(round(median(vFoo)))

})

## Collect results
NBsoftplusINGARCH_OOS = list(Filter_full=Filter_full, PIT_oos=PIT_oos, FE_oos=FE_oos, Mean_oos = vMean_oos,
                             Median_oos = Median_oos)

save(NBsoftplusINGARCH_OOS, file = paste(sPath, "output/NBsoftplusINGARCH_ps_OOS_1m.rdata", sep = ""))

### Organize results
load(paste(sPath, "output/lPred_4_3_8_ps.rdata", sep = ""))
load(paste(sPath, "output/lPred_INAR_ps.rdata", sep = ""))
load(paste(sPath, "output/lPred_MixINAR_ps.rdata", sep = ""))
load(paste(sPath, "output/lPred_HMMINAR_j1_ps.rdata", sep = ""))
load(paste(sPath, "output/lPred_HMMINAR_l1_ps.rdata", sep = ""))
load(paste(sPath, "output/lPred_HMMINAR_k1_ps.rdata", sep = ""))
load(paste(sPath, "output/NBsoftplusINGARCH_ps_OOS_1m.rdata", sep = ""))

iF = length(vN)/2

vN_Diff = diff(vN)
vN_Diff_78 = diff(vN, 390)

aForecastError_ps = array(0, dim = c(iF, 7, 2), dimnames = list(NULL,
                                                                c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1",
                                                                  "HMMINAR_j1", "HMMINAR_k1", "INGARCH"), c("Median", "Mean")))

vY = c(lPred$Run$dY0, lPred$Run$vY)

aForecastError_ps[, "HMMINAR", "Median"] = tail((vY - lPred$Filter$vMedian)[-1], iF)
aForecastError_ps[, "INAR", "Median"] = tail((vY - lPred_INAR$Filter$vMedian)[-1], iF)
aForecastError_ps[, "MixINAR", "Median"] = tail((vY - lPred_MixINAR$Filter$vMedian)[-1], iF)
aForecastError_ps[, "HMMINAR_j1", "Median"] = tail((vY - lPred_HMMINAR_j1$Filter$vMedian)[-1], iF)
aForecastError_ps[, "HMMINAR_l1", "Median"] = tail((vY - lPred_HMMINAR_l1$Filter$vMedian)[-1], iF)
aForecastError_ps[, "HMMINAR_k1", "Median"] = tail((vY - lPred_HMMINAR_k1$Filter$vMedian)[-1], iF)
aForecastError_ps[, "INGARCH", "Median"] = tail(vY, iF) - NBsoftplusINGARCH_OOS$Median_oos

aForecastError_ps[, "HMMINAR", "Mean"] = tail(lPred$Filter$vRes, iF)
aForecastError_ps[, "INAR", "Mean"] = tail(lPred_INAR$Filter$vRes, iF)
aForecastError_ps[, "MixINAR", "Mean"] = tail(lPred_MixINAR$Filter$vRes, iF)
aForecastError_ps[, "HMMINAR_j1", "Mean"] = tail(lPred_HMMINAR_j1$Filter$vRes, iF)
aForecastError_ps[, "HMMINAR_l1", "Mean"] = tail(lPred_HMMINAR_l1$Filter$vRes, iF)
aForecastError_ps[, "HMMINAR_k1", "Mean"] = tail(lPred_HMMINAR_k1$Filter$vRes, iF)
aForecastError_ps[, "INGARCH", "Mean"] = NBsoftplusINGARCH_OOS$FE_oos

save(aForecastError_ps, file = paste(sPath, "output/aForecastError_ps.rdata", sep = ""))
load(file = paste(sPath, "output/aForecastError.rdata", sep = ""))
load(file = paste(sPath, "output/aMSE.rdata", sep = ""))

vModelNames = dimnames(aForecastError_ps)[[2]]

vPeriods = c("Full", "Opening", "MidDay", "Closing") #34=150
# in terms of seasons
Opening = 1:10 #first half hour (note that first seasons are 1,2,3 minutes and then 4:5, 6:10, etc.)
MidDay  = 34:58 #from 2.5h to 4.5h
Closing = 76:81 #from 6h to 6.5h
mPeridos = cbind(TRUE, apply(mD_oos[, Opening] == 1, 1, any), apply(mD_oos[, MidDay] == 1, 1, any), apply(mD_oos[, Closing] == 1, 1, any))
colnames(mPeridos) = vPeriods

aMSE_ps = array(0, dim = c(length(vModelNames), length(vPeriods)), dimnames = list(vModelNames, vPeriods))
aDM_ps = array(9999, dim = c(length(vModelNames), length(vPeriods), 2), dimnames = list(vModelNames, vPeriods, c("Test", "Pval")))

for (sPeriod in vPeriods) {
  for (sModel in vModelNames) {
    vFoo = aForecastError_ps[, sModel, "Median"][mPeridos[, sPeriod]]^2
    #1% trimming
    vFoo = vFoo[vFoo<quantile(vFoo, 1-0.01)]
    aMSE_ps[sModel, sPeriod] = mean(vFoo)

    vFoo1 = aForecastError_ps[, sModel, "Median"][mPeridos[, sPeriod]]^2
    vFoo2 = aForecastError[, sModel, "Median"][mPeridos[, sPeriod]]^2

    vBaz = vFoo1<quantile(vFoo1, 1-0.01) & vFoo2<quantile(vFoo2, 1-0.01)
    vFoo1 = vFoo1[vBaz]
    vFoo2 = vFoo2[vBaz]

    aDM_ps[sModel, sPeriod,] = unlist(f.DMTest(vFoo1, vFoo2)[c("dTest", "dPvalue")])
  }
}

for (sPeriod in vPeriods) {
  aMSE_ps[, sPeriod] = aMSE_ps[, sPeriod]/aMSE[rownames(aMSE_ps), sPeriod]
}

mTab = NULL
mPval = NULL

for (sPeriod in vPeriods) {
  mTab = rbind(mTab,  aMSE_ps[, sPeriod])
  mPval = rbind(mPval, aDM_ps[, sPeriod, "Pval"] )
}

mTab_lower = mTab < 1
mTab_upper = mTab > 1

mTab = format(round(mTab, 2), scientific = FALSE)
mTab[mPval < 0.05 & mTab_lower] = paste("\\cellcolor{red}", mTab[mPval < 0.05 & mTab_lower], sep = "")
mTab[mPval < 0.05 & mTab_upper] = paste("\\cellcolor{green}", mTab[mPval < 0.05 & mTab_upper], sep = "")

mTab= cbind(mTab, vPeriods)

mTab[, ncol(mTab)] = paste(mTab[, ncol(mTab)], "\\\\")

write.table(mTab,
            file =  paste(sPath, "output/TradesPred_Median_psvsnps_1m.txt", sep = ""),
            sep = " & ",
            dec = ".",
            quote = FALSE,
            row.names = TRUE
)

############## ACF Forecast error ###################

iLag = 50

vModels = c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")

mACF_ps = matrix(0, iLag, length(vModels), dimnames = list(NULL, vModels))

for (sModel in vModels) {
  vRes = aForecastError_ps[, sModel, "Mean"]
  mACF_ps[,sModel] = acf(vRes, lag.max = iLag, plot = FALSE)$acf[-1,1,1]
}

for (sModel in vModels) {

  vACF = mACF_ps[,sModel]
  pdf(file = paste(sPath, "output/ACF_forecasterror_", sModel ,"_ps1m.pdf", sep = ""), width = 5, height = 5)
  par(mar = c(2,4,1,1))
  plot(1:length(vACF), vACF, ylim = range(c(mACF_ps, 1.96/sqrt(iF), -1.96/sqrt(iT))), las = 1, ylab = "", xlab = "", type = "n", cex.axis = 1.5)
  grid(10, 10, "gray80")

  lines(vACF, lwd = 2, lty = 1)

  abline(h = 1.96/sqrt(iF), col = "red")
  abline(h = -1.96/sqrt(iT), col = "red")

  dev.off()
}


############## OUT of Sample PIT ###################
library(hmminar)

aForecastPIT_ps = array(0, dim = c(iF, 7), dimnames = list(NULL, c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")))

aForecastPIT_ps[, "HMMINAR"] = tail(RandomizedPIT(lPred$Run)$vZ, iF)
aForecastPIT_ps[, "INAR"] = tail(RandomizedPIT(lPred_INAR$Run)$vZ, iF)
aForecastPIT_ps[, "MixINAR"] = tail(RandomizedPIT(lPred_MixINAR$Run)$vZ, iF)
aForecastPIT_ps[, "HMMINAR_j1"] = tail(RandomizedPIT(lPred_HMMINAR_j1$Run)$vZ, iF)
aForecastPIT_ps[, "HMMINAR_l1"] = tail(RandomizedPIT(lPred_HMMINAR_l1$Run)$vZ, iF)
aForecastPIT_ps[, "HMMINAR_k1"] = tail(RandomizedPIT(lPred_HMMINAR_k1$Run)$vZ, iF)
aForecastPIT_ps[, "INGARCH"] = NBsoftplusINGARCH_OOS$PIT_oos

save(aForecastPIT_ps, file = paste(sPath, "output/aForecastPIT_ps.rdata", sep = ""))

vModelNames = dimnames(aForecastPIT_ps)[[2]]
for (sModel in vModelNames) {

  vU = aForecastPIT_ps[,sModel]

  vU[vU > 0.9999] = 0.9999
  vU[vU < 1.0 - 0.9999] = 1.0 - 0.9999

  sPath_Save = paste(sPath, "output/", sep = "")

  pdf(file = paste(sPath_Save, "PIT_", sModel, "ps_oos_1m.pdf", sep = ""), width = 5, height = 5)
  par(mar = c(2,4,1,1))
  foo =  BinTest(vU, plot = TRUE, main = "", g = 10, mar = c(1.5, 1.5, 0, 1.5), las = 1, cex.axis = 1.5)

  dev.off()
}

rm(list = ls())



