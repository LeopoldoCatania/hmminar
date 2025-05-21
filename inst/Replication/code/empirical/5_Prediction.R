library(parallel)
library(hmminar)

# the main directory
sPath = "XXX"

#Load Data
load(file = paste(sPath, "data/vNT_1min.rdata", sep = ""))
#Load seasonal matrix
load(file = paste(sPath, "data/mD_1min.rdata", sep = ""))
# load seasonal indicator
load(file = paste(sPath, "data/vSeason_1min.rdata", sep = ""))
# load other functions
source(paste(sPath, "code/functions.R", sep = ""))

#Compute opening dummy + in sample and out of sample

iDays = length(vN)/390
iT = length(vN)

vO = numeric(iT)
vO[seq(1, iT, 390)] = 1

mD_is = mD[1:(iDays/2 * 390), ]
vN_is = vN[1:(iDays/2 * 390)]
vO_is = vO[1:(iDays/2 * 390)]

mD_oos = mD[(iDays/2 * 390 + 1):iT,]
vN_oos = vN[(iDays/2 * 390 + 1):iT]
vO_oos = vO[(iDays/2 * 390 + 1):iT]

bSeasonal = TRUE
iD = ncol(mD_is)

vSeason_is = vSeason[(iDays/2 * 390 + 1):iT]

## The following code performs filtering over the full sample using
## in sample estimated parameters

###############################################################################################################
########################################### HMM INAR ##########################################################

load(paste(sPath, "output/Fit4_3_8.rdata", sep = ""))

iJ = 4
iL = 3
iK = 8

iMax = max(vN)

vLogFactorial_K = lfactorial(0:(iMax))

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit$mGamma_Eta,
                                             Fit$mGamma_Alpha,
                                             Fit$mOmega,
                                             Fit$vLambda,
                                             Fit$vAlpha,
                                             Fit$vBeta,
                                             Fit$dVarPhi,
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

lPred = list(Filter = Filter, Run = Run, Fit_is = Fit)

save(lPred, file = paste(sPath, "output/lPred_4_3_8.rdata", sep = ""))

###############################################################################################################
########################################### INAR ##########################################################


load(paste(sPath, "output/Fit1_1_1.rdata", sep = ""))

iJ = 1
iL = 1
iK = 1

###

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit$mGamma_Eta,
                                             Fit$mGamma_Alpha,
                                             Fit$mOmega,
                                             Fit$vLambda,
                                             Fit$vAlpha,
                                             Fit$vBeta,
                                             Fit$dVarPhi,
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

lPred_INAR = list(Filter = Filter, Run = Run, Fit_is = Fit)

save(lPred_INAR, file = paste(sPath, "output/lPred_INAR.rdata", sep = ""))

gc()
###############################################################################################################
########################################### Mix INAR ##########################################################

load(paste(sPath, "output/Fit1_1_9.rdata", sep = ""))

iJ = 1
iL = 1
iK = 9

###
Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit$mGamma_Eta,
                                             Fit$mGamma_Alpha,
                                             Fit$mOmega,
                                             Fit$vLambda,
                                             Fit$vAlpha,
                                             Fit$vBeta,
                                             Fit$dVarPhi,
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

lPred_MixINAR = list(Filter = Filter, Run = Run, Fit_is = Fit)

save(lPred_MixINAR, file = paste(sPath, "output/lPred_MixINAR.rdata", sep = ""))

###############################################################################################################
########################################### HMM INAR L=1 ##########################################################

load(paste(sPath, "output/Fit5_1_6.rdata", sep = ""))

iJ = 5
iL = 1
iK = 6

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit$mGamma_Eta,
                                             Fit$mGamma_Alpha,
                                             Fit$mOmega,
                                             Fit$vLambda,
                                             Fit$vAlpha,
                                             Fit$vBeta,
                                             Fit$dVarPhi,
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

lPred_HMMINAR_l1 = list(Filter = Filter, Run = Run, Fit_is = Fit)

save(lPred_HMMINAR_l1, file = paste(sPath, "output/lPred_HMMINAR_l1.rdata", sep = ""))

gc()
###############################################################################################################
########################################### HMM INAR J=1 ##########################################################

load(paste(sPath, "output/Fit1_5_7.rdata", sep = ""))

iJ = 1
iL = 5
iK = 7

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit$mGamma_Eta,
                                             Fit$mGamma_Alpha,
                                             Fit$mOmega,
                                             Fit$vLambda,
                                             Fit$vAlpha,
                                             Fit$vBeta,
                                             Fit$dVarPhi,
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

lPred_HMMINAR_j1 = list(Filter = Filter, Run = Run, Fit_is = Fit)

save(lPred_HMMINAR_j1, file = paste(sPath, "output/lPred_HMMINAR_j1.rdata", sep = ""))

gc()

###############################################################################################################
########################################### HMM INAR K=1 ##########################################################

load(paste(sPath, "output/Fit5_1_1.rdata", sep = ""))

iJ = 5
iL = 1
iK = 1

Run = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vN, iJ, iK, iL,
                                             Fit$mGamma_Eta,
                                             Fit$mGamma_Alpha,
                                             Fit$mOmega,
                                             Fit$vLambda,
                                             Fit$vAlpha,
                                             Fit$vBeta,
                                             Fit$dVarPhi,
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

lPred_HMMINAR_k1 = list(Filter = Filter, Run = Run, Fit_is = Fit)

save(lPred_HMMINAR_k1, file = paste(sPath, "output/lPred_HMMINAR_k1.rdata", sep = ""))

gc()

###############################################################################################################
###################################### PIT and FE results #####################################################

load(paste(sPath, "output/NBsoftplusINGARCH_OOS_1m.rdata", sep = ""))
load(paste(sPath, "output/lPred_4_3_8.rdata", sep = ""))
load(paste(sPath, "output/lPred_INAR.rdata", sep = ""))
load(paste(sPath, "output/lPred_MixINAR.rdata", sep = ""))
load(paste(sPath, "output/lPred_HMMINAR_j1.rdata", sep = ""))
load(paste(sPath, "output/lPred_HMMINAR_l1.rdata", sep = ""))
load(paste(sPath, "output/lPred_HMMINAR_k1.rdata", sep = ""))

iF = length(vN_oos)

vN_Diff = diff(vN)
vN_Diff_78 = diff(vN, 390)

aForecastError = array(0, dim = c(iF, 9, 2), dimnames = list(NULL,
                                                          c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1",
                                                            "HMMINAR_j1", "HMMINAR_k1", "RW", "RW78", "INGARCH"), c("Median", "Mean")))

vY = c(lPred$Run$dY0, lPred$Run$vY)

aForecastError[, "HMMINAR", "Median"] = tail((vY - lPred$Filter$vMedian)[-1], iF)
aForecastError[, "INAR", "Median"] = tail((vY - lPred_INAR$Filter$vMedian)[-1], iF)
aForecastError[, "MixINAR", "Median"] = tail((vY - lPred_MixINAR$Filter$vMedian)[-1], iF)

aForecastError[, "HMMINAR_j1", "Median"] = tail((vY - lPred_HMMINAR_j1$Filter$vMedian)[-1], iF)
aForecastError[, "HMMINAR_l1", "Median"] = tail((vY - lPred_HMMINAR_l1$Filter$vMedian)[-1], iF)
aForecastError[, "HMMINAR_k1", "Median"] = tail((vY - lPred_HMMINAR_k1$Filter$vMedian)[-1], iF)

aForecastError[, "RW", "Median"] = tail(vN_Diff, iF)
aForecastError[, "RW78", "Median"] = tail(vN_Diff_78, iF)
aForecastError[, "INGARCH", "Median"] = tail((vY - NBsoftplusINGARCH_OOS$Median_oos)[-1], iF)

aForecastError[, "HMMINAR", "Mean"] = tail(lPred$Filter$vRes, iF)
aForecastError[, "INAR", "Mean"] = tail(lPred_INAR$Filter$vRes, iF)
aForecastError[, "MixINAR", "Mean"] = tail(lPred_MixINAR$Filter$vRes, iF)

aForecastError[, "HMMINAR_j1", "Mean"] = tail(lPred_HMMINAR_j1$Filter$vRes, iF)
aForecastError[, "HMMINAR_l1", "Mean"] = tail(lPred_HMMINAR_l1$Filter$vRes, iF)
aForecastError[, "HMMINAR_k1", "Mean"] = tail(lPred_HMMINAR_k1$Filter$vRes, iF)

aForecastError[, "RW", "Mean"] = tail(vN_Diff, iF)
aForecastError[, "RW78", "Mean"] = tail(vN_Diff_78, iF)
aForecastError[, "INGARCH", "Mean"] = NBsoftplusINGARCH_OOS$FE_oos

save(aForecastError, file = paste(sPath, "output/aForecastError.rdata", sep = ""))

vModelNames = dimnames(aForecastError)[[2]]

sBenc = "HMMINAR"

vPeriods = c("Full", "Opening", "MidDay", "Closing") #34=150
# in terms of seasons
Opening = 1:10 #first half hour (note that first seasons are 1,2,3 minutes and then 4:5, 6:10, etc.)
MidDay  = 34:58 #from 2.5h to 4.5h
Closing = 76:81 #from 6h to 6.5h
mPeridos = cbind(TRUE, apply(mD_oos[, Opening] == 1, 1, any), apply(mD_oos[, MidDay] == 1, 1, any), apply(mD_oos[, Closing] == 1, 1, any))
colnames(mPeridos) = vPeriods

aMSE = array(0, dim = c(length(vModelNames), length(vPeriods)), dimnames = list(vModelNames, vPeriods))
aDM = array(9999, dim = c(length(vModelNames), length(vPeriods), 2), dimnames = list(vModelNames, vPeriods, c("Test", "Pval")))

for (sPeriod in vPeriods) {

  for (sModel in vModelNames) {

    vFoo = aForecastError[, sModel, "Median"][mPeridos[, sPeriod]]^2
    #1% outlier trimming
    vFoo = vFoo[vFoo<quantile(vFoo, 1-0.01)]
    aMSE[sModel, sPeriod] = mean(vFoo)

    if (sModel != sBenc) {
      vFoo1 = aForecastError[, sBenc, "Median"][mPeridos[, sPeriod]]^2
      vFoo2 = aForecastError[, sModel, "Median"][mPeridos[, sPeriod]]^2

      ## 1% outlier trimming
      vBaz = vFoo1<quantile(vFoo1, 1-0.01) & vFoo2<quantile(vFoo2, 1-0.01)
      vFoo1 = vFoo1[vBaz]
      vFoo2 = vFoo2[vBaz]

      aDM[sModel, sPeriod,] = unlist(f.DMTest(vFoo1, vFoo2)[c("dTest", "dPvalue")])
    }
  }
}

save(aMSE, file = paste(sPath, "output/aMSE.rdata", sep = ""))
aMSE_stand = aMSE

for (sPeriod in vPeriods) {
  aMSE_stand[, sPeriod] = aMSE[, sPeriod]/aMSE[sBenc, sPeriod]
}

mTab = NULL
mPval = NULL

for (sPeriod in vPeriods) {
  mTab = rbind(mTab,  aMSE_stand[, sPeriod])
  mPval = rbind(mPval, aDM[, sPeriod, "Pval"] )
}

mTab_lower = mTab < 1
mTab_upper = mTab > 1

mTab = format(round(mTab, 2), scientific = FALSE)
mTab[mPval < 0.05 & mTab_lower] = paste("\\cellcolor{red}", mTab[mPval < 0.05 & mTab_lower], sep = "")
mTab[mPval < 0.05 & mTab_upper] = paste("\\cellcolor{green}", mTab[mPval < 0.05 & mTab_upper], sep = "")

mTab= cbind(mTab, vPeriods)

mTab[, ncol(mTab)] = paste(mTab[, ncol(mTab)], "\\\\")

write.table(mTab,
            file =  paste(sPath, "output/TradesPred_Median_1m.txt", sep = ""),
            sep = " & ",
            dec = ".",
            quote = FALSE,
            row.names = TRUE
)


############## OUT of Sample PIT ###################
library(hmminar)
aForecastPIT = array(0, dim = c(iF, 7), dimnames = list(NULL, c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")))

aForecastPIT[, "HMMINAR"] = tail(RandomizedPIT(lPred$Run)$vZ, iF)
aForecastPIT[, "INAR"] = tail(RandomizedPIT(lPred_INAR$Run)$vZ, iF)
aForecastPIT[, "MixINAR"] = tail(RandomizedPIT(lPred_MixINAR$Run)$vZ, iF)
aForecastPIT[, "HMMINAR_j1"] = tail(RandomizedPIT(lPred_HMMINAR_j1$Run)$vZ, iF)
aForecastPIT[, "HMMINAR_l1"] = tail(RandomizedPIT(lPred_HMMINAR_l1$Run)$vZ, iF)
aForecastPIT[, "HMMINAR_k1"] = tail(RandomizedPIT(lPred_HMMINAR_k1$Run)$vZ, iF)
aForecastPIT[, "INGARCH"] = tail(NBsoftplusINGARCH_OOS$PIT_oos, iF)

save(aForecastPIT, file = paste(sPath, "output/aForecastPIT.rdata", sep = ""))

vModelNames = dimnames(aForecastPIT)[[2]]
for (sModel in vModelNames) {

  vU = aForecastPIT[,sModel]

  vU[vU > 0.9999] = 0.9999
  vU[vU < 1.0 - 0.9999] = 1.0 - 0.9999

  pdf(file = paste(sPath, "output/PIT_", sModel, "_oos_1m.pdf", sep = ""), width = 5, height = 5)
  par(mar = c(2,4,1,1))
  foo =  BinTest(vU, plot = TRUE, main = "", g = 10, mar = c(1.5, 1.5, 0, 1.5), las = 1, cex.axis = 1.5)

  dev.off()
}

############## ACF Forecast error ###################

iLag = 50
vModels = c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")

mACF = matrix(0, iLag, length(vModels), dimnames = list(NULL, vModels))

for (sModel in vModels) {
  vRes = aForecastError[, sModel, "Mean"]
  mACF[,sModel] = acf(vRes, lag.max = iLag, plot = FALSE)$acf[-1,1,1]
}

for (sModel in vModels) {

  vACF = mACF[,sModel]
  pdf(file = paste(sPath, "output/ACF_forecasterror_", sModel ,"1m.pdf", sep = ""), width = 5, height = 5)
  par(mar = c(2,4,1,1))
  plot(1:length(vACF), vACF, ylim = range(mACF), las = 1, ylab = "", xlab = "", type = "n", cex.axis = 1.5)
  grid(10, 10, "gray80")

  lines(vACF, lwd = 2, lty = 1)

  abline(h = 1.96/sqrt(iF), col = "red")
  abline(h = -1.96/sqrt(iT), col = "red")

  dev.off()
}

rm(list = ls())
