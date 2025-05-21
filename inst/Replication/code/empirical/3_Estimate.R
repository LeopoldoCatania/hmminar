library(hmminar)
library(parallel)

sPath = "XXX"

#Load Data
load(file = paste(sPath, "data/vNT_1min.rdata", sep = ""))
#Load seasonal matrix
load(file = paste(sPath, "data/mD_1min.rdata", sep = ""))
#load extra functions
source(paste(sPath, "code/functions.R", sep = ""))

iW = 48 #Number of cores
## create clusters
cluster = makeCluster(iW)
clusterEvalQ(cluster, library(hmminar))

## Create the opening dummy
iT = length(vN)
vO = numeric(iT)
vO[seq(1, iT, 390)] = 1

# split the sample in two
vN_is = vN[1:(iT/2)]
vO_is = vO[1:(iT/2)]
mD_is = mD[1:(iT/2), ]
iT = length(vN_is)

## all possible combinations
vK = 1:10
vL = 1:7
vJ = 1:7

lComb = list()
iC = 1

for (iL in vL) {
  for (iK in vK) {
    for (iJ in vJ) {
      if(iK >= iL) {
        lComb[[iC]] = c(iJ, iK, iL)
        iC = iC + 1
      }
    }
  }
}

#Export objects
clusterExport(cluster, c("vN_is", "sPath", "mD_is", "vO_is"))

############################################
### THIS IS COMPUTATIONALLY EXPENSIVE !!!###
############################################
# Estimate all combinations
parLapplyLB(cluster, lComb, function(vComb) {

  iJ = vComb[1]
  iK = vComb[2]
  iL = vComb[3]

  sFile = paste(sPath, "output/Fit",
                iJ, "_", iL, "_", iK, ".rdata", sep = "")

  Fit = Estimate_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vN_is, iJ, iK, iL, mD_is, vO_is, dTol = 1e-6)

  save(Fit, file = sFile)

})

### Stop the clusters
stopCluster(cluster)

### Model Selection
vNames = numeric(length(lComb))
for(i in 1:length(lComb)) {
  vComb = lComb[[i]]
  iJ = vComb[1]
  iK = vComb[2]
  iL = vComb[3]
  vNames[i] = paste(iJ, "_", iL, "_", iK, sep = "")
}

## Extract all the information criteria
lIC = matrix(NA, length(vNames), 8,
             dimnames = list(vNames, c("AIC", "BIC", "HQC", "np", "llk", "iJ", "iL", "iK")))

for (i in 1:length(vNames)) {
  sFile = paste(sPath, "output/Fit",
                iJ, "_", iL, "_", iK, ".rdata", sep = "")
  load(sFile)
  lIC[i, c("AIC", "BIC", "HQC", "np", "llk")] = Fit$IC
  lIC[i, "iJ"] = Fit$iJ
  lIC[i, "iL"] = Fit$iL
  lIC[i, "iK"] = Fit$iK
  print(round(i/length(vFiles) * 100, 2))
}

save(lIC, file = paste(sPathDrop, "output/lIC_1m.rdata", sep = ""))

# Select models according to BIC
sCriteria = "BIC"
### Model (J,K,L)
which.min(lIC[, sCriteria]) #iJ = 4, iL = 3, iK = 8
### Model (1,K,L)
which.min(lIC[lIC[,"iJ"] == 1, sCriteria]) #iJ = 1, iL = 5, iK = 7
### Model (J,K,1)
which.min(lIC[lIC[,"iL"] == 1, sCriteria]) #iJ = 5, iL = 1, iK = 6
### Model (1,K,1)
which.min(lIC[lIC[,"iL"] == 1 & lIC[,"iJ"] == 1, sCriteria]) #iJ = 1, iL = 1, iK = 9
### Model (J,1,1)
which.min(lIC[lIC[,"iL"] == 1 & lIC[,"iK"] == 1, sCriteria]) #iJ = 5, iL = 1, iK = 1


# IN SAMPLE PIT PLOT

mPIT = matrix(0, iT-1, 7, dimnames = list(NULL, c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")))
mRes = matrix(0, iT-1, 7, dimnames = list(NULL, c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")))

########################################### HMM INAR ##########################################################
load(paste(sPath, "output/Fit4_3_8.rdata", sep = ""))
mPIT[, "HMMINAR"] = RandomizedPIT(Fit, seed = 123)$vZ
mRes[, "HMMINAR"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit)$vRes

########################################### INAR ##########################################################
load(paste(sPath, "output/Fit1_1_1.rdata", sep = ""))
mPIT[, "INAR"] = RandomizedPIT(Fit, seed = 123)$vZ
mRes[, "INAR"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit)$vRes

########################################### Mix INAR ##########################################################
load(paste(sPath, "output/Fit1_1_9.rdata", sep = ""))
mPIT[, "MixINAR"] = RandomizedPIT(Fit, seed = 123)$vZ
mRes[, "MixINAR"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit)$vRes

########################################### HMM INAR L=1 ##########################################################
load(paste(sPath, "output/Fit5_1_6.rdata", sep = ""))
mPIT[, "HMMINAR_l1"] = RandomizedPIT(Fit, seed = 123)$vZ
mRes[, "HMMINAR_l1"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit)$vRes

########################################### HMM INAR J=1 ##########################################################
load(paste(sPath, "output/Fit1_5_7.rdata", sep = ""))
mPIT[, "HMMINAR_j1"] = RandomizedPIT(Fit, seed = 123)$vZ
mRes[, "HMMINAR_j1"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit)$vRes

########################################### HMM INAR K=1 ##########################################################
load(paste(sPath, "output/Fit5_1_1.rdata", sep = ""))
mPIT[, "HMMINAR_k1"] = RandomizedPIT(Fit, seed = 123)$vZ
mRes[, "HMMINAR_k1"] = Filtering_MSMixInar_TwoChains_AlfaCor(Fit)$vRes

########################################### SoftPlusINGARCH ##########################################################
load(paste(sPath, "output/NBsoftplusINGARCH_IS.rdata", sep = ""))
mPIT[, "INGARCH"] = NBsoftplusINGARCH_IS$PIT_is[-1]
mRes[, "INGARCH"] = NBsoftplusINGARCH_IS$Res[-1]

save(mPIT, file = paste(sPath, "output/mPIT.rdata", sep = ""))
save(mRes, file = paste(sPath, "output/mRes.rdata", sep = ""))

vModelNames = colnames(mPIT)
sPath_Save = paste(sPath, "output/", sep = "")

for (sModel in vModelNames) {

  vU = mPIT[,sModel]

  vU[vU > 0.9999] = 0.9999
  vU[vU < 1.0 - 0.9999] = 1.0 - 0.9999

  pdf(file = paste(sPath_Save, "PIT_", sModel, "_is_1m.pdf", sep = ""), width = 5, height = 5)
  par(mar = c(2,4,1,1))
  foo =  BinTest(vU, plot = TRUE, main = "", g = 10, mar = c(1.5, 1.5, 0, 1.5), las = 1, cex.axis = 1.5)

  dev.off()
}

## Residuals ACF

iLag = 50
vModels = c("HMMINAR", "INAR", "MixINAR", "HMMINAR_l1", "HMMINAR_j1", "HMMINAR_k1", "INGARCH")
mACF = matrix(0, iLag, length(vModels), dimnames = list(NULL, vModels))

for (sModel in vModels) {
  vRes = mRes[, sModel]
  mACF[,sModel] = acf(vRes, lag.max = iLag, plot = FALSE)$acf[-1,1,1]
}

for (sModel in vModels) {

  vACF = mACF[,sModel]
  pdf(file = paste(sPath_Save, "ACF_residuals_", sModel ,"1m.pdf", sep = ""), width = 5, height = 5)
  par(mar = c(2,4,1,1))
  plot(1:length(vACF), vACF, ylim = range(mACF),
       las = 1, ylab = "", xlab = "", type = "n", cex.axis = 1.5)
  grid(10, 10, "gray80")

  lines(vACF, lwd = 2, lty = 1)

  abline(h = 1.96/sqrt(iT), col = "red")
  abline(h = -1.96/sqrt(iT), col = "red")

  dev.off()
}

# SEASONALITY PLOT
# for graphical purposes we scale the seasonality at the 1 minute frequency using
# linear interpolation via the make390seasons() function

mBeta = matrix(0, 390, 7, dimnames = list(NULL, c("HMMINAR", "INAR", "MixINAR", "HMMINARl1", "HMMINARj1", "HMMINARk1", "INGARCH")))

########################################### HMM INAR ##########################################################
load(paste(sPath, "output/Fit4_3_8.rdata", sep = ""))
mBeta[,"HMMINAR"] = make390seasons(Fit$vBeta)

########################################### INAR ##########################################################
load(paste(sPath, "output/Fit1_1_1.rdata", sep = ""))
mBeta[,"INAR"] = make390seasons(Fit$vBeta)

########################################### Mix INAR ##########################################################
load(paste(sPath, "output/Fit1_1_9.rdata", sep = ""))
mBeta[,"MixINAR"] = make390seasons(Fit$vBeta)

########################################### HMM INAR L=1 ##########################################################
load(paste(sPath, "output/Fit5_1_6.rdata", sep = ""))
mBeta[,"HMMINARl1"] = make390seasons(Fit$vBeta)

########################################### HMM INAR J=1 ##########################################################
load(paste(sPath, "output/Fit1_5_7.rdata", sep = ""))
mBeta[,"HMMINARj1"] = make390seasons(Fit$vBeta)

########################################### HMM INAR K=1 ##########################################################
load(paste(sPath, "output/Fit5_1_1.rdata", sep = ""))
mBeta[,"HMMINARk1"] = make390seasons(Fit$vBeta)

########################################### INGARCH ##########################################################
load(paste(sPath, "output/NBsoftplusINGARCH_IS.rdata", sep = ""))
mBeta[, "INGARCH"] = make390seasons(NBsoftplusINGARCH_IS$Fit_is$Filter$vBeta_s)

for (i in 1:ncol(mBeta)) {
  mBeta[, i] = mBeta[, i]/mBeta[1, i]
}

iD = 390

vSeq = seq(0, 390, 30)

vLabels = c("09:30", "10:00", "10:30", "11:00", "11:30",
            "12:00", "12:30", "13:00", "13:30", "14:00",
            "14:30", "15:00", "15:30", "16:00")

vModels = colnames(mBeta)

for(sModel in vModels){

  pdf(paste(sPath_Save, sModel, "_beta_coeff_1m.pdf", sep = ""), width = 5, height = 5)

  par(mar = c(2,4,1,1))
  plot(1:390, log(mBeta[, sModel]), type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 1.5)
  grid(10, 10, col = "gray80")
  lines(1:390, log(mBeta[, sModel]))

  axis(1, at = vSeq, labels = vLabels, cex.axis = 1.5)
  axis(1, at = seq(0, 390, 5), labels = FALSE, tcl = -0.25)
  axis(1, at = seq(0, 390, 1), labels = FALSE, tcl = -0.25/2)

  dev.off()

}

rm(list = ls())
