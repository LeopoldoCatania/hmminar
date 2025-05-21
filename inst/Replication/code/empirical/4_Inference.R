library(hmminar)

sPath = "XXX"

#Load Data
load(file = paste(sPath, "data/vNT_1min.rdata", sep = ""))
source(paste(sPath, "code/functions.R", sep = ""))

lInference = list()
lFit = list()

########################################### INAR ##########################################################
load(paste(sPath, "output/Fit1_1_1.rdata", sep = ""))

Fit = StandardizeFit(Fit)
lFit[["INAR"]] = Fit

lInference[["INAR"]] = Compute_SE_MSMixPoisInar_TwoChains_EM(Fit)

########################################### Mix INAR ##########################################################
load(paste(sPath, "output/Fit1_1_9.rdata", sep = ""))

Fit = StandardizeFit(Fit)
lFit[["MixINAR"]] = Fit
lInference[["MixINAR"]] = Compute_SE_MSMixPoisInar_TwoChains_EM(Fit)

########################################### HMM INAR L=1 ##########################################################
load(paste(sPath, "output/Fit5_1_6.rdata", sep = ""))

Fit = StandardizeFit(Fit)
lFit[["HMMINAR_L1"]] = Fit
lInference[["HMMINAR_L1"]] = Compute_SE_MSMixPoisInar_TwoChains_EM(Fit)

########################################### HMM INAR J=1 ##########################################################
load(paste(sPath, "output/Fit1_5_7.rdata", sep = ""))

Fit = StandardizeFit(Fit)
lFit[["HMMINAR_J1"]] = Fit
lInference[["HMMINAR_J1"]] = Compute_SE_MSMixPoisInar_TwoChains_EM(Fit)

########################################### HMM INAR K=1 ##########################################################
load(paste(sPath, "output/Fit5_1_1.rdata", sep = ""))

Fit = StandardizeFit(Fit)
lFit[["HMMINAR_K1"]] = Fit
lInference[["HMMINAR_K1"]] = Compute_SE_MSMixPoisInar_TwoChains_EM(Fit)

########################################### HMM INAR ##########################################################
load(paste(sPath, "output/Fit4_3_8.rdata", sep = ""))

Fit = StandardizeFit(Fit)
lFit[["HMMINAR"]] = Fit
lInference[["HMMINAR"]] = Compute_SE_MSMixPoisInar_TwoChains_EM(Fit)

## save results
save(lInference, file = paste(sPath, "output/lInference.rdata", sep = ""))

## Create tables for the tex files.
for (sModel in names(lFit)) {
  ComputeTables(Fit = lFit[[sModel]], vSD = lInference[[sModel]]$vSE, sPath_tab = paste(sPath, "output/", sep = ""), sModel)
}


