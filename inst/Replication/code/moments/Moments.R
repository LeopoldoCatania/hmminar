library(volume)
library(numDeriv)

sPath = "XXX"
sPathFig = paste(sPath, "output/", sep = "")
source(paste(sPath, "code/moments/ACF_Plot.R", sep = ""))

ilag_max = 20

######################################################################################

# iJ = 1, iK = 1, iL = 1

iJ = 1
vAlpha_111 = c(0.7)
mGamma_Alpha_111 = matrix(1)

iK = 1
vLambda_111 = c(3)

iL = 1

iM_111 = iJ*(iJ-1) + iL*(iL-1) + (iK-1)*iL + iJ + iK


mOmega_111 = matrix(1)
mGamma_Eta_111 = matrix(c(1), ncol = iL, byrow = TRUE)

lMoments_111 = Unconditional_Moments(mGamma_Eta_111, mGamma_Alpha_111,
                                     mOmega_111, vLambda_111, vAlpha_111, lag = ilag_max)

pdf(file = paste(sPathFig, "ACF_111_10.pdf", sep = ""), width = 5, height = 5)
par(mar = c(2,4,1,1))
ACF_Plot(lMoments_111, cov = FALSE, percentage = FALSE)
dev.off()
lMoments_111$VarDec

iT = 1000
set.seed(123)
lSim = fSim(iT, mOmega_111, vAlpha_111, mGamma_Alpha_111, mGamma_Eta_111, vLambda_111)
vY1 = lSim$vY[51:iT]
min_111=min(vY1);
max_111=max(vY1);

OD_111=lMoments_111$VarDec/lMoments_111$dMean_Y

######################################################################################

# iJ = 2, iK = 1, iL = 1
iJ = 2
vAlpha_211 = c(0.80, 0.565)
mGamma_Alpha_211 = matrix(c(0.85, 0.15,
                        0.15, 0.85), ncol = iJ, byrow = TRUE)

iK = 1
vLambda_211 = c(3)

iL = 1

iM_211 = iJ*(iJ-1) + iL*(iL-1) + (iK-1)*iL + iJ + iK

mOmega_211 = matrix(1)
mGamma_Eta_211 = matrix(c(1), ncol = iL, byrow = TRUE)

lMoments_211 = Unconditional_Moments(mGamma_Eta_211, mGamma_Alpha_211,
                                 mOmega_211, vLambda_211, vAlpha_211, lag = ilag_max)

pdf(file = paste(sPathFig, "ACF_211_10.pdf", sep = ""), width = 5, height = 5)
par(mar = c(2,4,1,1))
ACF_Plot(lMoments_211, cov = FALSE, percentage = FALSE)
dev.off()

lMoments_211$VarDec

OD_211=lMoments_211$VarDec/lMoments_211$dMean_Y

iT = 1000
set.seed(123)
lSim = fSim(iT, mOmega_211, vAlpha_211, mGamma_Alpha_211, mGamma_Eta_211, vLambda_211)
vY1 = lSim$vY[51:iT]
min_211=min(vY1);
max_211=max(vY1);

######################################################################################

# iJ = 1, iK = 2, iL = 1
iJ = 1
vAlpha_121 = c(0.7)
mGamma_Alpha_121 = matrix(1)

iK = 2
vLambda_121 = c(1, 5)

iL = 1

iM_121 = iJ*(iJ-1) + iL*(iL-1) + (iK-1)*iL + iJ + iK


mOmega_121 = matrix(c(0.5,
                  0.5), ncol = iL, byrow = TRUE)

mGamma_Eta_121 = matrix(c(1), ncol = iL, byrow = TRUE)

lMoments_121 = Unconditional_Moments(mGamma_Eta_121, mGamma_Alpha_121,
                                 mOmega_121, vLambda_121, vAlpha_121, lag = ilag_max)

pdf(file = paste(sPathFig, "ACF_121_10.pdf", sep = ""), width = 5, height = 5)
par(mar = c(2,4,1,1))
ACF_Plot(lMoments_121, cov = FALSE, percentage = FALSE)
dev.off()

lMoments_121$VarDec
OD_121=lMoments_121$VarDec/lMoments_121$dMean_Y

iT = 1000
set.seed(123)
lSim = fSim(iT, mOmega_121, vAlpha_121, mGamma_Alpha_121, mGamma_Eta_121, vLambda_121)
vY1 = lSim$vY[51:iT]
min_121=min(vY1);
max_121=max(vY1);

######################################################################################

# iJ = 2, iK = 2, iL = 1
iJ = 2
vAlpha_221 = c(0.80, 0.565)
mGamma_Alpha_221 = matrix(c(0.85, 0.15,
                            0.15, 0.85), ncol = iJ, byrow = TRUE)

iK = 2
vLambda_221 = c(1, 5)

iL = 1

iM_221 = iJ*(iJ-1) + iL*(iL-1) + (iK-1)*iL + iJ + iK


mOmega_221 = matrix(c(0.5,
                      0.5), ncol = iL, byrow = TRUE)

mGamma_Eta_221 = matrix(c(1), ncol = iL, byrow = TRUE)

lMoments_221 = Unconditional_Moments(mGamma_Eta_221, mGamma_Alpha_221,
                                     mOmega_221, vLambda_221, vAlpha_221, lag = ilag_max)
OD_221=lMoments_221$VarDec/lMoments_221$dMean_Y

pdf(file = paste(sPathFig, "ACF_221_10.pdf", sep = ""), width = 5, height = 5)
par(mar = c(2,4,1,1))
ACF_Plot(lMoments_221, cov = FALSE, percentage = FALSE)
dev.off()

lMoments_221$VarDec
iT = 10000
set.seed(123)
lSim = fSim(iT, mOmega_221, vAlpha_221, mGamma_Alpha_221, mGamma_Eta_221, vLambda_221)
vY1 = lSim$vY[51:iT]
min_221=min(vY1);
max_221=max(vY1);

######################################################################################

# iJ = 1, iK = 2, iL = 2
iJ = 1
vAlpha_122 = c(0.7)
mGamma_Alpha_122 = matrix(c(1), ncol = iJ, byrow = TRUE)

iK = 2
vLambda_122 = c(1, 4.65)

iL = 2

iM_122 = iJ*(iJ-1) + iL*(iL-1) + (iK-1)*iL + iJ + iK

mOmega_122 = matrix(c(0.2,0.7,
                      0.8,0.3), ncol = iL, byrow = TRUE)

mGamma_Eta_122 = matrix(c(0.95, 0.05,
                         0.05, 0.95), ncol = iL, byrow = TRUE)

lMoments_122 = Unconditional_Moments(mGamma_Eta_122, mGamma_Alpha_122,
                                     mOmega_122, vLambda_122, vAlpha_122, lag = ilag_max)

pdf(file = paste(sPathFig, "ACF_122_10.pdf", sep = ""), width = 5, height = 5)
par(mar = c(2,4,1,1))
ACF_Plot(lMoments_122, cov = FALSE, percentage = FALSE)
dev.off()

lMoments_122$VarDec
OD_122=lMoments_122$VarDec/lMoments_122$dMean_Y

iT = 1000
set.seed(123)
lSim = fSim(iT, mOmega_122, vAlpha_122, mGamma_Alpha_122, mGamma_Eta_122, vLambda_122)
vY1 = lSim$vY[51:iT]
min_122=min(vY1);
max_122=max(vY1);

######################################################################################

# iJ = 2, iK = 2, iL = 2
iJ = 2
vAlpha_222 = c(0.80, 0.565)
mGamma_Alpha_222 = matrix(c(0.85, 0.15,
                            0.15, 0.85), ncol = iJ, byrow = TRUE)

iK = 2
vLambda_222 = c(1, 5)

iL = 2

iM_222 = iJ*(iJ-1) + iL*(iL-1) + (iK-1)*iL + iJ + iK


mOmega_222 = matrix(c(0.2,0.7,
                      0.8,0.3), ncol = iL, byrow = TRUE)

mGamma_Eta_222 = matrix(c(0.95, 0.05,
                          0.05, 0.95), ncol = iL, byrow = TRUE)

lMoments_222 = Unconditional_Moments(mGamma_Eta_222, mGamma_Alpha_222,
                                     mOmega_222, vLambda_222, vAlpha_222, lag = ilag_max)

pdf(file = paste(sPathFig, "ACF_222_10.pdf", sep = ""), width = 5, height = 5)
par(mar = c(2,4,1,1))
ACF_Plot(lMoments_222, cov = FALSE, percentage = FALSE)
dev.off()

lMoments_222$VarDec
OD_222=lMoments_222$VarDec/lMoments_222$dMean_Y

iT = 10000
set.seed(123)
lSim = fSim(iT, mOmega_222, vAlpha_222, mGamma_Alpha_222, mGamma_Eta_222, vLambda_222)
vY1 = lSim$vY[51:iT]
min_222=min(vY1);
max_222=max(vY1);

######################################################################################

# iJ = 1, iK = 5, iL = 1
iJ = 1
vAlpha_151_c = c(0.7)
mGamma_Alpha_151_c = matrix(1)

iK = 5
vLambda_151_c = c(1, 2,3,15, 26.5)

iL = 1

iM_151_c = iJ*(iJ-1) + iL*(iL-1) + (iK-1)*iL + iJ + iK


mOmega_151_c = matrix(c(0.18,
                        0.22,0.57,0.01,0.02), ncol = iL, byrow = TRUE)

mGamma_Eta_151_c = matrix(c(1), ncol = iL, byrow = TRUE)

lMoments_151_c = Unconditional_Moments(mGamma_Eta_151_c, mGamma_Alpha_151_c,
                                       mOmega_151_c, vLambda_151_c, vAlpha_151_c, lag = ilag_max)

pdf(file = paste(sPathFig, "ACF_151_c_10.pdf", sep = ""), width = 5, height = 5)
par(mar = c(2,4,1,1))
ACF_Plot(lMoments_151_c, cov = FALSE, percentage = FALSE)
dev.off()

lMoments_151_c$VarDec
OD_151_c=lMoments_151_c$VarDec/lMoments_151_c$dMean_Y

iT = 1000
set.seed(123)
lSim = fSim(iT, mOmega_151_c, vAlpha_151_c, mGamma_Alpha_151_c, mGamma_Eta_151_c, vLambda_151_c)
vY1 = lSim$vY[51:iT]
min_151_b=min(vY1)
max_151_b=max(vY1)

