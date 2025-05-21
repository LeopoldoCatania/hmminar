library(volume)
library(parallel)

FixProb <- function(vP) {
  foo = vP < 1e-3
  if(any(foo)) {
    vP[foo] = 1e-3
  }
  foo = vP > 1 - 1e-3
  if(any(foo)) {
    vP[foo] = 1 - 1e-3
  }
  vP = vP/sum(vP)
  return(vP)
}

FixGamma <- function(mGamma) {
  for (i in 1:nrow(mGamma)) {
    mGamma[i, ] = FixProb(mGamma[i, ])
  }
  return(mGamma)
}

FixOmega <- function(mOmega) {
  for (i in 1:ncol(mOmega)) {
    mOmega[,i] = FixProb(mOmega[, i])
  }
  return(mOmega)
}

FixAlpha <- function(vAlpha) {
  foo = vAlpha < 1e-3
  if (any(foo)) {
    vAlpha[foo] = 1e-3
  }
  foo = vAlpha > 1 - 1e-3
  if (any(foo)) {
    vAlpha[foo] = 1 - 1e-3
  }
  return(vAlpha)
}

# number of replicates
iB = 10000

vToDo = 1:iB

set.seed(123)

sPath = "XXX"

#number of clusters
iG = 6

vSeed = sample(1:1e7, iB, replace = FALSE)

#
vT = c(250, 500, 1000, 5000)

# // iJ is the state space of S^alpha
# // iL is the State space of S^eta
# // iK is the State space of Z

iJ = 2
iK = 2
iL = 2

mD = matrix(1, max(vT), 1)
vBeta = 1

cluster = makeCluster(iG)
clusterExport(cluster, c("vT", "iJ", "iK", "iL", "mD", "vBeta", "sPath", "vSeed",
                         "FixGamma", "FixOmega", "FixAlpha", "FixProb"))

clusterEvalQ(cluster, {library(volume)})

parLapply(cluster, vToDo, function(b) {

  try({
    set.seed(vSeed[b])

    vLambda = c(1, 7)
    vAlpha = c(0.4, 0.9)
    mGamma_Alpha = diag(rep(0.9, 2))
    mGamma_Alpha[!diag(iJ)] = 0.1

    mGamma_Eta = mGamma_Alpha

    mOmega = matrix(NA, nrow = length(vLambda), iL)
    mOmega[, 1] = c(0.7, 0.3)
    mOmega[, 2] = c(0.3, 0.7)

    lSim = Sim_MSMixInarTwoChains(max(vT), mOmega, vLambda, vAlpha, mGamma_Alpha, mGamma_Eta, vBeta, mD)
    vY = lSim$vY

    lPn_TRUE = list(
      "mGamma_Eta" = mGamma_Eta,
      "mGamma_Alpha" = mGamma_Alpha,
      "mOmega" = mOmega,
      "vLambda" = vLambda,
      "vAlpha" = vAlpha,
      "vBeta" = vBeta
    )

    foo = pn2pw_MSMixPoisInar_TwoChains_EM(lPn_TRUE, iJ, iL, iK, iD, bSeasonal = FALSE, bAlphacor = FALSE)
    vPn_TRUE = pw2pn_MSMixPoisInar_TwoChains_EM(foo, iJ, iL, iK, iD, bSeasonal = FALSE, bAlphacor = FALSE, vecout = TRUE)

    for(iT in vT) {

      Fit = Estimate_MSMixPoisInar_TwoChains_Seasonal_EM(vY[1:iT], iJ, iK, iL, mD = NULL, dTol = 1e-7, lPn_Starting = lPn_TRUE)

      Fit$mGamma_Eta = FixGamma(Fit$mGamma_Eta)
      Fit$mGamma_Alpha = FixGamma(Fit$mGamma_Alpha)
      Fit$mOmega = FixOmega(Fit$mOmega)
      Fit$vAlpha = FixAlpha(Fit$vAlpha)

      lPn_Est = list(
        "mGamma_Eta" = (Fit$mGamma_Eta),
        "mGamma_Alpha" = (Fit$mGamma_Alpha),
        "mOmega" = (Fit$mOmega),
        "vLambda" = Fit$vLambda,
        "vAlpha" = (Fit$vAlpha)
      )

      foo = pn2pw_MSMixPoisInar_TwoChains_EM(lPn_Est, iJ, iL, iK, iD, bSeasonal = FALSE, bAlphacor = FALSE)
      vPn_Est = pw2pn_MSMixPoisInar_TwoChains_EM(foo, iJ, iL, iK, iD, bSeasonal = FALSE, bAlphacor = FALSE, vecout = TRUE)

      vEps = vPn_TRUE - vPn_Est

      vSE = Compute_SE_MSMixPoisInar_TwoChains_EM(Fit)
      vTest = (vPn_Est - vPn_TRUE)/vSE
      vP = 2 * pnorm(-abs(vTest))
      vReject = as.numeric(vP < 0.05)

      lOut = list(true = lPn_TRUE,
                  est = lPn_Est,
                  vSE = vSE,
                  vTest = vTest,
                  vP = vP,
                  vReject = vReject,
                  eps = vEps
      )

      save(lOut, file = paste(sPath, "output/T", iT, "_b", b, ".rdata", sep = ""))

    }

    rm(list = c("Fit", "lSim"))
    gc()

  })

})


stopCluster(cluster)
rm(list = ls())
