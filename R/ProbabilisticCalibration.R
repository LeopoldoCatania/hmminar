PIT_Test_MixINAR <- function(Fit, vY, iB, mD) {

  dAlpha = Fit$dAlpha
  vLambda = Fit$vLambda
  vOmega = Fit$vOmega
  vBeta = Fit$vBeta

  vLogFactorial_K = lfactorial(0:max(vY));


  ## warnings are related to log(pbinom) which can go to -Inf. This has no effect on the evaluation
  PIT = suppressWarnings(PIT_Bin_MixINAR(iB, vY,  dAlpha, vLambda, vOmega, vBeta, mD, vLogFactorial_K))

  return(PIT)

}

RandomizedPIT <- function(Fit, seed = NULL) {

  if (is.null(seed)) {
    seed = sample(1:1e5, 1)
  }

  set.seed(seed)

  vY  = Fit$vY
  dY0 = Fit$dY0
  iT = length(vY)
  vU = runif(iT)
  vZ = numeric(iT)

  mIndices = Fit$mIndices + 1

  vLogFactorial_K = lfactorial(0:(max(vY)))

  iJ = Fit$iJ
  iK = Fit$iK
  iL = Fit$iL

  iQ = iJ * iK * iL

  vSeason = Fit$vSeason + 1
  vO = Fit$vO
  if (is.null(vO)) vO = rep(0, iT)

  mOmega = Fit$mOmega
  vLambda = Fit$vLambda
  vAlpha = Fit$vAlpha
  dVarPhi = Fit$dVarPhi
  if(is.null(dVarPhi)) dVarPhi = 0.0
  vBeta = Fit$vBeta

  PredictedProb = Fit$PredictedProb
  PredictedProb[1, ] = getDelta(Fit$mGamma_OmegaZB, iQ)
  PredictedProb = PredictedProb[-(iT + 1), ]

  if (iQ == 1) {
    PredictedProb = matrix(PredictedProb)
  }

  # clusterEvalQ(cluster, library(volume))
  # clusterExport(cluster, c("vLambda", "vBeta", "vSeason", "iQ", "PredictedProb", "vY", "dY0", "vAlpha", "mIndices", "vO", "dVarPhi",
  #                          "vLogFactorial_K", "vU"))
  #
  # parSapply(cluster, 1:iT, function(t) {
  for(t in 1:iT) {
    dF_y = 0
    dF_ym1 = 0

    vLambda_s = vLambda * vBeta[vSeason[t]]

    for (q in 1:iQ) {
      if (t == 1) {


        dF_ym1 = dF_ym1 + PredictedProb[t, q] * pPB(vY[t] - 1, dY0, vAlpha[mIndices[q, 1]] * (1 - vO[t]) + vO[t] * dVarPhi, vLambda_s[mIndices[q, 2]],
                                                    vLogFactorial_K, FALSE)

      } else {

        dF_ym1 = dF_ym1 + PredictedProb[t, q] * pPB(vY[t] - 1, vY[t - 1], vAlpha[mIndices[q, 1]] * (1 - vO[t]) + vO[t] * dVarPhi, vLambda_s[mIndices[q, 2]],
                                                    vLogFactorial_K, FALSE)

      }

    }

    dF_y = dF_ym1

    for (q in 1:iQ) {
      if (t == 1) {


        dF_y = dF_y + PredictedProb[t, q] * dPB2(vY[t], dY0, vAlpha[mIndices[q, 1]] * (1 - vO[t]) + vO[t] * dVarPhi, vLambda_s[mIndices[q, 2]],
                                                 vLogFactorial_K, FALSE)

      } else {

        dF_y = dF_y + PredictedProb[t, q] * dPB2(vY[t], vY[t - 1], vAlpha[mIndices[q, 1]] * (1 - vO[t]) + vO[t] * dVarPhi, vLambda_s[mIndices[q, 2]],
                                                 vLogFactorial_K, FALSE)

      }

    }

    # dZ = (1.0 - vU[t]) * dF_ym1 + vU[t] * dF_y
    # return(dZ)
    # })

    vZ[t] = (1.0 - vU[t]) * dF_ym1 + vU[t] * dF_y
  }

  return(list("vZ" = vZ, "seed" = seed))

}

RandomizedPIT_NB <- function(vY, vDelta, dBeta, seed = NULL) {

  if (is.null(seed)) {
    seed = sample(1:1e5, 1)
  }

  set.seed(seed)

  iT = length(vY)
  vU = runif(iT)
  vZ = numeric(iT)

  for (t in 1:iT) {
    dF_y   = pNB(vY[t], vDelta[t], dBeta)
    dF_ym1 = pNB(vY[t] - 1, vDelta[t], dBeta)
    vZ[t] = (1.0 - vU[t]) * dF_ym1 + vU[t] * dF_y
  }

  return(list("vZ" = vZ, "seed" = seed))

}
