Estimate_SumZTNB <- function(vV, vN) {

  zeroes = which(vN == 0)
  if (length(zeroes) > 0) {
    vV = vV[-zeroes]
    vN = vN[-zeroes]
  }

  zeroes = which(vV == 0)
  if (length(zeroes) > 0) {
    vV = vV[-zeroes]
    vN = vN[-zeroes]
  }

  vPw = Starting_SumZTNB(vV, vN)
  iT  = length(vV)

  optimizer = optim(vPw, SumZTNB_Optimizer, vV = vV, vN = vN, iT = iT, method = "BFGS", hessian = TRUE)
  # optimizer = solnp(vPw, SumZTNB_Optimizer, vV = vV, vN = vN, iT = iT)

  vPn = pw2pn_SumZTNB(optimizer$par)

  mInference = InferenceFun_Vol(mHessian = optimizer$hessian * iT, optimizer$par)

  return(list(mInference = mInference, optimizer = optimizer, vPn = vPn))
}

Estimate_SumPois <- function(vV, vN) {

  zeros = unique(c(which(vV == 0), which(vN == 0)))

  vN = vN[-zeros]
  vV = vV[-zeros]

  iT = length(vV)

  dLambda = mean(vV/vN)

  dLLK = sum(dpois(vV, vN + dLambda, log = TRUE)/iT)

  return(dLambda)
}

## iJ number of mixture components
Estimate_MixPoisInar_Seasonal_EM <- function(vY, iJ, mD, dTol = 1e-4, iMaxiter = 1000) {

  Start = Sys.time()

  lPn = Starting_MixPoisInar2(vY, iJ)

  vOmega = lPn$vOmega
  dAlpha = lPn$dAlpha
  vLambda = lPn$vLambda
  vBeta = rep(1, ncol(mD))

  iMax = max(vY)

  vLogFactorial_K = lfactorial(0:(iMax))

  Fit = EMMixInar_Seasonal(vY, iJ, vOmega, dAlpha, vLambda, vBeta,
                           vLogFactorial_K, mD, dTol, iMaxiter)

  dLLK = Fit[["dLLK"]]
  np   = 2 * iJ + ncol(mD)
  iT   = length(vY)

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}

## iJ number of mixture components
Estimate_MixPoisInarP_Seasonal_EM <- function(vY, iJ, mD, iP = 1, dTol = 1e-4, iMaxiter = 1000) {

  Start = Sys.time()

  vPw = Starting_MixPoisInar(vY, iJ)
  lPn = pw2pn_MixPoissonINAR(vPw, iJ, iD = 1)

  vOmega = lPn$vOmega
  dAlpha = lPn$dAlpha
  vLambda = lPn$vLambda
  vBeta = rep(1, ncol(mD))

  iMax = max(vY)

  Fit = EMMixInarP_Seasonal(vY, iJ, iP, vOmega, dAlpha, vLambda, vBeta,
                            mD, dTol, iMaxiter)

  dLLK = Fit[["dLLK"]]
  np   = 2 * iJ + ncol(mD)
  iT   = length(vY)

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}

## iJ number of regimes
## iM number of mixture components
Estimate_MSAlphaMixPoisInar_Seasonal_EM <- function(vY, iJ, iM, mD, dTol = 1e-4, iMaxiter = 1000) {

  Start = Sys.time()

  iT   = length(vY)

  if (is.null(mD)) {
    mD = matrix(1, iT)
  }

  iD = ncol(mD)

  if (iD == 1 & sum(mD) == iT) {
    bSeasonal = FALSE
  } else {
    bSeasonal = TRUE
  }

  lPn = Starting_MSMixPoisInar_TwoChains(vY, iJ, iM, 1, mD, bSeasonal)

  vOmega = lPn$mOmega[, 1]
  vLambda = lPn$vLambda
  vBeta = lPn$vBeta
  mGamma = lPn$mGamma_Alpha
  vAlpha = lPn$vAlpha

  vDelta = getDelta(mGamma, iJ)

  iMax = max(vY)

  vLogFactorial_K = lfactorial(0:(iMax))

  Fit = EM_MSMixInar(vY, iJ, iM, vOmega, vLambda, vAlpha, mGamma, vDelta, vBeta, mD,
                     vLogFactorial_K, dTol, iMaxiter)

  dLLK = Fit[["dLLK"]]
  np   = 2 * iM + ncol(mD) + 2 * iJ + iJ^2 - 2
  iT   = length(vY)

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}

## iJ number of regimes
## iM number of mixture components
Estimate_MSOmegaMixPoisInar_Seasonal_EM <- function(vY, iJ, iM, mD, dTol = 1e-4, iMaxiter = 1000) {

  Start = Sys.time()

  iT   = length(vY)

  if (is.null(mD)) {
    mD = matrix(1, iT)
  }

  iD = ncol(mD)

  if (iD == 1 & sum(mD) == iT) {
    bSeasonal = FALSE
  } else {
    bSeasonal = TRUE
  }

  lPn = Starting_MSMixPoisInar_TwoChains(vY, 1, iM, iJ, mD, bSeasonal)

  mOmega = lPn$mOmega
  vLambda = lPn$vLambda
  vBeta = lPn$vBeta
  mGamma = lPn$mGamma_Eta
  dAlpha = lPn$vAlpha[1]

  vDelta = getDelta(mGamma, iJ)

  iMax = max(vY)

  vLogFactorial_K = lfactorial(0:(iMax))

  Fit = EM_MSMixInar_OmegaChain(vY, iJ, iM, mOmega, vLambda, dAlpha, mGamma, vDelta, vBeta, mD,
                     vLogFactorial_K, dTol, iMaxiter)

  dLLK = Fit[["dLLK"]]
  np   = 2 * iM + ncol(mD) + 2 * iJ + iJ^2 - 2
  iT   = length(vY)

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}

## iJ number of regimes
## iM number of mixture components
Estimate_MSMixPoisInar_v3_Seasonal_EM <- function(vY, iJ, iM, mD, dTol = 1e-4, iMaxiter = 1000) {

  Start = Sys.time()

  vPw = Starting_MixPoisInar(vY, iM)
  lPn = pw2pn_MixPoissonINAR(vPw, iM)

  vOmega = lPn$vOmega
  vLambda = lPn$vLambda
  vBeta = rep(1, ncol(mD))
  dAlpha = lPn$dAlpha

  if (iJ > 1) {
    mGamma = matrix(0.05/(iJ - 1), iJ, iJ)
    diag(mGamma) = 0.95
    mLambda = matrix(rep(vLambda, iJ), ncol = iM, byrow = TRUE)
    vSeq = seq(0.7, 1.3, length.out = iJ)
    for (j in 1:iJ) {
      mLambda[j, ] = mLambda[j, ] * vSeq[j]
    }
  } else {
    mGamma = matrix(1)
  }

  vDelta = rep(1/iJ, iJ)

  iMax = max(vY)

  Fit = EM_MSMixInar_v3(vY, iJ, iM, vOmega, mLambda, dAlpha, mGamma, vDelta, vBeta, mD,
                      dTol, iMaxiter)

  dLLK = Fit[["dLLK"]]
  np   = ncol(mD) + (iM - 1) + iJ * (iJ - 1) + iJ * iM + 1 + (iJ - 1)
  iT   = length(vY)

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}

## iJ number of regimes
## iM number of mixture components
Estimate_MSMixPoisInar_2_Seasonal_EM <- function(vY, iJ, iM, mD, dTol = 1e-4, iMaxiter = 1000) {

  Start = Sys.time()

  vPw = Starting_MixPoisInar(vY, iM)
  lPn = pw2pn_MixPoissonINAR(vPw, iM, iD = 1)

  vOmega = lPn$vOmega
  vLambda = lPn$vLambda
  vBeta = rep(1, ncol(mD))
  dAlpha = lPn$dAlpha

  if (iJ > 1) {
    vMu = seq(0.5, 1.5, length.out = iJ)
    vMu = vMu/sum(vMu)
    mGamma = matrix(0.05/(iJ - 1), iJ, iJ)
    diag(mGamma) = 0.95

  } else {
    mGamma = matrix(1)
    vMu = 1
  }

  vDelta = rep(1/iJ, iJ)

  iMax = max(vY)

  Fit = EM_MSMixInar_v2(vY, iJ, iM, vOmega, vLambda, vMu, dAlpha, mGamma, vDelta, vBeta, mD,
                        dTol, iMaxiter)

  dLLK = Fit[["dLLK"]]
  np   = 2 * iM + ncol(mD) + 2 * iJ + iJ^2 - 2
  iT   = length(vY)

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}


## iJ number of regimes
## iM number of mixture components
Estimate_MSMixPoisInar_2_Seasonal_Inv_EM <- function(vY, iJ, iM, mD, dTol = 1e-4, iMaxiter = 1000) {

  Start = Sys.time()

  vPw = Starting_MixPoisInar(vY, iM)
  lPn = pw2pn_MixPoissonINAR(vPw, iM, iD = 1)

  vOmega = lPn$vOmega
  vLambda = lPn$vLambda
  vBeta = rep(1, ncol(mD))
  dAlpha = lPn$dAlpha

  if (iJ > 1) {
    vMu = seq(0.5, 1.5, length.out = iJ)
    vMu = vMu/sum(vMu)
    mGamma = matrix(0.05/(iJ - 1), iJ, iJ)
    diag(mGamma) = 0.95

  } else {
    mGamma = matrix(1)
    vMu = 1
  }

  vDelta = rep(1/iJ, iJ)

  iMax = max(vY)

  Fit = EM_MSMixInar_v2_Inv(vY, iJ, iM, vOmega, vLambda, vMu, dAlpha, mGamma, vDelta, vBeta, mD,
                        dTol, iMaxiter)

  dLLK = Fit[["dLLK"]]
  np   = 2 * iM + ncol(mD) + 2 * iJ + iJ^2 - 2
  iT   = length(vY)

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}


## iJ number of regimes for alpha
## iK number of mixture components
## iL number of regimes for omega
Estimate_MSMixPoisInar_TwoChains_Seasonal_EM <- function(vY, iJ, iK, iL, mD = NULL,
                                                         dTol = 1e-5, iMaxiter = 1000, lPn_Starting = NULL) {

  Start = Sys.time()

  iT   = length(vY)

  if (is.null(mD)) {
    mD = matrix(1, iT)
    bSeasonal = FALSE
  }

  iD = ncol(mD)

  if (iD == 1 & sum(mD) == iT) {
    iFixSeason = 0
  } else {
    bSeasonal = TRUE
  }

  vSeason = numeric(iT)
  if (bSeasonal) {
    for (t in 1:iT) {
      for (d in 1:iD) {
        if ( abs(mD[t, d] - 1.0) < 1e-5) {
          vSeason[t] = d;
        }
      }
    }
    vSeason = vSeason - 1
  }

  iMax = max(vY)

  vLogFactorial_K = lfactorial(0:iMax)

  if(is.null(lPn_Starting)) {
    lPn = Starting_MSMixPoisInar_TwoChains(vY, iJ, iK, iL, mD, vSeason, bSeasonal, vLogFactorial_K)
  } else{
    lPn = lPn_Starting
  }


  Fit = EM_MSMixInarTwoChains_Dummy_C2(vY, iJ, iK, iL,
                                      lPn$mGamma_Eta,
                                      lPn$mGamma_Alpha,
                                      lPn$mOmega,
                                      lPn$vLambda,
                                      lPn$vAlpha,
                                      lPn$vBeta,
                                      vSeason,
                                      vLogFactorial_K,
                                      maxIter = iMaxiter, tol = dTol, bSeasonal)



  dLLK = Fit[["dLLK"]]
  np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC
  Fit[["iJ"]]   = iJ
  Fit[["iK"]]   = iK
  Fit[["iL"]]   = iL
  Fit[["iD"]]   = iD
  Fit[["mD"]]   = mD
  Fit[["bSeasonal"]]   = bSeasonal
  Fit[["bAlphacor"]]   = FALSE
  Fit[["bParSeason"]]  = FALSE

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}

Estimate_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor <- function(vY, iJ, iK, iL, mD = NULL, vO = NULL,
                                                                 dTol = 1e-5, iMaxiter = 1000, lPn = NULL) {

  Start = Sys.time()

  iT   = length(vY)

  if (is.null(mD)) {
    mD = matrix(1, iT)
    bSeasonal = FALSE
  }

  if (is.null(vO)) {
    vO = numeric(iT)
  }

  iD = ncol(mD)

  if (iD == 1 & sum(mD) == iT) {
    iFixSeason = 0
  } else {
    bSeasonal = TRUE
  }

  vSeason = numeric(iT)
  if (bSeasonal) {
    for (t in 1:iT) {
      for (d in 1:iD) {
        if ( abs(mD[t, d] - 1.0) < 1e-5) {
          vSeason[t] = d;
        }
      }
    }
    vSeason = vSeason - 1
  }

  iMax = max(vY)

  vLogFactorial_K = lfactorial(0:(iMax))

  if(is.null(lPn)) {
    lPn = Starting_MSMixPoisInar_TwoChains(vY, iJ, iK, iL, mD, vSeason, bSeasonal, vLogFactorial_K)
  }

  Fit = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vY, iJ, iK, iL,
                                       lPn$mGamma_Eta,
                                       lPn$mGamma_Alpha,
                                       lPn$mOmega,
                                       lPn$vLambda,
                                       lPn$vAlpha,
                                       lPn$vBeta,
                                       lPn$dVarPhi,
                                       vSeason,
                                       vLogFactorial_K,
                                       vO,
                                       maxIter = iMaxiter, tol = dTol, bSeasonal)

  lPn = list("vLambda" = as.numeric(Fit$vLambda),
             "vAlpha" = as.numeric(Fit$vAlpha),
             "dVarPhi" = Fit$dVarPhi,
             "mOmega" = Fit$mOmega,
             "mGamma_Eta" = Fit$mGamma_Eta,
             "mGamma_Alpha" = Fit$mGamma_Alpha,
             "vBeta" = Fit$vBeta)


  dLLK = Fit[["dLLK"]]
  np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD + 1

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC
  Fit[["iJ"]]   = iJ
  Fit[["iK"]]   = iK
  Fit[["iL"]]   = iL
  Fit[["iD"]]   = iD
  Fit[["mD"]]   = mD
  Fit[["lPn"]]  = lPn
  Fit[["bSeasonal"]]   = bSeasonal
  Fit[["bParSeason"]]  = FALSE
  Fit[["bAlphacor"]]   = TRUE

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}

## Omega = mI, iL = iK

Estimate_MSMixPoisInar_TwoChains_OmegaI_Seasonal_EM_AlfaCor <- function(vY, iJ, iK, mD = NULL, vO = NULL,
                                                                 dTol = 1e-5, iMaxiter = 1000, lPn = NULL) {

  Start = Sys.time()

  iT   = length(vY)

  if (is.null(mD)) {
    mD = matrix(1, iT)
    bSeasonal = FALSE
  }

  if (is.null(vO)) {
    vO = numeric(iT)
  }

  iD = ncol(mD)

  if (iD == 1 & sum(mD) == iT) {
    iFixSeason = 0
  } else {
    bSeasonal = TRUE
  }

  vSeason = numeric(iT)
  if (bSeasonal) {
    for (t in 1:iT) {
      for (d in 1:iD) {
        if ( abs(mD[t, d] - 1.0) < 1e-5) {
          vSeason[t] = d;
        }
      }
    }
    vSeason = vSeason - 1
  }

  iMax = max(vY)

  vLogFactorial_K = lfactorial(0:(iMax))

  if(is.null(lPn)) {
    lPn = Starting_MSMixPoisInar_TwoChains(vY, iJ, iK, iL = iK, mD, vSeason, bSeasonal, vLogFactorial_K)
  }

  lPn$mOmega = diag(iK)

  Fit = EM_MSMixInarTwoChains_Dummy_C_AlphaCor_OmegaIdentity(vY, iJ, iK, iL = iK,
                                       lPn$mGamma_Eta,
                                       lPn$mGamma_Alpha,
                                       lPn$mOmega,
                                       lPn$vLambda,
                                       lPn$vAlpha,
                                       lPn$vBeta,
                                       lPn$dVarPhi,
                                       vSeason,
                                       vLogFactorial_K,
                                       vO,
                                       maxIter = iMaxiter, tol = dTol, bSeasonal)

  lPn = list("vLambda" = as.numeric(Fit$vLambda),
             "vAlpha" = as.numeric(Fit$vAlpha),
             "dVarPhi" = Fit$dVarPhi,
             "mOmega" = Fit$mOmega,
             "mGamma_Eta" = Fit$mGamma_Eta,
             "mGamma_Alpha" = Fit$mGamma_Alpha,
             "vBeta" = Fit$vBeta)


  dLLK = Fit[["dLLK"]]
  np   = iJ + iJ * (iJ - 1) + iK + iK * (iK - 1) + iD + 1

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC
  Fit[["iJ"]]   = iJ
  Fit[["iK"]]   = iK
  Fit[["iD"]]   = iD
  Fit[["mD"]]   = mD
  Fit[["lPn"]]  = lPn
  Fit[["bSeasonal"]]   = bSeasonal
  Fit[["bParSeason"]]  = FALSE
  Fit[["bAlphacor"]]   = TRUE

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}




Estimate_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor_Direct <- function(vY, iJ, iK, iL, mD = NULL, vO = NULL,
                                                                        lPn = NULL) {

  Start = Sys.time()

  iT   = length(vY)

  if (is.null(mD)) {
    mD = matrix(1, iT)
    bSeasonal = FALSE
  }

  if (is.null(vO)) {
    vO = numeric(iT)
  }

  iD = ncol(mD)

  if (iD == 1 & sum(mD) == iT) {
    iFixSeason = 0
  } else {
    bSeasonal = TRUE
  }

  vSeason = numeric(iT)
  if (bSeasonal) {
    for (t in 1:iT) {
      for (d in 1:iD) {
        if ( abs(mD[t, d] - 1.0) < 1e-5) {
          vSeason[t] = d;
        }
      }
    }
    vSeason = vSeason - 1
  }

  iMax = max(vY)

  vLogFactorial_K = lfactorial(0:(iMax))

  if(is.null(lPn)) {
    lPn = Starting_MSMixPoisInar_TwoChains(vY, iJ, iK, iL, mD, vSeason, bSeasonal, vLogFactorial_K)
  }

  vPw = pn2pw_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(lPn, iJ, iL, iK, iD)

  optimizer = optim(vPw, function(vPw, vY, iJ, iK, iL, vSeason,
                                  vLogFactorial_K,
                                  vO) {

    lPn_foo = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw, iJ, iL, iK, iD, vecout = FALSE)

    dOut = -EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vY, iJ, iK, iL,
                                                   lPn_foo$mGamma_Eta,
                                                   lPn_foo$mGamma_Alpha,
                                                   lPn_foo$mOmega,
                                                   lPn_foo$vLambda,
                                                   lPn_foo$vAlpha,
                                                   lPn_foo$vBeta,
                                                   lPn_foo$dVarPhi,
                                                   vSeason,
                                                   vLogFactorial_K,
                                                   vO,
                                                   maxIter = 1, tol = 1e-4, bSeasonal)$dLLK

    return(dOut)
  }, method = "BFGS", hessian = FALSE, vY = vY, iJ = iJ, iK = iK, iL = iL,
  vSeason = vSeason,
  vLogFactorial_K = vLogFactorial_K,
  vO = vO)

  vPw = optimizer$par
  vPn = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw, iJ, iL, iK, iD, vecout = TRUE)
  lPn = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw, iJ, iL, iK, iD, vecout = FALSE)

  dLLK = -optimizer$value
  np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD + 1

  IC = ICfun(dLLK, np, iT)

  Fit = list()
  Fit[["optimizer"]] = optimizer
  Fit[["vPw"]] = vPw
  Fit[["vPn"]] = vPn
  Fit[["lPn"]] = lPn
  Fit[["dLLK"]] = dLLK
  Fit[["IC"]] = IC

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}
#
Compute_SE_MSMixPoisInar_TwoChains_EM <- function(Fit,
                                                  method = c("simple", "Richardson"),
                                                  method.args=list(eps=1e-1, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
                                                                   r=4, v=2, show.details=FALSE)) {

  lPn = Fit$lPn

  vY = c(Fit$dY0, Fit$vY)
  iT = length(vY)
  iJ = Fit$iJ
  iK = Fit$iK
  iL = Fit$iL
  iD = Fit$iD
  iH = Fit$iH

  vLogFactorial_K = lfactorial(0:(max(vY)))

  bSeasonal = Fit$bSeasonal
  bAlphacor = Fit$bAlphacor
  bParSeason = FALSE

  if(bSeasonal) {
    vSeason = c(Fit$iSeason0, Fit$vSeason)
  } else {
    vSeason = rep(0, iT)
  }

  if(bAlphacor) {
    vO = c(Fit$dO0, Fit$vO)
  } else {
    vO = rep(0, iT)
  }

  vPw = pn2pw_MSMixPoisInar_TwoChains_EM(lPn, iJ, iL, iK, iD, bSeasonal, bAlphacor, bParSeason)
  vPn = pw2pn_MSMixPoisInar_TwoChains_EM(vPw, iJ, iL, iK, iD, bSeasonal, bAlphacor, vecout = TRUE, bParSeason, iH)

  # if(!bParSeason) {

    vPw_fix = c("beta_1" = 0)
    vPn_fix = c("beta_1" = 1)
    vNames_Fixpw = names(vPw_fix)
    vNames_Fixpn = names(vPn_fix)

    vPw = vPw[-which(names(vPw) %in% vNames_Fixpw)]
    vPn = vPn[-which(names(vPn) %in% vNames_Fixpn)]

  # } else {
  #   vPw_fix = NULL
  #   vPn_fix = NULL
  #   vNames_Fixpw = NULL
  #   vNames_Fixpn = NULL
  # }

  mS = numDeriv::jacobian(function(vPw_foo, vY, iJ, iL, iK, iD, iH, vSeason,
                                   vLogFactorial_K, vO, bSeasonal, bAlphacor, vPw_fix) {

    if(!is.null(vPw_fix)) {
      vPw_foo = c(vPw_foo, vPw_fix)
    }

    lPn_foo = pw2pn_MSMixPoisInar_TwoChains_EM(vPw_foo, iJ, iL, iK, iD, bSeasonal, bAlphacor, vecout = FALSE, bParSeason = FALSE, iH)

    if(bAlphacor) {
      dVarPhi = lPn_foo$dVarPhi
    } else {
      dVarPhi = 999
    }
    if(bSeasonal) {
        vSeasonals = lPn_foo$vBeta
    } else {
      vSeasonals = 0
    }

    Filter = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vY, iJ, iK, iL,
                                                    lPn_foo$mGamma_Eta,
                                                    lPn_foo$mGamma_Alpha,
                                                    lPn_foo$mOmega,
                                                    as.numeric(lPn_foo$vLambda),
                                                    as.numeric(lPn_foo$vAlpha),
                                                    as.numeric(vSeasonals),
                                                    dVarPhi,
                                                    vSeason,
                                                    vLogFactorial_K,
                                                    vO,
                                                    maxIter = 1, tol = 1e-4, bSeasonal)

    mLLK = t(Filter$mLLK)
    mPredictedProb = Filter$PredictedProb[1:nrow(mLLK),]

    vLLK = apply(mLLK + log(mPredictedProb), 1, LogSumExp)

    return(vLLK)

  }, vPw,
  vY = vY, iJ = iJ, iL = iL, iK = iK, iD = iD,
  vSeason = vSeason, vPw_fix = vPw_fix,
  vLogFactorial_K = vLogFactorial_K, vO = vO,
  bSeasonal = bSeasonal, bAlphacor = bAlphacor,
  method = method[1], method.args = method.args)

  colnames(mS) = names(vPw)
  vToFix = which(apply(mS, 2, var) < 1e-7)

  if(length(vToFix) > 0) {
    vPw_fix = c(vPw_fix, vPw[vToFix])
    vPn_fix = c(vPn_fix, vPn[vToFix])
  }

  mS_full = mS

  if(length(vToFix) > 0) {
    mS = mS_full[, -vToFix]
    vPw = vPw[-which(names(vPw) %in% names(vPw_fix))]
    vPn = vPn[-which(names(vPn) %in% names(vPn_fix))]
  }

  mJ = numDeriv::jacobian(function(vPw_foo, vPw_fix, vPn_fix, iJ, iK, iL, iD) {

    vPw_foo = c(vPw_foo, vPw_fix)

    vPn_foo = pw2pn_MSMixPoisInar_TwoChains_EM(vPw_foo, iJ, iL, iK, iD, bSeasonal, bAlphacor, vecout = TRUE)

    vPn_foo = vPn_foo[-which(names(vPn_foo) %in% names(vPn_fix))]

    return(vPn_foo)

  }, vPw, iJ = iJ, iK = iK, iL = iL, iD = iD, vPw_fix = vPw_fix, vPn_fix = vPn_fix,
  method = "Richardson")
  # method = "simple")
  #

  mH = (t(mS) %*% mS)
  mI_inv = mJ %*% (solve(mH) * iT) %*% t(mJ)
  #
  vSE = sqrt(diag(mI_inv))/sqrt(iT)
  names(vSE) = names(vPn)

  lOut = list(vSE = vSE, mS = mS, mJ = mJ, vPw_fix = vPw_fix, vPn_fix = vPn_fix, mS_full = mS_full)

  return(lOut)

}

Compute_SE_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor <- function(Fit, eps = 1e-4) {

  lPn = getlPar(Fit)

  vY = Fit$vY
  iT = length(vY)
  iJ = Fit$iJ
  iK = Fit$iK
  iL = Fit$iL
  iD = Fit$iD
  vSeason = Fit$vSeason
  vLogFactorial_K = lfactorial(0:(max(vY)))
  vO = Fit$vO
  bSeasonal = TRUE

  vPw = pn2pw_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(lPn, iJ, iL, iK, iD)

  mS = matrix(NA, iT-1, length(vPw))

  for (i in 1:length(vPw)) {

    print(paste(round(i/length(vPw) * 100, 2), Sys.time()))

    vPw_plus = vPw
    vPw_minus = vPw

    vPw_plus[i] = vPw_plus[i] + eps
    vPw_minus[i] = vPw_minus[i] - eps

    lPn_plus = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw_plus, iJ, iL, iK, iD)
    lPn_minus = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw_minus, iJ, iL, iK, iD)

    Filter_plus = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vY, iJ, iK, iL,
                                                       lPn_plus$mGamma_Eta,
                                                       lPn_plus$mGamma_Alpha,
                                                       lPn_plus$mOmega,
                                                       as.numeric(lPn_plus$vLambda),
                                                       as.numeric(lPn_plus$vAlpha),
                                                       as.numeric(lPn_plus$vBeta),
                                                       lPn_plus$dVarPhi,
                                                       vSeason,
                                                       vLogFactorial_K,
                                                       vO,
                                                       maxIter = 1, tol = 1e-4, bSeasonal)

    mLLK_plus = t(Filter_plus$mLLK)
    mPredictedProb_plus = Filter_plus$PredictedProb[1:nrow(mLLK_plus),]

    vLLK_plus = apply(mLLK_plus + mPredictedProb_plus, 1, LogSumExp)


    Filter_minus = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vY, iJ, iK, iL,
                                                          lPn_minus$mGamma_Eta,
                                                          lPn_minus$mGamma_Alpha,
                                                          lPn_minus$mOmega,
                                                         as.numeric(lPn_minus$vLambda),
                                                         as.numeric(lPn_minus$vAlpha),
                                                         as.numeric(lPn_minus$vBeta),
                                                         lPn_minus$dVarPhi,
                                                         vSeason,
                                                         vLogFactorial_K,
                                                         vO,
                                                         maxIter = 1, tol = 1e-4, bSeasonal)

    mLLK_minus = t(Filter_minus$mLLK)
    mPredictedProb_minus = Filter_minus$PredictedProb[1:nrow(mLLK_minus),]

    vLLK_minus = apply(mLLK_minus + mPredictedProb_minus, 1, LogSumExp)

    mS[, i] = (vLLK_plus - vLLK_minus)/(2*eps)

  }

  #
  mJ = numDeriv::jacobian(function(vPw, iJ, iK, iL, iD) {

    vPn = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw, iJ, iL, iK, iD, vecout = TRUE)

    return(vPn)

  }, vPw, iJ = iJ, iK = iK, iL = iL, iD = iD,
  # method = "Richardson")
  method = "simple")
  #
  # mH = (t(mS) %*% mS)/iT
  # mI_inv = mJ %*% solve(mH) %*% t(mJ)
  # vSE = sqrt(diag(mI_inv))/sqrt(iT)

  lOut = list(mS = mS, mJ = mJ)

  return(lOut)

}


##### Parametric seasonality
# The seasonality is non parametric for the first H seasons
# After from season H+1 till season D we employ a parametric specification
# H = 0 means no parametric specification

Estimate_MSMixPoisInar_TwoChains_ParSeasonal_EM_AlfaCor <- function(vY, iJ, iK, iL, iD, iH = 0, vO = NULL,
                                                                    lPn = NULL) {

  Start = Sys.time()

  iT   = length(vY)

  if (is.null(vO)) {
    vO = numeric(iT)
  }

  vSeason = rep(1:iD, round(iT/iD) + 1)
  vSeason = vSeason[1:iT]
  vSeason = vSeason - 1

  iMax = max(vY)

  vLogFactorial_K = lfactorial(0:(iMax))

  if(is.null(lPn)) {
    lPn = Starting_MSMixPoisInar_TwoChains(vY, iJ, iK, iL, mD = NULL, vSeason, bSeasonal = TRUE,
                                           vLogFactorial_K, bParSeason = TRUE, iH = iH)
  }

  vPw = pn2pw_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(lPn, iJ, iL, iK, iD = NULL, bParSeason = TRUE)

  optimizer = optim(vPw, function(vPw, vY, iJ, iK, iL, iD, iH, vSeason,
                                  vLogFactorial_K,
                                  vO) {

    lPn_foo = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw, iJ, iL, iK, iD,
                                                                vecout = FALSE, bParSeason = TRUE)
    vLogDelta_S = numeric(iD);

    dOmega0 = vPw["omega0"]
    dOmega1 = vPw["omega1"]
    dOmega2 = vPw["omega2"]

    for(s in 1:iD) {
      vLogDelta_S[s] = dOmega0 + dOmega1 * cos(2.0 * pi * s/iD - dOmega2 * pi);
    }
    if(iH > 0) {
      vLogDelta_S[1:iH] = vPw[paste("w", 1:iH, sep = "")]
    }

    vDelta_S = exp(vLogDelta_S)

    dOut = 1e10
    try({
      dOut = -EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vY, iJ, iK, iL,
                                                     lPn_foo$mGamma_Eta,
                                                     lPn_foo$mGamma_Alpha,
                                                     lPn_foo$mOmega,
                                                     lPn_foo$vLambda,
                                                     lPn_foo$vAlpha,
                                                     vDelta_S,
                                                     lPn_foo$dVarPhi,
                                                     vSeason,
                                                     vLogFactorial_K,
                                                     vO,
                                                     maxIter = 1, tol = 1e-4, bSeasonal = TRUE)$dLLK
      if(is.nan(dOut)) {
        dOut = 1e10
      }
      if(!is.finite(dOut)) {
        dOut = 1e10
      }
    })

    return(dOut)
  }, method = "BFGS", hessian = FALSE, vY = vY, iJ = iJ, iK = iK, iL = iL, iD = iD,
  iH = iH,
  vSeason = vSeason,
  vLogFactorial_K = vLogFactorial_K,
  vO = vO)

  vPw = optimizer$par
  vPn = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw, iJ, iL, iK, iD, vecout = TRUE, bParSeason = TRUE, iH = iH)
  lPn = pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor(vPw, iJ, iL, iK, iD, vecout = FALSE, bParSeason = TRUE, iH = iH)

  vLogDelta_S = numeric(iD);
  dOmega0 = vPn["omega0"]
  dOmega1 = vPn["omega1"]
  dOmega2 = vPn["omega2"]

  for(s in 1:iD) {
    vLogDelta_S[s] = dOmega0 + dOmega1 * cos(2.0 * pi * s/iD - dOmega2 * pi);
  }
  if(iH > 0) {
    vLogDelta_S[1:iH] = vPw[paste("w", 1:iH, sep = "")]
  }

  vDelta_S = exp(vLogDelta_S)

  Fit = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vY, iJ, iK, iL,
                                               lPn$mGamma_Eta,
                                               lPn$mGamma_Alpha,
                                               lPn$mOmega,
                                               lPn$vLambda,
                                               lPn$vAlpha,
                                               vDelta_S,
                                               lPn$dVarPhi,
                                               vSeason,
                                               vLogFactorial_K,
                                               vO,
                                               maxIter = 1, tol = 1e-4, bSeasonal = TRUE)
  dLLK = Fit[["dLLK"]]
  np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + 3 + 1

  IC = ICfun(dLLK, np, iT)

  Fit[["iT"]]   = iT
  Fit[["IC"]]   = IC
  Fit[["iJ"]]   = iJ
  Fit[["iK"]]   = iK
  Fit[["iL"]]   = iL
  Fit[["iD"]]   = iD
  Fit[["iH"]]   = iH
  Fit[["lPn"]]  = lPn
  Fit[["bSeasonal"]]   = TRUE
  Fit[["bParSeason"]]        = TRUE
  Fit[["bAlphacor"]]   = TRUE
  Fit[["vLogDelta_S"]] = vLogDelta_S
  Fit[["vDelta_S"]] = vDelta_S

  dTime = Sys.time() - Start
  Fit[["Time"]] = dTime

  return(Fit)

}

BetaStep <- function(Fit, vBeta, dTol = 1e-5, iMaxiter = 1000) {

  vY = c(Fit$dY0, Fit$vY)
  vO = c(Fit$dO0, Fit$vO)

  vSeason = c(Fit$iSeason0, Fit$vSeason)

  if (max(vSeason > 0)) {
    bSeasonal = TRUE
  } else {
    bSeasonal = FALSE
  }

  vLogFactorial_K = lfactorial(0:max(vY))

  iJ = Fit$iJ
  iK = Fit$iK
  iL = Fit$iL
  iT = Fit$iT
  iD = Fit$iD

  Fit_New = EM_MSMixInarTwoChains_Dummy_C_AlphaCor(vY, Fit$iJ, Fit$iK, Fit$iL,
                                                   Fit$mGamma_Eta,
                                                   Fit$mGamma_Alpha,
                                                   Fit$mOmega,
                                                   Fit$vLambda,
                                                   Fit$vAlpha,
                                                   vBeta,
                                                   Fit$dVarPhi,
                                                   vSeason,
                                                   vLogFactorial_K,
                                                   vO,
                                                   maxIter = iMaxiter, tol = dTol, bSeasonal = bSeasonal)

  dLLK = Fit_New[["dLLK"]]
  np   = iJ + iJ * (iJ - 1) + iK + iK * iL + iL * (iL - 1) + iD

  IC = ICfun(dLLK, np, iT)

  Fit_New[["iT"]]   = iT
  Fit_New[["IC"]]   = IC
  Fit_New[["iJ"]]   = iJ
  Fit_New[["iK"]]   = iK
  Fit_New[["iL"]]   = iL
  Fit_New[["iD"]]   = iD
  Fit_New[["mD"]]   = Fit$mD

  return(Fit_New)


}

optimizerPNBAR <- function(vPar, vY, iD) {

  dDelta0 = vPar["delta0"]
  dBeta = vPar["beta"]
  dOmega0 = vPar["omega0"]
  dOmega1 = vPar["omega1"]
  dOmega2 = vPar["omega2"]

  dnLLK = try(-PNBAR(vY, dBeta, dDelta0, dOmega0, dOmega1, dOmega2, iD)$dLLK, silent = TRUE)

  if(is(dnLLK, "try-error")) {
    dnLLK = 1e4
  }

  return(dnLLK)
}

Estimate_PNBAR <- function(vY, iD) {

  vPar = c("delta0" = mean(vY) * 0.5, "beta" = 0.5, "omega0" = mean(vY) * 0.5, "omega1" = 0, "omega2" = 0.1)

  LB = c(0.00001, 0.00001, -50, -5, 0.000001)
  UB = c(1e5, 0.9999, 10, 5, 0.9999)

  optimizer = optim(vPar, optimizerPNBAR, lower = LB, upper = UB, vY = vY, iD = iD, method = "L-BFGS-B")

  dLLK = optimizer$value
  vPar = optimizer$par

  Filter = PNBAR(vY, vPar["beta"], vPar["delta0"], vPar["omega0"], vPar["omega1"], vPar["omega2"], iD)

  return(list(Filter = Filter, dLLK = dLLK, vPar = vPar))

}

optimizerNBAR <- function(vPar, vY) {

  dDelta0 = vPar["delta0"]
  dBeta = vPar["beta"]
  dOmega0 = vPar["omega0"]

  dnLLK = try(-NBAR(vY, dBeta, dDelta0, dOmega0)$dLLK, silent = TRUE)

  if(is(dnLLK, "try-error")) {
    dnLLK = 1e4
  }

  return(dnLLK)
}

Estimate_NBAR <- function(vY) {

  vPar = c("delta0" = mean(vY) * 0.7, "beta" = 0.7, "omega0" = mean(vY) * 0.7)

  LB = c(0.00001, 0.00001, -50)
  UB = c(1e3, 0.9999, 10)

  optimizer = optim(vPar, optimizerNBAR, lower = LB, upper = UB, method = "L-BFGS-B", vY = vY)

  dLLK = optimizer$value
  vPar = optimizer$par

  Filter = NBAR(vY, vPar["beta"], vPar["delta0"], vPar["omega0"])

  return(list(Filter = Filter, dLLK = dLLK, vPar = vPar))

}


## NB-softplus INGARCH with seasonal components
optimizerNBsoftplusINGARCH <- function(vPar, vY, iD, dC, iH) {

  dM0 = vPar["m0"]
  dAlpha0 = vPar["alpha0"]
  dAlpha1 = vPar["alpha1"]
  dBeta1 = vPar["beta1"]
  dP = vPar["p"]


  vLogDelta_S = numeric(iD);

  dOmega0 = vPar["omega0"]
  dOmega1 = vPar["omega1"]
  dOmega2 = vPar["omega2"]

  for(s in 1:iD) {
    vLogDelta_S[s] = dOmega0 + dOmega1 * cos(2.0 * pi * s/iD - dOmega2 * pi);
  }
  if(iH > 0) {
    vLogDelta_S[1:iH] = vPar[paste("w", 1:iH, sep = "")]
  }

  vDelta_S = exp(vLogDelta_S)

  dnLLK = try(-NBSoftplusSeasonalINGARCH_filter(vY, dM0,
                                                dAlpha0,
                                                dAlpha1,  dBeta1,
                                                vDelta_S, iD, dC,
                                                dP)$dLLK, silent = TRUE)

  # print(dnLLK)

  if(is(dnLLK, "try-error")) {
    dnLLK = 1e7
  }
  if(is.nan(dnLLK)){
    dnLLK = 1e7
  }

  return(dnLLK)
}

Estimate_NBsoftplusINGARCH <- function(vY, iD, dC = 1, vStarting = NULL, iH = 0) {

  if(is.null(vStarting)) {

  vPar = c("m0" = mean(vY) * 0.5,
           "alpha0" = mean(vY) * 0.5,
           "alpha1" = 0.4,
           "beta1" = 0.5,
           "p" = 0.7)

  vParS = StartingParSeasonal(vY, iD, iH = iH)

  vPar = c(vPar, vParS)

  } else {
    vPar = vStarting
  }

  if (iH == iD) {
    LB = c(1e-4, -10, -0.9999, -0.9999, 0.000001, rep(-10, iH))
    UB = c(10, 10, 0.9999, 0.9999, 0.9999, rep(10, iH))
  } else {
    LB = c(-10, -10, -0.9999, -0.9999, 0.000001, rep(-10, 3 + iH))
    UB = c(10, 10, 0.9999, 0.9999, 0.9999, rep(10, 3 + iH))
  }

  optimizer = optim(vPar, optimizerNBsoftplusINGARCH,
                    lower = LB, upper = UB, vY = vY, iD = iD, dC = dC,
                    iH = iH,
                    method = "L-BFGS-B")

  dLLK = optimizer$value
  vPar = optimizer$par

  dM0 = vPar["m0"]
  dAlpha0 = vPar["alpha0"]
  dAlpha1 = vPar["alpha1"]
  dBeta1 = vPar["beta1"]
  dP = vPar["p"]

  vLogDelta_S = numeric(iD);

  if(iH < iD) {
    dOmega0 = vPar["omega0"]
    dOmega1 = vPar["omega1"]
    dOmega2 = vPar["omega2"]

    for(s in 1:iD) {
      vLogDelta_S[s] = dOmega0 + dOmega1 * cos(2.0 * pi * s/iD - dOmega2 * pi);
    }
    if(iH > 0) {
      vLogDelta_S[1:iH] = vPar[paste("w", 1:iH, sep = "")]
    }
  } else {
    vLogDelta_S = vPar[paste("w", 1:iH, sep = "")]
  }

  vDelta_S = exp(vLogDelta_S)

  Filter = NBSoftplusSeasonalINGARCH_filter(vY, dM0,
                                                dAlpha0,
                                                dAlpha1,  dBeta1,
                                            vDelta_S, iD, dC,
                                                dP)

  iC = ICfun(Filter$dLLK, np = length(vPar), iT = length(vY))

  return(list(Filter = Filter, dLLK = dLLK, vPar = vPar, optimizer = optimizer, iC = iC, iH = iH))

}

## with same dummy spec as HMM-INAR
optimizerNBsoftplusINGARCH_dummy <- function(vPar, vY, vSeason, iD, dC) {

  dM0 = vPar["m0"]
  dAlpha0 = vPar["alpha0"]
  dAlpha1 = vPar["alpha1"]
  dBeta1 = vPar["beta1"]
  dP = vPar["p"]

  vBeta_s = vPar[paste("beta_s", 2:iD, sep = "")];
  vBeta_s = c("beta_s1" = 1, vBeta_s)

  dnLLK = try(-NBSoftplusSeasonalINGARCH_DummySeason_filter(vY, dM0,
                                                dAlpha0,
                                                dAlpha1,  dBeta1,
                                                vSeason, vBeta_s, dC,
                                                dP)$dLLK, silent = TRUE)

  # print(dnLLK)
  if(is(dnLLK, "try-error")) {
    dnLLK = 1e7
  }
  if(is.nan(dnLLK)){
    dnLLK = 1e7
  }

  return(dnLLK)
}

Estimate_NBsoftplusINGARCH_DummySeason <- function(vY, mD = NULL, dC = 1, vStarting = NULL) {

  iD = ncol(mD)
  iT = length(vY)

  vSeason = numeric(iT)
  for (t in 1:iT) {
    for (d in 1:iD) {
      if ( abs(mD[t, d] - 1.0) < 1e-5) {
        vSeason[t] = d;
      }
    }
  }

  if(is.null(vStarting)) {

    vPar = c("m0" = mean(vY) * 0.5,
             "alpha0" = mean(vY) * 0.5,
             "alpha1" = 0.4,
             "beta1" = 0.5,
             "p" = 0.7)

    vBeta_s = numeric(iD)
    for (d in 1:iD) {
      vBeta_s[d] = mean(vY[vSeason == d])
    }

    vBeta_s = vBeta_s/vBeta_s[1]
    vBeta_s = vBeta_s[-1]

    names(vBeta_s) = paste("beta_s", 2:iD, sep = "")

    vPar = c(vPar, vBeta_s)

  } else {
    vPar = vStarting
  }

  LB = c(-10, -10, -0.9999, -0.9999, 0.000001, rep(1e-4, iD-1))
  UB = c(10, 10, 0.9999, 0.9999, 0.9999, rep(1e4, iD - 1))

  vSeason = vSeason - 1 #C indexing
  optimizer = optim(vPar, optimizerNBsoftplusINGARCH_dummy,
                    lower = LB, upper = UB, vY = vY, vSeason = vSeason, iD = iD,
                    dC = dC,
                    method = "L-BFGS-B")

  dLLK = optimizer$value
  vPar = optimizer$par

  dM0 = vPar["m0"]
  dAlpha0 = vPar["alpha0"]
  dAlpha1 = vPar["alpha1"]
  dBeta1 = vPar["beta1"]
  dP = vPar["p"]

  vBeta_s = vPar[paste("beta_s", 2:iD, sep = "")];
  vBeta_s = c("beta_s1" = 1, vBeta_s)

  Filter = NBSoftplusSeasonalINGARCH_DummySeason_filter(vY, dM0,
                                                            dAlpha0,
                                                            dAlpha1,  dBeta1,
                                                            vSeason, vBeta_s, dC,
                                                            dP)

  iC = ICfun(Filter$dLLK, np = length(vPar), iT = length(vY))

  return(list(Filter = Filter, dLLK = dLLK, vPar = vPar, optimizer = optimizer, iC = iC))

}



Estimate_NBsoftplusINGARCH_DummySeason_DEoptim <- function(vY, mD = NULL, dC = 1,
                                                           control = DEoptim.control(itermax = 200, NP = 1000)) {

  iD = ncol(mD)
  iT = length(vY)

  vSeason = numeric(iT)
  for (t in 1:iT) {
    for (d in 1:iD) {
      if ( abs(mD[t, d] - 1.0) < 1e-5) {
        vSeason[t] = d;
      }
    }
  }

  vPar = c("m0" = mean(vY) * 0.5,
           "alpha0" = mean(vY) * 0.5,
           "alpha1" = 0.4,
           "beta1" = 0.5,
           "p" = 0.7)

  vBeta_s = numeric(iD)
  for (d in 1:iD) {
    vBeta_s[d] = mean(vY[vSeason == d])
  }

  vBeta_s = vBeta_s/vBeta_s[1]
  vBeta_s = vBeta_s[-1]

  names(vBeta_s) = paste("beta_s", 2:iD, sep = "")
  vPar = c(vPar, vBeta_s)
  vNames = names(vPar)

  LB = c(1e-4, -10, -0.9999, -0.9999, 0.000001, rep(1e-4, iD-1))
  UB = c(10, 10, 0.9999, 0.9999, 0.9999, rep(100, iD - 1))

  vSeason = vSeason - 1 #C indexing

  optimizer = DEoptim(fn = function(vPar, vY, vSeason, iD, dC, vNames) {

    names(vPar) = vNames
    optimizerNBsoftplusINGARCH_dummy(vPar, vY, vSeason, iD, dC)

  },
  lower = LB, upper = UB, control = control,
  vY = vY, vSeason = vSeason, iD = iD,
  dC = dC, vNames = vNames)

  dLLK = -optimizer$optim$bestval
  vPar = optimizer$optim$bestmem
  names(vPar) = vNames

  dM0 = vPar["m0"]
  dAlpha0 = vPar["alpha0"]
  dAlpha1 = vPar["alpha1"]
  dBeta1 = vPar["beta1"]
  dP = vPar["p"]

  vBeta_s = vPar[paste("beta_s", 2:iD, sep = "")];
  vBeta_s = c("beta_s1" = 1, vBeta_s)

  Filter = NBSoftplusSeasonalINGARCH_DummySeason_filter(vY, dM0,
                                                        dAlpha0,
                                                        dAlpha1,  dBeta1,
                                                        vSeason, vBeta_s, dC,
                                                        dP)

  iC = ICfun(Filter$dLLK, np = length(vPar), iT = length(vY))

  return(list(Filter = Filter, dLLK = dLLK, vPar = vPar, optimizer = optimizer, iC = iC))

}

