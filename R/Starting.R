Starting_MixPois <- function(vY, iJ) {

  # Starting Values
  dMean = mean(vY);
  dVar  = var(vY);

  if (iJ > 1) {
    vWeight = seq(0.9, 1.1, length.out = iJ)

  } else {
    vWeight = 1.0;
  }

  if (dVar > dMean) {
    dMu = dMean^2.0/(dVar - dMean);
    vNu = vWeight * rep((dVar/dMean - 1.0), iJ);
  } else {
    dMu = dMean;
    vNu = vWeight * dVar;
  }

  vOmega = vWeight/sum(vWeight);

  names(vOmega) = paste("omega", 1:iJ, sep = "")
  names(vNu) = paste("nu", 1:iJ, sep = "")
  names(dMu) = "mu"

  vPar = c(dMu, vNu, vOmega)

  optimizer = optim(vPar, fn = function(vPar, vY, iJ) {

    vOmega = vPar[paste("omega", 1:iJ, sep = "")]
    vOmega = vOmega/sum(vOmega)

    dMu = vPar["mu"]
    vNu = vPar[paste("nu", 1:iJ, sep = "")]

    -dLLK_SumMixPois(vY, vOmega, vNu, dMu)/length(vY)

  }, lower = c(rep(0.01, iJ + 1), rep(0.0001, iJ)), upper = c(rep(100000, iJ + 1), rep(0.9999, iJ)),
  vY = vY, iJ = iJ, method = "L-BFGS-B",  control = list(maxit = 2))

  vPar = optimizer$par

  vOmega = vPar[paste("omega", 1:iJ, sep = "")]
  vOmega = vOmega/sum(vOmega)

  dMu = vPar["mu"]
  vNu = vPar[paste("nu", 1:iJ, sep = "")]

  return(list(dMu = dMu, vNu = vNu, vOmega = vOmega))

}

Starting_MSSumMixPois <- function(vY, iJ, iL) {

  # Starting Values
  lPar_MixPois = Starting_MixPois(vY, iJ)

  if (iL > 1) {
    vWeight = seq(0.9, 1.1, length.out = iL)
  } else {
    vWeight = 1.0;
  }

  vMu = rep(lPar_MixPois$dMu, iL) * vWeight
  vNu = lPar_MixPois$vNu
  vOmega = lPar_MixPois$vOmega

  if (iL > 1) {
    mGamma = diag(iL)
    mGamma[as.logical(diag(iL))] = 0.99
    mGamma[!diag(iL)] = (1.0 - 0.99)/(iL - 1)
  } else {
    mGamma = matrix(1.0);
  }

  return(list(vMu = vMu, vNu = vNu, vOmega = vOmega, mGamma = mGamma))

}

Starting_SumZTNB <- function(vV, vN) {

  dB = 0.8
  dA = (1.0 - dB)/dB * median(vV)/median(vN)

  vPw = pn2pw_SumZTNB(c("a" = dA, "b" = dB))

  return(vPw)

}

Starting_PoissonSumMixThreeBinom <- function(vN) {

  vPi = c(0.2, 0.4, 0.9)
  vOmega = c(0.1, 0.4, 0.5)
  dLambda = mean(vN) * 4

  names(vPi) = paste("pi", 1:3, sep = "")
  names(vOmega) = paste("omega", 1:3, sep = "")

  vPn = c(lambda = dLambda, vPi, vOmega)

  vPw = pn2pw_PoissonSumMixThreeBinom(vPn)

  return(vPw)

}

Starting_PoissonSumMixTwoBinom <- function(vN) {

  vPi = c(0.2,  0.9)
  dOmega = 0.4
  dLambda = max(vN)

  names(vPi) = paste("pi", 1:2, sep = "")

  vPn = c(lambda = dLambda, vPi, omega = dOmega)

  vPw = pn2pw_PoissonSumMixTwoBinom(vPn)

  return(vPw)

}

Starting_PoissonSumMixTwoBinom_Seasonal <- function(vN) {

  vPi = c(0.2,  0.9)
  dOmega = 0.4
  vEta = c(eta1 = log(max(vN)), eta2 = 0.5, eta3 = 0.2)

  names(vPi) = paste("pi", 1:2, sep = "")

  vPn = c(vPi, omega = dOmega, vEta)

  vPw = pn2pw_PoissonSumMixTwoBinom_Seasonal(vPn)

  return(vPw)

}

Starting_MixPoisInar <- function(vN, iJ) {

  dAlpha = cor(vN[-1], vN[-length(vN)])
  if(dAlpha < 0) {
    dAlpha = 0.1
  }

  MixPois_Fit = EM_Mixture_Pois(vN, iJ)

  vLambda = c(MixPois_Fit$vLambda) * (1 - dAlpha)
  vOmega  = c(MixPois_Fit$vP)

  names(vOmega)  = paste("omega", 1:iJ, sep = "")
  names(vLambda) = paste("lambda", 1:iJ, sep = "")

  vPn = c("alpha" = dAlpha, vLambda, vOmega)
  vPw = pn2pw_MixPoissonINAR(vPn, iJ, iD = 1)

  return(vPw)

}

Starting_MixPoisInar2 <- function(vN, iJ) {

  dAlpha = cor(vN[-1], vN[-length(vN)])
  if(dAlpha < 0) {
    dAlpha = 0.1
  }

  MixPois_Fit = EM_Mixture_Pois(vN, iJ, maxIter = 50, tol = 1e-4)

  vLambda = c(MixPois_Fit$vLambda) * (1 - dAlpha)
  vOmega  = c(MixPois_Fit$vP)

  lPn = list(vLambda = vLambda,
             vOmega = vOmega,
             dAlpha = dAlpha)

  return(lPn)

}

Starting_MixPoisInarSeasonal <- function(vN, iJ) {

  vPw = Starting_MixPoisInar(vN, iJ)

  vPw["eta1"] = 0.0
  vPw["eta2"] = logit(0.5)

  return(vPw)
}

Starting_MSMixPoisInar_TwoChains <- function(vY, iJ, iK, iL, mD, vSeason,
                                             bSeasonal, vLogFactorial_K,
                                             bParSeason = FALSE, iH = 0) {

  iT = length(vY)
  dMean = mean(vY)

  ## Alpha
  if (iJ > 1) {
    iC = round(iT * 0.3)
    vCor = numeric(iT - iC)
    for (t in 1:(iT - iC)) {
      vCor[t] =  cor(vY[t:(iC + t)][-1], vY[t:(iC + t)][-iC])
    }
    vAlpha = sort(vCor)[round(seq(1, iT - iC, length.out = iJ))]

    mGamma_Alpha = diag(rep(0.99, iJ))
    mGamma_Alpha[!diag(iJ)] = 0.01/(iJ - 1)

  } else {
    vAlpha =  cor(vY[-1], vY[-iT])
    mGamma_Alpha = matrix(1)
  }

  if(any(vAlpha < 0)) {
    vAlpha[vAlpha < 0] = 0.1
  }

  ## lambda and omega
  if (iK > 1) {
    if (iL > 1) {
      vT = round(seq(1, iT, length.out = iL + 1))
      lFit_MixPois = list()
      for (l in 1:iL) {
        lFit_MixPois[[l]] = EM_Mixture_Pois(vY[vT[l]:vT[l + 1]], iK, maxIter = 50, tol = 1e-4)
      }
      mLambda = do.call(cbind, lapply(lFit_MixPois, function(x) x$vLambda))
      mOmega  = do.call(cbind, lapply(lFit_MixPois, function(x) x$vP))

      for (l in 1:iL) {
        for (k in 1:iK) {
          if(mOmega[k, l] > 0.8) {
            mOmega[k, l] = 0.8
          }
          if(mOmega[k, l] < 0.2) {
            mOmega[k, l] = 0.2
          }
        }
      }
      for (l in 1:iL) {
        mOmega[, l] = abs(mOmega[, l])/sum(abs(mOmega[, l]))
      }

      vLambda = rowMeans(mLambda)

      vSeq = seq(0.5, 1.5, length.out = iK)
      for (k in 1:iK) {
        if (min(abs(vLambda[k] - vLambda[-k])/vLambda[k]) < 0.1) {
          vLambda[k] = dMean * vSeq[k]
        }
      }

      mGamma_Eta = diag(rep(0.99, iL))
      mGamma_Eta[!diag(iL)] = 0.01/(iL - 1)
    } else {

      Fit_MixPois = EM_Mixture_Pois(vY, iK, maxIter = 50, tol = 1e-4)
      vLambda = Fit_MixPois$vLambda
      mOmega = matrix(Fit_MixPois$vP)
      mGamma_Eta = matrix(1)

    }
  } else {

    vLambda = dMean
    mOmega = matrix(1)
    mGamma_Eta = matrix(1)
  }

  lPn = list("vLambda" = vLambda * (1.0 - vAlpha[1]),
             "vAlpha" = vAlpha,
             "dVarPhi" = min(vAlpha)/2.0,
             "mOmega" = mOmega,
             "mGamma_Eta" = mGamma_Eta,
             "mGamma_Alpha" = mGamma_Alpha)

  ## seasonal
  if (bSeasonal) {
    if(bParSeason) {
      vParS = StartingParSeasonal(vY, iD, iH = iH)
      lPn[["vParS"]] = vParS
    } else {
      iD = ncol(mD)
      vBeta = numeric(iD)
      for (d in 1:iD) {
        vBeta[d] = mean(vY[mD[, d] == 1])/dMean
      }

      vBeta_start = vBeta

      Fit_Start = InitializeBeta_INAR(vY, vBeta_start, vSeason, vLogFactorial_K, iMaxiter = 20, dTol = 1e-4)

      vBeta = Fit_Start$vBeta

      if (any(is.nan(vBeta))) {
        vBeta[is.nan(vBeta)] = vBeta_start[is.nan(vBeta)]
      }
      lPn[["vBeta"]] = vBeta
    }
  } else {
    vBeta = 1
    lPn[["vBeta"]] = vBeta
  }

  return(lPn)

}

StartingParSeasonal <- function(vY, iD, vD = NULL, iH = 0) {

  if(is.null(vD)) {
    iT = length(vY)
    vD = numeric(iD)
    for (s in 1:iD) {
      vD[s] = mean(vY[seq(s, iT, iD)])
    }
  } else {
    iD = length(vD)
  }

  if(iH < iD) {
  vParS = c("omega0" = log(as.numeric(mean(vD))), "omega1" = 0, "omega2" = 0)

  minimizer = optim(vParS, function(vParS, vD, iD) {

    dOmega0 = vParS["omega0"]
    dOmega1 = vParS["omega1"]
    dOmega2 = vParS["omega2"]

    vLogDelta_S = numeric(iD);
    for(s in 1:iD) {
      vLogDelta_S[s] = dOmega0 + dOmega1 * cos(2.0 * pi * s/iD - dOmega2 * pi);
    }

    sum((exp(vLogDelta_S) - vD)^2)

  }, vD = vD, iD = iD, method = "BFGS")

  vParS = minimizer$par

  }
  if(iH > 0) {
    vW = log(vD[1:iH])
    names(vW) = paste("w", 1:iH, sep = "")
    if (iH == iD) {
      vParS = vW
    } else {
      vParS = c(vW, vParS)
    }
  }

  return(vParS)

}
