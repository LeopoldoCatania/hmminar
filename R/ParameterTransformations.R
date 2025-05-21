
logitInv <- function(dX) {
  dExp = exp(dX)
  if(any(dExp > 1e4)) {
    dExp[dExp > 1e4] = 1e4
  }
  if(any(dExp < 1e-4)) {
    dExp[dExp < 1e-4] = 1e-4
  }
  foo = dExp/(1 + dExp)
  return(foo)
}

logit <- function(dX) {
  log(dX) - log(1.0 - dX)
}

pw2pn_SumZTNB <- function(vPw) {

  dA = exp(vPw["a"])
  dB = logitInv(vPw["b"])

  if (dA < 1e-3) {
    dA = 1e-3
  }

  vPn = c(dA, dB)

  return(vPn)

}
pn2pw_SumZTNB <- function(vPn) {

  dA_tilde = log(vPn["a"])
  dB_tilde = logit(vPn["b"])

  vPw = c(dA_tilde, dB_tilde)

  return(vPw)

}

pn2pw_PoissonSumMixThreeBinom <- function(vPn) {

  vOmega = vPn[paste("omega", 1:3, sep = "")]
  dLambda = vPn["lambda"]
  vPi = vPn[paste("pi", 1:3, sep = "")]

  dLambda_tilde = log(dLambda)

  vPi_tilde = logit(vPi)

  vOmega_tilde = SimplexUnMapping(vOmega, 3)
  names(vOmega_tilde) = paste("omega", 1:2, sep = "")

  vPw = c(dLambda_tilde, vPi_tilde, vOmega_tilde)
  return(vPw)

}
pw2pn_PoissonSumMixThreeBinom <- function(vPw) {

  vOmega_tilde = vPw[paste("omega", 1:2, sep = "")]
  dLambda_tilde = vPw["lambda"]
  vPi_tilde = vPw[paste("pi", 1:3, sep = "")]

  dLambda = exp(dLambda_tilde)

  vPi = logitInv(vPi_tilde)

  vOmega = SimplexMapping(vOmega_tilde, 3)
  names(vOmega) = paste("omega", 1:3, sep = "")

  vPn = c(dLambda, vPi, vOmega)
  return(vPn)

}

pn2pw_PoissonSumMixTwoBinom <- function(vPn) {

  dOmega = vPn["omega"]
  dLambda = vPn["lambda"]
  vPi = vPn[paste("pi", 1:2, sep = "")]

  dLambda_tilde = log(dLambda)

  vPi_tilde = logit(vPi)

  dOmega_tilde = logit(dOmega)

  vPw = c(dLambda_tilde, vPi_tilde, dOmega_tilde)
  return(vPw)

}
pw2pn_PoissonSumMixTwoBinom <- function(vPw) {

  dOmega_tilde = vPw["omega"]
  dLambda_tilde = vPw["lambda"]
  vPi_tilde = vPw[paste("pi", 1:2, sep = "")]

  dLambda = exp(dLambda_tilde)

  vPi = logitInv(vPi_tilde)
  dOmega = logitInv(dOmega_tilde)

  vPn = c(dLambda, vPi, dOmega)
  return(vPn)

}

pn2pw_PoissonSumMixTwoBinom_Seasonal <- function(vPn) {

  dOmega = vPn["omega"]
  vPi = vPn[paste("pi", 1:2, sep = "")]
  vEta = vPn[paste("eta", 1:3, sep = "")]

  vEta_tilde = vEta

  vEta_tilde[1] = vEta[1]
  vEta_tilde[2] = vEta[2]
  vEta_tilde[3] = logit(vEta[3])

  vPi_tilde = logit(vPi)

  dOmega_tilde = logit(dOmega)

  vPw = c(vPi_tilde, dOmega_tilde, vEta_tilde)
  return(vPw)

}
pw2pn_PoissonSumMixTwoBinom_Seasonal <- function(vPw) {

  dOmega_tilde = vPw["omega"]
  vEta_tilde = vPw[paste("eta", 1:3, sep = "")]
  vPi_tilde = vPw[paste("pi", 1:2, sep = "")]

  vEta = vEta_tilde

  vEta[1] = vEta_tilde[1]
  vEta[2] = vEta_tilde[2]
  vEta[3] = logitInv(vEta_tilde[3])

  vPi = logitInv(vPi_tilde)
  dOmega = logitInv(dOmega_tilde)

  vPn = c(vPi, dOmega, vEta)
  return(vPn)

}
pw2pn_MixPoissonINAR <- function(vPw, iJ, iD = 1, Vector = FALSE) {

  vOmega_tilde  = vPw[paste("omega", 1:(iJ - 1), sep = "")]
  dAlpha_tilde  = vPw["alpha"]
  vLambda_tilde = vPw[paste("lambda", 1:iJ, sep = "")]
  vBeta_tilde = vPw[paste("beta", 1:iD, sep = "")]

  vLambda = exp(vLambda_tilde)
  vBeta   = exp(vBeta_tilde)
  dAlpha  = logitInv(dAlpha_tilde)

  if (iJ > 1) {
    vOmega  = SimplexMapping(vOmega_tilde, iJ)
  } else {
    vOmega = 1
  }
  names(vOmega) = paste("omega", 1:iJ, sep = "")

  if (Vector) {
    return(c(vOmega, dAlpha, vLambda, vBeta))
  } else {
    return(list(vOmega = vOmega,
                dAlpha = dAlpha,
                vLambda = vLambda,
                vBeta = vBeta))
  }

}

pn2pw_MixPoissonINAR <- function(vPn, iJ, iD) {

  vOmega  = vPn[paste("omega", 1:(iJ - 1), sep = "")]
  dAlpha  = vPn["alpha"]
  vLambda = vPn[paste("lambda", 1:iJ, sep = "")]
  vBeta = vPn[paste("beta", 1:iD, sep = "")]

  vLambda_tilde = log(vLambda)
  vBeta_tilde = log(vBeta)
  dAlpha_tilde  = logit(dAlpha)

  if(iJ > 1) {
    vOmega_tilde  = SimplexUnMapping(vOmega, iJ)
    names(vOmega_tilde) = paste("omega", 1:(iJ - 1), sep = "")
  } else {
    vOmega_tilde = NULL
  }

  return(c(vOmega_tilde, dAlpha_tilde, vLambda_tilde, vBeta_tilde))

}


UnmapGamma <- function(mGamma, name = NULL) {

  iF = ncol(mGamma)

  vPw = NULL
  for (i in 1:iF) {
    vPw_i = c(SimplexUnMapping(mGamma[i, ], iF))
    if(is.null(name)) {
      names(vPw_i) = paste("gamma_tilde", name, i, 1:(iF-1), sep = "_")
    } else {
      names(vPw_i) = paste("gamma_tilde", name, i, 1:(iF-1), sep = "_")
    }
    vPw = c(vPw, vPw_i)
  }

  return(vPw)
}

MapGamma <- function(vPw, iF, name = NULL) {

  mGamma = matrix(0, iF, iF)

  for (i in 1:iF) {
    if(is.null(name)) {
      vPw_i = vPw[paste("gamma_tilde", name, i, 1:(iF-1), sep = "_")]
    } else {
      vPw_i = vPw[paste("gamma_tilde", name, i, 1:(iF-1), sep = "_")]
    }
    mGamma[i, ] = SimplexMapping(vPw_i, iF)
  }
  return(mGamma)
}

UnmapOmega <- function(mOmega) {

  iK = nrow(mOmega)
  iL = ncol(mOmega)
  vPw = NULL
  for (i in 1:iL) {
    vPw_i = c(SimplexUnMapping(mOmega[,i], iK))
    names(vPw_i) = paste("omega_tilde", i, 1:(iK-1), sep = "_")
    vPw = c(vPw, vPw_i)
  }
  return(vPw)
}

MapOmega <- function(vPw, iK, iL) {
  mOmega = matrix(0, iK, iL)
  for (i in 1:iL) {
    vPw_i = vPw[paste("omega_tilde", i, 1:(iK-1), sep = "_")]
    mOmega[, i] = SimplexMapping(vPw_i, iK)
  }
  return(mOmega)
}

getlPar <- function(Fit) {

  mGamma_Eta = Fit$mGamma_Eta
  mGamma_Alpha = Fit$mGamma_Alpha
  mOmega = Fit$mOmega
  vAlpha = Fit$vAlpha
  vBeta = Fit$vBeta
  vLambda = Fit$vLambda

  if(!is.null(Fit$dVarPhi)) {

    dVarPhi = Fit$dVarPhi
    lPar = list(mGamma_Eta = mGamma_Eta,
                mGamma_Alpha = mGamma_Alpha,
                mOmega = mOmega,
                vAlpha = vAlpha,
                dVarPhi = dVarPhi,
                vBeta = vBeta,
                vLambda = vLambda)

  } else {

    lPar = list(mGamma_Eta = mGamma_Eta,
                mGamma_Alpha = mGamma_Alpha,
                mOmega = mOmega,
                vAlpha = vAlpha,
                vBeta = vBeta,
                vLambda = vLambda)

  }

  return(lPar)
}

###############
pn2pw_MSMixPoisInar_TwoChains_EM <- function(lPar, iJ, iL, iK, iD, bSeasonal = FALSE, bAlphacor = FALSE, bParSeason = FALSE) {

  mGamma_Eta = lPar$mGamma_Eta
  mGamma_Alpha = lPar$mGamma_Alpha
  mOmega = lPar$mOmega
  vAlpha = lPar$vAlpha

  if(bAlphacor) {
    dVarPhi = lPar$dVarPhi
  }

  if (bSeasonal) {
    if(bParSeason) {
      vParS = lPar$vParS
      vBeta = NULL
    } else {
      vBeta = lPar$vBeta
      vParS = NULL
    }
  }

  vLambda = lPar$vLambda

  if (iL > 1) {
    vPw_gamma_eta = UnmapGamma(mGamma_Eta, "eta")
  } else {
    vPw_gamma_eta = NULL
  }
  if (iJ > 1) {
    vPw_gamma_alpha = UnmapGamma(mGamma_Alpha, "alpha")
  } else {
    vPw_gamma_alpha = NULL
  }

  if(iK > 1) {
    vPw_omega = UnmapOmega(mOmega)
  } else {
    vPw_omega = NULL
  }

  vPw_alpha = as.numeric(sapply(vAlpha, logit))
  names(vPw_alpha) = paste("alpha_tilde", 1:iJ, sep = "_")

  if(bAlphacor) {
    vPw_varphi = as.numeric(logit(dVarPhi))
    names(vPw_varphi) = "varphi_tilde"
  } else {
    vPw_varphi = NULL
  }

  if (bSeasonal) {
    if(bParSeason) {
      vPw_pars = vParS
      vPw_beta = NULL
    } else {
      vPw_beta = as.numeric(log(vBeta))
      names(vPw_beta) = paste("beta", 1:iD, sep = "_")
      vPw_pars = NULL
    }
  } else {
    vPw_beta = NULL
    vPw_pars = NULL
  }

  vPw_lambda = as.numeric(log(vLambda))
  names(vPw_lambda) = paste("lambda", 1:iK, sep = "_")

  vPw = c(vPw_gamma_eta, vPw_gamma_alpha, vPw_omega, vPw_alpha, vPw_varphi, vPw_beta, vPw_pars, vPw_lambda)

  return(vPw)

}
pw2pn_MSMixPoisInar_TwoChains_EM <- function(vPw, iJ, iL, iK, iD, bSeasonal = FALSE, bAlphacor = FALSE, vecout = FALSE, bParSeason = FALSE, iH = 0) {

  lPar = list()
  vPn = NULL

  if(iL > 1) {
    mGamma_Eta = GammaCheck(MapGamma(vPw, iL, "eta"), iL)
    vPn_gamma_eta = c(t(mGamma_Eta[, 1:(iL - 1)]))
    names(vPn_gamma_eta) = c(sapply(1:iL, function(i) paste("gamma_eta", i, 1:(iL-1), sep = "_")))
  } else {
    mGamma_Eta = matrix(1)
    vPn_gamma_eta = NULL
  }
  if (iJ > 1) {
    mGamma_Alpha = GammaCheck(MapGamma(vPw, iJ, "alpha"), iJ)
    vPn_gamma_alpha = c(t(mGamma_Alpha[, 1:(iJ - 1)]))
    names(vPn_gamma_alpha) = c(sapply(1:iJ, function(i) paste("gamma_alpha", i, 1:(iJ-1), sep = "_")))
  } else {
    mGamma_Alpha = matrix(1)
    vPn_gamma_alpha = NULL
  }
  if (iK > 1) {
    mOmega = OmegaCheck(MapOmega(vPw, iK, iL), iK, iL)
    vPn_omega = c(mOmega[1:(iK-1),])
    names(vPn_omega) = c(sapply(1:iL, function(i) paste("omega", i, 1:(iK-1), sep = "_")))
  } else {
    mOmega = matrix(1)
    vPn_omega = NULL
  }

  vAlpha = as.matrix(as.numeric(logitInv(vPw[paste("alpha_tilde", 1:iJ, sep = "_")])))

  vLambda = as.matrix(as.numeric(exp(vPw[paste("lambda", 1:iK, sep = "_")])))
  if(any(vLambda < 1e-10)) {
    vLambda[vLambda<1e-10] = 1e-10
  }
  if(any(vLambda > 1e10)) {
    vLambda[vLambda > 1e10] = 1e10
  }

  lPar = list(mGamma_Eta = mGamma_Eta,
              mGamma_Alpha = mGamma_Alpha,
              mOmega = mOmega,
              vAlpha = vAlpha,
              vLambda = vLambda)

  if (bAlphacor) {
    dVarPhi = as.numeric(logitInv(vPw["varphi_tilde"]))
    lPar[["dVarPhi"]] = dVarPhi
  } else {
    dVarPhi = NULL
  }

  if (bSeasonal) {

    if(bParSeason) {
      if (iH == 0) {
        vParS = vPw[c("omega0", "omega1", "omega2")]
      } else {
        vParS = vPw[c(paste("w", 1:iH, sep = ""), "omega0", "omega1", "omega2")]
      }

      vBeta = NULL
      lPar[["vParS"]] = vParS
    } else {
      vParS = NULL
      vBeta = as.matrix(as.numeric(exp(vPw[paste("beta", 1:iD, sep = "_")])))
      lPar[["vBeta"]] = vBeta
    }
  } else {
    vBeta = NULL
    vParS = NULL
  }

  if(vecout) {
    vPn_alpha = c(vAlpha)
    names(vPn_alpha) = paste("alpha", 1:iJ, sep = "_")

    if (bAlphacor) {
      vPn_varphi = dVarPhi
      names(vPn_varphi) = "varphi"
    } else {
      vPn_varphi = NULL
    }

    if (bSeasonal) {
      if(bParSeason) {
        vPn_pars = vParS
        vPn_beta = NULL
      } else {
        vPn_beta = c(vBeta)
        names(vPn_beta) = paste("beta", 1:iD, sep = "_")
        vPn_pars = NULL
      }
    } else {
      vPn_beta = NULL
      vPn_pars = NULL
    }

    vPn_lambda = c(vLambda)
    names(vPn_lambda) = paste("lambda", 1:iK, sep = "_")

    vPn = c(vPn_gamma_eta, vPn_gamma_alpha, vPn_omega, vPn_alpha, vPn_varphi, vPn_beta, vPn_pars, vPn_lambda)

    return(vPn)
  } else {
    return(lPar)
  }
}

##############
pn2pw_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor <- function(lPar, iJ, iL, iK, iD, bParSeason = FALSE) {

  mGamma_Eta = lPar$mGamma_Eta
  mGamma_Alpha = lPar$mGamma_Alpha
  mOmega = lPar$mOmega
  vAlpha = lPar$vAlpha
  dVarPhi = lPar$dVarPhi

  if(bParSeason) {
    vParS = lPar$vParS
    vBeta = NULL
  } else {
    vBeta = lPar$vBeta
    vParS = NULL
  }

  vLambda = lPar$vLambda

  if (iL > 1) {
    vPw_gamma_eta = UnmapGamma(mGamma_Eta, "eta")
  } else {
    vPw_gamma_eta = NULL
  }
  if (iJ > 1) {
    vPw_gamma_alpha = UnmapGamma(mGamma_Alpha, "alpha")
  } else {
    vPw_gamma_alpha = NULL
  }

  if(iK > 1) {
    vPw_omega = UnmapOmega(mOmega)
  } else {
    vPw_omega = NULL
  }

  vPw_alpha = as.numeric(sapply(vAlpha, logit))
  names(vPw_alpha) = paste("alpha_tilde", 1:iJ, sep = "_")

  vPw_varphi = as.numeric(logit(dVarPhi))
  names(vPw_varphi) = "varphi_tilde"

  if(bParSeason) {
    vPw_pars = vParS
    vPw_beta = NULL
  } else {
    vPw_beta = as.numeric(log(vBeta))
    names(vPw_beta) = paste("beta", 1:iD, sep = "_")
    vPw_pars = NULL
  }

  vPw_lambda = as.numeric(log(vLambda))
  names(vPw_lambda) = paste("lambda", 1:iK, sep = "_")

  vPw = c(vPw_gamma_eta, vPw_gamma_alpha, vPw_omega, vPw_alpha, vPw_varphi, vPw_beta, vPw_pars, vPw_lambda)

  return(vPw)

}

pw2pn_MSMixPoisInar_TwoChains_Seasonal_EM_AlfaCor <- function(vPw, iJ, iL, iK, iD, vecout = FALSE, bParSeason = FALSE, iH = 0) {

  lPar = list()
  vPn = NULL

  if(iL > 1) {
    mGamma_Eta = GammaCheck(MapGamma(vPw, iL, "eta"), iL)
    vPn_gamma_eta = c(t(mGamma_Eta[, 1:(iL - 1)]))
    names(vPn_gamma_eta) = c(sapply(1:iL, function(i) paste("gamma_eta", i, 1:(iL-1), sep = "_")))
  } else {
    mGamma_Eta = matrix(1)
    vPn_gamma_eta = NULL
  }
  if (iJ > 1) {
    mGamma_Alpha = GammaCheck(MapGamma(vPw, iJ, "alpha"), iJ)
    vPn_gamma_alpha = c(t(mGamma_Alpha[, 1:(iJ - 1)]))
    names(vPn_gamma_alpha) = c(sapply(1:iJ, function(i) paste("gamma_alpha", i, 1:(iJ-1), sep = "_")))
  } else {
    mGamma_Alpha = matrix(1)
    vPn_gamma_alpha = NULL
  }
  if (iK > 1) {
    mOmega = OmegaCheck(MapOmega(vPw, iK, iL), iK, iL)
    vPn_omega = c(mOmega[1:(iK-1),])
    names(vPn_omega) = c(sapply(1:iL, function(i) paste("omega", i, 1:(iK-1), sep = "_")))
  } else {
    mOmega = matrix(1)
    vPn_omega = NULL
  }

  vAlpha = as.matrix(as.numeric(logitInv(vPw[paste("alpha_tilde", 1:iJ, sep = "_")])))
  dVarPhi = as.numeric(logitInv(vPw["varphi_tilde"]))

  vLambda = as.matrix(as.numeric(exp(vPw[paste("lambda", 1:iK, sep = "_")])))
  if(any(vLambda < 1e-10)) {
    vLambda[vLambda<1e-10] = 1e-10
  }
  if(any(vLambda > 1e10)) {
    vLambda[vLambda > 1e10] = 1e10
  }

  lPar = list(mGamma_Eta = mGamma_Eta,
              mGamma_Alpha = mGamma_Alpha,
              mOmega = mOmega,
              vAlpha = vAlpha,
              dVarPhi = dVarPhi,
              vLambda = vLambda)

  if(bParSeason) {
    if (iH == 0) {
      vParS = vPw[c("omega0", "omega1", "omega2")]
    } else {
      vParS = vPw[c(paste("w", 1:iH, sep = ""), "omega0", "omega1", "omega2")]
    }

    vBeta = NULL
    lPar[["vParS"]] = vParS
  } else {
    vParS = NULL
    vBeta = as.matrix(as.numeric(exp(vPw[paste("beta", 1:iD, sep = "_")])))
    lPar[["vBeta"]] = vBeta
  }

  if(vecout) {
    vPn_alpha = c(vAlpha)
    names(vPn_alpha) = paste("alpha", 1:iJ, sep = "_")

    vPn_varphi = dVarPhi
    names(vPn_varphi) = "varphi"

    if(bParSeason) {
      vPn_pars = vParS
      vPn_beta = NULL
    } else {
      vPn_beta = c(vBeta)
      names(vPn_beta) = paste("beta", 1:iD, sep = "_")
      vPn_pars = NULL
    }

    vPn_lambda = c(vLambda)
    names(vPn_lambda) = paste("lambda", 1:iK, sep = "_")

    vPn = c(vPn_gamma_eta, vPn_gamma_alpha, vPn_omega, vPn_alpha, vPn_varphi, vPn_beta, vPn_pars, vPn_lambda)

    return(vPn)
  } else {
    return(lPar)
  }
}

