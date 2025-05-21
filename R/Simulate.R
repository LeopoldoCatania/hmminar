Simulate <- function(Fit, iT = NULL) {

  vLambda = Fit$vLambda
  mOmega  = Fit$mOmega
  vAlpha = Fit$vAlpha
  vBeta = Fit$vBeta
  dVarPhi = Fit$dVarPhi
  mGamma_Alpha = Fit$mGamma_Alpha
  mGamma_Eta = Fit$mGamma_Eta
  mD = Fit$mD
  vO = c(Fit$dO0, Fit$vO)

  vSeason = c(Fit$iSeason0, Fit$vSeason)

  iJ = Fit$iJ
  iL = Fit$iL
  iK = Fit$iK
  iD = Fit$iD

  if (is.null(iT)) {
    iT = Fit$iT
  } else {

    if (iT > Fit$iT) {

      iTimes = round(iT/Fit$iT + 1)
      vO = rep(vO, iTimes)
      vSeason = rep(vSeason, iTimes)

      vO = vO[1:iT]
      vSeason = vSeason[1:iT]
    }

  }

  Sim = Sim_MSMixInarTwoChains_AlphaCor(iT, mOmega, vLambda, vAlpha,
                                        mGamma_Alpha, mGamma_Eta, vBeta,
                                        dVarPhi, vO, vSeason)

  Sim[["vSeason"]] = vSeason
  Sim[["vO"]] = vO
  Sim[["mD"]] = mD
  Sim[["vLambda"]] = vLambda
  Sim[["mOmega"]] = mOmega
  Sim[["vAlpha"]] = vAlpha
  Sim[["vBeta"]] = vBeta
  Sim[["dVarPhi"]] = dVarPhi
  Sim[["mGamma_Alpha"]] = mGamma_Alpha
  Sim[["mGamma_Eta"]] = mGamma_Eta
  Sim[["mD"]] = mD
  Sim[["iJ"]] = iJ
  Sim[["iL"]] = iL
  Sim[["iK"]] = iK
  Sim[["iD"]] = iD

  return(Sim)

}

fSim<- function(iT, mOmega, vAlpha, mGamma_Alpha, mGamma_Eta, vLambda,
                vO = NULL, vBeta = 1, vSeason = NULL, dVarPhi = 0) {

  iJ = length(vAlpha)
  iL = ncol(mOmega)
  iK = nrow(mOmega)

  mIndices = OmegaZB_Indices(iJ, iK, iL);

  mGamma = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);

  vDelta = getDelta(mGamma, iJ*iK*iL);
  vAlpha_bar = vAlpha[mIndices[, 1] + 1]
    # kronecker(rep(1, iK*iL), vAlpha)
  vLambda_bar = vLambda[mIndices[, 2] + 1]
    # kronecker(rep(1, iJ*iL), vLambda)

  if (is.null(vO)) vO = rep(0, iT)
  if (is.null(vSeason)) vSeason = rep(0, iT)

  lSim = Sim_MSMixInarOneChain_AlphaCor(iT, vLambda_bar, vAlpha_bar,
                                        mGamma, vBeta,
                                        dVarPhi, vO, vSeason)

  return(lSim)
}
