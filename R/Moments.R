Moments <- function(Fit) {

  iJ = Fit$iJ
  iL = Fit$iL
  mD = Fit$mD
  vMu = Fit$vMu
  vNu = Fit$vNu
  vBeta = Fit$vBeta
  vOmega = Fit$vOmega
  mGamma = Fit$mGamma
  PredictedProb = Fit$PredictedProb

  if(is.null(mD)) {
    mD = matrix(1, nrow = nrow(PredictedProb) - 1)
    vBeta = c(1)
  }

  lM = Moments_C(PredictedProb,
                 iJ, iL,
                 vMu, vNu,
                 vOmega,
                 mGamma,
                 vBeta,
                 mD)

  return(lM)

}
#
# mGamma_Eta = Fit$mGamma_Eta
# mGamma_Alpha = Fit$mGamma_Alpha
# mOmega = Fit$mOmega
# vLambda = Fit$vLambda
# vAlpha = Fit$vAlpha
# Fit$iJ
# Fit$iK
# Fit$iL

fCross_enlarged <- function(k, mA, mLambda, mDelta, mGamma) {

  iH = ncol(mLambda)
  mI = diag(iH)
  mGtilde = solve(mDelta) %*% t(mGamma) %*% mDelta
  vY_avg = solve(mI - mA %*% mGtilde) %*% diag(mLambda)
  vY2 = solve(mI - mA %*% mA %*% mGtilde) %*%( (mA %*% (mI - mA) + 2*mLambda %*% mA)%*% mGtilde %*%
                                                 solve(mI - mA%*%mGtilde)%*%diag(mLambda) + (mI + mLambda)%*%diag(mLambda))

  mFoo = matrix(0, iH, iH)
  if (k > 1) {
    for (i in 0:(k - 2)) {
      mFoo = mFoo + powmat(mA %*% t(mGamma), i) %*% mLambda %*% t(powmat(mGamma, k-i-1))
    }
  }

  dYtYtmk = t(diag(mA)) %*% t(mGamma) %*% mFoo %*% mDelta %*% vY_avg + t(diag(mA)) %*% t(mGamma) %*% powmat(mA %*% t(mGamma), k - 1) %*% mDelta %*% vY2 + t(diag(mLambda)) %*% t(powmat(mGamma, k)) %*% mDelta %*% vY_avg
  dAtAtmk = t(diag(mA)) %*% t(mGamma) %*% mFoo %*% mDelta %*% mA %*% mGtilde %*% vY_avg + t(diag(mA)) %*% t(mGamma) %*% powmat(mA %*% t(mGamma), k - 1) %*% mDelta %*% (vY2 - (as.numeric(mLambda%*%vY_avg) + diag(mLambda)))
  dEtatEtatmk = t(diag(mLambda)) %*% t(powmat(mGamma, k)) %*% mDelta %*% diag(mLambda)
  dAtEtatmk = t(diag(mA)) %*% t(mGamma) %*% mFoo %*% mDelta %*% diag(mLambda) + t(diag(mA)) %*% t(mGamma) %*% powmat(mA %*% t(mGamma), k - 1) %*% mDelta %*% ((as.numeric(mLambda%*%vY_avg) + diag(mLambda)))
  dEtatAtmk = t(diag(mLambda)) %*% t(powmat(mGamma, k)) %*% mDelta %*% mA %*% mGtilde %*% vY_avg

  return(
  c("dYtYtmk" = dYtYtmk,
    "dAtAtmk" = dAtAtmk,
    "dEtatEtatmk" = dEtatEtatmk,
    "dAtEtatmk" = dAtEtatmk,
    "dEtatAtmk" = dEtatAtmk)
  )

}
fCross <- function(k, mOmega, vAlpha, mGamma_Alpha, mGamma_Eta, vLambda) {

  iJ = length(vAlpha)
  iL = ncol(mOmega)
  iK = nrow(mOmega)

  mIndices = OmegaZB_Indices(iJ, iK, iL);

  mGamma = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);

  vDelta = getDelta(mGamma, iJ*iK*iL);

  vAlpha_bar =
    vAlpha[mIndices[, 1] + 1]
  # kronecker(rep(1, iK*iL), vAlpha)
  vLambda_bar =
    vLambda[mIndices[, 2] + 1]
  # kronecker(rep(1, iJ*iL), vLambda)

  mA = diag(as.numeric(vAlpha_bar), nrow = length(vAlpha_bar))
  mLambda = diag(as.numeric(vLambda_bar), nrow = length(vLambda_bar))
  mDelta = diag(as.numeric(vDelta), nrow = length(vDelta))

  fCross_enlarged(k, mA, mLambda, mDelta, mGamma)
}
powmat <- function(mX, j) {
  iH = ncol(mX)
  if (j == 0) {
    return(diag(iH))
  } else if (j == 1) {
    return(mX)
  } else {
    mX %*% powmat(mX, j - 1)
  }
}

Unconditional_Moments <- function(mGamma_Eta, mGamma_Alpha,
                                  mOmega, vLambda, vAlpha, lag = 200) {

  iL = ncol(mGamma_Eta)
  iJ = ncol(mGamma_Alpha)
  iK = nrow(mOmega)

  # //indices in the enlarged reparameterization
  # // First index for S_t^alpha
  # // Second index for Z_t
  # // Third index for S_t^eta
  mIndices = OmegaZB_Indices(iJ, iK, iL);

  mGamma = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);

  vDelta = getDelta(mGamma, iJ*iK*iL);

  vAlpha_bar =
    vAlpha[mIndices[, 1] + 1]
  # kronecker(rep(1, iK*iL), vAlpha)
  vLambda_bar =
    vLambda[mIndices[, 2] + 1]
  # kronecker(rep(1, iJ*iL), vLambda)

  mA = diag(as.numeric(vAlpha_bar), nrow = length(vAlpha_bar))
  mLambda = diag(as.numeric(vLambda_bar), nrow = length(vLambda_bar))
  mDelta = diag(as.numeric(vDelta), nrow = length(vDelta))
  iH = iJ*iK*iL
  mI = diag(iH)
  mGtilde = solve(mDelta) %*% t(mGamma) %*% mDelta
  vY_avg = solve(mI - mA %*% mGtilde) %*% diag(mLambda)

  vY2 = solve(mI - mA %*% mA %*% mGtilde) %*%( (mA %*% (mI - mA) + 2*mLambda %*% mA)%*% mGtilde %*%
                                                 solve(mI - mA%*%mGtilde)%*%diag(mLambda) + (mI + mLambda)%*%diag(mLambda))

  dEY = t(vDelta) %*% solve(mI - mA %*% solve(mDelta) %*% t(mGamma) %*% mDelta) %*% vLambda_bar
  dEY2 = t(vDelta) %*% solve(mI - mA %*% mA %*% solve(mDelta) %*% t(mGamma) %*% mDelta) %*% (
    (mI + mLambda) %*% vLambda_bar + (mA %*% (mI - mA) + 2 * mLambda %*% mA) %*% solve(mDelta) %*% t(mGamma) %*% mDelta %*%
      solve(mI - mA %*% solve(mDelta) %*% t(mGamma) %*% mDelta) %*% vLambda_bar
  )

  dEA = t(vDelta) %*% mA %*% solve(mDelta) %*% t(mGamma) %*% mDelta %*% vY_avg
  dEEta = t(vDelta) %*% vLambda_bar
  dEA_Eta = t(vDelta) %*% mLambda %*% mA %*% solve(mDelta) %*% t(mGamma) %*% mDelta %*% vY_avg

  dEA2 = t(vDelta) %*% mA %*% ((mI - mA) %*% solve(mDelta) %*% t(mGamma) %*% mDelta %*% vY_avg + mA %*% solve(mDelta) %*% t(mGamma) %*% mDelta %*% vY2)
  dEEta2 = t(vDelta) %*% (mI + mLambda) %*% vLambda_bar

  dVar_y = dEY2 - dEY^2
  dVar_a = dEA2 - dEA^2
  dVar_eta = dEEta2 - dEEta^2

  vCross_y_y = numeric(lag)
  vCross_a_a = numeric(lag)
  vCross_eta_eta = numeric(lag)
  vCross_a_eta = numeric(lag)
  vCross_eta_a = numeric(lag)

  for (i in 1:lag) {
    foo = fCross(i, mOmega, vAlpha, mGamma_Alpha, mGamma_Eta, vLambda)
    vCross_y_y[i] = foo["dYtYtmk"]
    vCross_a_a[i] = foo["dAtAtmk"]
    vCross_eta_eta[i] = foo["dEtatEtatmk"]
    vCross_a_eta[i] = foo["dAtEtatmk"]
    vCross_eta_a[i] = foo["dEtatAtmk"]
  }

  #VarianceDecomposition
  VarDec = c("varY" = dVar_y,
             "varA" = dVar_a,
             "varEta" = dVar_eta,
             "covA_Eta" = dEA_Eta - dEA*dEEta)

  # vCross_y_y = c(dEY2, vCross_y_y)
  # vCross_a_a = c(dEA2, vCross_a_a)
  # vCross_eta_eta = c(dEEta2, vCross_eta_eta)
  # vCross_a_eta = c(dEA_Eta, vCross_a_eta)
  # vCross_eta_a = c(dEA_Eta, vCross_eta_a)

  vAutoCov_y = vCross_y_y - as.numeric(dEY)^2
  vAutoCov_a = vCross_a_a - as.numeric(dEA)^2
  vAutoCov_eta = vCross_eta_eta - as.numeric(dEEta)^2

  vAutoCov_a_eta = vCross_a_eta - as.numeric(dEEta) * as.numeric(dEA)
  vAutoCov_eta_a = vCross_eta_a - as.numeric(dEEta) * as.numeric(dEA)

  vAutoCor_y = vAutoCov_y/as.numeric(dVar_y)
  vAutoCor_a = vAutoCov_a/as.numeric(dVar_a)
  vAutoCor_eta = vAutoCov_eta/as.numeric(dVar_eta)

  vAutoCor_eta_a = vAutoCov_eta_a/as.numeric(sqrt(dVar_eta*dVar_a))
  vAutoCor_a_eta = vAutoCov_a_eta/as.numeric(sqrt(dVar_eta*dVar_a))

  lOut = list(dMean_Y = dEY,
              dVar_Y = dVar_y,
              dEY2 = dEY2,
              dMean_A = dEA,
              dVar_A = dVar_a,
              dEA2 = dEA2,
              dMean_Eta = dEEta,
              dVar_Eta = dVar_eta,
              dEEta2 = dEEta2,
              vCross_Y = vCross_y_y,
              vAutoCov_Y = vAutoCov_y,
              vAutoCor_Y = vAutoCor_y,
              vCross_A = vCross_a_a,
              vAutoCov_A = vAutoCov_a,
              vAutoCor_A = vAutoCor_a,
              vCross_Eta = vCross_eta_eta,
              vAutoCov_Eta = vAutoCov_eta,
              vAutoCor_Eta = vAutoCor_eta,
              vCross_Eta_A = vCross_eta_a,
              vAutoCov_Eta_A = vAutoCov_eta_a,
              vAutoCor_Eta_A = vAutoCor_eta_a,
              vCross_A_Eta = vCross_a_eta,
              vAutoCov_A_Eta = vAutoCov_a_eta,
              vAutoCor_A_Eta = vAutoCor_a_eta,
              VarDec = VarDec
              )
  return(lOut)

}


Conditional_Moments <- function(Fit) {

  iJ = Fit$iJ
  iK = Fit$iK
  iL = Fit$iL
  iH = iJ * iK * iL
  iT = Fit$iT
  vY = c(dY0, Fit$vY)
  vO = Fit$vO
  vSeason = c(iSeason0, Fit$vSeason) + 1
  vBeta = Fit$vBeta
  vAlpha = Fit$vAlpha
  vLambda = Fit$vLambda
  dVarPhi = Fit$dVarPhi

  mPredProbs = Fit$PredictedProb

  mPredProbs_S_t_Alpha = matrix(0, iT, iJ)
  mPredProbs_S_t_Eta = matrix(0, iT, iL)
  mPredProbs_Z_t = matrix(0, iT, iK)

  vE_A_t = numeric(iT)
  vE_Eta_t = numeric(iT)
  vE_A_t_Eta_t = numeric(iT)

  vE_A2_t = numeric(iT)
  vE_Eta2_t = numeric(iT)

  # //indices in the enlarged reparameterization
  # // First index for S_t^alpha
  # // Second index for Z_t
  # // Third index for S_t^eta
  mIndices = OmegaZB_Indices(iJ, iK, iL) + 1;

  for (t in 2:iT) {

    for (h in 1:iH) {
      mPredProbs_S_t_Alpha[t, mIndices[h, 1]] = mPredProbs_S_t_Alpha[t, mIndices[h, 1]] + mPredProbs[t, h]
      mPredProbs_S_t_Eta[t, mIndices[h, 3]] = mPredProbs_S_t_Eta[t, mIndices[h, 3]] + mPredProbs[t, h]
      mPredProbs_Z_t[t, mIndices[h, 2]] = mPredProbs_Z_t[t, mIndices[h, 2]] + mPredProbs[t, h]

      vE_A_t[t] = vE_A_t[t] + mPredProbs[t, h] * vY[t-1] * (vAlpha[mIndices[h, 1]]*(1 - vO[t]) + dVarPhi * vO[t])

      foo_a = 0
      if(vY[t-1]>0) {
        for (s in 1:vY[t-1]) {
          for (m in 1:vY[t-1]) {
            if (s == m) {
              foo_a = foo_a + mPredProbs[t, h] * (vAlpha[mIndices[h, 1]]*(1 - vO[t]) + dVarPhi * vO[t])
            } else {
              foo_a = foo_a + mPredProbs[t, h] * (vAlpha[mIndices[h, 1]]*(1 - vO[t]) + dVarPhi * vO[t])^2
            }
          }
        }
      }

      vE_A2_t[t] = vE_A2_t[t] + foo_a

      vE_Eta_t[t] = vE_Eta_t[t] + mPredProbs[t, h] * vLambda[mIndices[h, 2]] * vBeta[vSeason[t]]
      vE_Eta2_t[t] = vE_Eta2_t[t] + mPredProbs[t, h] * vLambda[mIndices[h, 2]] * vBeta[vSeason[t]]*(vLambda[mIndices[h, 2]] * vBeta[vSeason[t]] + 1)

      vE_A_t_Eta_t[t] = vE_A_t_Eta_t[t] + mPredProbs[t, h] * vY[t-1] * (vAlpha[mIndices[h, 1]]*(1 - vO[t]) + dVarPhi * vO[t]) * vLambda[mIndices[h, 2]] * vBeta[vSeason[t]]

    }
  }

  vVar_A_t = vE_A2_t - vE_A_t^2
  vVar_Eta_t = vE_Eta2_t - vE_Eta_t^2
  vCov_A_t_Eta_t = vE_A_t_Eta_t - vE_A_t * vE_Eta_t

  vE_Y_t = vE_A_t + vE_Eta_t
  vE_Y2_t = vE_A2_t + vE_Eta2_t + 2*vE_A_t_Eta_t
  vVar_Y_t = vVar_A_t + vVar_Eta_t + 2*vCov_A_t_Eta_t

  lOut = list(
    vE_Y_t = vE_Y_t,
    vE_Y2_t = vE_Y2_t,
    vVar_Y_t = vVar_Y_t,
    vE_A_t = vE_A_t,
    vE_A2_t = vE_A2_t,
    vVar_A_t = vVar_A_t,
    vE_Eta_t = vE_Eta_t,
    vE_Eta2_t = vE_Eta2_t,
    vVar_Eta_t = vVar_Eta_t,
    vCov_A_t_Eta_t = vCov_A_t_Eta_t
  )

  return(lOut)

}

