Filtering_MixInar <- function(Fit, vY, mD) {
  
  iD = ncol(mD)
  iT = length(vY)
  np = Fit$IC["np"]
  iJ = (np - iD)/2.0
  
  vLambda = Fit$vLambda
  vOmega  = Fit$vOmega
  dAlpha = Fit$dAlpha
  vBeta = Fit$vBeta
  
  lOut = Filtering_MixInar_C(vY, iJ, vOmega, dAlpha, vLambda, vBeta, mD)
  
  vRes = (vY[-1] - lOut$vMu)/sqrt(lOut$vSigma2)
  
  lOut[["vRes"]] = vRes
  
  return(lOut)
}

Filtering_MSMixInar_v2 <- function(Fit, vY, mD) {
  
  iD = ncol(mD)
  iT = length(vY)
  np = Fit$IC["np"]
  
  vLambda = Fit$vLambda
  vOmega  = Fit$vOmega
  dAlpha = Fit$dAlpha
  vBeta = Fit$vBeta
  vMu   = Fit$vMu
  mGamma = Fit$mGamma
  vDelta = Fit$vDelta
  
  iJ = length(vDelta)
  iM = length(vLambda)
  
  lOut = Filtering_MSMixInar_v2_Seasonal_Dummy(vY, iJ, iM, vOmega, vLambda, vMu,
                                               dAlpha, mGamma, vDelta, vBeta, mD) 
  
  vRes = (vY[-1] - lOut$vMu[-c(1, iT + 1)])/sqrt(lOut$vSigma2[-c(1, iT + 1)])
  
  lOut[["vRes"]] = vRes
  
  return(lOut)
}

Filtering_MSMixInar_v2_Inv <- function(Fit, vY, mD) {
  
  iD = ncol(mD)
  iT = length(vY)
  np = Fit$IC["np"]
  
  vLambda = Fit$vLambda
  vOmega  = Fit$vOmega
  dAlpha = Fit$dAlpha
  vBeta = Fit$vBeta
  vMu   = Fit$vMu
  mGamma = Fit$mGamma
  vDelta = Fit$vDelta
  
  iJ = length(vDelta)
  iM = length(vLambda)
  
  lOut = Filtering_MSMixInar_v2_Seasonal_Dummy_Inv(vY, iJ, iM, vOmega, vLambda, vMu,
                                               dAlpha, mGamma, vDelta, vBeta, mD) 
  
  vRes = (vY[-1] - lOut$vMu[-c(1, iT + 1)])/sqrt(lOut$vSigma2[-c(1, iT + 1)])
  
  lOut[["vRes"]] = vRes
  
  return(lOut)
}

Filtering_MSMixInar <- function(Fit, vY, mD) {
  
  iD = ncol(mD)
  iT = length(vY)
  np = Fit$IC["np"]
  
  vLambda = Fit$vLambda
  vOmega  = Fit$vOmega
  vAlpha = Fit$vAlpha
  vBeta = Fit$vBeta
  mGamma = Fit$mGamma
  vDelta = Fit$vDelta
  
  iJ = length(vDelta)
  iM = length(vLambda)
  
  lOut = Filtering_MSMixInar_Seasonal_Dummy(vY, iJ, iM, vOmega, vLambda,  
                                            vAlpha, mGamma, vDelta, vBeta, mD) 
  
  vRes = (vY[-1] - lOut$vMu[-c(1, iT + 1)])/sqrt(lOut$vSigma2[-c(1, iT + 1)])
  
  lOut[["vRes"]] = vRes
  
  return(lOut)
}

Filtering_MSMixInar_v3 <- function(Fit, vY, mD) {
  
  iD = ncol(mD)
  iT = length(vY)
  np = Fit$IC["np"]
  
  mLambda = Fit$mLambda
  vOmega  = Fit$vOmega
  dAlpha = Fit$vAlpha
  vBeta = Fit$vBeta
  mGamma = Fit$mGamma
  vDelta = Fit$vDelta
  
  iJ = length(vDelta)
  iM = length(vOmega)
  
  lOut = Filtering_MSMixInar_v3_Seasonal_Dummy(vY, iJ, iM, vOmega, mLambda,  
                                              dAlpha, mGamma, vDelta, vBeta, mD) 
  
  vRes = (vY[-1] - lOut$vMu[-c(1, iT + 1)])/sqrt(lOut$vSigma2[-c(1, iT + 1)])
  
  lOut[["vRes"]] = vRes
  
  return(lOut)
}

Filtering_MSMixInar_TwoChains <- function(Fit) {
 
  vLambda = Fit$vLambda
  mOmega  = Fit$mOmega
  vAlpha = Fit$vAlpha
  vBeta = Fit$vBeta
  mGamma_Alpha = Fit$mGamma_Alpha
  mGamma_Eta = Fit$mGamma_Eta
  mD = Fit$mD
  vY = Fit$vY
  
  iJ = Fit$iJ
  iL = Fit$iL
  iK = Fit$iK
  iD = Fit$iD
  iT = length(vY)
  vSeason = Fit$vSeason + 1
  
  PredictedProb = Fit$PredictedProb
  mIndices = Fit$mIndices + 1
  
  vMu = numeric(iT)
  vE2 = numeric(iT)
  
  iQ = iJ * iL * iK
  for (t in 2:iT) {
    for (q in 1:iQ) {
      vMu[t] = vMu[t] + PredictedProb[t, q] * (vAlpha[mIndices[q, 1]] * vY[t - 1] + vLambda[mIndices[q, 2]] * vBeta[vSeason[t]])
      vE2[t] = vE2[t] + PredictedProb[t, q] * (
        # //2nd moment B
        vAlpha[mIndices[q, 1]] * (1.0 - vAlpha[mIndices[q, 1]]) * vY[t - 1]  + 
          (vAlpha[mIndices[q, 1]] * vY[t - 1])^2.0 + 
          # //2nd moment eps
        (vLambda[mIndices[q, 2]] * vBeta[vSeason[t]]) * ((vLambda[mIndices[q, 2]] * vBeta[vSeason[t]]) + 1.0) + 
          # // 2 times prod of first moments
        2.0 * (vAlpha[mIndices[q, 1]] * vY[t - 1] * vLambda[mIndices[q, 2]] * vBeta[vSeason[t]]) 
        );
    }
  }
  
  vSigma2 = vE2 - vMu^2
  
  vRes = ((vY - vMu)/sqrt(vSigma2))[-1]
  
  lOut = list()
  lOut[["vRes"]] = vRes
  lOut[["vMu"]] = vMu
  lOut[["vSigma2"]] = vSigma2
  
  return(lOut)
}

Filtering_MSMixInar_TwoChains_AlfaCor <- function(Fit, ComputeMedian = FALSE, iB = 1e4-1) {
  
  vLambda = Fit$vLambda
  mOmega  = Fit$mOmega
  vAlpha = Fit$vAlpha
  vBeta = Fit$vBeta
  dVarPhi = Fit$dVarPhi
  mGamma_Alpha = Fit$mGamma_Alpha
  mGamma_Eta = Fit$mGamma_Eta
  mD = Fit$mD
  vY = Fit$vY
  vO = Fit$vO
  
  vY = c(Fit$dY0, vY)
  vO = c(Fit$dO0, vO)
  vSeason = c(Fit$iSeason0, Fit$vSeason) + 1
  
  iJ = Fit$iJ
  iL = Fit$iL
  iK = Fit$iK
  iD = Fit$iD
  iT = Fit$iT
  
  PredictedProb = Fit$PredictedProb
  mIndices = Fit$mIndices + 1
  
  vMu = numeric(iT)
  vE2 = numeric(iT)
  vMedian = numeric(iT)
  
  iQ = iJ * iL * iK
  for (t in 2:iT) {
    for (q in 1:iQ) {
      vMu[t] = vMu[t] + PredictedProb[t - 1, q] * ((vAlpha[mIndices[q, 1]] * (1.0 - vO[t]) + vO[t] * dVarPhi) * vY[t - 1] + vLambda[mIndices[q, 2]] * vBeta[vSeason[t]])
      vE2[t] = vE2[t] + PredictedProb[t - 1, q] * (
        # //2nd moment B
        (vAlpha[mIndices[q, 1]] * (1.0 - vO[t]) + vO[t] * dVarPhi) * (1.0 - (vAlpha[mIndices[q, 1]] * (1.0 - vO[t]) + vO[t] * dVarPhi)) * vY[t - 1]  + 
          ((vAlpha[mIndices[q, 1]] * (1.0 - vO[t]) + vO[t] * dVarPhi) * vY[t - 1])^2.0 + 
          # //2nd moment eps
          (vLambda[mIndices[q, 2]] * vBeta[vSeason[t]]) * ((vLambda[mIndices[q, 2]] * vBeta[vSeason[t]]) + 1.0) + 
          # // 2 times prod of first moments
          2.0 * ((vAlpha[mIndices[q, 1]] * (1.0 - vO[t]) + vO[t] * dVarPhi) * vY[t - 1] * vLambda[mIndices[q, 2]] * vBeta[vSeason[t]]) 
      );
    }
    if (ComputeMedian) {
      
      vS = sample(1:iQ, iB, replace = TRUE, prob = PredictedProb[t - 1, ])
      vSim = sapply(vS, function(iS) {
        
        dAlpha = vAlpha[mIndices[iS, 1]]
        dLambda = vLambda[mIndices[iS, 2]]
        
        iA = rbinom(1, vY[t-1], (dAlpha * (1.0 - vO[t]) + vO[t] * dVarPhi))
        iEta = rpois(1, dLambda * vBeta[vSeason[t]])
        
        return(iA + iEta)
        
      })
      vMedian[t] = median(vSim)
    }
  }
  
  vSigma2 = vE2 - vMu^2
  
  vRes = ((vY - vMu)/sqrt(vSigma2))[-1]
  
  lOut = list()
  lOut[["vRes"]] = vRes
  lOut[["vMu"]] = vMu
  lOut[["vSigma2"]] = vSigma2
  lOut[["vMedian"]] = vMedian
  
  return(lOut)
}

Filtering_Volumes <- function(Fit_V, vN, vV) {
  
  dA = Fit_V$vPn["a"]
  dB = Fit_V$vPn["b"]
  
  dMean  = dB * dA/((1 - (1 - dB)^dA) * (1 - dB)) 
  dE2    = dB * dA * (dB * dA + 1)/((1 - (1 - dB)^dA) * (1 - dB)^2) 
  dVar   = dE2 - dMean^2
  
  vMean_Given_N = vN * dMean
  vVar_Given_N  =  vN * dVar
  
  vRes = (vV - vMean_Given_N)/sqrt(vVar_Given_N)
  
  lOut = list()
  
  lOut[["vMean_Given_N"]] = vMean_Given_N
  lOut[["vVar_Given_N"]] = vVar_Given_N
  lOut[["vRes"]] = vRes
  
  return(lOut)
  
}

Predicting_Volumes <- function(Fit_V, Fit_N) {
  
  dA = Fit_V$vPn["a"]
  dB = Fit_V$vPn["b"]
  
  dMean  = dB * dA/((1 - (1 - dB)^dA) * (1 - dB)) 
  dE2    = dB * dA * (dB * dA + 1)/((1 - (1 - dB)^dA) * (1 - dB)^2) 
  dVar   = dE2 - dMean^2
  
  Filtering_N = Filtering_MSMixInar_TwoChains_AlfaCor(Fit)
  
  vMean = Filtering_N$vMu * dMean
  vVar = Filtering_N$vMu * dVar + dMean^2 * Filtering_N$vSigma2
  
  vN = Fit_N$vY
  
  vMean_GivenN = vN * dMean

  return(list(vMean_GivenN = vMean_GivenN,
    vMean = vMean,
              vVar  = vVar, 
              Filtering_N = Filtering_N))
  
}

