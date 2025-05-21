SumZTNB_Optimizer <- function(vPw, vV, vN, iT, vLogfactorial) {
  
  vPn = pw2pn_SumZTNB(vPw)
  
  # negative average llk
  # dnaLLK = try(-dLLK_SumZTNB_par(vV, vN, vPn["a"], vPn["b"], iCores = 6)/iT, silent = TRUE)
  dnaLLK = try(-dLLK_SumZTNB(vV, vN, vPn["a"], vPn["b"])/iT, silent = TRUE)
  if (is(dnaLLK, "try-error")) {
    dnaLLK = 1e30
  }
  
  if (!is.finite(dnaLLK)) {
    dnaLLK = 1e30
  }
  
  
  return(dnaLLK)
  
}


PoissonSumMixThreeBinom_Optimizer <- function(vPw, vY, iT, iMax, vLogFactorial_K, vNames) {
  
  names(vPw) = vNames
  
  vPn = pw2pn_PoissonSumMixThreeBinom(vPw)
  
  vOmega = vPn[paste("omega", 1:3, sep = "")]
  dLambda = vPn["lambda"]
  vPi = vPn[paste("pi", 1:3, sep = "")]
  
  dnaLLK = try(-dLLKPoissonSumMixThreeBinom_fast(vY, dLambda, vPi,
                              vOmega, vLogFactorial_K, iT, iMax), silent = TRUE)
  
  if (is(dnaLLK, "try-error")) {
    dnaLLK = 1e30
  }
  
  if (!is.finite(dnaLLK)) {
    dnaLLK = 1e30
  }
  
  return(dnaLLK)
  
}

PoissonSumMixTwoBinom_Optimizer <- function(vPw, vY, iT, vLogFactorial_K, iMax, vNames, ...) {
  
  names(vPw) = vNames
  
  vPn = pw2pn_PoissonSumMixTwoBinom(vPw)
  
  dOmega = vPn["omega"]
  dLambda = vPn["lambda"]
  vPi = vPn[paste("pi", 1:2, sep = "")]
  
  dnaLLK = 
    try(
      -dLLKPoissonSumMixTwoBinom(vY, dLambda, vPi,
                                                 dOmega, vLogFactorial_K, iT, iMax, 7)
      , silent = TRUE)
  
  if (is(dnaLLK, "try-error")) {
    dnaLLK = 1e30
  }
  
  if (!is.finite(dnaLLK)) {
    dnaLLK = 1e30
  }
  
  return(dnaLLK)
  
}

PoissonSumMixTwoBinomSeasonal_Optimizer <- function(vPw, vY, iT, iS, vLogFactorial_K, iMax, vNames, ...) {
  
  names(vPw) = vNames
  
  vPn = pw2pn_PoissonSumMixTwoBinom_Seasonal(vPw)
  
  dOmega = vPn["omega"]
  vPi = vPn[paste("pi", 1:2, sep = "")]
  vEta = vPn[paste("eta", 1:3, sep = "")]
  
  dnaLLK = 
    try(
      -dLLKPoissonSumMixTwoBinom_Seasonal(vY, vPi,
                                 dOmega, vEta, vLogFactorial_K, iT, iS, iMax, 7)$dLLK
      , silent = TRUE)
  
  if (is(dnaLLK, "try-error")) {
    dnaLLK = 1e30
  }
  
  if (!is.finite(dnaLLK)) {
    dnaLLK = 1e30
  }
  
  return(dnaLLK)
  
}

SumMixThreeBinom_Optimizer <- function(vPw, vY, iT, dK, vLogFactorial_K, vNames) {
  
  names(vPw) = vNames
  
  vPn = pw2pn_PoissonSumMixThreeBinom(vPw)[-1]
  
  vOmega = vPn[paste("omega", 1:3, sep = "")]
  vPi = vPn[paste("pi", 1:3, sep = "")]
  
  dnaLLK = try(-dLLKSumMixThreeBinom(vY, dK, vPi,
                                     vOmega, vLogFactorial_K, iT), silent = TRUE)
  
  if (is(dnaLLK, "try-error")) {
    dnaLLK = 1e30
  }
  if (!is.finite(dnaLLK)) {
    dnaLLK = 1e30
  }
  return(dnaLLK)
  
}


MixPoissonINAR_Optimizer <- function(vPw, vY, iT, iJ, vLogFactorial_K, vNames) {
  
  names(vPw) = vNames
  
  lPn = pw2pn_MixPoissonINAR(vPw, iJ)
  
  dLLK = try(-dLLK_MixPB(vY, lPn$vLambda, lPn$vOmega, vLogFactorial_K, lPn$dAlpha, iT, iJ), silent = TRUE)
  
  if (is(dLLK, "try-error")){
    dLLK = 1e10
  }
  
  if (!is.finite(dLLK)){
    dLLK = 1e10
  }
  
  return(dLLK)
  
}


MixPoissonINARSeasonal_Optimizer <- function(vPw, vY, iT, iJ, iS, vLogFactorial_K, vNames) {
  
  names(vPw) = vNames
  
  lPn = pw2pn_MixPoissonINARSeasonal(vPw, iJ)
  
  dLLK = try(-dLLK_MixPB_seasonal(vY, lPn$vLambda, lPn$vOmega, lPn$vEta, vLogFactorial_K, lPn$dAlpha, iT, iJ, iS)$dLLK, silent = TRUE)
  
  if (is(dLLK, "try-error")){
    dLLK = 1e10
  }
  
  if (!is.finite(dLLK)){
    dLLK = 1e10
  }
  
  return(dLLK)
  
}

