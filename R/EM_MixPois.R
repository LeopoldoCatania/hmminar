SumMixPoisEM <- function(vY, iJ, dTol = 1e-6, iMaxiter = 4000) {
  
  # Starting Values
  lPar = Starting_MixPois(vY, iJ) 
  
  Fit = EM_SumMixPois(vY = vY, iJ = iJ, 
                      dMu = lPar$dMu, vNu = lPar$vNu, 
                      vOmega = lPar$vOmega,
                      dTol = dTol, iMaxiter = iMaxiter)
  
  return(Fit)
  
}

EM <- function(vY, iJ, iL, mD = NULL, dTol = 1e-7, iMaxiter = 4000) {
  
  dStartTime = Sys.time()
  
  # Starting Values
  lPar = Starting_MSSumMixPois(vY, iL, iJ)
  
  if (is.null(mD)) {
    
    iD = 0
  
  Fit = EM_MSSumMixPois(vY = vY, iJ = iJ, iL = iL, 
                         vMu = lPar$vMu, vNu = lPar$vNu, 
                         vOmega = lPar$vOmega, mGamma = lPar$mGamma,
                         dTol = dTol, iMaxiter = iMaxiter)
  } else {
    iD = ncol(mD)
    
    vBeta = rep(1, iD)
    
    Fit = EM_MSSumMixPois_Dummy(vY = vY, iJ = iJ, iL = iL, 
                                vMu = lPar$vMu, vNu = lPar$vNu, 
                                vOmega = lPar$vOmega, mGamma = lPar$mGamma,
                                vBeta = vBeta, mD = mD,
                                dTol = dTol, iMaxiter = iMaxiter)
    
  }
  
  dTime = Sys.time() - dStartTime
  
  iT = length(vY)
  
  Fit[["iT"]] = iT
  Fit[["iJ"]] = iJ
  Fit[["iL"]] = iL
  Fit[["iD"]] = iL
  Fit[["mD"]] = mD
  
  np = iJ^2 +2 * iL + iD - 1
  
  IC = ICfun(Fit$dLLK, np, iT)
  
  Fit[["IC"]] = IC
  
  Fit[["Time"]] = dTime
  
  return(Fit)
  
}
