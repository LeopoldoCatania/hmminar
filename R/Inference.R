ICfun <- function(llk, np, iT) {
  
  AIC <- -2 * (llk - np)
  BIC <- -2 * llk + np * log(iT)
  HQC <- -2 * llk + 2 * np * log(log(iT))
  
  IC = c(AIC = AIC,
         BIC = BIC,
         HQC = HQC,
         np  = np,
         llk = llk)
  
  return(IC)
}

InferenceFun_Vol <- function(mHessian, vPw) {
  
  vPn = pw2pn_SumZTNB(vPw)
  iK_s = length(vPn)
  
  out = matrix(NA, iK_s, 4L, dimnames = list(names(vPn), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
  
  out[, "Estimate"] = vPn
  
  mJacob = jacobian(pw2pn_SumZTNB, vPw)
  mInvHessian = ginv(mHessian)
  mSandwitch = t(mJacob) %*% mInvHessian %*% mJacob
  
  vSE = sqrt(diag(mSandwitch))
  vTest = vPn/vSE
  vPvalues = 1 - pnorm(abs(vTest))
  
  
  out[, "Std. Error"] = vSE
  out[, "t value"] = vTest
  out[, "Pr(>|t|)"] = vPvalues
  
  return(out)
}


InferenceFun_MixInar <- function(Fit, vY, mD) {
  
  iD = ncol(mD)
  iT = length(vY)
  np = Fit$IC["np"]
  iJ = (np - iD)/2.0
  
  vLambda = Fit$vLambda
  vOmega  = Fit$vOmega
  dAlpha = Fit$dAlpha
  vBeta = Fit$vBeta
  
  names(vLambda) = paste("lambda", 1:iJ, sep = "")
  names(vOmega) = paste("omega", 1:iJ, sep = "")
  names(vBeta) = paste("beta", 1:iD, sep = "")
  
  vPn = c(vLambda, vOmega, vBeta, "alpha" = dAlpha)
  
  vPw = pn2pw_MixPoissonINAR(vPn, iJ, iD)
  
  vLogFactorial_K = lfactorial(0:max(vY))
  
  mS = matrix(NA, iT - 1, np)
  
  for (t in 2:iT) {
    
    mS[t - 1, ] = grad(function(vPw, iJ, iD, dY, dX, vLogFactorial_K){
      
      lPn = pw2pn_MixPoissonINAR(vPw, iJ, iD)
      
      -dMixPB(dY, dX, lPn$dAlpha, lPn$vLambda * lPn$vBeta[mD[t, ] == 1],
             lPn$vOmega, vLogFactorial_K, TRUE)
      
    }, vPw, iJ = iJ, iD = iD, dY = vY[t], dX = vY[t - 1], vLogFactorial_K = vLogFactorial_K)
    
  }
  
  mHessian = t(mS) %*% mS
  
  out = matrix(NA, length(vPn), 4L, dimnames = list(names(vPn), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
  
  out[, "Estimate"] = vPn
  
  mJacob = jacobian(pw2pn_MixPoissonINAR, vPw, iJ = iJ, iD = iD, Vector = TRUE)
  mInvHessian = ginv(mHessian)
  mSandwitch = mJacob %*% mInvHessian %*% t(mJacob)
  
  vSE = sqrt(diag(mSandwitch)/iT)
  vTest = vPn/vSE
  vPvalues = 1 - pnorm(abs(vTest))
  
  out[, "Std. Error"] = vSE
  out[, "t value"] = vTest
  out[, "Pr(>|t|)"] = vPvalues
  
  return(out)
}


