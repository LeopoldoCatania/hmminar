RandomizedPIT_NB_Softplus <- function(vY, vN, dP, seed = NULL) {

  if (is.null(seed)) {
    seed = sample(1:1e5, 1)
  }

  set.seed(seed)

  iT = length(vY)
  vU = runif(iT)
  vZ = numeric(iT)

  for (t in 1:iT) {
    dF_y   = pNBin(vY[t], vN[t], dP)
    dF_ym1 = pNBin(vY[t] - 1, vN[t], dP)
    vZ[t] = (1.0 - vU[t]) * dF_ym1 + vU[t] * dF_y
  }

  return(list("vZ" = vZ, "seed" = seed))

}

f.freq.2 = function(vPvals, dAlphaTest, bDoFdr = FALSE) {
  vPvals = apply(as.matrix(vPvals),
        MARGIN = 1,
        function(x) min(x, na.rm = TRUE))
  if (isTRUE(bDoFdr)) {
    vPvals.adj = p.adjust(vPvals, method = "fdr")
    sFreq  = sum(vPvals.adj < dAlphaTest) / length(vPvals.adj)
  } else {
    sFreq = sum(vPvals < dAlphaTest) / length(vPvals)
  }
  return(list(sFreq= sFreq, vPvals = vPvals))
}
.f.roll = function(x) {

  tmp = matrix(data = NA, nrow = length(251:2000), ncol = 3)
  for (i in 1:1750) {
    x_sub = x[i:(250 + i - 1)]
    tmp[i,1] = sd(x_sub)
    tmp[i,2] = moments::skewness(x_sub)
    tmp[i,3] = moments::kurtosis(x_sub)
  }
  out = apply(tmp, 2, mean)
  return(out)
}
f.roll = compiler::cmpfun(.f.roll)

## Additional Function
TickLoss <- function(vVaR, vY, dAlpha) {
  vHit  = as.numeric(vY < vVaR)
  vLoss = (dAlpha - vHit)*(vY - vVaR)
  return(vLoss)
}

FZLoss <- function(vY, vVaR, dTau, vES) {
  vHit  = as.numeric(vY < vVaR)
  vLoss = -1/(dTau*vES) * vHit * (vVaR - vY) + vVaR/vES + log(-vES) - 1
  return(vLoss)
}

# Taylor (2017)
ALLoss <- function(vY, vVaR, dTau, vES) {
  vHit  = as.numeric(vY < vVaR)
  vLoss = -log((dTau - 1)/vES) - ((vY - vVaR) * (dTau - vHit))/(dTau*vES)
  return(vLoss)
}

CleanVaR <- function(vVaR, dThreshold = 1e3) {
  require(zoo)

  if (is.na(vVaR[1])) vVaR[1] = vVaR[which(!is.na(vVaR))[1]]

  vAboveThreshold = abs(vVaR) > dThreshold

  if (any(vAboveThreshold)) {
    vVaR[vAboveThreshold] = NA
    warning(paste("There are", length(which(vAboveThreshold)), "abs(VaR) values above the threshold", dThreshold))
  }

  vNAS = is.na(vVaR)

  if (any(vNAS)) {
    vVaR = na.locf(vVaR)
    warning(paste("There are", length(which(vNAS)), "VaR NAs"))
  }

  return(vVaR)
}

f.ttest = function(x, y, dAlphaTest, bDoFdr) {
  require("sandwich")
  require("lmtest")

  m.x = f.freq(x, dAlphaTest, bDoFdr)
  m.y = f.freq(y, dAlphaTest, bDoFdr)

  n.x   = length(x)
  if(isTRUE(bDoFdr)){
	  x = p.adjust(x, method = "fdr")
	  y = p.adjust(y, method = "fdr")
  }
  tmp.x = as.numeric(x < dAlphaTest)
  tmp.y = as.numeric(y < dAlphaTest)
  tmp.diff = tmp.x - tmp.y
  out = 0
  if (sd(tmp.diff) > 0) {
    fit = lm(tmp.diff ~ 1)
    tmp = lmtest::coeftest(fit)
    try({
      tmp = lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit, prewhite = TRUE))
    })
    tstat = mean(tmp.diff) / tmp[,2]
    #pval = 2 * (1.0 - stats::pnorm(abs(tstat)))
    pval = 2 * (1.0 - stats::pt(abs(tstat), df = n.x - 1))

    if (pval <= dAlphaTest) {
      if (m.x <= m.y) {
        out = -1
      } else {
        out = 1
      }
    }
  }
  return(out)
}

f.PIT.tests <- function(x, alpha = 1){

  f.LR <- function(x) {
    tmp = rugarch::BerkowitzTest(x, lags = 1, significance = 0.05, tail.test = FALSE, alpha = 0.1)
    out = list(pval = tmp$LRp, stat = tmp$LR)
    return(out)
  }

  id        <- x < alpha  # <=> obs. below predicted VaR
  x         <- x[id]
  x         <- x / alpha
  x[x == 0] <- 1e-10
  x         <- qnorm(p = x)

  s.LR   <- f.LR(x)$stat
  p.LR   <- f.LR(x)$pval

  out  <- list(Stat = s.LR, Pvalue = p.LR)
  return(out)
}

## rejection frequency
f.freq = function(vPvals, dAlphaTest, bDoFdr = FALSE) {

  if (isTRUE(bDoFdr)) {
    vPvals = p.adjust(vPvals, method = "fdr")
    sFreq  = sum(vPvals < dAlphaTest) / length(vPvals)
  } else {
    sFreq = sum(vPvals < dAlphaTest) / length(vPvals)
  }

  return(sFreq)
}

## HAC standard errors
f.std_mean <- function(vY, bDoHac = TRUE){

  fit = lm(vY ~ 1)

  if (isTRUE(bDoHac)) {
    library("sandwich")
    library("lmtest")

    sumfit = lmtest::coeftest(fit, vcov. = sandwich::vcovHAC(fit, prewhite = TRUE))
    std_mean = sumfit[1,2]
  } else {
    sumfit = summary(fit)
    std_mean = sumfit$coefficients[2]
  }

  return(std_mean)
}

f.pvalue <- function(dTest) {
  dPvalue = 2 * (1.0 - stats::pnorm(abs(dTest)))
  return(dPvalue)
}

# Diebold Mariano 1995 test
f.DMTest <- function(vLoss1, vLoss2, vD = NULL, bDoHac = TRUE){

  if (is.null(vD)) {
    vD = vLoss1 - vLoss2
  }

  uql = quantile(vD, 0.99)
  vD[vD >= uql] = uql

  dMu  = mean(vD, na.rm = TRUE)
  #dMu  = median(vD)

  dSig = f.std_mean(vD, bDoHac = bDoHac)
  #dSig = (sd(vD) / sqrt(lengh(vD))
  #dSig = diff(quantile(vD, probs = c(0.25, 0.75))) / 1.35

  # test and pvalue
  dTest = dMu / dSig
  dPvalue = f.pvalue(dTest)

  out = list(dMu = dMu, dTest = dTest, dPvalue = dPvalue)
  return(out)
}

f.lDMTest <- function(lLoss, bDoHac = TRUE){
  out = f.DMTest(vLoss1 = lLoss$vLoss1, vLoss2 = lLoss$vLoss2, vD = NULL, bDoHac = bDoHac)
  return(out)
}

# Empirical cdf
f.ECDF <- function(vObs){

  iK = length(vObs)

  vX = sort(vObs)
  vY = (1:iK)/iK

  return(list(x = vX, y = vY))
}

BinTest <- function(pit, g = 20L, alpha = 0.05, plot = FALSE, main, ...) {

  h = hist(pit, nclass = g, plot = FALSE)

  n_i = h$counts
  test = sum((n_i - mean(n_i))^2/mean(n_i))
  crit = qchisq(1.0 - alpha, g - 1L)
  pvalue = 1 - pchisq(test, g - 1L)
  confidence = mean(n_i) + c(-qnorm(1 - alpha) * sqrt(mean(n_i)), +qnorm(1 - alpha) * sqrt(mean(n_i)))

  confidence = confidence/mean(n_i)

  h = hist(pit, nclass = g, plot = plot, freq = FALSE, main = main, ylab = "", xlab = "", col = "blue",
           ylim = c(0, max(h$density) * 1.1), ... )

  if (plot) {
    abline(h = confidence, col = "red", lwd = 2L, xlim = c(0, 1))
  }

  out = list(test = test, crit = crit, pvalue = pvalue, hist = h, confidence = confidence)

  return(out)

}


IIDTest <- function(pit, alpha = 0.05, lag = 20) {

  N = length(pit)

  m1 = as.numeric(pit - mean(pit))
  m2 = as.numeric((pit - mean(pit))^2)
  m3 = as.numeric((pit - mean(pit))^3)
  m4 = as.numeric((pit - mean(pit))^4)

  data1 = do.call(cbind, lapply(1:lag, function(i) c(m1[-(1:i)], rep(NA, i))))
  data1 = data.frame(head(data1, N - lag))
  data2 = do.call(cbind, lapply(1:lag, function(i) c(m2[-(1:i)], rep(NA, i))))
  data2 = data.frame(head(data2, N - lag))
  data3 = do.call(cbind, lapply(1:lag, function(i) c(m3[-(1:i)], rep(NA, i))))
  data3 = data.frame(head(data3, N - lag))
  data4 = do.call(cbind, lapply(1:lag, function(i) c(m4[-(1:i)], rep(NA, i))))
  data4 = data.frame(head(data4, N - lag))

  m1 = head(m1, N - lag)
  m2 = head(m2, N - lag)
  m3 = head(m3, N - lag)
  m4 = head(m4, N - lag)

  fit1 = lm(m1 ~ ., data = data1)
  fit2 = lm(m2 ~ ., data = data2)
  fit3 = lm(m3 ~ ., data = data3)
  fit4 = lm(m4 ~ ., data = data4)

  test1 = (N - lag) * summary(fit1)$r.squared
  test2 = (N - lag) * summary(fit2)$r.squared
  test3 = (N - lag) * summary(fit3)$r.squared
  test4 = (N - lag) * summary(fit4)$r.squared

  crit = qchisq(1 - alpha, lag)

  pvalue1 = 1 - pchisq(test1, lag)
  pvalue2 = 1 - pchisq(test2, lag)
  pvalue3 = 1 - pchisq(test3, lag)
  pvalue4 = 1 - pchisq(test4, lag)


  out = list(test = c(test1 = test1, test2 = test2, test3 = test3, test4 = test4), crit = crit, pvalue = c(pvalue1,
                                                                                                           pvalue2, pvalue3, pvalue4),
             lag = lag)
  return(out)
}


StandardizeFit <- function(Fit) {

  vBeta = Fit$lPn$vBeta
  Fit$lPn$vBeta = vBeta/vBeta[1]
  Fit$lPn$vLambda = Fit$vLambda*vBeta[1]
  Fit$lPn$mGamma_Eta = GammaCheck(Fit$lPn$mGamma_Eta, Fit$iL)
  Fit$lPn$mGamma_Alpha = GammaCheck(Fit$lPn$mGamma_Alpha, Fit$iJ)
  Fit$lPn$mOmega = OmegaCheck(Fit$lPn$mOmega, Fit$iK, Fit$iL)

  return(Fit)
}


### create tables for estimated parameters

addBar <- function(x) {
  x[, ncol(x)] = paste(x[, ncol(x)], "\\\\")
  return(x)
}
AddSE <- function(mX, mX_SE, digits = 4) {
  NAs = is.na(mX_SE)
  mX_SE = format(round(mX_SE, digits), scientific = FALSE)
  if(any(NAs)) {
    mX_SE[NAs] = "--"
  }
  mX[] = paste("\\underset{(", mX_SE, ")}{", mX, "}", sep = "")
  return(mX)
}

ComputeTables <- function(Fit, vSD, sPath_tab, sModel) {

  iJ = Fit$iJ
  iK = Fit$iK
  iL = Fit$iL
  iD = Fit$iD
  lPn = Fit$lPn

  mGamma_Alpha = lPn$mGamma_Alpha
  mGamma_Eta = lPn$mGamma_Eta
  mOmega = lPn$mOmega
  vLambda = lPn$vLambda
  vAlpha = as.matrix(lPn$vAlpha)
  dVarPhi = lPn$dVarPhi
  vBeta = lPn$vBeta

  mGamma_Alpha_SD = mGamma_Alpha
  mGamma_Eta_SD = mGamma_Eta
  mOmega_SD = mOmega
  vLambda_SD = vLambda
  vAlpha_SD = vAlpha
  dVarPhi_SD = dVarPhi
  vBeta_SD = vBeta

  mGamma_Alpha_SD[] = NA
  mGamma_Eta_SD[] = NA
  mOmega_SD[] = NA
  vLambda_SD[] = NA
  vBeta_SD[] = NA
  vAlpha_SD[] = NA
  dVarPhi_SD[] = NA

  if (iJ > 1) {
    for (j in 1:iJ) {
      foo = paste("gamma_alpha_", j, "_", 1:(iJ-1), sep = "")
      mGamma_Alpha_SD[j, 1:(iJ-1)] = vSD[foo]
    }
  }

  if(iL > 1) {
    for (l in 1:iL) {
      foo = paste("gamma_eta_", l, "_", 1:(iL-1), sep = "")
      mGamma_Eta_SD[l, 1:(iL-1)] = vSD[foo]

      foo = paste("omega_", l, "_", 1:(iK-1), sep = "")
      mOmega_SD[1:(iK-1), l] = vSD[foo]
    }
  } else {
    if (iK > 1) {
      foo = paste("omega_", 1, "_", 1:(iK-1), sep = "")
      mOmega_SD[1:(iK-1), 1] = vSD[foo]
    }
  }

  foo = paste("alpha_", 1:iJ, sep = "")
  vAlpha_SD = vSD[foo]

  foo = paste("lambda_", 1:iK, sep = "")
  vLambda_SD = vSD[foo]

  foo = paste("beta_", 1:iD, sep = "")
  vBeta_SD = vSD[foo]

  mBeta = matrix(vBeta, 9, 9)
  mBeta_SD = matrix(vBeta_SD, 9, 9)

  dVarPhi_SD = vSD["varphi"]

  mGamma_Alpha = format(round(mGamma_Alpha, 3), scientific = FALSE)
  mGamma_Eta = format(round(mGamma_Eta, 3), scientific = FALSE)
  mOmega = format(round(mOmega, 3), scientific = FALSE)
  vAlpha = format(round(vAlpha, 3), scientific = FALSE)
  vLambda = format(round(vLambda, 3), scientific = FALSE)
  dVarPhi = format(round(dVarPhi, 3), scientific = FALSE)
  mBeta = format(round(mBeta, 3), scientific = FALSE)

  mGamma_Alpha = AddSE(mGamma_Alpha, mGamma_Alpha_SD)
  mGamma_Eta = AddSE(mGamma_Eta, mGamma_Eta_SD)
  mOmega = AddSE(mOmega, mOmega_SD)
  vAlpha = AddSE(vAlpha, vAlpha_SD)
  vLambda = AddSE(vLambda, vLambda_SD)

  mBeta = AddSE(mBeta, mBeta_SD)
  dVarPhi = AddSE(dVarPhi, dVarPhi_SD)

  mGamma_Alpha = addBar(mGamma_Alpha)
  mGamma_Eta = addBar(mGamma_Eta)
  mOmega = addBar(mOmega)
  vAlpha = addBar(vAlpha)
  vLambda = addBar(vLambda)
  mBeta = addBar(mBeta)

  write.table(mGamma_Alpha, file = paste(sPath_tab, sModel, "_mGamma_Alpha_1m.txt", sep = ""), quote = FALSE, sep = " & ", row.names = FALSE, col.names = FALSE)
  write.table(mGamma_Eta, file = paste(sPath_tab, sModel, "_mGamma_Eta_1m.txt", sep = ""), quote = FALSE, sep = " & ", row.names = FALSE, col.names = FALSE)
  write.table(mOmega, file = paste(sPath_tab, sModel, "_mOmega_1m.txt", sep = ""), quote = FALSE, sep = " & ", row.names = FALSE, col.names = FALSE)
  write.table(vAlpha, file = paste(sPath_tab, sModel, "_vAlpha_1m.txt", sep = ""), quote = FALSE, sep = " & ", row.names = FALSE, col.names = FALSE)
  write.table(vLambda, file = paste(sPath_tab, sModel, "_vLambda_1m.txt", sep = ""), quote = FALSE, sep = " & ", row.names = FALSE, col.names = FALSE)
  write.table(mBeta, file = paste(sPath_tab, sModel, "_mBeta_1m.txt", sep = ""), quote = FALSE, sep = " & ", row.names = FALSE, col.names = FALSE)
  write.table(dVarPhi, file = paste(sPath_tab, sModel, "_vVarpi_1m.txt", sep = ""), quote = FALSE, sep = " & ", row.names = FALSE, col.names = FALSE)


}

f.to.date <- function(x){
  year = substr(x, 1, 4)
  month = substr(x, 5, 6)
  day = substr(x, 7, 8)
  paste(year, month, sep = "/")
}

make390seasons <- function(vBeta) {
  iD = 390

  #Seasonality in the NPS specification
  mD_1 = matrix(0, iD, 81)
  #opening
  mD_1[1,1] = 1
  mD_1[2,2] = 1
  mD_1[3,3] = 1
  mD_1[4:5,4] = 1

  vSeasons_range = seq(5, iD, 5)
  for (j in 1:(length(vSeasons_range)-1)) {
    mD_1[(vSeasons_range[j]+1):vSeasons_range[j+1], j+4] = 1
  }

  vD = rep(NA, iD)
  vD[1] = vBeta[1]
  vD[2] = vBeta[2]
  vD[3] = vBeta[3]
  vD[5] = vBeta[4]
  vD[seq(10, iD, 5)] = vBeta[-c(1:4)]

  #Interpolate
  vD[4] = (vD[5] + vD[3])/2

  vEps = (vD[seq(10, iD, 5)] - vD[seq(5, iD-5, 5)])/5

  i = 0
  for (l in seq(5, iD-5, 5)) {
    i = i + 1
    for (q in 1:4) {
      vD[l + q] = vD[l + q - 1] + vEps[i]
    }
  }
  return(vD)
}


AdjustFitForFiltering <- function(Fit_PS) {
  Fit_PS$vBeta = Fit_PS$vDelta_S
  Fit_PS$mGamma_Eta = Fit_PS$lPn$mGamma_Eta
  Fit_PS$mGamma_Alpha = Fit_PS$lPn$mGamma_Alpha
  Fit_PS$mOmega = Fit_PS$lPn$mOmega
  Fit_PS$vLambda = Fit_PS$lPn$vLambda
  Fit_PS$vAlpha = Fit_PS$lPn$vAlpha
  Fit_PS$dVarPhi = Fit_PS$lPn$dVarPhi
  return(Fit_PS)
}

FromNPStoPS <- function(Fit_NPS, iH = 0) {

  iJ = Fit_NPS$iJ
  iL = Fit_NPS$iL
  iK = Fit_NPS$iK

  vBeta = Fit_NPS$vBeta

  iD = 390

  #Seasonality in the NPS specification
  mD_1 = matrix(0, iD, 81)
  #opening
  mD_1[1,1] = 1
  mD_1[2,2] = 1
  mD_1[3,3] = 1
  mD_1[4:5,4] = 1

  vSeasons_range = seq(5, iD, 5)
  for (j in 1:(length(vSeasons_range)-1)) {
    mD_1[(vSeasons_range[j]+1):vSeasons_range[j+1], j+4] = 1
  }

  vD = rep(NA, iD)
  vD[1] = vBeta[1]
  vD[2] = vBeta[2]
  vD[3] = vBeta[3]
  vD[5] = vBeta[4]
  vD[seq(10, iD, 5)] = vBeta[-c(1:4)]

  #Interpolate
  vD[4] = (vD[5] + vD[3])/2

  vEps = (vD[seq(10, iD, 5)] - vD[seq(5, iD-5, 5)])/5

  i = 0
  for (l in seq(5, iD-5, 5)) {
    i = i + 1
    for (q in 1:4) {
      vD[l + q] = vD[l + q - 1] + vEps[i]
    }
  }

  vParS = StartingParSeasonal(vD = vD, iH = iH)

  lPn = Fit_NPS$lPn
  lPn = lPn[-which(names(lPn) == "vBeta")]
  lPn[["vParS"]] = vParS

  mOmega = lPn$mOmega
  for(l in 1:iL) {
    for (k in 1:iK) {
      if(mOmega[k, l] < 1e-2) {
        mOmega[k, l] = 1e-2
      }
      if(mOmega[k, l] > 1 - 1e-2) {
        mOmega[k, l] = 1 - 1e-2
      }
    }
    mOmega[, l] = mOmega[, l]/sum(mOmega[, l])
  }

  lPn$mOmega = mOmega

  if(iJ > 1) {
    mGamma_Alpha = lPn$mGamma_Alpha
    for(j1 in 1:iJ) {
      for(j2 in 1:iJ) {
        if(mGamma_Alpha[j1, j2] < 1e-2) {
          mGamma_Alpha[j1, j2] = 1e-2
        }
        if(mGamma_Alpha[j1, j2] > 1 - 1e-2) {
          mGamma_Alpha[j1, j2] = 1 - 1e-2
        }
      }
      mGamma_Alpha[j1, ] = mGamma_Alpha[j1, ]/sum(mGamma_Alpha[j1, ])
    }
    lPn$mGamma_Alpha = mGamma_Alpha
  }

  if(iL > 1) {
    mGamma_Eta = lPn$mGamma_Eta
    for(l1 in 1:iL) {
      for(l2 in 1:iL) {
        if(mGamma_Eta[l1, l2] < 1e-2) {
          mGamma_Eta[l1, l2] = 1e-2
        }
        if(mGamma_Eta[l1, l2] > 1 - 1e-2) {
          mGamma_Eta[l1, l2] = 1 - 1e-2
        }
      }
      mGamma_Eta[l1, ] = mGamma_Eta[l1, ]/sum(mGamma_Eta[l1, ])
    }
    lPn$mGamma_Eta = mGamma_Eta
  }

  for(j in 1:iJ) {
    if(lPn$vAlpha[j] < 1e-2) {
      lPn$vAlpha[j] = 1e-2
    }
    if(lPn$vAlpha[j] > 1 - 1e-2) {
      lPn$vAlpha[j] = 1 - 1e-2
    }
  }
  for(k in 1:iK) {
    if(lPn$vLambda[k] < 1e-2) {
      lPn$vLambda[k] = 1e-2
    }
    if(lPn$vLambda[k] > 1e3) {
      lPn$vLambda[k] = 1e3
    }
  }

  vY = c(Fit_NPS$dY0, Fit_NPS$vY)
  vO = c(Fit_NPS$dO0, Fit_NPS$vO)

  Fit_PS = Estimate_MSMixPoisInar_TwoChains_ParSeasonal_EM_AlfaCor(vY, iJ, iK, iL, iD, vO = vO,
                                                                   lPn = lPn, iH = iH)

  return(Fit_PS)

}



