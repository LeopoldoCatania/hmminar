library(volume)
library(parallel)

sPath = "XXX"

# number of replicates
iB = 10000
#
vT = c(250, 500, 1000, 5000)

aOut = array(NA, c(iB, 10, length(vT), 5), dimnames = list(1:iB, NULL, vT, c("eps", "se", "ttest", "pval", "reject")))
aEst = array(NA, c(iB, 10, length(vT)), dimnames = list(1:iB, NULL, vT))

vFiles = list.files(paste(sPath, "output", sep = ""), full.names = TRUE)

for(b in 1:iB) {
  for (iT in vT) {

    sFile = paste(sPath, "output/T", iT, "_b", b, ".rdata", sep = "")
    if(sFile %in% vFiles) {
      load(paste(sPath, sFile, sep = ""))

      aOut[ib,,iT, "eps"] = c(lOut$eps)
      aOut[ib,,iT, "se"] = c(lOut$vSE)
      aOut[ib,,iT, "ttest"] = c(lOut$vTest)
      aOut[ib,,iT, "pval"] = c(lOut$vP)
      aOut[ib,,iT, "reject"] = c(lOut$vReject)

      foo = pn2pw_MSMixPoisInar_TwoChains_EM(lOut$est, 2, 2, 2, 1, bSeasonal = FALSE, bAlphacor = FALSE)
      vPn_Est = pw2pn_MSMixPoisInar_TwoChains_EM(foo, 2, 2, 2, 1, bSeasonal = FALSE, bAlphacor = FALSE, vecout = TRUE)

      aEst[ib,,iT] = vPn_Est
    }
  }
}

save(aOut, file = paste(sPath, "output/montecarlo_fixed.rdata", sep = ""))

nas = unique(do.call(c, lapply(1:length(vT), function(i) which(is.na(aOut[,1,i,1])))))
aOut_foo = aOut[-nas,,,]
aEst_foo = aEst[-nas,,]

mBias = matrix(NA, dim(aOut_foo)[2], length(vT))
mSE_num = matrix(NA, dim(aOut_foo)[2], length(vT))
mRMSE = matrix(NA, dim(aOut_foo)[2], length(vT))
mReject = matrix(NA, dim(aOut_foo)[2], length(vT))

for (i in 1:length(vT)) {
  mBias[,i] = apply(aOut_foo[,,i,"eps"], 2, mean, na.rm = TRUE)
  mSE_num[,i] = apply(aOut_foo[,,i,"se"], 2, mean, na.rm = TRUE)
  mRMSE[,i] = sqrt(apply(aOut_foo[,,i, "eps"]^2, 2, mean, na.rm = TRUE))
  mReject[,i] = apply(aOut_foo[,,i, "reject"], 2, mean, na.rm = TRUE)
}

vNames = paste("$", c("\\gamma_{11}^\\eta", "\\gamma_{22}^\\eta",
                      "\\gamma_{11}^\\alpha", "\\gamma_{22}^\\alpha",
                      "\\omega_{11}",
                      "\\omega_{21}",
                      "\\alpha_1", "\\alpha_2",
                      "\\lambda_1", "\\lambda_2"), "$", sep = "")

rownames(mBias) = vNames
rownames(mRMSE) = vNames
rownames(mReject) = vNames

colnames(mBias) = vT
colnames(mRMSE) = vT
colnames(mReject) = vT

mBias = format(round(mBias, 3), scientific = FALSE)
mRMSE = format(round(mRMSE, 3), scientific = FALSE)
mReject = format(round(mReject, 3), scientific = FALSE)

mTab = cbind(mBias, mRMSE, mReject)
mTab[, ncol(mTab)] = paste(mTab[, ncol(mTab)], "\\\\")

write.table(mTab, file = paste(sPath, "output/montecarlo.txt"),
            sep = " & ", quote = FALSE, row.names = TRUE, col.names = TRUE)
