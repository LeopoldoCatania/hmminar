library(e1071)

#Set path
sPath = "XXX"

#Load Data
load(file = paste(sPath, "data/vNT_1min.rdata", sep = ""))
load(file = paste(sPath, "data/mD_1min.rdata", sep = ""))
load(file = paste(sPath, "data/vDates_1min.rdata", sep = ""))

source(paste(sPath, "code/functions.R", sep = ""))

vPeriods = c("Full", "Opening", "MidDay", "Closing")

Opening = 1:10 #first half hour (note that first seasons are 1,2,3 minutes and then 4:5, 6:10, etc.)
MidDay  = 34:58 #from 2.5h to 4.5h
Closing = 76:81 #from 6h to 6.5h

mStat = matrix(NA, length(vPeriods), 9, dimnames = list(vPeriods, c("Mean", "Median", "Mode", "Max", "Min", "SD", "Kurt", "Skew", "id")))

mPeridos = cbind(TRUE, apply(mD[, Opening] == 1, 1, any), apply(mD[, MidDay] == 1, 1, any), apply(mD[, Closing] == 1, 1, any))

colnames(mPeridos) = vPeriods

ID_fun <- function(x) var(x)/mean(x)
Mode_fun <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

vFun = c(mean, median, Mode_fun, max, min, sd, kurtosis, skewness, ID_fun)

for (sPeriod in vPeriods) {
  for (i in 1:length(vFun)) {
    f = vFun[[i]]
    mStat[sPeriod, i] = f(vN[mPeridos[, sPeriod]])
  }
}

mStat = format(round(mStat, 2), scientific = FALSE)

mStat= cbind(mStat, vPeriods)

mStat[, ncol(mStat)] = paste(mStat[, ncol(mStat)], "\\\\")

write.table(mStat,
            file =  paste(sPath, "output/DescStat_1m.txt", sep = ""),
            sep = " & ",
            dec = ".",
            quote = FALSE,
            row.names = TRUE
)


## Seasonality
vMean = sapply(1:390, function(s) mean(vN[seq(s, length(vN), 390)]))
vSeq = seq(0, 390, 30)
vLabels = c("09:30", "10:00", "10:30", "11:00", "11:30",
            "12:00", "12:30", "13:00", "13:30", "14:00",
            "14:30", "15:00", "15:30", "16:00")

pdf(file = paste(sPath, "output/vSeasonality_Empirical.pdf", sep = ""), width = 6, height = 4)

par(mar = c(2,4,1,1))
plot(1:390, vMean, type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, cex.axis = 1.5)
grid(10, 10, col = "gray80")
lines(1:390, vMean)

axis(1, at = vSeq, labels = vLabels, cex.axis = 1.5)
axis(1, at = seq(0, 390, 5), labels = FALSE, tcl = -0.25)
axis(1, at = seq(0, 390, 1), labels = FALSE, tcl = -0.25/2)

dev.off()

## AUTOCORRELATION
vACF = acf(vN, plot = FALSE, lag.max = 390 * 5)$acf[-1,,]

pdf(file = paste(sPath, "output/ACF_1m.pdf", sep = ""), width = 6, height = 4)

par(mar = c(2,4,1,1))
plot(1:(390*5), vACF, type = "n", xlab = "", ylab = "", las = 1, ylim = range(vACF), xaxt = "n", cex.axis = 1.5)
grid(10, 10, "gray80")

lines(1:(390*5), vACF, lwd = 2)

abline(h = 1.96/sqrt(length(vN)), col = "red")
abline(h = -1.96/sqrt(length(vN)), col = "red")

axis(1, at = c(1, seq(390, 390*5, 390)), labels = paste(1:6, "D", sep = ""), cex.axis = 1.5)

dev.off()

vDates = sapply(vDates, f.to.date)
vDates = c(sapply(vDates, rep, 390))

pdf(file = paste(sPath_Save, "Trades_1m.pdf", sep = ""), width = 6, height = 4)

par(mar = c(2,4,1,1))
plot(1:(length(vN)), vN, type = "n", xlab = "", ylab = "", las = 1, ylim = range(vN), xaxt = "n", cex.axis = 1.5)
grid(10, 10, "gray80")

lines(1:(length(vN)), vN, lwd = 2)

vSeq = c(1, seq(390*5*5, length(vN), 390*5*5))

axis(1, at = vSeq, labels = vDates[vSeq], cex.axis = 1.5)

dev.off()

vDates_rep = c(sapply(vDates, rep, 390))

vN_day = vN[vDates_rep == "20010509"]

plot.ts(vN_day)

pdf(file = paste(sPath, "output/Trades_20010509_1m.pdf", sep = ""), width = 6, height = 4)

par(mar = c(2,4,1,1))

plot(1:(length(vN_day)), vN_day, type = "n", xlab = "", ylab = "", las = 1, ylim = range(vN_day), xaxt = "n", cex.axis = 1.5)
grid(10, 10, "gray80")

lines(1:(length(vN_day)), vN_day, lwd = 2)

vSeq = seq(0, 390, 30)

vLabels = c("09:30", "10:00", "10:30", "11:00", "11:30",
            "12:00", "12:30", "13:00", "13:30", "14:00",
            "14:30", "15:00", "15:30", "16:00")

axis(1, at = vSeq, labels = vLabels, cex.axis = 1.5)
axis(1, at = seq(0, 390, 5), labels = FALSE, tcl = -0.25)
axis(1, at = seq(0, 390, 1), labels = FALSE, tcl = -0.25/2)

dev.off()

### Test by Harris, D., & McCabe, B. (2019). Semiparametric independence testing for time series of counts and the role of the support. Econometric Theory, 35(6), 1111-1145.

f_Test <- function(vY) {

  iT = length(vY)
  vSupport = 0:max(vY)

  vPi = sapply(vSupport, function(x) mean(vY==x))

  dMu = mean(vY)
  dS_T = 0
  vG = numeric(iT-1)
  for (t in 2:iT) {
    iY_t = vY[t]
    if (iY_t == 0) {
      dG = 0
    } else {
      dG = vPi[(iY_t-1) + 1]/vPi[iY_t + 1] # it includes 0
    }
    vG[t-1] = dG
    dS_T = dS_T + (vY[t-1] - dMu)*(dG - 1)
  }

  dXi_T = dS_T*sqrt(iT)/(var(vY)*var(vG))

  return(dXi_T)

}

vSeq = seq(0, length(vN), 390)

iD = length(vN)/390
vTest = numeric(iD)

for (d in 1:iD) {
  vY = vN[(vSeq[d]+1):vSeq[d+1]]
  vTest[d] = f_Test(vY)
}

## plot
pdf(file = paste(sPath, "output/HM_Test.pdf", sep = ""), width = 6, height = 4)

par(mar = c(2,4,1,1))

plot(1:(iD), vTest, type = "n", xlab = "", ylab = "", las = 1, ylim = range(vTest), cex.axis = 1.5)
grid(10, 10, "gray80")

lines(1:(iD), vTest, lwd = 2)

abline(h = qnorm(0.99), col = "red")

dev.off()





