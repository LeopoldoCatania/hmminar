ACF_Plot <- function(lMoments, cov = TRUE, percentage = FALSE) {



  ilag_max = length(lMoments$vAutoCor_Y)

  dVar_Y = as.numeric(lMoments$dVar_Y)

  zero = lMoments$vAutoCov_Y
  first = lMoments$vAutoCov_A_Eta
  second = (lMoments$vAutoCov_A_Eta + lMoments$vAutoCov_Eta_A)
  third = (lMoments$vAutoCov_A_Eta +
             lMoments$vAutoCov_Eta_A +
             lMoments$vAutoCov_Eta)
  fourth = (lMoments$vAutoCov_A_Eta +
              lMoments$vAutoCov_Eta_A +
              lMoments$vAutoCov_Eta +
              lMoments$vAutoCov_A)

  if (!cov) {
    zero = zero/dVar_Y
    first = first/dVar_Y
    second = second/dVar_Y
    third = third/dVar_Y
    fourth = fourth/dVar_Y
  }

  if (percentage) {
    zero = zero/fourth
    first = first/fourth
    second = second/fourth
    third = third/fourth
    fourth = fourth/fourth
#    plot(1:ilag_max, zero, type = "s",
#         las = 0, ylab = "", xlab = "", ylim = c(0,1),panel.first=grid())
    plot(1:ilag_max, zero, lwd=4, type='l',
         las = 0, ylab = "", xlab = "", ylim = c(0,1),panel.first=grid(), cex.axis = 1.5)
  } else {
    if (!cov) {
      plot(1:ilag_max, zero, lwd=4, type='l', ylim = c(0,1),
           las = 1, ylab = "", xlab = "", yaxt = "n",panel.first=grid(), cex.axis = 1.5)
    } else {
      plot(1:ilag_max, zero, lwd=4, type='l',
         las = 1, ylab = "", xlab = "",panel.first=grid(), cex.axis = 1.5)
    }
  }

  lines(1:ilag_max, first, col = "orange")
  lines(1:ilag_max, second, col = "purple")
  lines(1:ilag_max, third, col = "red")
  lines(1:ilag_max, fourth, col = "blue")
  # lines(1:ilag_max, zero, col = "black", lwd = 2)

  x <- c(1:ilag_max, ilag_max:1)
  y = c(rep(0, ilag_max), rev(first))
  polygon(x, y,
          col = "orange", lwd = 1, border = NA)

  y = c(first, rev(second))
  polygon(x, y,
          col = "purple", lwd = 1, border = NA)

  y = c(second, rev(third))
  polygon(x, y,
          col = "red", lwd = 1, border = NA)

  y = c(third, rev(fourth))
  polygon(x, y,
          col = "blue", lwd = 1, border = NA)
  #if (cov) {
   # legend("topright", legend = c(expression(paste(Cov,"(",A[t],",", eta[t-k], ")", sep = "")),
  #                               expression(paste(Cov,"(",eta[t],",", A[t-k], ")", sep = "")),
  #                               expression(paste(Cov,"(",eta[t],",", eta[t-k], ")", sep = "")),
  #                               expression(paste(Cov,"(",A[t],",", A[t-k], ")", sep = ""))),
  #        col = c("orange", "purple", "red", "blue"), lty = 1, lwd = 6)

  #} else {
  # legend("topright", legend = c(expression(paste("(",A[t],",", eta[t-k], ")", sep = "")),
  #                               expression(paste("(",eta[t],",", A[t-k], ")", sep = "")),
  #                               expression(paste("(",eta[t],",", eta[t-k], ")", sep = "")),
  #                               expression(paste("(",A[t],",", A[t-k], ")", sep = ""))),
  #        col = c("orange", "purple", "red", "blue"), lty = 1, lwd = 6)
  #}

  axis(1, at = 1:ilag_max, tcl = -0.2, labels = FALSE)
  if (!cov) {
    axis(2, at = seq(0,1,0.2), labels = format(seq(0,1,0.2),1), las = 2, cex.axis = 1.5)
    axis(2, at = seq(0,1,0.1), tcl = -0.2, labels = FALSE, cex.axis = 1.5)
  }
  points(1:ilag_max, zero, bg='tomato2', pch=21, cex=2, lwd=3)
}
