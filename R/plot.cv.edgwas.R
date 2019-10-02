#' Plot the cross-validation curve produced by cv.edgwas
#'
#' Plots the cross-validation curve, and empirical 95% confidence band, as a function of the rho values used.
#'
#' @param x Fitted "cv.edgwas" object.
#'
#' @details A plot is produced, and nothing is returned.
#' @export plot.cv.edgwas
#'
plot.cv.edgwas <- function(x, ...) {
  cvobj <- x

  if(all(abs(diff(cvobj$rho)) == 1)) {
    xvals <- cvobj$rho
    xlab <- bquote(log(rho))
  } else {
    xvals <- log(cvobj$rho)
    xlab <- bquote(log(rho))
  }

  plot.default(x = xvals, y = cvobj$cvm,
       ylim = range(cvobj$cvup, cvobj$cvlo),
       xlab = xlab, ylab = cvobj$name,
       type = "n", bty = "n")
  segments(x0 = xvals, x1 = xvals,
           y0 = cvobj$cvlo, y1 = cvobj$cvup,
           col = "grey")
  points(xvals, cvobj$cvm, pch = 20, col = "orangered")
  abline(v = log(cvobj$rho.1se), lty = 3)
  abline(v = log(cvobj$rho.min), lty = 3)

  invisible()
}
