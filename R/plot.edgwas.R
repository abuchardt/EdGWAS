#' Plot diagnostics for an edgwas object
#'
#' Two plots (selectable by which) are currently available: a plot of the mean standard error curve produced by edgwas i.e., the MSE and upper and lower standard deviation curves plotted against the values of rho used in the fits, and an adjacency matrix plot for rho.min.
#'
#' @param x Fitted "edgwas" object.
#' @param which If a subset of the plots is required, specify a subset of the numbers \code{1:5}.
#' @param main An overall title for the plot: see \code{\link{title}}
#' @param xlab A title for the x axis: see \code{\link{title}}
#' @param ylab A title for the y axis: see \code{\link{title}}
#' @param ask Logical; if \code{TRUE}, the user is \emph{ask}ed before each plot, see \code{\link{par}(ask=.)}.
#' @param pointsCol The color for points corresponding to rho.min and rho.1se. Default is \code{pointsCol = "#D95F02"}.
#' @param col The color to be used for the curve. Default is \code{col = "grey40"}.
#' @param ylim The y limits of the plot.
#' @param xlim The x limits of the plot.
#' @param ... Other paramters to be passed through to plotting functions.
#'
#' @details Plots are produced, and nothing is returned.
#'
#' @seealso{\code{\link{edgwas}}}
#'
#' @examples
#' N <- 500
#' q <- 10
#' set.seed(1)
#' x <- matrix(sample(0:2, N*q, replace=TRUE), nrow=N, ncol=q)
#' B <- matrix(0, nrow = q, ncol = q)
#' B[1, 1:2] <- 10
#' y <- x %*% B + matrix(rnorm(N*q), nrow = N, ncol = q)
#' ###
#' pc <- edgwas(x, y)
#' plot(pc, 1)
#' plot(pc, 2)
#'
#' @export
#'

plot.edgwas <- function (x, which = c(1L:2L),
                            main = "", xlab = NULL, ylab = NULL,
                            ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                            pointsCol = "#D95F02", ylim = NULL, xlim = NULL, col = NULL,
                            ...) {

  if (!inherits(x, "edgwas"))
    stop("use only with \"edgwas\" objects")
  if(!is.numeric(which) || any(which < 1) || any(which > 2))
    stop("'which' must be in 1:2")
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  rho <- x$rho[(!sapply(x$clusters, is.null))]
  stderr <- x$betaSD[(!sapply(x$clusters, is.null))]
  unr <- length(rho)

  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  ##---------- Do the individual plots : ----------
  if (show[1L]) {

    if(is.null(xlim)) xlim <- range(rho[-unr]) else xlim = xlim

    if(x$logrho == 1) {
      xvals <- log(rho[-unr])
      xlim <- log(xlim)
      if(is.null(xlab))
        xlab <- bquote(log(rho))
    } else {
      xvals <- rho[-unr]
      if(is.null(xlab))
        xlab <- bquote(rho)
    }

    if(is.null(ylim)) ylim <- range(stderr) else ylim <- ylim

    ptCol <- rep(ifelse(is.null(col), "grey40", col), length(rho)-1)
    ptCol[which.min(stderr)] <- pointsCol

    plot(x = xvals, y = x$betaSD[1:(unr-1)], main = main,
         ylim = ylim, xlim = xlim,
         xlab = xlab, ylab = ifelse(is.null(ylab), "Std. error", ylab),
         type = "n", bty = "n", ...)
    abline(h = x$betaSD[unr], col = "grey", lty = 2)
    points(xvals, x$betaSD[1:(unr-1)], pch = 20,
           col = ptCol)

  }
  if (show[2L]) {

    matplot.edgwas(x, s = "rho.min")

  }

  invisible()
}


