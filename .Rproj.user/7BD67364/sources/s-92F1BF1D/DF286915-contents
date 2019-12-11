#' Plot diagnostics for a cv.edgwas object
#'
#' Four plots (selectable by which) are currently available: a plot of the cross-validation curve produced by cv.edgwas i.e., the mean cross-validated error and upper and lower standard deviation curves plotted against the values of rho used in the fits, and three types of adjacency matrix plot; static plots for each value of rho used in the fits, interactive plots for the values of rho used in the fits, and a single plot for the optimal value of rho.
#'
#' @param x Fitted "cv.edgwas" object.
#' @param which If a subset of the plots is required, specify a subset of the numbers \code{1:5}.
#' @param zoom A non-negative integer determining the number of points closest to rho.min to show in the plot of the cross-validation curve (\code{which = 1}). Default is \code{zoom = 0L} for plotting the full curve.
#' @param main An overall title for the plot: see \code{\link{title}}
#' @param xlab A title for the x axis: see \code{\link{title}}
#' @param ylab A title for the y axis: see \code{\link{title}}
#' @param ask Logical; if \code{TRUE}, the user is \emph{ask}ed before each plot, see \code{\link{par}(ask=.)}.
#' @param pointsCol The color for points corresponding to rho.min and rho.1se. Default is pointsCol = "#D95F02".
#' @param ylim The y limits of the plot.
#' @param xlim The x limits of the plot.
#' @param ... Other paramters to be passed through to plotting functions.
#'
#' @details Plots are produced, and nothing is returned.
#'
#' @seealso{\code{\link{edgwas}} and \code{\link{cv.edgwas}}}
#'
#' @examples
#' N <- 500 #
#' q <- 10 #
#' p <- 20 #
#' set.seed(1)
#' x <- matrix(sample(0:2, N*p, replace=TRUE), nrow=N, ncol=p)
#' B <- matrix(0, nrow = p, ncol = q)
#' B[1, 1:2] <- 10
#' y <- x %*% B + matrix(rnorm(N*q), nrow = N, ncol = q)
#' ###
#' pc <- cv.edgwas(x, y, scores = FALSE)
#' plot(pc, 1)
#' plot(pc, 2)
#'
#' @export
#'

plot.cv.edgwas <- function (x, which = c(1L:4L), zoom = 0L,
                            main = "", xlab = NULL, ylab = NULL,
                            ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                            pointsCol = "#D95F02", ylim = NULL, xlim = NULL,
                            ...) {

    if (!inherits(x, "cv.edgwas"))
      stop("use only with \"cv.edgwas\" objects")
    if(!is.numeric(which) || any(which < 1) || any(which > 4))
      stop("'which' must be in 1:4")
    show <- rep(FALSE, 4)
    show[which] <- TRUE

    rho <- x$rho

    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    ##---------- Do the individual plots : ----------
    if (show[1L]) {

      if(zoom > 0 && is.numeric(zoom)) {
        idx <- (which(x$rho == x$rho.min) - zoom) : (which(x$rho == x$rho.min) + zoom)
        idx <- idx[idx %in% seq_along(x$rho)]

      } else if (zoom == 0L) {
        idx <- seq_along(x$rho)
      } else {
        stop("zoom should be a non-negative integer")
      }

      if(is.null(xlim)) xlim <- range(x$rho[idx]) else xlim = xlim

      if(x$edgwas.fit$logrho == 1) {
        xvals <- log(x$rho[idx])
        xlim <- log(xlim)
        if(is.null(xlab))
          xlab <- bquote(log(rho))
      } else {
        xvals <- x$rho[idx]
        if(is.null(xlab))
          xlab <- bquote(rho)
      }

      if(is.null(ylim)) ylim <- range(x$cvup[idx], x$cvlo[idx]) else ylim <- ylim

      ptCol <- rep("grey40", length(x$rho[idx]))
      ptCol[which(x$rho[idx] %in% c(x$rho.min, x$rho.1se))] <- pointsCol

      plot(x = xvals, y = x$cvm[idx], main = main,
                   ylim = ylim, xlim = xlim,
                   xlab = xlab, ylab = ifelse(is.null(ylab), x$name, ylab),
                   type = "n", bty = "n", ...)
      segments(x0 = xvals, x1 = xvals,
               y0 = x$cvlo[idx], y1 = x$cvup[idx],
               col = "grey")
      points(xvals, x$cvm[idx], pch = 20,
             col = ptCol)

    }
    if (show[2L]) {

      matplot.edgwas(x, s = "rho.min")

    }
    if (show[3L]) {

      matplot.edgwas(x, s = "rho.1se")

    }
    if (show[4L]) {

      .interactiveHeatmap(x)

    }
    invisible()
  }


.interactiveHeatmap <- function(x) {

  data_heatmap_1 <- lapply(seq_along(x$edgwas.fit$A), function(a) {
    AA <- x$edgwas.fit$A[[a]]
    colnames(AA) <- rownames(AA) <- paste0("V", seq(ncol(AA)))
    cormat <- AA #reorderA(AA, x$edgwas.fit$P[[a]])
    cbind(melt(cormat), rho = x$rho[a])
  })

  data_heatmap <- do.call(rbind, data_heatmap_1)


  gp <- ggplot(data_heatmap, aes(data_heatmap$Var1,
                                 data_heatmap$Var2,
                                 frame = data_heatmap$rho)) +
    geom_tile(aes(fill = data_heatmap$value)) +
    scale_fill_gradientn(colors = c("#ffffff", "#16161d")) +
    scale_y_discrete(limits = rev(unique(data_heatmap$Var2))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks = element_blank()) +
    geom_vline(xintercept=seq(1.5, 10-0.5, 1),
               lwd=.5, colour="white") +
    geom_hline(yintercept=seq(1.5, 10-0.5, 1),
               lwd=.5, colour="white") +
    coord_fixed() + theme(legend.position = "none")

  p <- animation_slider(ggplotly(gp),
                        currentvalue = list(prefix = "rho ", font = list(color="black")))

  print(p)

}

reorderA <- function(A, P = NULL){
  # Use correlation between variables as distance
  if(is.null(P)) dd <- stats::as.dist(A)
  else dd <- stats::as.dist(P)
  hc <- stats::hclust(dd)
  A <- A[hc$order, hc$order]
  A
}
