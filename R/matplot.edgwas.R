#' Adjacency matrix plot for a cv.edgwas object
#'
#' @param x Fitted "cv.edgwas" object.
#' @param s Value(s) of the penalty parameters rho for which the corresponding adjacency matrix is plotted. Default is rho.min.
#' @param A True adjacency matrix (optional).
#' @param col Color(s) to fill or shade the non-zero entries with. The default is \code{col = "lightgrey"}.
#' @param bty A character string which determine the type of box which is drawn about plots. Default is "n" which suppresses the box.
#' @param axes A logical value indicating whether axes should be drawn on the plot.
#' @param frame.plot A logical value indicating whether a box should be drawn around the plot.
#' @param ... Other paramters to be passed through to plotting functions.
#'
#' @details Plots are produced, and nothing is returned.
#'
#' @export matplot.edgwas
#'

matplot.edgwas <- function (x, s = c("rho.min", "rho.1se"), A = NULL,
                            col = "darkgrey", bty = "n", axes = TRUE, frame.plot = FALSE, ...) {

  if (is.numeric(s)) rho <- s
  else
    if (is.character(s)) {
      s <- match.arg(s)
      rho <- x[[s]]
    }
  else stop("Invalid form for s")

  for (i in 1:length(rho)) {
    rhoi <- rho[i]
    Ahat <- x$edgwas.fit$A[[which(x$rho == rhoi)]]

    q <- ncol(Ahat)
    plot(c(0,q), c(0,q), type = "n", ylim = c(0, q+2),
         xlab="", ylab="", xaxt='n', yaxt='n', bty = bty, asp = 1)
    # create the matrix
    for (i in seq(q)) {
      for (j in seq(q)) {
        graphics::rect(i-1,q+1-j,i,q-j,
             col = ifelse(Ahat[i,j] == 0, NA, col),
             density = ifelse(Ahat[i,j] == 0, 0, 50),
             border = ifelse(!is.null(A) && A[i,j] == 1, "#D95F02", NA))
      }
    }
    rect(-.5,-.5,q+.5,q+.5,
         border = ifelse(frame.plot, "darkgrey", NA))
    if (isTRUE(axes)) {
      text(x = 0, y = (1:q) - .5, labels = paste0(q:1), pos = 2)
      text(x = (1:q) - .5, y = q, labels = paste0(1:q), pos = 3)
    }
  }

  invisible()
}
