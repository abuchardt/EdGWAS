#' Plot diagnostics for a cv.edgwas object
#'
#' Three plots (selectable by which) are currently available: a plot of the cross-validation curve produced by cv.edgwas i.e., the mean cross-validated error and upper and lower standard deviation curves plotted against the values of rho used in the fits, and an adjacency matrix plot for each value of rho used in the fits (static or interactive).
#'
#' @param x Fitted "cv.edgwas" object.
#'
#' @details Plots are produced, and nothing is returned.
#' @export plot.cv.edgwas
#'
plot.cv.edgwas <- function (x, which = c(1L:3L),
                            main = "",
                            ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                            ...) {

    if (!inherits(x, "cv.edgwas"))
      stop("use only with \"cv.edgwas\" objects")
    if(!is.numeric(which) || any(which < 1) || any(which > 3))
      stop("'which' must be in 1:3")
    show <- rep(FALSE, 3)
    show[which] <- TRUE

    rho <- x$rho

    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    ##---------- Do the individual plots : ----------
    if (show[1L]) {
      if(all(abs(diff(x$rho)) == 1)) {
        xvals <- x$rho
        xlab <- bquote(log(rho))
      } else {
        xvals <- log(x$rho)
        xlab <- bquote(log(rho))
      }

      dev.hold()
      plot(x = xvals, y = x$cvm, main = main,
                   ylim = range(x$cvup, x$cvlo),
                   xlab = xlab, ylab = x$name,
                   type = "n", bty = "n", ...)
      segments(x0 = xvals, x1 = xvals,
               y0 = x$cvlo, y1 = x$cvup,
               col = "grey")
      points(xvals, x$cvm, pch = 20, col = "orangered")
      abline(v = log(x$rho.1se), lty = 3)
      abline(v = log(x$rho.min), lty = 3)

      if (one.fig)
        title(sub = sub.caption, ...)
      mtext(getCaption(1), 3, 0.25, cex = cex.caption)

      dev.flush()
    }
    if (show[2L]) {

      for(i in seq_along(rho)) {
        dev.hold()
        A <- x$edgwas.fit$A[[i]]
        colnames(A) <- rownames(A) <- paste0("V", seq(ncol(A)))

        reorderA <- function(A){
          # Use correlation between variables as distance
          dd <- as.dist(A)
          hc <- hclust(dd)
          A <-A[hc$order, hc$order]
        }
        # Reorder the correlation matrix
        cormat <- reorderA(A)
        # Melt the correlation matrix
        data_heatmap <- melt(cormat)

        gtitle1 <- bquote("Adjacency matrix for" ~ rho == .(round(rho[i], 3)) ~ "(rho.1se)")
        gtitle2 <- bquote("Adjacency matrix for" ~ rho == .(round(rho[i], 3)) ~ "(rho.min)")
        gtitle3 <- bquote("Adjacency matrix for" ~ rho == .(round(rho[i], 3)))
        gtitle <- if (rho[i] == x$rho.1se) gtitle1 else if
                     (rho[i] == x$rho.min) gtitle2 else gtitle3

        heatmap <- ggplot(data_heatmap, aes(Var1, Var2)) +
          geom_tile(aes(fill = value)) +
          scale_fill_gradientn(colors = c("#ffffff", "#16161d")) +
          scale_y_discrete(limits = rev(unique(data_heatmap$Var2))) +
          ggtitle(gtitle) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), #rect(linetype = "dashed", fill = NA),
                panel.background = element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.ticks = element_blank()) +
          geom_vline(xintercept=seq(1.5, 10-0.5, 1),
                     lwd=.5, colour="white") +
          geom_hline(yintercept=seq(1.5, 10-0.5, 1),
                     lwd=.5, colour="white") +
          coord_fixed() + theme(legend.position = "none")

        print(heatmap)

        dev.flush()
      }
    }
    if (show[3L]) {

      reorderA <- function(A, P){
        # Use correlation between variables as distance
        dd <- as.dist(P)
        hc <- hclust(dd)
        A <- A[hc$order, hc$order]
        A
      }

      data_heatmap_1 <- lapply(seq_along(cvfit$edgwas.fit$A), function(a) {
        AA <- cvfit$edgwas.fit$A[[a]]
        colnames(AA) <- rownames(AA) <- paste0("V", seq(ncol(AA)))
        cormat <- AA #reorderA(AA, cvfit$edgwas.fit$P[[a]])
        cbind(melt(cormat), rho = cvfit$rho[a])
      })

      data_heatmap <- do.call(rbind, data_heatmap_1)

        dev.hold()
        gp <- ggplot(data_heatmap, aes(Var1, Var2, frame = rho)) +
          geom_tile(aes(fill = value)) +
          scale_fill_gradientn(colors = c("#ffffff", "#16161d")) +
          scale_y_discrete(limits = unique(data_heatmap$Var2)) +
          #ggtitle(gtitle) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(), #rect(linetype = "dashed", fill = NA),
                panel.background = element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.ticks = element_blank()) +
          geom_vline(xintercept=seq(1.5, 10-0.5, 1),
                     lwd=.5, colour="white") +
          geom_hline(yintercept=seq(1.5, 10-0.5, 1),
                     lwd=.5, colour="white") +
          coord_fixed() + theme(legend.position = "none")

        p <- layout(animation_slider(ggplotly(gp),
                                     currentvalue = list(prefix = "rho ", font = list(color="black"))),
                    margin = list(pad = 0, b = 90, l = 0, r = 90))

        print(p)
        dev.flush()
    }
    #if (!one.fig && par("oma")[3L] >= 1)
    #  mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
    invisible()
  }
