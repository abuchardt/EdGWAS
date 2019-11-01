#' Plot diagnostics for a cv.edgwas object
#'
#' Four plots (selectable by which) are currently available: a plot of the cross-validation curve produced by cv.edgwas i.e., the mean cross-validated error and upper and lower standard deviation curves plotted against the values of rho used in the fits, and three types of adjacency matrix plot; static plots for each value of rho used in the fits, interactive plots for the values of rho used in the fits, and a single plot for the optimal value of rho.
#'
#' @param x Fitted "cv.edgwas" object.
#'
#' @details Plots are produced, and nothing is returned.
#' @export plot.cv.edgwas
#'
plot.cv.edgwas <- function (x, which = c(1L:5L),
                            main = "", xlab = NULL, ylab = NULL,
                            ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                            pointsCol = NULL, ylim = NULL,
                            ...) {

    if (!inherits(x, "cv.edgwas"))
      stop("use only with \"cv.edgwas\" objects")
    if(!is.numeric(which) || any(which < 1) || any(which > 5))
      stop("'which' must be in 1:5")
    show <- rep(FALSE, 5)
    show[which] <- TRUE

    rho <- x$rho

    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    ##---------- Do the individual plots : ----------
    if (show[1L]) {
      if(x$edgwas.fit$logrho == 1) {
        xvals <- log(x$rho)
        if(is.null(xlab))
          xlab <- bquote(log(rho))
        rho1se <- log(x$rho.1se)
        rhomin <- log(x$rho.min)
      } else {
        xvals <- x$rho
        if(is.null(xlab))
          xlab <- bquote(rho)
        rho1se <- x$rho.1se
        rhomin <- x$rho.min
      }
      if(is.null(ylim)) ylim <- range(x$cvup, x$cvlo) else ylim <- ylim
      dev.hold()
      plot(x = xvals, y = x$cvm, main = main,
                   ylim = ylim,
                   xlab = xlab, ylab = ifelse(is.null(ylab), x$name, ylab),
                   type = "n", bty = "n", ...)
      abline(v = rho1se, lty = 2)
      abline(v = rhomin, lty = 2)
      segments(x0 = xvals, x1 = xvals,
               y0 = x$cvlo, y1 = x$cvup,
               col = "grey", ...)
      points(xvals, x$cvm, pch = 20,
             col = ifelse(is.null(pointsCol), "#D95F02", pointsCol), ...)

      dev.flush()
    }
    if (show[2L]) {

      for(i in seq_along(rho)) {
        dev.hold()
        A <- x$edgwas.fit$A[[i]]
        colnames(A) <- rownames(A) <- paste0("V", seq(ncol(A)))

        # Reorder the correlation matrix
        cormat <- A #reorderA(A)
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

      data_heatmap_1 <- lapply(seq_along(x$edgwas.fit$A), function(a) {
        AA <- x$edgwas.fit$A[[a]]
        colnames(AA) <- rownames(AA) <- paste0("V", seq(ncol(AA)))
        cormat <- AA #reorderA(AA, x$edgwas.fit$P[[a]])
        cbind(melt(cormat), rho = x$rho[a])
      })

      data_heatmap <- do.call(rbind, data_heatmap_1)

        dev.hold()
        gp <- ggplot(data_heatmap, aes(Var1, Var2, frame = rho)) +
          geom_tile(aes(fill = value)) +
          scale_fill_gradientn(colors = c("#ffffff", "#16161d")) +
          scale_y_discrete(limits = rev(unique(data_heatmap$Var2))) +
          scale_x_discrete(position = 'top') +
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
        dev.flush()
    }
    if (show[4L]) {

      AA <- x$edgwas.fit$A[[which(x$rho == x$rho.1se)]]
      q <- ncol(AA)

      dev.hold()
      plot(c(0,q), c(0,q), type = "n", ylim = c(0, q+2),
           xlab="", ylab="", xaxt='n', yaxt='n', bty = "n", asp = 1)
      # create the matrix
      for(i in seq(q)) {
        for(j in seq(q)) {
          rect(i-1,q+1-j,i,q-j,
               col = ifelse(AA[i,j] == 0,
                            "white", "gray"),
               border = "white")
        }
      }
      text(x = 0, y = (1:q) - .5, labels = paste0(q:1), pos = 2)
      text(x = (1:q) - .5, y = q, labels = paste0(1:q), pos = 3)
      mtext(text = "rho.1se", side = 3)
      dev.flush()
    }
    if (show[5L]) {

      AA <- x$edgwas.fit$A[[which(x$rho == x$rho.min)]]
      q <- ncol(AA)

      dev.hold()
      plot(c(0,q), c(0,q), type = "n", ylim = c(0, q+2),
           xlab="", ylab="", xaxt='n', yaxt='n', bty = "n", asp = 1)
      # create the matrix
      for(i in seq(q)) {
        for(j in seq(q)) {
          rect(i-1,q+1-j,i,q-j,
               col = ifelse(AA[i,j] == 0,
                            "white", "gray"),
               border = "white")
        }
      }
      text(x = 0, y = (1:q) - .5, labels = paste0(q:1), pos = 2)
      text(x = (1:q) - .5, y = q, labels = paste0(1:q), pos = 3)
      mtext(text = "rho.min", side = 3)
      dev.flush()
    }
    #if (!one.fig && par("oma")[3L] >= 1)
    #  mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
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

  dev.hold()
  gp <- ggplot(data_heatmap, aes(Var1, Var2, frame = rho)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colors = c("#ffffff", "#16161d")) +
    scale_y_discrete(limits = unique(data_heatmap$Var2)) +
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

  p
}

reorderA <- function(A, P = NULL){
  # Use correlation between variables as distance
  if(is.null(P)) dd <- as.dist(A)
  else dd <- as.dist(P)
  hc <- hclust(dd)
  A <- A[hc$order, hc$order]
  A
}
