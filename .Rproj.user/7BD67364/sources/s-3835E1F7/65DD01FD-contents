#' Bootstrap for sensitivity analysis of EdGwas
#'
#' Does bootstrap for EdGwas, produces a plot, and returns a value for rho.
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.
#' @param y response matrix, of dimension nobs x nouts. Quantitative for family="gaussian". For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For "binomial" if y is presented as a vector, it will be coerced into a factor.
#' @param rho (Non-negative) optional user-supplied rho sequence; default is NULL, and EdGwas chooses its own sequence.
#' @param nfolds number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3.
#' @param type.measure loss to use for cross-validation. Currently one option; the default is type.measure="rand", which uses the rand index.
#'
#' @return Cluster associations. \item{clust}{returns a vector with group memberships}
#' @export boot.edgwas
#'
boot.edgwas <- function(x, y, rho = NULL, nfolds = 10, type.measure = "rand",...) {

  if (missing(type.measure)) {
    type.measure <- "default"
  } else type.measure <- match.arg(type.measure)

  if (!is.null(rho) && length(rho) < 2)
    stop("Need more than one value of rho for boot.edgwas")
  N <- nrow(x)

  edgwas.call <- match.call(expand.dots = TRUE)
  edgwas.call[[1]] <- as.name("edgwas")
  edgwas.object <- edgwas(x, y, rho = rho, ...)
  edgwas.object$call <- edgwas.call

  #nz <- sapply(predict(edgwas.object, type = "nonzero"),
  #             length)

  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")

  foldid <- sample(rep(seq(nfolds), length = N))
  outlist <- as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    which <- foldid == i
    if (is.matrix(y))
      y_sub <- y[!which, ]
    else y_sub <- y[!which]
    outlist[[i]] <- edgwas(x[!which, , drop = FALSE],
                           y_sub, rho = rho, ...)
  }

  fun <- paste("boot", type.measure, sep = ".")
  rho <- edgwas.object$rho
  bootstuff <- do.call(fun, list(outlist, rho, x, y,
                               foldid, type.measure))
  bootRI <- bootstuff$bootRI
  bootup <- bootstuff$bootup
  bootlo <- bootstuff$bootlo
  #bootsd <- bootstuff$bootsd
  #nas <- is.na(bootsd)
  #if (any(nas)) {
  #  rho <- rho[!nas]
  #  bootRI <- bootRI[!nas]
  #  bootsd <- bootsd[!nas]
  #  nz <- nz[!nas]
  #}
  bootname <- names(bootstuff$type.measure)
  rhoOpt <- rho[which.max(bootRI)]

  out <- list(rho = rho, bootRI = bootRI, #bootsd = bootsd,
              bootup = bootup, #bootRI + bootsd,
              bootlo = bootlo, #bootRI - bootsd,
              name = bootname, edgwas.fit = edgwas.object,
              rhoOpt = rhoOpt)
  class(obj) <- "boot.edgwas"
  obj
}


# Cross-validation with Rand index
boot.rand <- function(outlist, rho, x, y, foldid, type.measure) {

  RI <- list(NULL)
  nfolds <- nlevels(as.factor(foldid))

  for(i in seq(nfolds)) {
    which <- foldid == i
    if (is.matrix(y))
      y_sub <- y[which, ]
    else y_sub <- y[which]

    testClus <- edgwas(x = x[which, , drop = FALSE], y = y_sub, rho = rho)
    adjMat <- testClus$A

    RI[[i]] <- sapply(seq_along(adjMat), FUN = function(j) randIndex(outlist[[i]]$A[[j]], adjMat[[j]])$RI)
  }

  meanRI <- colMeans(do.call(rbind, RI))
  lRI <- apply(do.call(rbind, RI), 2, quantile, prob = .025)
  uRI <- apply(do.call(rbind, RI), 2, quantile, prob = .975)

  names(type.measure) <- "Rand index"

  list(bootRI = meanRI, bootup = uRI, bootlo = lRI, type.measure = type.measure)

}

# Rand index for adjecency matrix
randIndex <- function(clus1, clus2) {
  a <- sum(clus1[upper.tri(clus1)] == 0 & clus2[upper.tri(clus2)] == 0, na.rm = TRUE)
  b <- sum(clus1[upper.tri(clus1)] == 1 & clus2[upper.tri(clus2)] == 1, na.rm = TRUE)
  c <- sum(clus1[upper.tri(clus1)] == 0 & clus2[upper.tri(clus2)] == 1, na.rm = TRUE)
  d <- sum(clus1[upper.tri(clus1)] == 1 & clus2[upper.tri(clus2)] == 0, na.rm = TRUE)

  RI <- (a + b)/(a + b + c + d)

  list(RI = RI)
}
