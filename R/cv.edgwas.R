#' Cross-validation for EdGwas
#'
#' Does k-fold cross-validation for EdGwas, produces a plot, and returns a value for rho.
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.
#' @param y response matrix, of dimension nobs x nouts. Quantitative for family="gaussian". For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For "binomial" if y is presented as a vector, it will be coerced into a factor.
#' @param rho (Non-negative) optional user-supplied rho sequence; default is NULL, and EdGwas chooses its own sequence.
#' @param nfolds number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3.
#' @param type.measure loss to use for cross-validation. Currently one option; the default is type.measure="mse", which uses the mean-squared error.
#' @param ... Other arguments that can be passed to edgwas.
#'
#' @return Cluster associations. \item{clust}{returns a vector with group memberships}#'
#'
#' @examples
#' # Gaussian
#' N <- 100 #
#' q <- 9
#' p <- 1000 #
#' set.seed(1)
#' X <- matrix(sample(0:2, N*p, replace=TRUE), nrow=N, ncol=p)
#' B <- matrix(0, nrow = p, ncol = q)
#' B[1:2, 1:5] <- 1
#' Y <- X %*% B + matrix(rnorm(N*q), nrow = N, ncol = q)
#' ###
#' pc <- cv.edgwas(x = X, y = Y, nfolds = 5)
#'
#' @export cv.edgwas
#'
cv.edgwas <- function(x, y, rho = NULL, nfolds = 10, type.measure = "mse", ...) {

  if (missing(type.measure)) {
    type.measure <- "default"
  } else type.measure <- match.arg(type.measure)

  if (!is.null(rho) && length(rho) < 2)
    stop("Need more than one value of rho for cv.edgwas")

  edgwas.call <- match.call(expand.dots = TRUE)
  edgwas.call[[1]] <- as.name("edgwas")
  edgwas.object <- edgwas(x, y, rho = rho)#, ...)
  edgwas.object$call <- edgwas.call

  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")

  lengthRho <- length(edgwas.object$rho)

  fun <- paste("cv", type.measure, sep = ".")
  rho <- edgwas.object$rho
  cvstuff <- do.call(fun, list(rho, x, y, lengthRho, nfolds, type.measure))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- names(cvstuff$type.measure)

  rhoMin <- rho[which.min(cvsd)]
  idx <- max(which(rev(cvm[seq(which.min(cvm))]) < cvm[which.min(cvm)] + cvsd[which.min(cvm)]))
  rho1se <- rev(rho[seq(which.min(cvm))])[idx]

  out <- list(rho = rho, cvm = cvm, cvsd = cvsd,
              cvup = cvm + cvsd,
              cvlo = cvm - cvsd,
              name = cvname,
              rho.min = rhoMin,
              rho.1se = rho1se,
              edgwas.fit = edgwas.object)
  class(out) <- "cv.edgwas"
  out
}


# Cross-validation
cv.default <- function(rho, x, y, lengthRho, nfolds, type.measure, ...) {

  foldid <- sample(rep(seq(nfolds), length = nrow(x)))

  mse <- list(NULL)
  for (i in seq(nfolds)) {

    cat("i: ", i, ", ")
    fold <- foldid == i

    if (is.matrix(y)) {
      yTrain <- y[!fold, ]
      yTest <- y[fold, ]
    } else {
      yTrain <- y[!fold]
      yTest <- y[fold]
    }
    xTrain <- x[!fold, , drop = FALSE]
    xTest <- x[fold, , drop = FALSE]

    outlist <- edgwas(xTrain, yTrain, rho = NULL)#, ...)

    mse[[i]] <- vector("numeric", lengthRho)
    for (j in seq(lengthRho)) {
      cat(".")

      w <- expm::sqrtm(outlist$P[[j]]) ## qxq
      wy <- yTrain %*% w ## nrow(xTrain)xq
      xVex <- matrix(rep(c(xTrain), ncol(w)), ncol = ncol(w)) ## nrow(xTrain)*pxq
      wx <- xVex %*% w ## nrow(xTrain)*pxq

      fit <- list(NULL)
      for (l in seq(ncol(y))) {
        xsubsub <- matrix(wx[,l], ncol = ncol(x)) ## 900*10000
        lasFit <- glmnet::cv.glmnet(x = xsubsub, y = c(wy[, l]))
        nz <- which(coef(lasFit, s="lambda.min")[-1] != 0)
        if (length(nz) < 1) {
          fit[[l]] <- matrix(mean(wy[,l]), nrow = nrow(yTest))
        } else {
          trainData <- data.frame(y = wy[,l], x = I(xsubsub[,nz]))
          lmFit <- lm(y ~ x, data = trainData)
          fit[[l]] <- cbind(1, xTest[,nz]) %*% coef(lmFit)
        }
      }

      yHat <- do.call(cbind, fit)
      mse[[i]][j] <- mean((yTest - yHat)^2, na.rm = TRUE)
    }
    cat("\n")
  }

  cvm <- rowMeans(do.call(cbind, mse))
  cvsd <- apply(do.call(cbind, mse), 1,  sd)

  names(type.measure) <- "Mean-Squared Error"


  list(cvm = cvm, cvsd = cvsd, type.measure = type.measure)

}


