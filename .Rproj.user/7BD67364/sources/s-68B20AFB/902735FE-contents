#' Cross-validation for edgwas
#'
#' Does k-fold cross-validation for edgwas, produces a plot, and returns a value for rho.
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.
#' @param y response matrix, of dimension nobs x nouts. Quantitative for family="gaussian". For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For "binomial" if y is presented as a vector, it will be coerced into a factor.
#' @param rho (Non-negative) optional user-supplied rho sequence; default is NULL, and EdGwas chooses its own sequence.
#' @param nfolds number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3.
#' @param type.measure loss to use for cross-validation. Currently one option; the default is type.measure="mse", which uses the mean-squared error.
#' @param ... Other arguments that can be passed to edgwas.
#'
#' @details The function runs edgwas nfolds+1 times; the first to get the rho sequence, and then the remainder to compute the fit with each of the folds omitted. The error is accumulated, and the average error and standard deviation over the folds is computed. Note that the results of cv.edgwas are random, since the folds are selected at random. Users can reduce this randomness by running cv.edgwas many times, and averaging the error curves.
#'
#' @return Cluster associations. \item{clust}{returns a vector with group memberships}#'
#'
#' @examples
#' # Gaussian
#' N <- 100 #
#' q <- 9 #
#' p <- 1000 #
#' set.seed(1)
#' X <- matrix(sample(0:2, N*p, replace=TRUE), nrow=N, ncol=p)
#' B <- matrix(0, nrow = p, ncol = q)
#' B[1, 1:2] <- 1
#' Y <- X %*% B + matrix(rnorm(N*q), nrow = N, ncol = q)
#' ###
#' pc <- cv.edgwas(x = X, y = Y, nfolds = 5)
#'
#' @export cv.edgwas
#'
cv.edgwas <- function(x, y, rho = NULL, nfolds = 10, type.measure = "mse",
                      penalty = c("ridge", "lasso"), trace = NULL, logrho = TRUE, ...) {

  if (missing(type.measure)) {
    type.measure <- "default"
  } else type.measure <- match.arg(type.measure)

  if (!is.null(rho) && length(rho) < 2)
    stop("Need more than one value of rho for cv.edgwas")

  edgwas.call <- match.call(expand.dots = TRUE)
  edgwas.call[[1]] <- as.name("edgwas")
  edgwas.object <- edgwas(x, y, rho = rho, logrho = logrho, ...)
  edgwas.object$call <- edgwas.call

  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")

  rho <- edgwas.object$rho

  fun <- paste("cv", type.measure, sep = ".")
  cvstuff <- do.call(fun, list(rho, x, y, nfolds, penalty, type.measure, trace, logrho))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- names(cvstuff$type.measure)

  minCrit <- cvm
  rhoMin <- rho[which.min(minCrit)]
  idx <- max(which(rev(minCrit[seq(which.min(minCrit))]) < minCrit[which.min(minCrit)] + cvsd[which.min(minCrit)]))
  rho1se <- rev(rho[seq(which.min(minCrit))])[idx]

  out <- list(rho = rho, cvm = cvm, cvsd = cvsd, cvsd2 = cvstuff$cvsd2,
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
cv.default <- function(rho, x, y, nfolds, penalty, type.measure, trace, logrho, ...) {

  nrho <- length(rho)
  foldid <- sample(rep(seq(nfolds), length = nrow(x)))

  #mse <- list(NULL)
  predmat <- vector(mode = "list", length = nrho)
  predmat <- lapply(predmat, FUN = function(l) matrix(NA, nrow(y), ncol(y)))
  for (i in seq(nfolds)) {

    if(!is.null(trace)) cat("i: ", i, ", ")

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

    outlist <- edgwas(xTrain, yTrain, rho = NULL, nrho = nrho, logrho)#, ...)

    mse[[i]] <- vector("numeric", length(rho))
    for (j in seq(length(rho))) {

      if(!is.null(trace)) cat(".")

      w <- expm::sqrtm(outlist$P[[j]]) ## qxq
      wy <- yTrain %*% w ## nrow(xTrain)xq
      xVex <- matrix(rep(c(xTrain), ncol(w)), ncol = ncol(w)) ## nrow(xTrain)*pxq
      wx <- xVex %*% w ## nrow(xTrain)*pxq

      #fit <- list(NULL)
      for (l in seq(ncol(y))) {
        xsubsub <- matrix(wx[,l], ncol = ncol(x)) ## 900*10000

        if (penalty == "lasso") { # LASSO
          lasFit <- glmnet::cv.glmnet(x = xsubsub, y = c(wy[, l]))
          nz <- which(coef(lasFit, s="lambda.min")[-1] != 0)
          if (length(nz) < 1) {
            preds <- matrix(mean(wy[,l]), nrow = nrow(yTest))
          } else {
            trainData <- data.frame(y = wy[,l], x = I(xsubsub[,nz]))
            lmFit <- lm(y ~ x, data = trainData)
            preds <- cbind(1, xTest[,nz]) %*% coef(lmFit)
          }
        } else if (penalty == "ridge") { # RIDGE
          ridgeFit <- glmnet::cv.glmnet(x = xsubsub, y = c(wy[, l]), alpha = 0)
          preds <- cbind(1, xTest) %*% matrix(coef(ridgeFit), ncol = 1)
        }
        predmat[[j]][fold, l] <- preds

      }

      #yHat <- do.call(cbind, fit)
      #mse[[i]][j] <- mean((yTest - yHat)^2, na.rm = TRUE)
    }

    if(!is.null(trace)) cat("\n")
  }

  N <- nrow(y)
  #cvraw = switch(type.measure, mse = (y - predmat)^2,
  #               deviance = (y - predmat)^2, mae = abs(y - predmat))
  cvraw <- lapply(predmat, FUN = function(l) (y - l)^2)

  #cvraw <- do.call(rbind, mse)
  cvm <- sapply(cvraw, mean, na.rm = TRUE)
  cvsd <- sapply(cvraw, sd, na.rm = TRUE)
  #cvsd <- sqrt(sapply(scale(cvraw, center = cvm, scale = FALSE)^2, 2, mean, na.rm = TRUE)/(nfolds - 1))
  #cvsd <- sqrt(sapply(cvraw, FUN = function(l) mean((l - cvm)^2, na.rm = TRUE))/(N-1))

  #cvm <- rowMeans(cvraw)
  #cvsd1 <- sqrt(apply((cvraw - cvm)^2, 2, sum, na.rm = TRUE)/(nfolds - 1))
  #cvsd2 <- apply(cvraw, 2,  sd, na.rm = TRUE)

  names(type.measure) <- "Mean-Squared Error"


  list(cvm = cvm, cvsd = cvsd, cvsd2 = cvsd2, type.measure = type.measure)

}


