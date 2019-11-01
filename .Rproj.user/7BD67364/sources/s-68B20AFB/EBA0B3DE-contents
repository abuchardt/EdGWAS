# Cross-validation
cv.precma <- function(fit, rho, x, y, foldid, type.measure) {

  foldid <- sample(rep(seq(nfolds), length = N))

  mse <- list(NULL)
  yHat <- list(NULL)
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



  }

  cvm <- rowMeans(do.call(cbind, mse))
  cvsd <- apply(do.call(cbind, mse), 1,  sd)

  names(type.measure) <- "Precision matrix approximation"

  list(mse = mse, cvm = cvm, cvsd = cvsd, type.measure = type.measure)

}


