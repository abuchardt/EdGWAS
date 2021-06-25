#' Generate polygenic scores
#'
#' This function generates polygenic scores (PSs) by fitting a univariate simple linear regression model for each feature x on each outcome component y.
#'
#' For each outcome component \eqn{Y_l} we fit a univariate simple linear regression on the form \deqn{Y_l = X_j B_{jl} + E_l,} where the scalar \eqn{B_{jl}} is a regression coefficient \eqn{E_l} is is a vector of length nobs of independent Gaussian random errors. For a multivariate outcome \eqn{Y} we define the PS for each outcome component l = 1,...,q and each individual i = 1,...,N by \deqn{PS_{il} = \sum_{j=1}^p X_{ij}\hat{B}_{jl},} where \eqn{\hat{B}_{jl}} is the maximum likelihood estimate of \eqn{B_{jl}}.
#'
#' @param x Input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.
#' @param y Quantitative response matrix, of dimension nobs x nouts.
#'
#' @return An object of class "ps.edgwas" is returned. \item{PS}{A matrix of dimension nobs x nouts of polygenic scores.}
#'
#' @examples
#' N <- 500 #
#' q <- 10 #
#' p <- 20 #
#' set.seed(1)
#' x <- matrix(rbinom(n = N*p, size = 2, prob = 0.3), nrow=N, ncol=p)
#' B <- matrix(0, nrow = p, ncol = q)
#' B[1, 1:2] <- 2
#' y <- x %*% B + matrix(rnorm(N*q), nrow = N, ncol = q)
#' ###
#' ps <- ps.edgwas(x, y)
#'
#' @export
#'

ps.edgwas <- function(x, y) {

  # Compute q simple GWASs
  fit <- MESS::plr(y = y, x = x)
  betaList <- lapply(fit, FUN = function(X) X[, 1])
  beta <- do.call(cbind, betaList)

  # Compute Nxq PSs
  PS <- x %*% beta

  # Return
  obj <- list(PS = PS, beta = beta)
  class(obj) <- "ps.edgwas"
  obj
}
