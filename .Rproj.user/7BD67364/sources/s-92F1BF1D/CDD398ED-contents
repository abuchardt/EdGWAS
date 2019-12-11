#' Cluster multiple traits via polygenic scores
#'
#' This function clusters traits that share some genetic component via polygenic scores (PSs). It fits a sparse precision matrix via graphical lasso. The regularisation path is computed for the lasso penalty at a grid of values for the regularisation parameter rho.
#'
#' ...
#'
#' @param x Input matrix, of dimension nobs x nvars or nobs x nouts; each row is an observation vector. A matrix of PSs if \code{scores = TRUE} (default) and nvars = nouts. Can be in sparse matrix format.
#' @param y Quantitative response matrix, of dimension nobs x nouts.
#' @param scores Are PSs provided (default is TRUE) or should PSs be generated from x variable (FALSE).
#' @param rho (Non-negative) regularisation parameter for lasso passed to glasso. \code{rho=0} means no regularisation. Can be a scalar (usual) or a symmetric nouts by nouts matrix, or a vector of length nouts. In the latter case, the penalty matrix has jkth element sqrt(rho[j]*rho[k]).
#' @param nrho The number of rho values - default is 40.
#' @param logrho Logical flag for log transformation of the rho sequence. Default is \code{logrho = FALSE}.
#' @param rho.min.ratio Smallest value for rho, as a fraction of rho.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero) - default is 10e-04.
#'
#' @return An object of class "edgwas" is returned. \item{call}{The call that produced this object.} \item{alpha}{A matrix of intercepts of dimension nouts x length(rho)} \item{beta}{A matrix of coefficients for the PSs of dimension nouts x length(rho)} \item{A}{A length(rho) list of estimated adjacency matrices A of 0s and 1s, where A_{ij} is equal to 1 iff edges i and j are adjacent and A_{ii} is 0.} \item{P}{A length(rho) list of estimated precision matrices (matrix inverse of correlation matrices).} \item{Sigma}{A length(rho) list of estimated correlation matrices.} \item{rho}{The actual sequence of rho values used.} \item{PS}{Polygenic scores used. If  \code{scores = FALSE} they are computed by \code{\link{ps.edgwas}}} \item{logrho}{Logical flag for log transformation of the rho sequence. Default is \code{logrho = FALSE}.}
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
#' pc <- edgwas(x, y, scores = FALSE)
#'
#' @export
#'

edgwas <- function(x, y, scores = TRUE, rho = NULL,
                   nrho = ifelse(is.null(rho), 40, length(rho)),
                   logrho = FALSE,
                   rho.min.ratio = 10e-04
                   ) {

  this.call <- match.call()

  x <- as.matrix(x)
  y <- as.matrix(y)

  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  nouts <- ifelse(is.null(dimy), 1, dimy[2])
  if (nrowy != nobs)
    stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (",
               nobs, ")", sep = ""))
  if (isTRUE(scores) && nvars != nouts)
    stop(paste("number of scores (", nvars, ") not equal to number of outcome components (", nouts, ")", sep = ""))
  vnames <- colnames(x)
  if (is.null(vnames))
    vnames <- paste("V", seq(nvars), sep = "")

  if (is.null(rho)) {
    if (rho.min.ratio >= 1)
      stop("rho.min.ratio should be less than 1")
    flmin <- as.double(rho.min.ratio)
    urho <- double(1)
    #rholist <- rho
    #nrho <- 10
  } else {
    nrho = as.integer(nrho)
    flmin <- as.double(1)
    if (any(rho < 0))
      stop("rhos should be non-negative")
    urho <- as.double(rev(sort(rho)))
    nrho <- as.integer(length(rho))
    rholist <- as.double(rho)
  }


  if (isFALSE(scores)) {
    PS0 <- ps.edgwas(x, y)$PS
    PS <- scale(PS0, center = TRUE, scale = TRUE)
  } else {
    PS = scale(x, center = TRUE, scale = TRUE)
  }


  ####################
  # STEP 2
  ####################
  # Covariance matrix (qxq) of PSs
  SigmaPS <- cov(PS)

  # Calculate rho path (first get rho_max):
  if (is.null(rho)) {
    rho_max <- max(colSums(SigmaPS)) #max(abs(colSums(SigmaPS)))
    if(logrho == 0) {
      rholist <- round(seq(rho_max*flmin, rho_max,
                           length.out = nrho), digits = 10)
    } else if (logrho == 1) {
      rholist <- round(exp(seq(log(rho_max*flmin), log(rho_max),
                               length.out = nrho)), digits = 10)
    } else if (logrho == 2) {
      rholist <- round(log(seq(exp(rho_max*flmin), exp(rho_max),
                               length.out = nrho)), digits = 10)
    }
  } else {
    rholist <- rho
  }


  # Graphical lasso
  # Estimates a sparse inverse covariance matrix using a lasso (L1) penalty
  #glPS <- glasso::glassopath(s = SigmaPS, rho = rholist,
  #                            penalize.diagonal = FALSE, trace = 0)
  glPS2 <- lapply(seq_along(rholist), FUN = function(r) {
    glasso::glasso(s = SigmaPS, rho = rholist[r],
                   penalize.diagonal = FALSE,
                   trace = 0)
  })

  rho <- rholist


  A <- list(NULL)

  P <- lapply(seq_along(rho), function(r) glPS2[[r]]$wi)
  Sigma <- lapply(seq_along(rho), function(r) glPS2[[r]]$w)
  A <- lapply(seq_along(rho), function(r) {
    AA <- P[[r]]
    AA[abs(AA) > 1.5e-8] <- 1
    diag(AA) <- 0
    AA
  })

  #beta <- vector(mode = "list", length = nrho)
  #beta <- lapply(beta, FUN = function(l) matrix(NA, 2, nrho))
  alpha <- matrix(NA, nouts, nrho)
  beta <- matrix(NA, nouts, nrho)
  for (j in seq(nrho)) {

    # Rotate Y and PSs (to obtain independent Y)
    w <- expm::sqrtm(P[[j]]) ## qxq
    yIn <- y %*% w
    xIn <- PS %*% w

    for (l in seq(ncol(y))) {

      Sigma12 <- Sigma[[j]][l, -l, drop = FALSE] # 1 x (q-1)
      Sigma21 <- Sigma[[j]][-l, l, drop = FALSE] # (q-1) x 1
      Sigma22I <- solve(Sigma[[j]][-l,-l]) # (q-1) x (q-1)#

      # Update PSs
      xUp <- drop(Sigma12 %*% tcrossprod(Sigma22I, xIn[, -l])) + xIn[, l]

      fit <- lm(yIn[,l] ~ xUp)
      alpha[l, j] <- coef(fit)[1]
      beta[l, j] <- coef(fit)[2] #matrix(coef(fit), ncol = 1)
    }

  }

  # Return
  fit <- list(call = this.call,
              alpha = alpha, beta = beta,
              A = A, P = P, Sigma = Sigma,
              rho = rho, PS = PS, logrho = logrho)
  class(fit) <- "edgwas"
  fit
}
