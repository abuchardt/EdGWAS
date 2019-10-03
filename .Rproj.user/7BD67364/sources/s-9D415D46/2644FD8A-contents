#' Cluster multiple traits via PRS
#'
#' This function clusters traits that share some genetic component via polygenic risk scores (PRSs). It fits a sparse precision matrix via graphical lasso. The regularisation path is computed for the lasso penalty at a grid of values for the regularisation parameter rho.
#'
#' ...
#'
#' @param x Input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.
#' @param y Quantitative response matrix, of dimension nobs x nouts.
#' @param rho (Non-negative) regularisation parameter for lasso passed to glasso. rho=0 means no regularisation. Can be a scalar (usual) or a symmetric nouts by nouts matrix, or a vector of length nouts. In the latter case, the penalty matrix has jkth element sqrt(rho[j]*rho[k]).
#'
#' @return Cluster associations. \item{A}{A length(rho) list of matrices A of 0s and 1s, where A_{ij} is equal to 1 iff edges i and j are adjacent and A_{ii} is 0.} \item{rho}{The actual sequence of rho values used.} \item{P}{A length(rho) list of matrices P, the estimated precision matrix (matrix inverse of correlation matrix).}
#'
#'
#' @examples
#' # Gaussian
#' N <- 1000
#' q <- 9
#' p <- 10000
#' set.seed(1)
#' X <- matrix(sample(0:2, N*p, replace=TRUE), nrow=N, ncol=p)
#' B <- matrix(0, nrow = p, ncol = q)
#' B[1, 1:2] <- 1
#' Y <- X %*% B + matrix(rnorm(N*q), nrow = N, ncol = q)
#' ###
#' pc <- edgwas(x = X, y = Y, rho = 0.1)
#'
#' @export edgwas
#'

edgwas <- function(x, y, rho = NULL, nrho = ifelse(is.null(rho), 20, length(rho)),
                   rho.min.ratio = 10e-04
                   #thresh = 1e-07,
                   #dfmax = nvars + 1, pmax = min(dfmax * 2 + 20, nvars),
                   #exclude, penalty.factor = rep(1, nvars),
                   #lower.limits = -Inf, upper.limits = Inf,
                   #maxit = 1e+05,
                   #standardize.response = FALSE,
                   ) {

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


  #####################
  # STEP 1
  #####################
  # Compute q simple GWASs
  gwasResults <- lapply(1:nouts, FUN = function(k) MESS::mfastLmCpp(y = y[,k], x = x)$coefficients)
  # Save all pxq regression coefficients
  betaHatMat <- do.call(cbind, gwasResults)

  # Compute Nxq PRSs
  PRS <- mapply(1:nouts, FUN = function(k) x %*% betaHatMat[,k])

  ####################
  # STEP 2
  ####################
  # Covariance matrix (qxq) of PRSs
  SigmaPRS <- cor(PRS) # cov(PRS) #

  # Calculate rho path (first get rho_max):
  if (is.null(rho)) {
    rho_max <- mean(colSums(SigmaPRS)) #max(abs(colSums(SigmaPRS)))
    rholist <- round(exp(seq(log(rho_max*flmin), log(rho_max),
                             length.out = nrho)), digits = 10)
    #rholist <- round(seq(rho_max*flmin, rho_max,
    #                     length.out = nrho), digits = 10)
  } else {
    rholist <- rho
  }


  # Graphical lasso
  # Estimates a sparse inverse covariance matrix using a lasso (L1) penalty

  glPRS <- glasso::glassopath(s = SigmaPRS, rho = rholist, trace = FALSE)
  rho <- glPRS$rho

  A <- list(NULL)
  P <- list(NULL)

  P <- lapply(seq_along(rho), function(rhok) glPRS$wi[,,rhok])
  A <- lapply(seq_along(rho), function(rhok) {
    AA <- P[[rhok]]
    AA[abs(AA) > 1.5e-8] <- 1
    diag(AA) <- 0
    AA
  })

  # Return
  fit <- list(A = A, rho = rho, P = P) #clust = outcl, dist = d)
  class(fit) <- "edgwas"
  fit
}
