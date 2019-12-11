#' EdGwas:
#'
#' A package implementing methods for clustering traits that share some genetic component via polygenic scores (PSs). It fits a sparse precision matrix of the PSs via graphical lasso. The regularisation path is computed for the lasso penalty rho at a grid of values for the regularisation parameter rho. A cross-validation method for tuning the penalty rho is provided.
#'
#' The EdGwas package provides three categories of important functions:
#' simple fit, cross-validation and plots.
#'
#' @section fit:
#' The function \code{\link{edgwas}} clusters outcome components (traits) that share some feature (genetic component) via polygenic scores (PSs). It fits a sparse precision matrix via graphical lasso. The regularisation path is computed for the lasso penalty at a grid of values for the regularisation parameter rho.
#'
#' @section cross-validation:
#' The cross-validation function \code{\link{cv.edgwas}} does k-fold cross-validation for edgwas and returns a value for the penalty rho.
#'
#' @section bootstrap:
#' The bootstrap function ...
#'
#' @docType package
#' @name EdGwas-package
#' @author Ann-Sophie Buchardt \email{asbu@@sund.ku.dk}\cr Maintainer: Ann-Sophie Buchardt
#' \email{asbu@@sund.ku.dk}
#' @references glasso, glmnet
NULL

#' @importFrom stats cor cov lm coef predict as.dist hclust
#' @importFrom MESS mfastLmCpp
#' @importFrom expm sqrtm
#' @importFrom plotly layout animation_slider ggplotly
#' @importFrom graphics rect par plot points segments text
#' @importFrom grDevices dev.interactive devAskNewPage
NULL

#' @import glasso
#' @import glmnet
#' @import Matrix
#' @import reshape2
#' @import ggplot2
NULL
