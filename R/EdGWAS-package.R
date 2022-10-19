#' EdGWAS:
#'
#' A package implementing methods for clustering traits that share some genetic component via polygenic scores (PSs). It fits a sparse precision matrix of the PSs via graphical lasso. The regularisation path is computed for the lasso penalty rho at a grid of values for the regularisation parameter rho. A cross-validation method for tuning the penalty rho is provided.
#'
#' The EdGWAS package provides three categories of important functions:
#' simple fit, cross-validation and plots.
#'
#' @section fit:
#' The function \code{\link{edgwas}} clusters outcome components (traits) that share some feature (genetic component) via polygenic scores (PSs). It fits a sparse precision matrix via graphical lasso. The regularisation path is computed for the lasso penalty at a grid of values for the regularisation parameter rho. The clusters are used to specify the structure of the error covariance matrix of a GLS model and the feasible GLS estimator is used for estimating the unknown parameters in a linear regression model with a certain unknown degree of correlation between the residuals.
#'
#' @docType package
#' @name EdGWAS-package
#' @author Ann-Sophie Buchardt \email{asbu@@sund.ku.dk}\cr Maintainer: Ann-Sophie Buchardt
#' \email{asbu@@sund.ku.dk}
#' @references glasso, glmnet
NULL

#' @importFrom stats cor cov lm coef predict as.dist hclust sd resid var
#' @importFrom MESS mfastLmCpp plr
#' @importFrom expm sqrtm
#' @importFrom igraph graph.adjacency get.data.frame components
#' @importFrom lme4 lmer
#' @importFrom plotly layout animation_slider ggplotly
#' @importFrom graphics rect par plot points segments text abline
#' @importFrom grDevices dev.interactive devAskNewPage
NULL

#' @import glasso
#' @import Matrix
#' @import reshape2
NULL
