#' @useDynLib SAM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom splines ns
#' @importFrom stats median
#' @importFrom graphics matplot
NULL


#' Sparse Additive Modelling
#'
#' The package SAM targets at high dimensional predictive modeling (regression and classification) for complex data analysis. SAM is short for sparse additive modeling, and adopts the computationally efficient basis spline technique. We solve the optimization problems by various computational algorithms including the block coordinate descent algorithm, fast iterative soft-thresholding algorithm, and newton method. The computation is further accelerated by warm-start and active-set tricks.
#'
#' \tabular{ll}{
#'   Package: \tab SAM\cr
#'   Type: \tab Package\cr
#'   Version: \tab 1.0.5\cr
#'   Date: \tab 2014-02-11\cr
#'   License: \tab GPL-2 \cr
#' }
#' @docType package
#' @author Tuo Zhao, Xingguo Li, Haoming Jiang, Han Liu, and Kathryn Roeder\cr
#' Maintainers: Haoming Jiang<hjiang98@gatech.edu>;
#' @references 
#' P. Ravikumar, J. Lafferty, H.Liu and L. Wasserman. "Sparse Additive Models", \emph{Journal of Royal Statistical Society: Series B}, 2009.\cr
#' T. Zhao and H.Liu. "Sparse Additive Machine", \emph{International Conference on Artificial Intelligence and Statistics}, 2012.\cr
#' @seealso \code{\link{samQL}},\code{\link{samHL}},\code{\link{samLL}},\code{\link{samEL}}
"_PACKAGE"
#> [1] "_PACKAGE"