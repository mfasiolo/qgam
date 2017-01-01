##########################
#' Generic checking function
#' 
#' @description Generic function for checking R objects which produces, for instance, convergence tests or diagnostic plots. 
#'              For \code{qgam} objects \code{check.qgam()} will be used. 
#' @param obj the object to be checked. 
#' @param ... extra arguments, mainly used by graphic functions. 
#' @return Reports the results of convergence tests and/or produces diagnostic plots.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' #######
#' # Using check.qgam
#' #######
#' library(qgam)
#' set.seed(0)
#' dat <- gamSim(1, n=200)
#' b<-qgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, qu = 0.5)
#' plot(b, pages=1)
#' check(b, pch=19, cex=.3) 
#' @exportMethod check
#' @docType methods
#'
check <- function (obj, ...) {
  UseMethod("check", obj)
}