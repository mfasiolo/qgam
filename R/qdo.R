##########################
#' Manipulating the output of \code{mqgam}
#' 
#' @description Contrary to \code{qgam}, \code{mqgam} does not output a standard \code{gamObject}, hence
#'              methods such as \code{predict.gam} or \code{plot.gam} cannot be used directly. \code{qdo}
#'              provides a simple wrapper for such methods.
#'  
#' @param obj A list which is the output of a \code{mqgam} call. 
#' @param qu A scalar in (0, 1) representing the quantile of interest, which should be an element of \code{names(obj$fit)}.
#' @param fun The method or function that we want to use on the \code{gamObject} corresponding to quantile \code{qu}. For instance
#'            \code{predict}, \code{plot} or \code{summary}.
#' @param ... Additional arguments to be passed to \code{fun}.
#' @return The output of \code{fun}, whatever that is.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' library(qgam); library(MASS)
#' 
#' quSeq <- c(0.4, 0.6)
#' set.seed(737)
#' fit <- mqgam(accel~s(times, k=20, bs="ad"), data = mcycle, err = 0.05, qu = quSeq, 
#'              control = list("tol" = 0.01)) # <- semi-sloppy tolerance to speed-up calibration 
#' 
#' qdo(fit, 0.4, summary)
#' invisible(qdo(fit, 0.4, plot, pages = 1))
#' @export qdo
#'
qdo <- function(obj, qu, fun, ...){
  
  if( !(qu %in% names(obj[["fit"]])) ) stop("qu is not in obj[[\"qu\"]].")
  
  tmpObj <- obj[["fit"]][[ which(names(obj[["fit"]]) == qu) ]]
  
  tmpObj[["model"]] <- obj[["model"]]
  tmpObj[["smooth"]] <- obj[["smooth"]]
  
  out <- fun(tmpObj, ...)
  
  return( out )
}
