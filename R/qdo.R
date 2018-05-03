##########################
#' Manipulating the output of \code{mqgam}
#' 
#' @description Contrary to \code{qgam}, \code{mqgam} does not output a standard \code{gamObject}, hence
#'              methods such as \code{predict.gam} or \code{plot.gam} cannot be used directly. \code{qdo}
#'              provides a simple wrapper for such methods.
#'  
#' @param obj the output of a \code{mqgam} call. 
#' @param qu A vector whose elements must be in (0, 1). Each element indicates a quantile of interest, 
#'           which should be an element of \code{names(obj$fit)}. If left to \code{NULL} the function
#'           \code{fun} will be applied to each of the quantile fits in \code{obj}.
#' @param fun The method or function that we want to use on the \code{gamObject} corresponding to quantile \code{qu}. For instance
#'            \code{predict}, \code{plot} or \code{summary}. By default this is the identity function (\code{I}), which
#'            means that the fitted model for quantile \code{qu} is returned.
#' @param ... Additional arguments to be passed to \code{fun}.
#' @return A list where the i-th entry is the output of \code{fun} (whatever that is) corresponding to quantile \code{qu[i]}.
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
#' 
#' # Return the object for qu = 0.6 and then plot it
#' tmp <- qdo(fit, 0.6)
#' plot(tmp)
#'
qdo <- function(obj, qu=NULL, fun=I, ...){
  
  if( is.null(qu) ) { qu <- names(obj$fit) } 
  
  if( length(qu)>1 ){
    
   out <- lapply(qu, function(.q) qdo(obj, .q, fun, ...)) 
    
  } else {
  
  if( !(qu %in% names(obj[["fit"]])) ) stop("qu is not in obj[[\"qu\"]].")
  
  tmpObj <- obj[["fit"]][[ which(names(obj[["fit"]]) == qu) ]]
  
  tmpObj[["model"]] <- obj[["model"]]
  tmpObj[["smooth"]] <- obj[["smooth"]]
  tmpObj[["call"]][["data"]] <- obj[["data"]]
  
  out <- fun(tmpObj, ...)
  
  }
  
  return( out )
}
