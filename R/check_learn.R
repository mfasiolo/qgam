##########################
#' Visual checks for the output of tuneLearn()
#' 
#' @description Provides some visual plots showing how the calibration criterion and the effective degrees of 
#'              freedom of each smooth component vary with the learning rate.  
#'  
#' @param obj the output of a call to \code{tuneLearn}.
#' @param ... currently not used, here only for compatibility reasons.
#' @return It produces several plots. 
#' @details The first plot shows how the calibrations loss, which we are trying to minimize, varies with the 
#'          log learning rate. This function should look quite smooth, if it doesn't then try to increase
#'          \code{err} or \code{control$K} (the number of bootstrap samples) in the original call to 
#'          \code{tuneLearn}. The second plot shows how the effective degrees of freedom of each smooth term
#'          vary with log(sigma). Generally as log(sigma) increases the complexity of the fit decreases, hence
#'          the slope is negative.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2017). Fast calibrated additive quantile regression. 
#'             Available at \url{https://github.com/mfasiolo/qgam/blob/master/draft_qgam.pdf}.
#' @examples
#' library(qgam)
#' set.seed(525)
#' dat <- gamSim(1, n=200)
#' b <- tuneLearn(lsig = seq(-0.5, 1, length.out = 10), 
#'                y~s(x0)+s(x1)+s(x2)+s(x3), 
#'                data=dat, qu = 0.5)
#' check(b) 
#'
check.learn <- function(obj, ...)
{  
  sig <- as.numeric( names( obj$loss ) )
  
  # readline(prompt = "Press <Enter> to see the next plot...")
  plot(sig, obj$loss, type = "b", ylab = "Calibration Loss", xlab = expression("log(" * sigma * ")"))
  rug(sig[obj$convProb], side = 3, col = 2, lwd = 2)
  
  if( !is.null(obj$edf) )
  {
    # readline(prompt = "Press <Enter> to see the next plot...")
    nc <- ncol(obj$edf)
    matplot(obj$edf[ , 1], obj$edf[ , 2:nc], type = 'b', ylab = "Penalized EDF", xlab = expression("log(" * sigma * ")"), 
            pch = 1:nc, col = 1:nc)
    legend("topright", colnames(obj$edf)[2:nc], pch = 1:nc, col = 1:nc, bg="transparent")
    rug(sig[obj$convProb], side = 3, col = 2, lwd = 2)
  }
  
  return( invisible(NULL) )
}
