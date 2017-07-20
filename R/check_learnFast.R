##########################
#' Visual checks for the output of tuneLearnFast()
#' 
#' @description Provides some visual checks to verify whether the Brent optimizer used by \code{tuneLearnFast()} worked correctly.
#' @param obj the output of a call to \code{tuneLearnFast}. 
#' @param sel integer vector determining which of the plots will be produced. For instance if \code{sel = c(1, 3)} only
#'            the 1st and 3rd plots are showed. No entry of \code{sel} can be bigger than the number of quantiles considered
#'            in the original \code{tuneLearnFast()} call. That is, if estimated the learning rate for \code{qu = c(0.1, 0.4)},
#'            then \code{max(sel)} must be <= 3.
#' @param ... currently not used, here only for compatibility reasons.
#' @return It produces several plots. 
#' @details The top plot in the first page shows the bracket used to estimate log(sigma) for each quantile.
#'          The brackets are delimited by the crosses and the red dots are the estimates. If a dot falls very close to one of the crosses, 
#'          that might indicate problems. The bottom plot shows, for each quantile, the value of parameter \code{err} used. Sometimes the algorithm
#'          needs to increase \code{err} above its user-defined value to achieve convergence. Subsequent plots show, for each quantile, the value
#'          of the loss function corresponding to each value of log(sigma) explored by Brent algorithm.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2017). Fast calibrated additive quantile regression. 
#'             Available at \url{https://github.com/mfasiolo/qgam/blob/master/draft_qgam.pdf}.
#' @examples
#' library(qgam)
#' set.seed(525)
#' dat <- gamSim(1, n=200)
#' b <- tuneLearnFast(y ~ s(x0)+s(x1)+s(x2)+s(x3), 
#'                    data = dat, qu = c(0.4, 0.5), 
#'                    control = list("tol" = 0.05)) # <- sloppy tolerance to speed-up calibration 
#' check(b) 
#' check(b, 3) # Produces only third plot
#'
check.learnFast <- function(obj, sel = NULL, ...)
{  
  est <- obj$store
  brac <- obj$ranges
  lsig <- obj$lsig
  errors <- obj$err
  qu <- as.numeric(names(obj$lsig))
  nq <- length(qu)
  sel <- if(is.null(sel)){ 1:(nq+1) } else { sort(sel) }
  
  oldPar <- par(no.readonly = TRUE)
  if( 1%in%sel ){
    layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(2, 1))
    par(mai = c(1, 1, 0.1, 0.1))
    plot(qu, lsig, ylim = range(as.vector(brac)), xlim = range(qu)+c(-1e-5,+1e-5), col = 2, 
         ylab = expression("Log(" * sigma * ")"), xlab = "Quantile")
    points(qu, brac[ , 1], pch = 3)
    points(qu, brac[ , 2], pch = 3)
    points(qu, rowMeans(brac), pch = 3)
    for(zz in 1:nq) segments(qu[zz], mean(brac[zz, ]) - abs(diff(brac[zz, ]))/4, 
                             qu[zz], mean(brac[zz, ]) + abs(diff(brac[zz, ]))/4, col = 1)
    plot(qu, errors, xlab = "Quantile")
    par(oldPar)
  }
  
  if(any(sel > 1))
  {
    selQ <- sel[sel>1] - 1
    # readline(prompt = "Press <Enter> to see the next plot...")
    pDim <- min( ceiling(sqrt(length(selQ))), 2 )
    par(mfrow = c(pDim, pDim))
    for( ii in selQ )
    {
      plot(sort(est[[ii]][1, ]), est[[ii]][2, order(est[[ii]][1, ])], 
           main = substitute(Quantile == x, list(x = round(qu[ii], 3))), 
           ylab = "loss", xlab = expression(log(sigma)), type = 'b')
      abline(v = est[[ii]][1, which.min(est[[ii]][2, ])], col = 2)
      #if((ii %% (pDim^2) == 0) && (ii!=nq)) readline(prompt = "Press <Enter> to see the next plot...")
    }
  }
  par(oldPar)
  
  return( invisible(NULL) )
}
