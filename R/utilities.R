###########
# Utilities
##########

####### Check loss function
checkloss <- function(y, mu, qu){
  
  tau <- 1 - qu
  
  d <- y - mu
  
  l <- - sum( tau*d[d<0] ) - sum( (tau-1)*d[d>0] )
  
  return( l )
  
} 

#### Vettorize check Loss function
checklossVett <- function(y, mu, p){
  
  n <- length( p )
  
  out <- sapply(1:n,
                function(ii){
                  return( checkloss(y, mu[ , ii], p[ii]) )
                })
  
  return( out )
}


#### Vettorized empirical cdf
qqVett <- function(y, mu){
  
  nq <- ncol( mu )
  nobs <- length( y )
  
  out <- sapply(1:nq,
                function(ii){
                  return( sum( (y - mu[ , ii]) < 0 ) / nobs )
                })
  
  return( out )
}

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
#' @param ... Addinal arguments to be passed to \code{fun}.
#' @return The output of \code{fun}, whatever that is.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' library(qgam); library(MASS)
#' 
#' quSeq <- c(0.4, 0.6)
#' set.seed(737)
#' fit <- mqgam(accel~s(times, k=20, bs="ad"), data = mcycle, err = 0.01, qu = quSeq, 
#'              control = list("tol" = 0.01)) # <- semi-sloppy tolerance to speed-up calibration 
#' 
#' qdo(fit, 0.4, summary)
#' qdo(fit, 0.4, plot)
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



##########################
#' Visual checks for the output of \code{tuneLearnFast}
#' 
#' @description Provides so visual checks to verify whether the Brent optimizer used by \code{tuneLearnFast} worked correctly.
#'  
#' @param cal The output of a call to \code{tuneLearnFast}. This is also contained in the \code{calibr} slot of the output of
#'            a call to \code{mqgam}.
#' @return It provides several plots. The top plot in first page shows the bracket used to estimate log(sigma) for each quantile.
#'         The brackets are delimited by the crosses and the red dots are the estimates. If a dot falls very close to one of the crosses, 
#'         that might indicate problems. The bottom plot shows, for each quantile, the value of parameter \code{err} used. Sometimes the algorithm
#'         needs to increase \code{err} above its user-defined value to achieve convergence. Subsequent plots show, for each quantile, the value
#'         of the loss function corresponding to each value of log(sigma) explored by Brent algorithm.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' library(qgam); library(MASS)
#' 
#' quSeq <- c(0.4, 0.5, 0.6)
#' set.seed(737)
#' fit <- mqgam(accel~s(times, k=20, bs="ad"), data = mcycle, err = 0.01, qu = quSeq, 
#'              control = list("tol" = 0.01)) # <- semi-sloppy tolerance to speed-up calibration 
#' 
#' checkLearn(fit$calibr)
#' @export qdo
#'
checkLearn <- function(cal)
{  
  est <- cal$store
  brac <- cal$ranges
  lsig <- cal$lsig
  errors <- cal$err
  
  qu <- as.numeric(names(cal$lsig))
  nq <- length(qu)
  
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), 
         heights=c(2, 1))
  oldPar <- par(mai = c(1, 1, 0.1, 0.1))
  plot(qu, lsig, ylim = range(as.vector(brac)), xlim = range(qu)+c(-1e-5,+1e-5), col = 2, 
       ylab = expression("Log(" * sigma * ")"), xlab = "Quantile")
  points(qu, brac[ , 1], pch = 3)
  points(qu, brac[ , 2], pch = 3)
  points(qu, rowMeans(brac), pch = 3)
  for(zz in 1:nq) segments(qu[zz], mean(brac[zz, ]) - abs(diff(brac[zz, ]))/4, 
                           qu[zz], mean(brac[zz, ]) + abs(diff(brac[zz, ]))/4, col = 1)
  plot(qu, errors, xlab = "Quantile")
  
  readline(prompt = "Press <Enter> to see the next plot...")
  
  par(oldPar)
  
  pDim <- min( ceiling(sqrt(nq)), 2 )
  par(mfrow = c(pDim, pDim))
  for( ii in 1:nq )
  {
    plot(sort(est[[ii]][1, ]), est[[ii]][2, order(est[[ii]][1, ])], 
         main = substitute(Quantile == x, list(x = round(qu[ii], 3))), 
         ylab = "loss", xlab = expression(log(sigma)), type = 'b')
    abline(v = est[[ii]][1, which.min(est[[ii]][2, ])], col = 2)
    
    if(ii %% (pDim^2) == 0) readline(prompt = "Press <Enter> to see the next plot...")
  }
  
  par(oldPar)
  
  return( invisible(NULL) )
  
}