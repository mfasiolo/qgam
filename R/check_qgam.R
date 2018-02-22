##########################
#' Some diagnostics for a fitted qgam model
#' 
#' @description Takes a fitted gam object produced by \code{qgam()} and produces some diagnostic information 
#'              about the fitting procedure and results. It is partially based on \code{mgcv::gam.check}.
#'  
#' @param obj the output of a \code{qgam()} call. 
#' @param nbin number of bins used in the internal call to \code{cqcheck()}. 
#' @param lev the significance levels used by \code{cqcheck()}, which determines the width of the confidence 
#'            intervals.
#' @param ... extra arguments to be passed to \code{plot()}
#' @return Simply produces some plots and prints out some diagnostics.
#' @details This function provides two plots. The first shows how the number of responses falling below the fitted
#'          quantile (y-axis) changes with the fitted quantile (x-axis). To be clear: if the quantile is fixed to, say, 0.5
#'          we expect 50\% of the responses to fall below the fit. See \code{?cqcheck()} for details. The second plot related
#'          to \code{|F(hat(mu)) - F(mu0)|}, which is the absolute bias attributable to the fact that qgam is using 
#'          a smoothed version of the pinball-loss. The absolute bias is evaluated at each observation, and an histogram
#'          is produced. See Fasiolo et al. (2017) for details. The function also prints out the average absolute bias,
#'          and the proportion of observations lying below the regression line. It also provides some convergence 
#'          diagnostics (regarding the optimization), which are the same as in \code{mgcv::gam.check}. 
#'          It reports also the maximum (k') and the selected degrees of freedom of each smooth term.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>, Simon N. Wood. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2017). Fast calibrated additive quantile regression. 
#'             Available at \url{https://arxiv.org/abs/1707.03307}.
#' @examples
#' library(qgam)
#' set.seed(0)
#' dat <- gamSim(1, n=200)
#' b<-qgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, qu = 0.5)
#' plot(b, pages=1)
#' check.qgam(b, pch=19, cex=.3)
#'
check.qgam <- function(obj,
                       nbin = 10,
                       lev = 0.05,
                       ...)
  ## takes a fitted gam object and produces some standard diagnostic plots
{
  
  svpar <- par(no.readonly = TRUE) 
  cqcheck(obj = obj, v = as.matrix(obj$fitted.values)[ , 1], nbin = nbin, xlab = "Fitted values", ...)
  par(svpar) 
  
  ## Checking bias induced by having smoothed the loss
  # Here we are estimating E( Phi(y, mu, lam*sig) - I(y > mu) | x ) using a Gaussian GAM
  co <- obj$family$getCo()
  sig <- exp( obj$family$getTheta() )
  if( is.list(obj$formula) ) { sig <- sig * obj$fitted.values[ , 2] }
  lam <- co / sig
  
  dat <- obj$model
  form <- if( is.list(obj$formula) ) { obj$formula[[1]] } else { obj$formula }
  res <- dat[[ form[[2]] ]] - as.matrix(obj$fitted.values)[ , 1]
  form[[2]] <- as.symbol( "bias" )
  
  dat$bias <- plogis(res, 0, sig*lam) - as.numeric(res > 0)
  
  fitBias <- gam(form, data = dat) 
  hist(fitBias$fitted.values, xlab = expression(F(hat(mu)) - F(mu[0])), main = "Bias due to smoothed loss")
  
  cat("Theor. proportion of neg. resid.:", obj$family$getQu(), "  Actual proportion:", mean(res<0))
  cat("\nMean absolute bias |F(mu) - F(mu0)| =", mean(abs(fitBias$fitted.values)))

  ## now summarize convergence information
  cat("\nMethod:",obj$method,"  Optimizer:",obj$optimizer)
  if (!is.null(obj$outer.info)) { ## summarize convergence information
    if (obj$optimizer[2]%in%c("newton","bfgs"))
    { boi <- obj$outer.info
    cat("\n",boi$conv," after ",boi$iter," iteration",sep="")
    if (boi$iter==1) cat(".") else cat("s.")
    cat("\nGradient range [",min(boi$grad),",",max(boi$grad),"]",sep="")
    cat("\n(score ",obj$gcv.ubre," & scale ",obj$sig2,").",sep="")
    ev <- eigen(boi$hess)$values
    if (min(ev)>0) cat("\nHessian positive definite, ") else cat("\n")
    cat("eigenvalue range [",min(ev),",",max(ev),"].\n",sep="")
    } else { ## just default print of information ..
      cat("\n"); print(obj$outer.info)
    }
  } else { ## no sp, perf iter or AM case
    if (length(obj$sp)==0) ## no sp's estimated
      cat("\nModel required no smoothing parameter selection")
    else {
      cat("\nSmoothing parameter selection converged after",obj$mgcv.conv$iter,"iteration")
      if (obj$mgcv.conv$iter>1) cat("s")

      if (!obj$mgcv.conv$fully.converged)
        cat(" by steepest\ndescent step failure.\n") else cat(".\n")
      cat("The RMS",obj$method,"score gradient at convergence was",obj$mgcv.conv$rms.grad,".\n")
      if (obj$mgcv.conv$hess.pos.def)
        cat("The Hessian was positive definite.\n") else cat("The Hessian was not positive definite.\n")
      #cat("The estimated model rank was ",obj$mgcv.conv$rank,
      #           " (maximum possible: ",obj$mgcv.conv$full.rank,")\n",sep="")
    }
  }
  if (!is.null(obj$rank)) {
    cat("Model rank = ",obj$rank,"/",length(obj$coefficients),"\n")
  }

  cat("\n")
  ## now check k
  kchck <- .kcheck(obj)
  if (!is.null(kchck)) {
    cat("Basis dimension (k) check: if edf is close too k\' (maximum possible edf) \n")
    cat("it might be worth increasing k. \n\n")
    printCoefmat(kchck,digits=3);
  }
  
  return( invisible(NULL) )
  
}