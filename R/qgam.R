##########################
#' Fit a smooth additive quantile regression model
#' 
#' @description This function fits a smooth additive regression model for a single quantile.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
#' @param qu The quantile of interest. Should be in (0, 1).
#' @param lsig The value of the log learning rate used to create the Gibbs posterior. By defauls \code{lsig=NULL} and this
#'             parameter is estimated by posterior calibration described in Fasiolo et al. (2016). Obviously, the function is much faster
#'             if the user provides a value. 
#' @param err An upper bound on the error of the estimated quantile curve. Should be in (0, 1). See Fasiolo et al. (2016) for details.
#' @param multicore If TRUE the calibration will happen in parallel.
#' @param ncores Number of cores used. Relevant if \code{multicore == TRUE}.
#' @param cluster An object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster.
#' @param paropts a list of additional options passed into the foreach function when parallel computation is enabled. 
#'                This is important if (for example) your code relies on external data or packages: 
#'                use the .export and .packages arguments to supply them so that all cluster nodes 
#'                have the correct environment set up for computing. 
#' @param control A list of control parameters to be used by \code{tuneLearnFast}. See \code{?tuneLearnFast} for details.
#' @param argGam A list of parameters to be passed to \code{mgcv::gam}. This list can potentially include all the arguments listed
#'               in \code{?gam}, with the exception of \code{formula}, \code{family} and \code{data}.
#' @param ... additional arguments passed to \code{mgcv::gam}.
#' @return A \code{gamObject}. See \code{?gamObject}.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2016). Fast calibrated additive quantile regression. Available at
#'             \url{https://github.com/mfasiolo/qgam/draft_qgam.pdf}.
#' @examples
#
#' #####
#' # Univariate "car" example
#' ####
#' library(qgam); library(MASS)
#' 
#' # Fit for quantile 0.8 using the best sigma
#' set.seed(6436)
#' fit <- qgam(accel~s(times, k=20, bs="ad"), data = mcycle, err = 0.05, qu = 0.8)
#' 
#' # Plot the fit
#' xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
#' pred <- predict(fit, newdata = xSeq, se=TRUE)
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
#' lines(xSeq$times, pred$fit, lwd = 1)
#' lines(xSeq$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
#' lines(xSeq$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)   
#' 
#' #####
#' # Multivariate Gaussian example
#' ####
#' library(qgam)
#' set.seed(2)
#' dat <- gamSim(1,n=400,dist="normal",scale=2)
#' 
#' fit <- qgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, err = 0.05, qu = 0.8)
#' plot(fit, scale = F, pages = 1)                          
#' @export qgam  
#'
qgam <- function(form, data, qu, lsig = NULL, err = 0.01, 
                 multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                 control = list(), argGam = NULL)
{
  if( length(qu) > 1 ) stop("length(qu) > 1, so you should use mqgam()")
  
  # Setting up control parameter (mostly used by tuneLearnFast)
  ctrl <- list( "gausFit" = NULL, "verbose" = FALSE, "b" = 0)
  
  # Checking if the control list contains unknown names entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control, verbose = FALSE)
  
  # Gaussian fit, used for initializations 
  if( is.formula(form) ) {
    fam <- "logF"
    if( is.null(ctrl[["gausFit"]]) ) { ctrl$gausFit <- do.call("gam", c(list("formula" = form, "data" = data), argGam)) }
    varHat <- ctrl$gausFit$sig2
  } else {
    fam <- "logFlss"
    if( is.null(ctrl[["gausFit"]]) ) { ctrl$gausFit <- do.call("gam", c(list("formula" = form, "data" = data, "family" = gaulss(b=ctrl[["b"]])), argGam)) }
    varHat <- 1/ctrl$gausFit$fit[ , 2]^2
  }  # Start = NULL in gamlss because it's not to clear how to deal with model for sigma 
  
  # Selecting the learning rate sigma
  if( is.null(lsig) ) {  
    learn <- tuneLearnFast(form = form, data = data, err = err, qu = qu, ncores = ncores, control = ctrl, argGam = argGam)
    lsig <- learn$lsig
  }
  
  # Fit model
  lam <- err * sqrt(2*pi*varHat) / (2*log(2)*exp(lsig))
  fit <- do.call("gam", c(list("formula" = form, "family" = get(fam)(qu = qu, lam = lam, theta = lsig), "data" = data), argGam))
  
  return( fit )
}


