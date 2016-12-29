##########################
#' Fit multiple smooth additive quantile regression models
#' 
#' @description This function fits a smooth additive regression model to several quantiles.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
#' @param qu A vectors of quantiles of interest. Each entry should be in (0, 1).
#' @param lsig The value of the log learning rate used to create the Gibbs posterior. By defauls \code{lsig=NULL} and this
#'             parameter is estimated by posterior calibration described in Fasiolo et al. (2016). Obviously, the function is much faster
#'             if the user provides a value. 
#' @param err An upper bound on the error of the estimated quantile curve. Should be in (0, 1). If it is a vector, it should be of the 
#'            same length of \code{qu}. See Fasiolo et al. (2016) for details.
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
#' @return A list with entries: \itemize{
#'                   \item{\code{fit} = a \code{gamObject}, one for each entry of \code{qu}.  Notice that the
#'                                      slots \code{model} and \code{smooth} of each object has been removed to save memory. 
#'                                      See \code{?gamObject}. }
#'                   \item{\code{model} = the \code{model} slot of the \code{gamObject}s in the \code{fit} slot. This is the same for every
#'                                        fit, hence only one copy is stored.}
#'                   \item{\code{smooth} = the \code{smooth} slot of the \code{gamObject}s in the \code{fit} slot. This is the same for every
#'                                        fit, hence only one copy is stored.}
#'                   \item{\code{calibr} = a list which is the output of an internal call to \code{tuneLearnFast}, which is used for calibrating
#'                                         the learning rate. See \code{?tuneLearnFast} for details.}
#' }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2016). Fast calibrated additive quantile regression. Available at
#'             \url{https://github.com/mfasiolo/qgam/blob/master/draft_qgam.pdf}.
#' @examples
#' #####
#' # Univariate "car" example
#' ####
#' library(qgam); library(MASS)
#' 
#' # Fit for quantile 0.8 using the best sigma
#' quSeq <- c(0.2, 0.4, 0.6, 0.8)
#' set.seed(6436)
#' fit <- mqgam(accel~s(times, k=20, bs="ad"), data = mcycle, err = 0.05, qu = quSeq)
#' 
#' # Plot the fit
#' xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
#' for(iq in quSeq){
#'   pred <- qdo(fit, iq, predict, newdata = xSeq)
#'   lines(xSeq$times, pred, col = 2)
#' }
#' 
#' #####
#' # Multivariate Gaussian example
#' ####
#' library(qgam)
#' set.seed(2)
#' dat <- gamSim(1, n=400, dist="normal", scale=2)
#' 
#' fit <- mqgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, err = 0.05, qu = c(0.2, 0.8), 
#'              control = list("tol" = 0.01)) # <- semi-sloppy tolerance to speed-up calibration
#' 
#' invisible( qdo(fit, 0.2, plot, pages = 1) )
#' @export mqgam  
#'
mqgam <- function(form, data, qu, lsig = NULL, err = 0.05, 
                  multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                  control = list(), argGam = NULL)
{
  nq <- length(qu)
  
  if( length(err) != nq ){
    if(length(err) == 1) { 
      err <- rep(err, nq) 
    } else {
      stop("\"err\" should either be a scalar or a vector of the same length as \"qu\".")
    }
  }
  
  # Setting up control parameter (mostly used by tuneLearnFast)
  ctrl <- list( "gausFit" = NULL, "verbose" = FALSE, "b" = 0)
  
  # Checking if the control list contains unknown names entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control, verbose = FALSE)
  
  # Initial Gaussian fit
  if( is.null(ctrl[["gausFit"]]) )
  {
    if( is.formula(form) ){
      gausFit <- do.call("gam", c(list("formula" = form, "data" = data), argGam))
    } else {
      gausFit <- do.call("gam", c(list("formula" = form, "data" = data, "family" = gaulss(b=ctrl[["b"]])), argGam))
    }
    ctrl[["gausFit"]] <- gausFit
  }
  
  # Output list
  out <- list()
  
  if( is.null(lsig) ) { # Selecting the learning rate sigma OR ....
    learn <- tuneLearnFast(form = form, data = data, err = err, qu = qu,
                           multicore = multicore, cluster = cluster, ncores = ncores, paropts = paropts,
                           control = ctrl, argGam = argGam)
    lsig <- learn$lsig
    out[["calibr"]] <- learn
  } else { # ... use the one provided by the user
    if( length(lsig) == 1 ) {
      lsig <- rep(lsig, nq)
    } else {
      if( length(lsig) != nq ) stop("lsig should either be scalar or a vector of length(qu) ")
    } }
  
  # Fitting a quantile model for each qu
  out[["fit"]] <- lapply(1:nq, function(ii){
    
    .out <- qgam(form, data, qu[ii], lsig = lsig[ii], err = err[ii], multicore = FALSE, control = ctrl, argGam = argGam)
    
    # Removing data and smooth matrix to reduce memory requirements. There quantities
    # are kept only inside the first fit ( qfit[[1]] )
    if(ii > 1){
      .out$model  <- NULL
      .out$smooth <- NULL 
    } 
    
    return( .out )
  })
  
  # Storing output list
  names(out[["fit"]]) <- qu
  out[["model"]] <- out[["fit"]][[1]][["model"]]
  out[["smooth"]] <- out[["fit"]][[1]][["smooth"]]
  out[["fit"]][[1]][["model"]] <- NULL
  out[["fit"]][[1]][["smooth"]] <- NULL
  
#   out[["qu"]] <- qu
#   out[["lambda"]] <- lam
#   out[["lsig"]] <- lsig
  
  return( out )
}

