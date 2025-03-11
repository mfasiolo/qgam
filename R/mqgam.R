##########################
#' Fit multiple smooth additive quantile regression models
#' 
#' @description This function fits a smooth additive regression model to several quantiles.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
#' @param qu A vectors of quantiles of interest. Each entry should be in (0, 1).
#' @param discrete If TRUE then covariate discretisation is used for faster model fitting. See \code{mgcv::}\link[mgcv]{bam} for details.
#' @param lsig The value of the log learning rate used to create the Gibbs posterior. By defauls \code{lsig=NULL} and this
#'             parameter is estimated by posterior calibration described in Fasiolo et al. (2017). Obviously, the function is much faster
#'             if the user provides a value. 
#' @param err An upper bound on the error of the estimated quantile curve. Should be in (0, 1). If it is a vector, it should be of the 
#'            same length of \code{qu}. Since qgam v1.3 it is selected automatically, using the methods of Fasiolo et al. (2017).
#'            The old default was \code{err=0.05}.
#' @param multicore If TRUE the calibration will happen in parallel.
#' @param ncores Number of cores used. Relevant if \code{multicore == TRUE}.
#' @param cluster An object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster.
#' @param paropts a list of additional options passed into the foreach function when parallel computation is enabled. 
#'                This is important if (for example) your code relies on external data or packages: 
#'                use the .export and .packages arguments to supply them so that all cluster nodes 
#'                have the correct environment set up for computing. 
#' @param control A list of control parameters. The only one relevant here is \code{link}, which is the link function
#'                used (see \code{?elf} and \code{?elflss} for defaults). All other control parameters are used by 
#'                \code{tuneLearnFast}. See \code{?tuneLearnFast} for details.
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
#' @references Fasiolo, M., Wood, S.N., Zaffran, M., Nedellec, R. and Goude, Y., 2020. 
#'             Fast calibrated additive quantile regression. 
#'             Journal of the American Statistical Association (to appear).
#'             \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1725521}.
#' @references Fasiolo, M., Wood, S.N., Zaffran, M., Nedellec, R. and Goude, Y., 2021. 
#'             qgam: Bayesian Nonparametric Quantile Regression Modeling in R. 
#'             Journal of Statistical Software, 100(9), 1-31, \doi{10.18637/jss.v100.i09}.
#' @examples
#' 
#' #####
#' # Multivariate Gaussian example
#' ####
#' library(qgam)
#' set.seed(2)
#' dat <- gamSim(1, n=300, dist="normal", scale=2)
#' 
#' fit <- mqgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, qu = c(0.2, 0.8))
#' 
#' invisible( qdo(fit, 0.2, plot, pages = 1) )
#' 
#' #####
#' # Univariate "car" example
#' ####
#' library(qgam); library(MASS)
#' 
#' # Fit for quantile 0.8 using the best sigma
#' quSeq <- c(0.2, 0.4, 0.6, 0.8)
#' set.seed(6436)
#' fit <- mqgam(accel~s(times, k=20, bs="ad"), data = mcycle, qu = quSeq)
#' 
#' # Plot the fit
#' xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
#' for(iq in quSeq){
#'   pred <- qdo(fit, iq, predict, newdata = xSeq)
#'   lines(xSeq$times, pred, col = 2)
#' }
#'
mqgam <- function(form, data, qu, discrete = FALSE, lsig = NULL, err = NULL, 
                  multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                  control = list(), argGam = NULL)
{
  nq <- length(qu)
  
  discrete <- .should_we_use_discrete(form = form, discrete = discrete)
  
  # Removing all NAs, unused variables and factor levels from data
  data <- .cleanData(.dat = data, .form = form, .drop = argGam$drop.unused.levels)
  
  if( !is.null(err) && length(err) != nq ){
    if(length(err) == 1) { 
      err <- rep(err, nq) 
    } else {
      stop("\"err\" should either be a scalar or a vector of the same length as \"qu\".")
    }
  }
  
  # Setting up control parameter (mostly used by tuneLearnFast)
  ctrl <- list("verbose" = FALSE, "b" = 0, "link" = "identity")
  
  # Checking if the control list contains unknown names entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control, verbose = FALSE)
  
  ctrl$init_qgam <- .init_gauss_fit(form = form, data = data, ctrl = ctrl, argGam = argGam, qu = qu, discrete = discrete)

  # Output list
  out <- list()
  
  if( is.null(lsig) ) { # Selecting the learning rate sigma OR ....
    learn <- tuneLearnFast(form = form, data = data, err = err, qu = qu, discrete = discrete,
                           multicore = multicore, cluster = cluster, ncores = ncores, paropts = paropts,
                           control = ctrl, argGam = argGam)
    lsig <- learn$lsig
    err <- learn$err # Over-writing err parameters!
    out[["calibr"]] <- learn
  } else { # ... use the one provided by the user
    if( length(lsig) == 1 ) {
      lsig <- rep(lsig, nq)
    } else {
      if( length(lsig) != nq ) stop("lsig should either be scalar or a vector of length(qu) ")
    } }
  
  # Fitting a quantile model for each qu
  out[["fit"]] <- lapply(1:nq, function(ii){
    
    if( !is.null(out$calibr) ){
      # Annoyingly, initial coeffs are supplied via "coef" argument in bam() and "start" in gam()
      argGam[[ ifelse(discrete, "coef", "start") ]] <- learn$final_fit[[ii]]$coefstart 
      argGam$mustart <- learn$final_fit[[ii]]$mustart
      argGam$in.out <- learn$final_fit[[ii]]$in.out
      ctrl$init_qgam$initM <- NA # initM should NOT be used by qgam, hence we want to get an error if it gets used.
    } else {
      ctrl$initM <- ctrl$initM[[ii]]
    }
      
    .out <- qgam(form, data, qu[ii], lsig = lsig[ii], err = err[ii], discrete = discrete, multicore = FALSE, control = ctrl, argGam = argGam)
    
    # Removing data and smooth matrix to reduce memory requirements. There quantities
    # are kept only inside the first fit ( qfit[[1]] )
    if(ii > 1){
      .out$model  <- NULL
      .out$smooth <- NULL 
      .out$call$data <- NULL
    } 
    
    return( .out )
  })
  
  # Storing output list
  names(out[["fit"]]) <- qu
  out[["model"]] <- out[["fit"]][[1]][["model"]]
  out[["smooth"]] <- out[["fit"]][[1]][["smooth"]]
  out[["data"]] <- out[["fit"]][[1]][["call"]][["data"]]
  out[["fit"]][[1]][["model"]] <- NULL
  out[["fit"]][[1]][["smooth"]] <- NULL
  out[["fit"]][[1]][["call"]][["data"]] <- NULL
  
  class(out) <- "mqgam"
  
  return( out )
}

