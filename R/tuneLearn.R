##########################
#' Tuning the learning rate for Gibbs posterior
#' 
#' @description The learning rate (sigma) of the Gibbs posterior is tuned either by calibrating the credible intervals for the fitted
#'              curve, or by minimizing the pinball loss on out-of-sample data. This is done by bootrapping or by k-fold cross-validation. 
#'              Here the calibration loss function is evaluated on a grid of values provided by the user.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
#' @param lsig A vector of value of the log learning rate (log(sigma)) over which the calibration loss function is evaluated.
#' @param qu The quantile of interest. Should be in (0, 1).
#' @param err An upper bound on the error of the estimated quantile curve. Should be in (0, 1). 
#'            Since qgam v1.3 it is selected automatically, using the methods of Fasiolo et al. (2017).
#'            The old default was \code{err=0.05}.
#' @param multicore If TRUE the calibration will happen in parallel.
#' @param ncores Number of cores used. Relevant if \code{multicore == TRUE}.
#' @param cluster An object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster.
#' @param paropts a list of additional options passed into the foreach function when parallel computation is enabled. 
#'                This is important if (for example) your code relies on external data or packages: 
#'                use the .export and .packages arguments to supply them so that all cluster nodes 
#'                have the correct environment set up for computing. 
#' @param control A list of control parameters for \code{tuneLearn} with entries: \itemize{
#'                   \item{\code{loss} = loss function use to tune log(sigma). If \code{loss=="cal"} is chosen, then log(sigma) is chosen so that
#'                                       credible intervals for the fitted curve are calibrated. See Fasiolo et al. (2017) for details.
#'                                       If \code{loss=="pin"} then log(sigma) approximately minimizes the pinball loss on the out-of-sample
#'                                       data.}
#'                   \item{\code{sam} = sampling scheme use: \code{sam=="boot"} corresponds to bootstrapping and \code{sam=="kfold"} to k-fold
#'                                      cross-validation. The second option can be used only if \code{ctrl$loss=="pin"}.}
#'                   \item{\code{K} = if \code{sam=="boot"} this is the number of boostrap datasets, while if \code{sam=="kfold"} this is the 
#'                                    number of folds. By default \code{K=50}.}
#'                   \item{\code{b} = offset parameter used by the mgcv::gauslss. By default \code{b=0}.}
#'                   \item{\code{vtype} = type of variance estimator used to standardize the deviation from the main fit in the calibration.
#'                                        If set to \code{"m"} the variance estimate obtained by the full data fit is used, if set to \code{"b"}
#'                                        than the variance estimated produced by the bootstrap fits are used. By default \code{vtype="m"}.}
#'                   \item{\code{epsB} = positive tolerance used to assess convergence when fitting the regression coefficients on bootstrap data.  
#'                                       In particular, if \code{|dev-dev_old|/(|dev|+0.1)<epsB} then convergence is achieved. 
#'                                       Default is \code{epsB=1e-5}.}
#'                   \item{\code{verbose} = if TRUE some more details are given. By default \code{verbose=FALSE}.}
#'                   \item{\code{link} = link function to be used. See \code{?elf} and \code{?elflss} for defaults.}
#'                   \item{\code{progress} = argument passed to plyr::llply. By default \code{progress="text"} so that progress
#'                                           is reported. Set it to \code{"none"} to avoid it.}
#' }
#' @param argGam A list of parameters to be passed to \code{mgcv::gam}. This list can potentially include all the arguments listed
#'               in \code{?gam}, with the exception of \code{formula}, \code{family} and \code{data}.
#' @return A list with entries: \itemize{
#'                   \item{\code{lsig} = the value of log(sigma) resulting in the lowest loss.}
#'                   \item{\code{loss} = vector containing the value of the calibration loss function corresponding 
#'                                       to each value of log(sigma).}
#'                   \item{\code{edf} = a matrix where the first colums contain the log(sigma) sequence, and the remaining
#'                                      columns contain the corresponding effective degrees of freedom of each smooth.}
#'                   \item{\code{convProb} = a logical vector indicating, for each value of log(sigma), whether the outer
#'                                           optimization which estimates the smoothing parameters has encountered convergence issues.
#'                                           \code{FALSE} means no problem.}
#'                                           }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2017). Fast calibrated additive quantile regression. Available at
#'             \url{https://arxiv.org/abs/1707.03307}.
#' @examples
#' library(qgam); library(MASS)
#' 
#' # Calibrate learning rate on a grid
#' set.seed(41444)
#' sigSeq <- seq(1.5, 5, length.out = 10)
#' closs <- tuneLearn(form = accel~s(times,k=20,bs="ad"), 
#'                    data = mcycle, 
#'                    lsig = sigSeq, 
#'                    qu = 0.5)
#' 
#' plot(sigSeq, closs$loss, type = "b", ylab = "Calibration Loss", xlab = "log(sigma)")
#' 
#' # Pick best log-sigma
#' best <- sigSeq[ which.min(closs$loss) ]
#' abline(v = best, lty = 2)
#' 
#' # Fit using the best sigma
#' fit <- qgam(accel~s(times,k=20,bs="ad"), data = mcycle, qu = 0.5, lsig = best)
#' summary(fit)
#' 
#' pred <- predict(fit, se=TRUE)
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", 
#'      ylim = c(-150, 80))
#' lines(mcycle$times, pred$fit, lwd = 1)
#' lines(mcycle$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
#' lines(mcycle$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)                        
#'
tuneLearn <- function(form, data, lsig, qu, err = NULL, 
                      multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                      control = list(), argGam = NULL)
{ 
  if( length(qu) > 1 ) stop("length(qu) > 1, but this method works only for scalar qu")
  
  # Removing all NAs, unused variables and factor levels from data
  data <- .cleanData(.dat = data, .form = form, .drop = argGam$drop.unused.levels)
  
  lsig <- sort( lsig )
  
  # Setting up control parameter
  ctrl <- list( "loss" = "calFast", "sam" = "boot", "K" = 50, "b" = 0, "vtype" = "m", "epsB" = 1e-5, "verbose" = FALSE, 
                "link" = "identity", 
                "progress" = ifelse(multicore, "none", "text") )
  
  # Checking if the control list contains unknown names. Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  if( ctrl$progress == FALSE ) { ctrl$progress <- "none" }
  if( !(ctrl$vtype%in%c("m", "b")) ) stop("control$vtype should be either \"m\" or \"b\" ")
  if( !(ctrl$loss%in%c("calFast", "cal", "pin")) ) stop("control$loss should be either \"cal\", \"pin\" or \"calFast\" ")
  if( !(ctrl$sam%in%c("boot", "kfold")) ) stop("control$sam should be either \"boot\" or \"kfold\" ")
  if( (ctrl$loss=="cal") && (ctrl$sam=="kfold")  ) stop("You can't use control$sam == \"kfold\" when ctrl$loss==\"cal\" ")
  
  n <- nrow(data)
  nt <- length(lsig)
  
  # Gaussian fit, used for initializations 
  # NB Initializing smoothing parameters using gausFit is a very BAD idea
  if( is.formula(form) ) {
    gausFit <- do.call("gam", c(list("formula" = form, "data" = quote(data), 
                                     "family" = gaussian(link=ctrl[["link"]]))), argGam)
    varHat <- gausFit$sig2
    initM <- list("start" = coef(gausFit) + c(quantile(gausFit$residuals, qu), rep(0, length(coef(gausFit))-1)), 
                  "in.out" = NULL) # let gam() initialize sp via initial.spg() 
    formL <- form
  } else {
    gausFit <- do.call("gam", c(list("formula" = form, "data" = quote(data), 
                                     "family" = gaulss(link=list(ctrl[["link"]], "logb"), b=ctrl[["b"]])), argGam))
    varHat <- 1/gausFit$fit[ , 2]^2
    initM <- list("start" = NULL, "in.out" = NULL) # Have no cluse
    formL <- form[[1]]
  }  
  
  # Get loss smoothness
  if( is.null(err) ){ err <- .getErrParam(qu = qu, gFit = gausFit) }
  
  # For each value of 'lsig' fit on full data
  main <- .tuneLearnFullFits(lsig = lsig, form = formL, fam = "elf", qu = qu, err = err,
                             ctrl = ctrl, data = data, argGam = argGam, gausFit = gausFit, 
                             varHat = varHat, initM = initM)
  
  # Get score for each value of 'lsig'
  outLoss <- if( ctrl$loss == "calFast" ){ # Fast calibration (loss already calculated) OR ...
    sapply(main[["store"]], "[[", "loss")
  } else { # ... bootstrapping or cross-validation
    .tuneLearnBootstrapping(lsig = lsig, form = formL, fam = "elf", qu = qu, ctrl = ctrl, 
                            data = data, store = main[["store"]], pMat = main[["pMat"]], 
                            argGam = argGam, multicore = multicore, cluster = cluster, 
                            ncores = ncores, paropts = paropts)
  }
  names( outLoss ) <- lsig
  
  # convProb indicates whether there have been convergence problems during smoothing parameter estimation
  convProb <- sapply(main[["store"]], "[[", "convProb")
  names(convProb) <- lsig
    
  out <- list("lsig" = lsig[which.min(outLoss)], "loss" = outLoss, 
              "edf" = main[["edfStore"]], "convProb" = convProb)
  attr(out, "class") <- "learn"
  
  return( out )
  
}






