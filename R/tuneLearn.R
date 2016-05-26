##########################
#' Tuning the learning rate for Gibbs posterior
#' 
#' @description The learning rate (sigma) of the Gibbs posterior is tuned using a calibration approach,
#'              based on boostrapping. Here the calibration loss function is evaluated on a grid of values
#'              provided by the user.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
#' @param lsig A vector of value of the log learning rate (log(sigma)) over which the calibration loss function is evaluated.
#' @param qu The quantile of interest. Should be in (0, 1).
#' @param err An upper bound on the error of the estimated quantile curve. Should be in (0, 1). See XXX for details.
#' @param multicore If TRUE the calibration will happen in parallel.
#' @param ncores Number of cores used. Relevant if \code{multicore == TRUE}.
#' @param cluster An object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster.
#' @param paropts a list of additional options passed into the foreach function when parallel computation is enabled. 
#'                This is important if (for example) your code relies on external data or packages: 
#'                use the .export and .packages arguments to supply them so that all cluster nodes 
#'                have the correct environment set up for computing. 
#' @param control A list of control parameters for \code{tuneLearn} with entries: \itemize{
#'                   \item{\code{K} = number of boostrap datasets used for calibration. By default \code{K=50}.}
#'                   \item{\code{b} = offset parameter used by the mgcv::gauslss. By default \code{b=0}.}
#'                   \item{\code{verbose} = if TRUE some more details are given. By default \code{verbose=FALSE}.}
#'                   \item{\code{progress} = argument passed to plyr::llply. By default \code{progress="text"} so that progress
#'                                           is reported. Set it to \code{"none"} to avoid it.}
#' }
#' @param controlGam A list of control parameters to be passed to the internal \code{mgcv::gam} calls. 
#'                   See the \code{control} argument in \code{?mgcv::gam}.
#' @return A vector containing the value of the calibration loss function corresponding to each value of log(sigma).
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' library(qgam); library(MASS)
#' 
#' # Calibrate learning rate on a grid
#' set.seed(41444)
#' sigSeq <- seq(1.5, 5, length.out = 10)
#' closs <- tuneLearn(form = accel~s(times,k=20,bs="ad"), 
#'                    data = mcycle, 
#'                    err = 0.01, 
#'                    lsig = sigSeq, 
#'                    qu = 0.5)
#' 
#' plot(sigSeq, closs, type = "b", ylab = "Calibration Loss", xlab = "log(sigma)")
#' 
#' # Pick best log-sigma
#' best <- sigSeq[ which.min(closs) ]
#' abline(v = best, lty = 2)
#' 
#' # Fit using the best sigma
#' fit <- qgam(accel~s(times,k=20,bs="ad"), data = mcycle, qu = 0.5, err = 0.01, lsig = best)
#' summary(fit)
#' 
#' pred <- predict(fit, se=TRUE)
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", 
#'      ylim = c(-150, 80))
#' lines(mcycle$times, pred$fit, lwd = 1)
#' lines(mcycle$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
#' lines(mcycle$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)                        
#' @export tuneLearn  
#'
tuneLearn <- function(form, data, lsig, qu, err = 0.01, 
                      multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                      control = list(), controlGam = list() )
{ 
  if( length(qu) > 1 ) stop("length(qu) > 1, but this method works only for scalar qu")
  
  lsig <- sort( lsig )
  
  # Setting up control parameter
  ctrl <- list( "K" = 50, "b" = 0, "verbose" = FALSE, "progress" = ifelse(multicore, "none", "text") )
   
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  n <- nrow(data)
  nt <- length(lsig)
  
  if( multicore ){ 
    # Making sure "qgam" is loaded on cluser
    paropts[[".packages"]] <- unique( c("qgam", paropts[[".packages"]]) )
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores) #, exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
  }
  
  # Create K boostrap dataset
  ind <- lapply(1:ctrl[["K"]], function(nouse) sample(1:n, n, replace = TRUE))
  boot <- lapply(ind, function(ff) data[ff, ] )

  # Gaussian fit, used for initializations 
  if( is.formula(form) ) {
    fam <- "logF"
    gausFit <- gam(form, data = data, control = controlGam)
    varHat <- gausFit$sig2
    initM <- list("start" = coef(gausFit) + c(qnorm(qu, 0, sqrt(gausFit$sig2)), rep(0, length(coef(gausFit))-1)), 
                  "in.out" = list("sp" = gausFit$sp, "scale" = 1)) 
  } else {
    fam <- "logFlss"
    gausFit <- gam(form, data = data, family = gaulss(b=ctrl[["b"]]), control = controlGam)
    varHat <- 1/gausFit$fit[ , 2]^2
    initM <- list("start" = NULL, "in.out" = list("sp" = gausFit$sp, "scale" = 1)) 
  }  # Start = NULL in gamlss because it's not to clear how to deal with model for sigma 
  
  # Create gam object for full data fits
  mainObj <- gam(form, family = get(fam)(qu = qu, lam = NA, theta = NA), data = data, control = controlGam, fit = FALSE)
  
  # FULL data fits, used to estimate the smoothing parameters 
  mainFit <- vector("list", nt)
  for( ii in nt:1 ) # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harder)
  {
    mainObj$family$putLam( err * sqrt(2*pi*varHat) / (2*log(2)*exp(lsig[ii])) )
    mainObj$family$putTheta( lsig[ii] )
    
    withCallingHandlers({
    fit <- gam(G = mainObj, in.out = initM[["in.out"]], start = initM[["start"]]) 
    }, warning = function(w) {
      if (length(grep("Fitting terminated with step failure", conditionMessage(w))) ||
          length(grep("Iteration limit reached without full convergence", conditionMessage(w))))
      {
        message( paste("log(sigma) = ", round(lsig[ii], 6), " : outer Newton did not converge fully.", sep = "") )
        invokeRestart("muffleWarning")
      }
    })
    
    initM <- list("start" = coef(fit), "in.out" = list("sp" = fit$sp, "scale" = 1))
    
    # Main fit will be used when fitting the bootstrap datasets
    mainFit[[ii]] <- list("sp" = fit$sp, "fit" = fit$fitted, "lam" = fit$family$getLam())
  }
  
  # Loop over bootstrap datasets to get standardized deviations from full data fit
  withCallingHandlers({
    z <- llply( .data = boot, 
                .fun = .getBootDev,
                .parallel = multicore,
                .progress = ctrl[["progress"]],
                .inform = ctrl[["verbose"]],
                .paropts = paropts,
                # ... from here
                data = data, lsig = lsig, form = form, fam = fam, qu = qu, mainFit = mainFit, controlGam = controlGam)
  }, warning = function(w) {
    # There is a bug in plyr concerning a useless warning about "..."
    if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
      invokeRestart("muffleWarning")
  })

  # Get stardardized deviations and calculate Anderson-Darling normality statistic 
  z <- do.call("cbind", z)
  loss <- apply(z, 1, function(.x) .adTest(as.vector(.x)))
  names(loss) <- lsig
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return( loss )
}


###########
# Internal function that fits the bootstrapped datasets and returns standardized deviations from full data fit. 
# To be run in parallel.
###########
.getBootDev <- function(bdat, data, lsig, form, fam, qu, mainFit, controlGam)
{ 
  nt <- length( lsig )
  n <- nrow( bdat )
  
  init <- NULL
  
  # Create gam object
  bObj <- gam(form, family = get(fam)(qu = qu, lam = NA, theta = NA), data = bdat, 
              sp = mainFit[[1]]$sp, control = controlGam, fit = FALSE)
  
  .z <- matrix(NA, nt, n)
  for( ii in nt:1 )  # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harded)
  {   
    bObj$lsp0 <- log( mainFit[[ii]]$sp )
    bObj$family$putLam( mainFit[[ii]]$lam )
    bObj$family$putTheta( lsig[ii] )
    
    fit <- gam(G = bObj, start = init)
    init <- betas <- coef(fit)
    Vp <- fit$Vp
    
    # Create prediction design matrix (only in first iteration)
    if(ii == nt) { 
      pMat <- predict.gam(fit, newdata = data, type = "lpmatrix") 
      lpi <- attr(pMat, "lpi")
      if( !is.null(lpi) ){ pMat <- pMat[ , lpi[[1]]] }
    }
    
    # In the gamlss case, we are interested only in the calibrating the location mode
    if( !is.null(lpi) ){
      betas <- betas[lpi[[1]]]
      Vp <- fit$Vp[lpi[[1]], lpi[[1]]]
    }
    
    mu <- pMat %*% betas
    sdev <- sqrt( diag( pMat%*%Vp%*%t(pMat) ) )
    
    .z[ii, ] <- (mu - as.matrix(mainFit[[ii]]$fit)[ , 1]) / sdev
  }
  
  return( .z )
}
