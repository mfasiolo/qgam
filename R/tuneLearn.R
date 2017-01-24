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
#' @param err An upper bound on the error of the estimated quantile curve. Should be in (0, 1). See Fasiolo et al. (2016) for details.
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
#'                                       credible intervals for the fitted curve are calibrated. See Fasiolo et al. (2016) for details.
#'                                       If \code{loss=="pin"} then log(sigma) approximately minimizes the pinball loss on the out-of-sample
#'                                       data.}
#'                   \item{\code{sam} = sampling scheme use: \code{sam=="boot"} corresponds to bootstrapping and \code{sam=="kfold"} to k-fold
#'                                      cross-validation. The second option can be used only if \code{ctrl$loss=="pin"}.}
#'                   \item{\code{K} = if \code{sam=="boot"} this is the number of boostrap datasets, while if \code{sam=="kfold"} this is the 
#'                                    number of folds. By default \code{K=50}.}
#'                   \item{\code{b} = offset parameter used by the mgcv::gauslss. By default \code{b=0}.}
#'                   \item{\code{verbose} = if TRUE some more details are given. By default \code{verbose=FALSE}.}
#'                  \item{\code{link} = Link function to be used. See \code{?elf} and \code{?elflss} for defaults.}
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
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2016). Fast calibrated additive quantile regression. Available at
#'             \url{https://github.com/mfasiolo/qgam/blob/master/draft_qgam.pdf}.
#' @examples
#' library(qgam); library(MASS)
#' 
#' # Calibrate learning rate on a grid
#' set.seed(41444)
#' sigSeq <- seq(1.5, 5, length.out = 10)
#' closs <- tuneLearn(form = accel~s(times,k=20,bs="ad"), 
#'                    data = mcycle, 
#'                    err = 0.05, 
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
#' fit <- qgam(accel~s(times,k=20,bs="ad"), data = mcycle, qu = 0.5, err = 0.05, lsig = best)
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
tuneLearn <- function(form, data, lsig, qu, err = 0.05, 
                      multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                      control = list(), argGam = NULL)
{ 
  if( length(qu) > 1 ) stop("length(qu) > 1, but this method works only for scalar qu")
  
  # Removing all NAs from data
  data <- na.omit( data )

  lsig <- sort( lsig )
  
  # Setting up control parameter
  ctrl <- list( "loss" = "cal", "sam" = "boot", "K" = 50, "b" = 0, "verbose" = FALSE, 
                "link" = if(is.formula(form)){"identity"}else{list("identity", "log")}, 
                "progress" = ifelse(multicore, "none", "text") )
   
  # Checking if the control list contains unknown names. Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  if( !(ctrl$loss%in%c("cal", "pin")) ) stop("control$loss should be either \"cal\" or \"pin\" ")
  if( !(ctrl$sam%in%c("boot", "kfold")) ) stop("control$sam should be either \"boot\" or \"kfold\" ")
  if( (ctrl$loss=="cal") && (ctrl$sam=="kfold")  ) stop("You can't use control$sam == \"kfold\" when ctrl$loss==\"cal\" ")
  
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
  
  if( ctrl$sam == "boot" ){ # Create indexes of K boostrap dataset OR...
   bootInd <- lapply(1:ctrl[["K"]], function(nouse) sample(1:n, n, replace = TRUE))
  } else { # ... OR for K training sets for CV 
   tmp <- sample(rep(1:ctrl[["K"]], length.out = n), n, replace = FALSE)
   bootInd <- lapply(1:ctrl[["K"]], function(ii) which(tmp != ii)) 
  }

  # Gaussian fit, used for initializations 
  if( is.formula(form) ) {
    fam <- "elf"
    gausFit <- do.call("gam", c(list("formula" = form, "data" = data), argGam))
    varHat <- gausFit$sig2
    initM <- list("start" = coef(gausFit) + c(qnorm(qu, 0, sqrt(gausFit$sig2)), rep(0, length(coef(gausFit))-1)), 
                  "in.out" = list("sp" = gausFit$sp, "scale" = 1)) 
  } else {
    fam <- "elflss"
    gausFit <- do.call("gam", c(list("formula" = form, "data" = data, "family" = gaulss(b=ctrl[["b"]])), argGam))
    varHat <- 1/gausFit$fit[ , 2]^2
    initM <- list("start" = NULL, "in.out" = list("sp" = gausFit$sp, "scale" = 1)) 
  }  # Start = NULL in gamlss because it's not to clear how to deal with model for sigma 
  
  # Create gam object for full data fits
  mainObj <- do.call("gam", c(list("formula" = form, "family" = get(fam)(qu = qu, lam = NA, theta = NA, link = ctrl$link), "data" = data, fit = FALSE), argGam))
  
  # Store degrees of freedom for each value of lsig
  tmp <- pen.edf(gausFit)
  if( length(tmp) )
  {
    edfStore <- matrix(NA, nt, length(tmp) + 1)
    colnames(edfStore) <- c("lsig", names( tmp ))
  } else {
    edfStore <- NULL
  } 
  
  # Vector indicating convergence problems
  convProb <- rep(FALSE, nt)
  names(convProb) <- lsig

  # FULL data fits, used to estimate the smoothing parameters 
  mainFit <- vector("list", nt)
  for( ii in nt:1 ) # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harder)
  {
    mainObj$family$putLam( err * sqrt(2*pi*varHat) / (2*log(2)*exp(lsig[ii])) )
    mainObj$family$putTheta( lsig[ii] )
    
    withCallingHandlers({
    fit <- do.call("gam", c(list("G" = mainObj, "in.out" = initM[["in.out"]], "start" = initM[["start"]]), argGam)) 
    }, warning = function(w) {
      if (length(grep("Fitting terminated with step failure", conditionMessage(w))) ||
          length(grep("Iteration limit reached without full convergence", conditionMessage(w))))
      {
        message( paste("log(sigma) = ", round(lsig[ii], 3), " : outer Newton did not converge fully.", sep = "") )
        convProb[ii] <<- TRUE
        invokeRestart("muffleWarning")
      }
    })
    
    if( !is.null(edfStore) ) { edfStore[ii, ] <- c(lsig[ii], pen.edf(fit)) }
    
    initM <- list("start" = coef(fit), "in.out" = list("sp" = fit$sp, "scale" = 1))
    
    # Main fit will be used when fitting the bootstrap datasets
    mainFit[[ii]] <- list("sp" = fit$sp, "fit" = fit$fitted, "lam" = fit$family$getLam())
  }
  
  # Loop over bootstrap datasets to get standardized deviations from full data fit
  withCallingHandlers({
    z <- llply( .data = bootInd, 
                .fun = .getBootDev,
                .parallel = multicore,
                .progress = ctrl[["progress"]],
                .inform = ctrl[["verbose"]],
                .paropts = paropts,
                # ... from here
                data = data, lsig = lsig, form = form, fam = fam, link = ctrl$link, qu = qu,
                loss = ctrl$loss, mainFit = mainFit, argGam = argGam)
  }, warning = function(w) {
    # There is a bug in plyr concerning a useless warning about "..."
    if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
      invokeRestart("muffleWarning")
  })

  # Get stardardized deviations and ... 
  z <- do.call("cbind", z)
  if( ctrl$loss == "cal"){ # ... calculate Anderson-Darling normality statistic OR ...
   outLoss <- apply(z, 1, function(.x) .adTest(as.vector(.x)))
  } else { # ... pinball loss
   outLoss <- apply(z, 1, function(.x) .checkloss(as.vector(.x), 0, qu = qu))
  }
  names(outLoss) <- lsig
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  out <- list("lsig" = lsig[which.min(outLoss)], "loss" = outLoss, "edf" = edfStore, "convProb" = convProb)
  attr(out, "class") <- "learn"
  
  return( out )
}


###########
# Internal function that fits the bootstrapped datasets and returns standardized deviations from full data fit. 
# To be run in parallel.
###########
.getBootDev <- function(bindex, data, lsig, form, fam, link, qu, loss, mainFit, argGam)
{ 
  nt <- length( lsig )
  n <- length( bindex )
  
  # Creating bootstrapped dataset. Redundant levels of factor variable are dropped
  bdat <- droplevels( data[bindex, , drop = FALSE] )
  
  if( loss == "cal" ){ # Testing data: full data if we calibrate OR ...
    odat <- data 
  } else { # ... and out-of-sample data and responses if we cross-validate on pinball los
    odat <- data[-bindex, ]  
    testObs <- odat[[ as.character(if(is.formula(form)){ as.list(form)[[2]] }else{ as.list(form[[1]])[[2]] }) ]]
  }
  
  init <- NULL
  
  # Create gam object
  bObj <- do.call("gam", c(list("formula" = form, "family" = get(fam)(qu = qu, lam = NA, theta = NA, link = link), "data" = bdat, 
                                "sp" = if(length(mainFit[[1]]$sp)){mainFit[[1]]$sp}else{NULL}, fit = FALSE), argGam))

  .z <- matrix(NA, nt, nrow(odat))
  for( ii in nt:1 )  # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harded)
  {   
    # In the gamlss case lambda is a vector, and we have to take only those values that are in the boostrapped dataset.
    lambda <- mainFit[[ii]]$lam[ if(is.formula(form)){ 1 } else { bindex } ]
    
    bObj$lsp0 <- log( mainFit[[ii]]$sp )
    bObj$family$putLam( lambda )
    bObj$family$putTheta( lsig[ii] )
    
    fit <- do.call("gam", c(list("G" = bObj, "start" = init), argGam))
    init <- betas <- coef(fit)
    Vp <- fit$Vp
    
    # Create prediction design matrix (only in first iteration)
    if(ii == nt) { 
      pMat <- predict.gam(fit, newdata = odat, type = "lpmatrix") 
      lpi <- attr(pMat, "lpi")
      if( !is.null(lpi) ){ pMat <- pMat[ , lpi[[1]]] }
    }
    
    # In the gamlss case, we are interested only in the calibrating the location mode
    if( !is.null(lpi) ){
      betas <- betas[lpi[[1]]]
      Vp <- fit$Vp[lpi[[1]], lpi[[1]]]
    }
    
    mu <- pMat %*% betas
    
    if( loss == "cal" ){ # Return stardardized deviations from full data fit OR ...
     sdev <- sqrt(rowSums((pMat %*% Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
     .z[ii, ] <- (mu - as.matrix(mainFit[[ii]]$fit)[ , 1]) / sdev
    } else { # ... out of sample observations minus their fitted values 
     .z[ii, ] <- testObs - mu
    }
  }
  
  return( .z )
}
