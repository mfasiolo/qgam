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
  ctrl <- list( "loss" = "cal", "sam" = "boot", "K" = 50, "b" = 0, "vtype" = "m", "epsB" = 1e-5, "verbose" = FALSE, 
                "link" = if(is.formula(form)){"identity"}else{list("identity", "log")}, 
                "progress" = ifelse(multicore, "none", "text") )
  
  # Checking if the control list contains unknown names. Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  if( !(ctrl$vtype%in%c("m", "b")) ) stop("control$vtype should be either \"m\" or \"b\" ")
  if( !(ctrl$loss%in%c("cal", "pin")) ) stop("control$loss should be either \"cal\" or \"pin\" ")
  if( !(ctrl$sam%in%c("boot", "kfold")) ) stop("control$sam should be either \"boot\" or \"kfold\" ")
  if( (ctrl$loss=="cal") && (ctrl$sam=="kfold")  ) stop("You can't use control$sam == \"kfold\" when ctrl$loss==\"cal\" ")
  
  n <- nrow(data)
  nt <- length(lsig)
  
  if( ctrl$sam == "boot" ){ # Create weights for K boostrap dataset OR...
    wb <- lapply(1:ctrl[["K"]], function(nouse) tabulate(sample(1:n, n, replace = TRUE), n))
  } else { # ... OR for K training sets for CV 
    tmp <- sample(rep(1:ctrl[["K"]], length.out = n), n, replace = FALSE)
    wb <- lapply(1:ctrl[["K"]], function(ii) tabulate(which(tmp != ii), n)) 
  }
  
  # Gaussian fit, used for initializations 
  if( is.formula(form) ) {
    fam <- "elf"
    gausFit <- do.call("gam", c(list("formula" = form, "data" = data), argGam))
    varHat <- gausFit$sig2
  } else {
    fam <- "elflss"
    gausFit <- do.call("gam", c(list("formula" = form, "data" = data, "family" = gaulss(b=ctrl[["b"]])), argGam))
    varHat <- 1/gausFit$fit[ , 2]^2
  }  
  
  # Initializing regression coefficient using gausFit is not good idea, especially for extreme values of lsig
  initM <- list("start" = NULL, "in.out" = list("sp" = gausFit$sp, "scale" = 1))
  
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
  
  ############## START FITTING ON FULL DATA
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
    
    # Create prediction matrix (only in the first iteration)
    if( ii == nt ){
      pMat <- predict.gam(fit, type = "lpmatrix") 
      lpi <- attr(pMat, "lpi")
      if( !is.null(lpi) ){ 
        pMat <- pMat[ , lpi[[1]]] # "lpi" attribute lost here, re-inserted in next line 
        attr(pMat, "lpi") <- lpi }
    }
    
    sdev <- NULL 
    if(ctrl$loss == "cal" && ctrl$vtype == "m"){
      Vp <- fit$Vp
      # In the gamlss case, we are interested only in the calibrating the location mode
      if( !is.null(lpi) ){  Vp <- fit$Vp[lpi[[1]], lpi[[1]]]  }
      sdev <- sqrt(rowSums((pMat %*% Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
    }
    
    initM <- list("start" = coef(fit), "in.out" = list("sp" = fit$sp, "scale" = 1))
    
    # Main fit will be used when fitting the bootstrap datasets
    mainFit[[ii]] <- list("sp" = fit$sp, "fit" = fit$fitted, "lam" = fit$family$getLam(), "init" = initM$start, "sdev" = sdev)
  } ############## END FITTING ON FULL DATA
  
  # Create gam object for bootstrap fits
  bObj <- do.call("gam", c(list("formula" = form, "family" = get(fam)(qu = qu, lam = NA, theta = NA, link = ctrl$link), "data" = data, 
                                "sp" = if(length(mainFit[[1]]$sp)){mainFit[[1]]$sp}else{NULL}, fit = FALSE), argGam))
  
  # Preparing bootstrap object for gam.fit3
  bObj <- .prepBootObj(obj = bObj, eps = ctrl$epsB, control = argGam$control)
  
  # Internal function that fits the bootstrapped datasets and returns standardized deviations from full data fit. To be run in parallel.
  # GLOBALS: lsig, ctrl, mainFit, pMat, bObj, argGam
  .getBootDev <- function(.wb)
  {   # # # # # # # # # .getBootDev START # # # # # # # # #
    y <- bObj$y 
    ns <- length(lsig); n  <- length(y)
    
    # Number of test observations
    nt <- ifelse(ctrl$loss == "cal", n, sum(!.wb)) 
    
    # Creating boot weights from boot indexes 
    bObj$w <- .wb
    
    lpi <- attr(pMat, "lpi")
    glss <- inherits(bObj$family, "general.family")
    
    init <- NULL
    .z <- vector("list", ns)
    for( ii in ns:1 )  # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harded)
    {   
      # In the gamlss case lambda is a vector, and we have to take only those values that are in the boostrapped dataset.
      lam <- mainFit[[ii]]$lam
      
      bObj$lsp0 <- log( mainFit[[ii]]$sp )
      bObj$family$putLam( lam )
      bObj$family$putTheta( lsig[ii] )
      
      init <- if(is.null(init)){ list(mainFit[[ii]]$init) } else { list(init, mainFit[[ii]]$init) }
      
      if( glss ){
        init <- lapply(init, function(inp) mgcv:::Sl.initial.repara(bObj$Sl, inp, inverse=FALSE, both.sides=FALSE))
        fit <- .gamlssFit(x=bObj$X, y=bObj$y, lsp=as.matrix(bObj$lsp0), Sl=bObj$Sl, weights=bObj$w, 
                          offset=bObj$offset, family=bObj$family, control=bObj$control, 
                          Mp=bObj$Mp, start=init, needVb=(ctrl$loss=="cal" && ctrl$vtype=="b"))
        # In gamlss, we want to calibrate only the location and we need to reparametrize the coefficients
        init <- betas <- mgcv:::Sl.initial.repara(bObj$Sl, fit$coef, inverse=TRUE, both.sides=FALSE)
        betas <- betas[lpi[[1]]] 
      } else {
        bObj$null.coef <- bObj$family$get.null.coef(bObj)$null.coef
        fit <- .egamFit(x=bObj$X, y=bObj$y, sp=as.matrix(bObj$lsp0), Eb=bObj$Eb, UrS=bObj$UrS,
                        offset=bObj$offset, U1=bObj$U1, Mp=bObj$Mp, family = bObj$family, weights=bObj$w,
                        control=bObj$control, null.coef=bObj$null.coef, 
                        start=init, needVb=(ctrl$loss == "cal" && ctrl$vtype == "b"))
        init <- betas <- fit$coef
      }
      
      mu <- pMat %*% betas
      
      if( ctrl$loss == "cal" ){ # (1) Return standardized deviations from full data fit OR ... 
        if( ctrl$vtype == "b" ){ # (2) Use variance of bootstrap fit OR ...
          Vp <- .getVp(fit, bObj, bObj$lsp0, lpi)
          sdev <- sqrt(rowSums((pMat %*% Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
        } else { # (2)  ... variance of the main fit
          sdev <- mainFit[[ii]]$sdev
        }
        .z[[ii]] <- drop((mu - as.matrix(mainFit[[ii]]$fit)[ , 1]) / sdev)
      } else { # (1) ... out of sample observations minus their fitted values 
        .z[[ii]] <- drop(y - mu)[ !.wb ]
      }
    }
    .z <- do.call("rbind", .z)
    return( .z )
  }  # # # # # # # # # .getBootDev END # # # # # # # # #
  
  if( multicore ){ 
    # Making sure "qgam" is loaded on cluser
    paropts[[".packages"]] <- unique( c("qgam", paropts[[".packages"]]) )
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores) #, exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
    
    # Exporting stuff. To about all environment being exported all the time, use .GlobalEnv  
    clusterExport(cluster, c("pMat", "bObj", "lsig", "ctrl", "mainFit", "argGam", ".getVp", ".egamFit", ".gamlssFit"), 
                  envir = environment())
    environment(.getBootDev) <- .GlobalEnv
  }
  
  # Loop over bootstrap datasets to get standardized deviations from full data fit
  withCallingHandlers({
    z <- llply( .data = wb, 
                .fun = .getBootDev,
                .parallel = multicore,
                .progress = ctrl[["progress"]],
                .inform = ctrl[["verbose"]],
                .paropts = paropts)
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