##########################
#' Fast learning rate calibration for the Gibbs posterior
#' 
#' @description The learning rate (sigma) of the Gibbs posterior is tuned using a calibration approach,
#'              based on boostrapping. Here the loss function is minimized, for each quantile, using a Brent search.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
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
#'                   \item{\code{init} = an initial value for the log learning rate (log(sigma)). 
#'                                       By default \code{init=NULL} and the optimization is initialized by other means.}
#'                   \item{\code{brac} = initial bracket for Brent method. By default \code{brac=c(0.5, 2)}, so the initial 
#'                                       search range is \code{(init - 0.5, init + 2)}.}
#'                   \item{\code{redWd} = parameter which determines when the bracket size needs to be reduced.
#'                                        If \code{redWd==10} then the bracket is halved if the nearest solution
#'                                        falls within the central 10\% of its width. By default \code{redWd = 10}.}
#'                   \item{\code{K} = number of boostrap datasets used for calibration. By default \code{K=50}.}
#'                   \item{\code{b} = offset parameter used by the mgcv::gauslss. By default \code{b=0}.}
#'                   \item{\code{tol} = tolerance used in the Brent search. By default \code{tol=.Machine$double.eps^0.25}.
#'                                      See \code{?optimize} for details.}
#'                   \item{\code{verbose} = if TRUE some more details are given. By default \code{verbose=FALSE}.}
#' }
#' @param argGam A list of parameters to be passed to \code{mgcv::gam}. This list can potentially include all the arguments listed
#'               in \code{?gam}, with the exception of \code{formula}, \code{family} and \code{data}.
#' @return A list with entries: \itemize{
#'                   \item{\code{lsig} = a vector containing the values of log(sigma) that minimize the loss function, 
#'                                       for each quantile.}
#'                   \item{\code{err} = the error bound used for each quantile. Generally each entry is identical to the
#'                                      argument \code{err}, but in some cases the function increases it to enhance stabily.}
#'                   \item{\code{ranges} = the search ranges by the Brent algorithm to find log-sigma, for each quantile. }
#'                   \item{\code{store} = a list, where the i-th entry is a matrix containing all the locations (1st row) at which
#'                                        the loss function has been evaluated and its value (2nd row), for the i-th quantile.}
#' }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2016). Fast calibrated additive quantile regression. Available at
#'             \url{https://github.com/mfasiolo/qgam/draft_qgam.pdf}.
#' @examples
#' library(qgam); library(MASS)
#' 
#' ###
#' # Single quantile fit
#' ###
#' # Calibrate learning rate on a grid
#' set.seed(5235)
#' tun <- tuneLearnFast(form = accel~s(times,k=20,bs="ad"), 
#'                      data = mcycle, 
#'                      err = 0.05, 
#'                      qu = 0.2)
#' 
#' # Fit for quantile 0.2 using the best sigma
#' fit <- qgam(accel~s(times, k=20, bs="ad"), data = mcycle, qu = 0.2,
#'             err = 0.05, lsig = tun$lsig)
#' 
#' pred <- predict(fit, se=TRUE)
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", 
#'      ylim = c(-150, 80))
#' lines(mcycle$times, pred$fit, lwd = 1)
#' lines(mcycle$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
#' lines(mcycle$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2) 
#' 
#' ###
#' # Multiple quantile fits
#' ###
#' # Calibrate learning rate on a grid
#' quSeq <- c(0.25, 0.5, 0.75)
#' set.seed(5235)
#' tun <- tuneLearnFast(form = accel~s(times, k=20, bs="ad"), 
#'                      data = mcycle, 
#'                      err = 0.05, 
#'                      qu = quSeq)
#' 
#' # Fit using estimated sigmas
#' fit <- mqgam(accel~s(times, k=20, bs="ad"), data = mcycle, qu = quSeq,
#'              err = 0.05, lsig = tun$lsig)
#' 
#' # Plot fitted quantiles
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", 
#'      ylim = c(-150, 80))
#' for(iq in quSeq){
#'   pred <- qdo(fit, iq, predict)
#'   lines(mcycle$times, pred, col = 2)
#' }                   
#' @export tuneLearnFast
#'
tuneLearnFast <- function(form, data, qu, err = 0.01,
                          multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                          control = list(), argGam = NULL)
{ 
  n <- nrow(data)
  nq <- length(qu)
  
  # Setting up control parameter
  ctrl <- list( "init" = NULL,
                "brac" = log( c(1/2, 2) ), 
                "K" = 50,
                "redWd" = 10,
                "tol" = .Machine$double.eps^0.25,
                "b" = 0,
                "gausFit" = NULL,
                "verbose" = FALSE )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  tol <- ctrl[["tol"]]
  brac <- ctrl[["brac"]]
  
  # Sanity check
  if( tol > 0.1 * abs(diff(brac)) ) stop("tol > bracket_widths/10, choose smaller tolerance or larger bracket")
  
  # Create indexes of K boostrap dataset
  bootInd <- lapply(1:ctrl[["K"]], function(nouse) sample(1:n, n, replace = TRUE))
  
  # Gaussian fit, used for initialization
  if( is.formula(form) ) {
    fam <- "logF"                      
    if( is.null(ctrl[["gausFit"]]) ) { gausFit <- do.call("gam", c(list("formula" = form, "data" = data), argGam)) } else { gausFit <- ctrl$gausFit }
    varHat <- gausFit$sig2
  } else {
    fam <- "logFlss"                 
    if( is.null(ctrl[["gausFit"]]) ) { gausFit <- do.call("gam", c(list("formula" = form, "data" = data, "family" = gaulss(b=ctrl[["b"]])), argGam)) } else { gausFit <- ctrl$gausFit }
    varHat <- 1/gausFit$fit[ , 2]^2
  }
  
  # Order quantiles so that those close to the median are dealt with first
  oQu <- order( abs(qu-0.5) )
  
  # (Optional) Initializing the search range for sigma
  if( is.null(ctrl[["init"]]) ){
    # We assume lam~0 and we match (2 times) the variance of a symmetric (median) Laplace density with that of the Gaussian fit.
    # This is an over-estimate for extreme quantiles, but experience suggests that it's better erring on the upper side.
    tmp <- 0.5 #qu[ oQu[1] ]
    if( !is.list(form) ){
      isig <- log(sqrt( 5 * gausFit$sig2 * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
    } else {
      isig <- log(sqrt( 5 * (ctrl[["b"]]+exp(coef(gausFit)["(Intercept).1"]))^2 * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
    }
  } else {
    isig <- ctrl[["init"]]
  }
  
  # Create gam object for full data fits
  mObj <- do.call("gam", c(list("formula" = form, "family" = get(fam)(qu = NA, lam = NA, theta = NA), 
                                "data" = data, "fit" = FALSE), argGam))
  
  # Create a gam object for each bootstrap sample
  bObj <- lapply(bootInd, function(.ind){
    out <- do.call("gam", c(list("formula" = form, "family" = get(fam)(qu = NA, lam = NA, theta = NA), 
                                 "data" = data[.ind, ], "sp" = gausFit$sp, "fit" = FALSE), argGam))
    # Save boostrap indexes to be used later 
    out$bootInd <- .ind
    return( out )
  })
  
  # Create prediction design matrices for each bootstrap sample
  pMat <- lapply(bObj, function(.obj){
    class( .obj ) <- c("gam", "glm", "lm") 
    .obj$coefficients <- rep(0, ncol(.obj$X)) # Needed to fool predict
    out <- predict.gam(.obj, newdata = data, type = "lpmatrix")
    
    # Calibration uses the linear predictor for the quantile location, we discard the rest 
    lpi <- attr(gausFit$formula, "lpi")
    if( !is.null(lpi) ){  
      out <- out[ , lpi[[1]]] #"lpi" attribute lost here, re-inserted in next line 
      attr(out, "lpi") <- lpi 
    }
    
    return( out )
  })
  
  if( multicore ){ 
    # Create cluster
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores) #, exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
    
    # Load "qgam" and user-specified packages
    tmp <- unique( c("qgam", paropts[[".packages"]]) )
    clusterExport(cluster, "tmp", envir = environment())
    clusterEvalQ(cluster, { lapply(tmp, library, character.only = TRUE) })
    paropts[[".packages"]] <- NULL
    
    # Export bootstrap objects, prediction matrix and user-defined stuff
    tmp <- unique( c("bObj", "pMat", paropts[[".export"]]) )
    clusterExport(cluster, tmp, envir = environment())
    paropts[[".export"]] <- NULL
  }
  
  # Estimated learning rates, num of bracket expansions, error rates and bracket ranges used in bisection
  sigs <- efacts <- errors <- numeric(nq)
  rans <- matrix(NA, nq, 2)
  store <- vector("list", nq)
  names(sigs) <- names(errors) <- rownames(rans) <- qu
  
  # Here we need bTol > aTol otherwise we the new bracket will be too close to the probable solution
  aTol <- 0.05
  bTol <- 0.2
  
  for(ii in 1:nq)
  {
    oi <- oQu[ii]
    
    ef <- 1
    
    repeat{
      
      # Compute bracket
      srange <- isig + ef * brac
      
      # Estimate log(sigma) using brent methods with current bracket (srange)
      res  <- .tuneLearnFast(mObj = mObj, bObj = bObj, pMat = pMat, qu = qu[oi], err = err, srange = srange, 
                             gausFit = gausFit, varHat = varHat,
                             multicore = multicore, cluster = cluster, ncores = ncores, paropts = paropts,  
                             control = ctrl, argGam = argGam)  
      
      # Store loss function evaluations
      store[[oi]] <- cbind(store[[oi]], res[["store"]])
      lsig <- res$minimum
      
      # If solution not too close to boundary store results and determine bracket for next iteration
      if( all(abs(lsig-srange) > aTol * abs(diff(srange))) ){ 
        
        sigs[oi] <- lsig
        rans[oi, ] <- srange
        efacts[oi] <- ef
        errors[oi] <- res$err
        
        # Determine what quantile needs to be dealt with next, then choose bracket and initialization using old results
        if(ii < nq)
        {
          kk <- oQu[ which.min(abs(qu[oQu[ii+1]] - qu[oQu[1:ii]])) ]
          isig <- sigs[kk] 
          wd <- abs(diff(rans[kk, ]))
          brac <- c(-1, 1) * wd / 2
          
          # If kk solution close to center of kk bracket, halve the bracket size 
          # (unless the size of the bracket is < 10*tol or the bracket has been expanded in the old iteration)
          if( (abs(isig - mean(rans[kk, ])) < wd/ctrl$redWd) && (wd > 10*tol) && (efacts[kk] == 1)) brac <- brac / 2
        }
        
        break
      }
      
      # If solution is close to bracket boundaries, we shift bracket and expand it
      # This (- wd + bTol*wd)/2 is divided by 2 to make the algorithm more reluctant to reduce lsig
      wd <- abs( diff(brac) )
      isig <- lsig + ifelse(lsig-srange[1] < aTol*wd, (- wd + bTol*wd)/2, wd - bTol*wd)
      ef <- 2*ef
    }
    
    if( ctrl$verbose && (nq>1) )
    {
      tseq <- oQu[1:ii]
      tmp <- rans[tseq, , drop = FALSE]
      layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
             heights=c(2, 1))
      par(mai = c(1, 1, 0.1, 0.1))
      plot(qu[tseq], sigs[oQu[1:ii]], ylim = range(as.vector(tmp)), xlim = range(qu), col = 2, 
           ylab = expression("Log(" * sigma * ")"), xlab = "qu")
      points(qu[tseq], tmp[ , 1], pch = 3)
      points(qu[tseq], tmp[ , 2], pch = 3)
      points(qu[tseq], rowMeans(tmp), pch = 3)
      for(zz in 1:ii) segments(qu[oQu[zz]], rowMeans(rans)[oQu[zz]] - abs(diff(tmp[zz, ]))/ctrl$redWd, 
                               qu[oQu[zz]], rowMeans(rans)[oQu[zz]] + abs(diff(tmp[zz, ]))/ctrl$redWd, col = 1)
      plot(qu, efacts, xlab = "qu", "ylab" = "Bracket expansions")  
      plot(qu, errors)
    }
  }
  
  names(sigs) <- qu
  
  out <- list("lsig" = sigs, "err" = errors, "ranges" = rans, "store" = store)
  attr(out, "class") <- "tuneLearnFast"
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return( out )
}

##########################################################################
### Internal version, which works for a single quantile qu
########################################################################## 
.tuneLearnFast <- function(mObj, bObj, pMat, qu, err, srange,
                           gausFit, varHat,
                           multicore, cluster, ncores, paropts, 
                           control, argGam)
{
  nbo <- length( bObj )
  
  # Initialization of regression and smoothing coefficients for full data (not bootstrap) fit, using Gaussian fit 
  if( is.formula(mObj$formula) ) { # Extended Gam OR ...
    initM <- list("start" = coef(gausFit) + c(qnorm(qu, 0, sqrt(gausFit$sig2)), rep(0, length(coef(gausFit))-1)), 
                  "in.out" = list("sp" = gausFit$sp, "scale" = 1)) 
    
  } else { # ... GAMLSS
    initM <- list("start" = NULL, "in.out" = list("sp" = gausFit$sp, "scale" = 1)) 
  }  # Start = NULL in gamlss because it's not to clear how to deal with model for sigma 
  
  # Loss function to be minimized using Brent method
  objFun <- function(lsig, mObj, bObj, initM, initB, pMat, qu, varHat, cluster)
  {
    nbo <- length( bObj )
    lam <- err * sqrt(2*pi*varHat) / (2*log(2)*exp(lsig))
    
    mObj$family$putQu( qu )
    mObj$family$putLam( lam )
    mObj$family$putTheta( lsig )
    
    # Full data fit
    withCallingHandlers({
      mFit <- do.call("gam", c(list("G" = mObj, "in.out" = initM[["in.out"]], "start" = initM[["start"]]), argGam))}, warning = function(w) {
        if (length(grep("Fitting terminated with step failure", conditionMessage(w))) ||
            length(grep("Iteration limit reached without full convergence", conditionMessage(w))))
        {
          message( paste("qu = ", qu, ", log(sigma) = ", round(lsig, 6), " : outer Newton did not converge fully.", sep = "") )
          invokeRestart("muffleWarning")
        }
      })
    
    mMU <- as.matrix(mFit$fit)[ , 1]
    mSP <- mFit$sp
    
    initM <- list("start" = coef(mFit), "in.out" = list("sp" = mSP, "scale" = 1))
    
    ## Function to be run in parallel (over boostrapped datasets)
    # GLOBAL VARS: bObj, pMat, initB, mSP, mMU, lam, lsig, qu
    .funToApply <- function(ind)
    {
      z <- init <- vector("list", length(ind))
      for( kk in ind){
        .obj <- bObj[[kk]]; .pMat <- pMat[[kk]]; .init <- initB[[kk]];
        
        # In gamlss case lambda is not constant, so we pick the 
        # right values for the k-th bootstrap dataset
        .lambda <- if(is.formula(.obj$formula)){ lam[1] } else { lam[.obj$bootInd] }
        
        .obj$lsp0 <- log( mSP )
        .obj$family$putQu( qu )
        .obj$family$putLam( .lambda )
        .obj$family$putTheta( lsig )
        
        fit <- do.call("gam", c(list("G" = .obj, "start" = .init), argGam))
        
        .init <- betas <- coef(fit)
        Vp <- fit$Vp
        
        # In the gamlss case, we are interested only in the calibrating the location model
        # so we drop the coefficients related to the scale model
        lpi <- attr(.pMat, "lpi")
        if( !is.null(lpi) ){
          betas <- betas[lpi[[1]]]
          Vp <- fit$Vp[lpi[[1]], lpi[[1]]]
        }
        
        # Calculating stardardized deviations from full data fit
        mu <- .pMat %*% betas
        sdev <- sqrt(rowSums((.pMat %*% Vp) * .pMat)) # same as sqrt(diag(.pMat%*%Vp%*%t(.pMat))) but (WAY) faster
        
        z[[kk]] <- (mu - mMU) / sdev
        init[[kk]] <- .init
      }
      
      z <- do.call("c", z)
      
      return( list("z" = z, "init" = init) )
    }
    
    if( !is.null(cluster) ){
      nc <- length(cluster)
      environment(.funToApply) <- .GlobalEnv
      clusterExport(cluster, c("initB", "mSP", "mMU", "lam", "lsig", "qu", "argGam"), envir = environment())
    } else {
      nc <- 1
    }
    
    # Divide work (boostrap datasets) between cluster workers
    sched <- mapply(function(a, b) rep(a, each = b), 1:nc, 
                    c(rep(floor(nbo / nc), nc - 1), floor(nbo / nc) + nbo %% nc), SIMPLIFY = FALSE ) 
    sched <- split(1:nbo, do.call("c", sched))
    
    # Loop over bootstrap datasets to get standardized deviations from full data fit
    withCallingHandlers({
      out <- llply(.data = sched,
                   .fun = .funToApply,
                   .parallel = multicore,
                   .inform = control[["verbose"]],
                   .paropts = paropts#,
                   ### ... arguments start here
      ) 
    }, warning = function(w) {
      # There is a bug in plyr concerning a useless warning about "..."
      if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
        invokeRestart("muffleWarning")
    })
    
    loss <- .adTest( do.call("c", lapply(out, "[[", "z")) )
    initB <- unlist(lapply(out, "[[", "init"), recursive=FALSE)
    
    return( list("loss" = loss, "initM" = initM, "initB" = initB) )
  }
  
  init <- list("initM" = initM, "initB" = vector("list", nbo))
  
  # If we get convergence error, we increase "err" up to 0.2. If the error persists (or if the 
  # error is of another nature) we throw an error
  repeat{
    res <- tryCatch(.brent(brac=srange, f=objFun, mObj = mObj, bObj = bObj, init = init, 
                           pMat = pMat, qu = qu, varHat = varHat, cluster = cluster, t = control$tol), 
                    error = function(e) e)
    
    if("error" %in% class(res)){
      if( grepl("can't correct step size", res) ) {
        if(err < 0.2){
          err <- min(2*err, 0.2)
          if(control$verbose) message( paste("Increase \"err\" to ", err, " to get convergence") )
        } else {
          stop("I increased \"err\" up to 0.2, but still didn't get convergence.")
        }
      } else {
        stop( res )
      }
    } else { break } }
  
  res[["err"]] <- err
  
  return( res )
}





