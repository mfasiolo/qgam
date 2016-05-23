####### Tuning the learning rate for Gibbs posterior

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
    # Force evaluation of everything in the environment, so it will available to .getBootDev on cluster
    #.forceEval(ALL = TRUE)
    
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
  for( ii in nt:1 ) # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harded)
  {
    mainObj$family$putLam( err * sqrt(2*pi*varHat) / (2*log(2)*exp(lsig[ii])) )
    mainObj$family$putTheta( lsig[ii] )
    
    withCallingHandlers({
    fit <- gam(G = mainObj, in.out = initM[["in.out"]], start = initM[["start"]])}, warning = function(w) {
      # There is a bug in plyr concerning a useless warning about "..."
      if (length(grep("Fitting terminated with step failure", conditionMessage(w))) ||
          length(grep("Iteration limit reached without full convergence", conditionMessage(w))))
      {
        message( paste("log(sigma) = ", round(lsig[ii], 6), " : outer Newton did not converge.", sep = "") )
        invokeRestart("muffleWarning")
      }
    })
    
    initM <- list("start" = coef(fit), "in.out" = list("sp" = fit$sp, "scale" = 1))
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

  z <- do.call("cbind", z)

  loss <- apply(z, 1, function(.x) .adTest(as.vector(.x)))
  names(loss) <- lsig
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return( loss )
}



###########
# Function that fits the bootstrapped datasets and returns standardized deviations from full data fit. 
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
