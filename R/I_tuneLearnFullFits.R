###############
#### Internal function that does the full-data fits
###############
.tuneLearnFullFits <- function(lsig, form, fam, qu, err, ctrl, data, argGam, gausFit, varHat, initM){
  
  n <- nrow(data)
  nt <- length(lsig)
  
  # Create gam object for full data fits
  mainObj <- do.call("gam", c(list("formula" = form, 
                                   "family" = quote(elf(qu = qu, co = NA, theta = NA, link = ctrl$link)), 
                                   "data" = quote(data), "fit" = FALSE), 
                              argGam))
  
  # Remove "sp" as it is already been fixed
  argGam <- argGam[ names(argGam) != "sp" ]
  
  # Create reparametrization list needed for sandwich calibration
  repar <- .prepBootObj(obj = mainObj, eps = NULL, control = argGam$control)[ c("UrS", "Mp", "U1") ]

  # Store degrees of freedom for each value of lsig
  tmp <- pen.edf( gausFit )
  if( length(tmp) ) {
    edfStore <- list( )
  } else {
    edfStore <- NULL
  } 
  
  # FULL data fits, used to estimate the smoothing parameters 
  store <- vector("list", nt)
  for( ii in 1:nt ) # START lsigma loop, from smallest to largest (because when lsig is large the smooth params diverge)
  {
    mainObj$family$putCo( err * sqrt(2*pi*varHat) / (2*log(2)) )
    mainObj$family$putTheta( lsig[ii] )
    
    convProb <- FALSE # Variable indicating convergence problems
    withCallingHandlers({
      fit <- do.call("gam", c(list("G" = quote(mainObj), "in.out" = initM[["in.out"]], "mustart" = initM[["mustart"]]), argGam)) 
    }, warning = function(w) {
      if (length(grep("Fitting terminated with step failure", conditionMessage(w))) ||
          length(grep("Iteration limit reached without full convergence", conditionMessage(w))))
      {
        message( paste("log(sigma) = ", round(lsig[ii], 3), " : outer Newton did not converge fully.", sep = "") )
        convProb <<- TRUE
        invokeRestart("muffleWarning")
      }
    })
  
    if( !is.null(edfStore) ) { 
      edfStore[[ii]] <- c(lsig[ii], pen.edf(fit)) 
    }
    
    # Create prediction matrix (only in the first iteration)
    if( ii == 1 ){
      pMat <-  predict.gam(fit, type = "lpmatrix") 
    }
    
    sdev <- NULL 
    if(ctrl$loss %in% c("cal", "calFast") && ctrl$vtype == "m"){
      Vp <- fit$Vp
      sdev <- sqrt(rowSums((pMat %*% Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
    }
    
    initM <- list("mustart" = fit$fitted.values, "in.out" = list("sp" = fit$sp, "scale" = 1))
    
    if( ctrl$loss == "calFast" ){ # Fast calibration OR ...
      if( ii == 1 ){
        EXXT <- crossprod(pMat, pMat) / n                       # E(xx^T)
        EXEXT <- tcrossprod( colMeans(pMat), colMeans(pMat) )   # E(x)E(x)^T
      }
      Vbias <- .biasedCov(fit = fit, X = pMat, EXXT = EXXT, EXEXT = EXEXT, mObj = mainObj)
      
      store[[ii]] <- list("loss" = .sandwichLoss(mFit = fit, X = pMat, sdev = sdev, repar = repar, 
                                                 alpha = Vbias$alpha, VSim = Vbias$V, mObj = mainObj), 
                          "convProb" = convProb)
    } else { # Bootstrapping or cross-validation: full data fit will be used when fitting the bootstrap datasets
      store[[ii]] <- list("sp" = fit$sp, "fit" = fit$fitted, "co" = fit$family$getCo(), 
                          "init" = initM$start, "sdev" = sdev, "weights" = fit$working.weights, 
                          "res" = fit$residuals, "convProb" = convProb)
    }
    
  } 
  
  if( !is.null(edfStore) ){
    edfStore <- do.call("rbind", edfStore)
    colnames(edfStore) <- c("lsig", names( pen.edf(fit) ))
  }
  
  return( list("store" = store, 
               "edfStore" = edfStore, 
               "pMat" = if(ctrl$loss != "calFast") { pMat } else { NULL } ) )

}
