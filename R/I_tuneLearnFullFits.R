###############
#### Internal function that does the full-data fits
###############
.tuneLearnFullFits <- function(lsig, form, fam, qu, err, ctrl, data, argGam, gausFit, varHat, initM){
  
  n <- nrow(data)
  nt <- length(lsig)
  
  # Create gam object for full data fits
  mainObj <- do.call("gam", c(list("formula" = form, 
                                   "family" = get(fam)(qu = qu, co = NA, theta = NA, link = ctrl$link), 
                                   "data" = data, "fit" = FALSE), 
                              argGam))
  
  # Create reparametrization list for... 
  repar <- if( is.list(form) ){ # ... GAMLSS case OR...
    Sl.setup( mainObj )
  } else { # ... extended GAM case
    .prepBootObj(obj = mainObj, eps = NULL, control = argGam$control)[ c("UrS", "Mp", "U1") ]
  } # these are needed for sandwich calibration
  
  # Store degrees of freedom for each value of lsig
  tmp <- pen.edf(gausFit)
  if( length(tmp) )
  {
    edfStore <- matrix(NA, nt, length(tmp) + 1)
    colnames(edfStore) <- c("lsig", names( tmp ))
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
      fit <- do.call("gam", c(list("G" = mainObj, "in.out" = initM[["in.out"]], "start" = initM[["start"]]), argGam)) 
    }, warning = function(w) {
      if (length(grep("Fitting terminated with step failure", conditionMessage(w))) ||
          length(grep("Iteration limit reached without full convergence", conditionMessage(w))))
      {
        message( paste("log(sigma) = ", round(lsig[ii], 3), " : outer Newton did not converge fully.", sep = "") )
        convProb <<- TRUE
        invokeRestart("muffleWarning")
      }
    })
    
    if( !is.null(edfStore) ) { edfStore[ii, ] <- c(lsig[ii], pen.edf(fit)) }
    
    # Create prediction matrix (only in the first iteration)
    if( ii == 1 ){
      pMat <- pMatFull <- predict.gam(fit, type = "lpmatrix") 
      lpi <- attr(pMat, "lpi")
      if( !is.null(lpi) ){ 
        pMat <- pMat[ , lpi[[1]]] # "lpi" attribute lost here
        attr(pMat, "lpi") <- lpi
      }
    }
    
    sdev <- NULL 
    if(ctrl$loss %in% c("cal", "calFast") && ctrl$vtype == "m"){
      Vp <- fit$Vp
      # In the gamlss case, we are interested only in the calibrating the location mode
      if( !is.null(lpi) ){  Vp <- fit$Vp[lpi[[1]], lpi[[1]]]  }
      sdev <- sqrt(rowSums((pMat %*% Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
    }
    
    initM <- list("start" = coef(fit), "in.out" = list("sp" = fit$sp, "scale" = 1))
    
    if( ctrl$loss == "calFast" ){ # Fast calibration OR ...
      store[[ii]] <- list("loss" = .sandwichLoss(mFit = fit, X = pMat, XFull = pMatFull, sdev = sdev, repar = repar), 
                          "convProb" = convProb)
    } else { # Bootstrapping or cross-validation: full data fit will be used when fitting the bootstrap datasets
      store[[ii]] <- list("sp" = fit$sp, "fit" = fit$fitted, "co" = fit$family$getCo(), 
                          "init" = initM$start, "sdev" = sdev, "weights" = fit$working.weights, 
                          "res" = fit$residuals, "convProb" = convProb)
    }
    
  } 
  
  return( list("store" = store, 
               "edfStore" = edfStore, 
               "pMat" = if(ctrl$loss != "calFast") { pMat } else { NULL } ) )

}
