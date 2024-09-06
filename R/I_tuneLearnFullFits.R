###############
#### Internal function that does the full-data fits
###############
.tuneLearnFullFits <- function(lsig, form, fam, qu, err, ctrl, data, discrete, argGam, gausFit, varHat, initM){
  
  n <- nrow(data)
  nt <- length(lsig)
  
  gam_name <- ifelse(discrete, "bam", "gam")
  
  # Create gam object for full data fits
  mObj <- do.call(gam_name, c(list("formula" = form, 
                                      "family" = quote(elf(qu = qu, co = NA, theta = NA, link = ctrl$link)), 
                                      "data" = quote(data), "discrete" = discrete, "fit" = FALSE), 
                                 argGam))
  
  # Remove "sp" as it is already been fixed
  argGam <- argGam[ names(argGam) != "sp" ]
  
  # Preparing reparametrization list and hide it within mObj. This will be needed by the sandwich calibration
  if( ctrl$loss == "calFast" ){
    mObj$hidRepara <- .prepBootObj(obj = mObj, eps = NULL, control = argGam$control)[ c("UrS", "Mp", "U1") ]
  }

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
    mObj$family$putCo( err * sqrt(2*pi*varHat) / (2*log(2)) )
    mObj$family$putTheta( lsig[ii] )
    
    convProb <- FALSE # Variable indicating convergence problems
    withCallingHandlers({
      call_list <- c(list("G" = quote(mObj), "in.out" = initM[["in.out"]], "mustart" = initM[["mustart"]], "discrete" = discrete), argGam)
      # Annoyingly, initial coeffs are supplied via "coef" argument in bam() and "start" in gam()
      call_list[[ ifelse(discrete, "coef", "start") ]] <- initM$coefstart 
      mFit <- do.call(gam_name, call_list) 
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
      edfStore[[ii]] <- c(lsig[ii], pen.edf(mFit)) 
    }
    
    # Create prediction matrix (only in the first iteration)
    if( ii == 1 ){
      # Stuff needed for the sandwich estimator
      if(discrete){
        pMat <- NULL
        colmX <- XWyd(X=mObj$Xd,w=rep(1,n),y=rep(1,n),k=mObj$kd,ks=mObj$ks,
                      ts=mObj$ts,dt=mObj$dt,v=mObj$v,qc=mObj$qc,drop=mObj$drop)/n
        sandStuff <- list("XFull" = pMat,
                          "EXXT" = XWXd(X=mObj$Xd,w=rep(1,n),k=mObj$kd,ks=mObj$ks,
                                        ts=mObj$ts,dt=mObj$dt,v=mObj$v,qc=mObj$qc,
                                        drop=mObj$drop)/n,  # E(xx^T)
                          "EXEXT" = tcrossprod( colmX, colmX)) # E(x)E(x)^T
      }else{
        pMat <- predict.gam(mFit, type = "lpmatrix")
        sandStuff <- list("XFull" = pMat,
                          "EXXT" = crossprod(pMat, pMat) / n,                    # E(xx^T)
                          "EXEXT" = tcrossprod( colMeans(pMat), colMeans(pMat))) # E(x)E(x)^T
      }
    }
    
    # Standard deviation of fitted quantile using full data
    sdev <- NULL 
    if(ctrl$loss %in% c("cal", "calFast") && ctrl$vtype == "m"){
      Vp <- mFit$Vp
      if(discrete){
        sdev <- diagXVXd(X=mObj$Xd,V=Vp,k=mObj$kd,ks=mObj$ks,ts=mObj$ts,
                         dt=mObj$dt,v=mObj$v,qc=mObj$qc,drop=mObj$drop)^.5
      }else{
        sdev <- sqrt(rowSums((pMat %*% Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
      }
    }
    
    initM <- list("mustart" = mFit$fitted.values, 
                  "coefstart" = coef(mFit),
                  "in.out" = list("sp" = if(gam_name == "bam" & !is.null(mFit$full.sp)){ mFit$full.sp } else { mFit$sp }, "scale" = 1))
    
    if( ctrl$loss == "calFast" ){ # Fast calibration OR ...
      Vbias <- .biasedCov(fit = mFit, X = sandStuff$XFull, EXXT = sandStuff$EXXT, EXEXT = sandStuff$EXEXT, mObj = mObj)
      outLoss <- .sandwichLoss(mFit = mFit, X = pMat, sdev = sdev, repar = mObj$hidRepara, 
                               alpha = Vbias$alpha, VSim = Vbias$V, mObj = mObj)
      store[[ii]] <- list("loss" = outLoss, "convProb" = convProb)
    } else { # Bootstrapping or cross-validation: full data fit will be used when fitting the bootstrap datasets
      store[[ii]] <- list("sp" = mFit$sp, "fit" = mFit$fitted, "co" = mFit$family$getCo(), 
                          "init" = initM$coefstart, "sdev" = sdev, "weights" = mFit$working.weights, 
                          "res" = mFit$residuals, "convProb" = convProb)
    }
    
  } 
  
  if( !is.null(edfStore) ){
    edfStore <- do.call("rbind", edfStore)
    colnames(edfStore) <- c("lsig", names( pen.edf(mFit) ))
  }
  
  return( list("store" = store, 
               "edfStore" = edfStore, 
               "pMat" = if(ctrl$loss != "calFast") { pMat } else { NULL } ) )

}
