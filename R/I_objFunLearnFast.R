#######
# Internal loss function to be minimized using Brent method
#
.objFunLearnFast <- function(lsig, mObj, bObj, wb, initM, initB, pMat, SStuff, qu, ctrl, varHat, 
                             err, argGam, cluster, multicore, paropts)
{ 
  if(ctrl$progress){ cat(".")}
  
  discrete <- !is.null(mObj$Xd)
  gam_name <- ifelse(discrete, "bam", "gam")

  co <- err * sqrt(2*pi*varHat) / (2*log(2))

  mObj$family$putQu( qu )
  mObj$family$putCo( co )
  mObj$family$putTheta( lsig )
  
  
  call_list <- c(list("G" = quote(mObj), "in.out" = initM[["in.out"]], "mustart" = initM[["mustart"]], "discrete" = discrete), argGam)
  
  # Annoyingly, initial coeffs are supplied via "coef" argument in bam() and "start" in gam()
  call_list[[ ifelse(discrete, "coef", "start") ]] <- initM$coefstart 

  withCallingHandlers({ mFit <- do.call(gam_name, call_list) }, warning = function(w) {
      if (length(grep("Fitting terminated with step failure", conditionMessage(w))) ||
          length(grep("Iteration limit reached without full convergence", conditionMessage(w))))
      {
        message( paste("qu = ", qu, ", log(sigma) = ", round(lsig, 6), " : outer Newton did not converge fully.", sep = "") )
        invokeRestart("muffleWarning")
      }
    })
  
  mMU <- mFit$fit
  initM <- list("mustart" = mFit$fitted.values, 
                "coefstart" = coef(mFit), 
                "in.out" = list("sp" = if(gam_name == "bam" && !is.null(mFit$full.sp)){ mFit$full.sp } else { mFit$sp }, "scale" = 1))
  
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
  
  if(ctrl$loss == "calFast"){ # Fast calibration OR ...
   
    Vbias <- .biasedCov(fit = mFit, X = SStuff$XFull, EXXT = SStuff$EXXT, EXEXT = SStuff$EXEXT, mObj = mObj)
    outLoss <- .sandwichLoss(mFit = mFit, X = pMat, sdev = sdev, repar = mObj$hidRepara, 
                             alpha = Vbias$alpha, VSim = Vbias$V, mObj = mObj)
    initB <- NULL
    
  } else { # ... bootstrapping or cross-validation
    if(discrete){ stop("discrete = TRUE allowed only when \"control$loss == calFast\"")}
    
    ## Function to be run in parallel (over boostrapped datasets)  
    # It has two sets of GLOBAL VARS
    # Set 1: bObj, pMat, wb, argGam, ctrl    (Exported by tuneLearnFast)
    # If multicore=F, .funToApply() will look for these inside the objFun call. That's why objFun need them as arguments.
    # If multicore=T, .funToApply() will look for them in .GlobalEnv. That's why we export them to cluster nodes in tuneLearnFast.
    # Set 2:  initB, initM, mMU, co, lsig, qu, sdev   (Exported by .tuneLearnFast)
    # As before but, if multicore=T, these are exported directly by objFun because they change from one call of objFun to another.
    .funToApply <- function(ind)
    {
      bObj$lsp0 <- log( initM$in.out$sp )
      bObj$family$putQu( qu )
      bObj$family$putCo( co )
      bObj$family$putTheta( lsig )
      
      z <- init <- vector("list", length(ind))
      for( ii in 1:length(ind) ){
        # Creating boot weights from boot indexes 
        kk <- ind[ ii ]
        .wb <- wb[[ kk ]]
        bObj$w <- .wb
        
        # Recycle boot initialization, but at first iteration this is NULL... 
        .init <- if(is.null(initB[[kk]])){ list(initM$coefstart) } else { list(initB[[kk]], initM$coefstart) }
        
        # I need to get null coefficients.
        bObj$null.coef <- bObj$family$get.null.coef(bObj)$null.coef
        .fit <- .egamFit(x=bObj$X, y=bObj$y, sp=as.matrix(bObj$lsp0), Eb=bObj$Eb, UrS=bObj$UrS,
                         offset=bObj$offset, U1=bObj$U1, Mp=bObj$Mp, family = bObj$family, weights=bObj$w,
                         control=bObj$control, null.coef=bObj$null.coef, 
                         start=.init, needVb=(ctrl$loss == "cal" && ctrl$vtype == "b"))
        .init <- .betas <- .fit$coef
        
        
        .mu <- pMat %*% .betas
        
        if( ctrl$loss == "cal" ){ # (1) Return standardized deviations from full data fit OR ... 
          if( ctrl$vtype == "b" ){ # (2) Use variance of bootstrap fit OR ...
            .Vp <- .fit$Vp
            .sdev <- sqrt(rowSums((pMat %*% .Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
          } else { # (2)  ... variance of the main fit
            .sdev <- sdev
          }
          z[[ii]] <- (.mu - mMU) / .sdev
        } else { # (1) ... out of sample observations minus their fitted values 
          z[[ii]] <- (bObj$y - .mu)[ !.wb ]
        }
        
        init[[ii]] <- .init
      }
      
      return( list("z" = z, "init" = init) )
    } 
    
    if( !is.null(cluster) ){
      nc <- length(cluster)
      environment(.funToApply) <- .GlobalEnv
      clusterExport(cluster, c("initB", "initM", "mMU", "co", "lsig", "qu", "sdev"), envir = environment())
    } else {
      nc <- 1
    }
    
    # Divide work (boostrap datasets) between cluster workers
    nbo <- ctrl$K
    sched <- mapply(function(a, b) rep(a, each = b), 1:nc, 
                    c(rep(floor(nbo / nc), nc - 1), floor(nbo / nc) + nbo %% nc), SIMPLIFY = FALSE ) 
    sched <- split(1:nbo, do.call("c", sched))
    
    # Loop over bootstrap datasets to get standardized deviations from full data fit
    withCallingHandlers({
      out <- llply(.data = sched,
                   .fun = .funToApply,
                   .parallel = multicore,
                   .inform = ctrl[["verbose"]],
                   .paropts = paropts#,
                   ### ... arguments start here
      ) 
    }, warning = function(w) {
      # There is a bug in plyr concerning a useless warning about "..."
      if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
        invokeRestart("muffleWarning")
    })
    
    # Get stardardized deviations and ... 
    .bindFun <- if( ctrl$loss == "cal" ) { "rbind" } else { "c" }
    z <- do.call(.bindFun, do.call(.bindFun, lapply(out, "[[", "z")))
    if( ctrl$loss == "cal"){ # ... calculate KL distance OR ...
      vrT <- .colVars(z)
      outLoss <- mean(vrT + colMeans(z)^2 - log(vrT))
    } else { # ... pinball loss
      outLoss <- .checkloss(as.vector(z), 0, qu = qu)
    }
    names(outLoss) <- lsig
    
    initB <- unlist(lapply(out, "[[", "init"), recursive=FALSE)
    
  } 
  
  return( list("outLoss" = outLoss, "initM" = initM, "initB" = initB) )
} 