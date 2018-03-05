#######
# Internal loss function to be minimized using Brent method
#
.objFunLearnFast <- function(lsig, mObj, bObj, wb, initM, initB, pMat, SStuff, qu, ctrl, varHat, 
                             err, argGam, cluster, multicore, paropts)
{ 
  if(ctrl$progress){ cat(".")}
  
  co <- err * sqrt(2*pi*varHat) / (2*log(2))
  lpi <- attr(pMat, "lpi")
  
  mObj$family$putQu( qu )
  mObj$family$putCo( co )
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
  initM <- list("start" = coef(mFit), "in.out" = list("sp" = mFit$sp, "scale" = 1))
  
  # Standard deviation of fitted quantile using full data
  sdev <- NULL 
  if(ctrl$loss %in% c("cal", "calFast") && ctrl$vtype == "m"){
    Vp <- mFit$Vp
    # In the gamlss case, we are interested only in the calibrating the location mode
    if( !is.null(lpi) ){  Vp <- mFit$Vp[lpi[[1]], lpi[[1]]]  }
    sdev <- sqrt(rowSums((pMat %*% Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
  }
  
  if(ctrl$loss == "calFast"){ # Fast calibration OR ...
   
    Vbias <- .biasedCov(fit = mFit, X = SStuff$XFull, EXXT = SStuff$EXXT, EXEXT = SStuff$EXEXT, lpi = lpi)
    outLoss <- .sandwichLoss(mFit = mFit, X = pMat, XFull = SStuff$XFull, sdev = sdev, repar = mObj$hidRepara, 
                             alpha = Vbias$alpha, VSim = Vbias$V)
    initB <- NULL
    
  } else { # ... bootstrapping or cross-validation
    ## Function to be run in parallel (over boostrapped datasets)  
    # It has two sets of GLOBAL VARS
    # Set 1: bObj, pMat, wb, argGam, ctrl    (Exported by tuneLearnFast)
    # If multicore=F, .funToApply() will look for these inside the objFun call. That's why objFun need them as arguments.
    # If multicore=T, .funToApply() will look for them in .GlobalEnv. That's why we export them to cluster nodes in tuneLearnFast.
    # Set 2:  initB, initM, mMU, co, lsig, qu, sdev   (Exported by .tuneLearnFast)
    # As before but, if multicore=T, these are exported directly by objFun because they change from one call of objFun to another.
    .funToApply <- function(ind)
    {
      .lpi <- attr(pMat, "lpi")
      glss <- inherits(bObj$family, "general.family")
      
      bObj$lsp0 <- log( initM$in.out$sp )
      bObj$family$putQu( qu )
      bObj$family$putCo( co )
      bObj$family$putTheta( lsig )
      
      z <- init <- vector("list", length(ind))
      for( kk in ind){
        # Creating boot weights from boot indexes 
        .wb <- wb[[kk]]
        bObj$w <- .wb
        
        # Recycle boot initialization, but at first iteration this is NULL... 
        .init <- if(is.null(initB[[kk]])){ list(initM$start) } else { list(initB[[kk]], initM$start) }
        
        if( glss ){ # In gamlss I need to reparametrize initialization and in Ex GAM I need to get null coefficients.
          .init <- lapply(.init, function(inp) Sl.initial.repara(bObj$Sl, inp, inverse=FALSE, both.sides=FALSE))
          .fit <- .gamlssFit(x=bObj$X, y=bObj$y, lsp=as.matrix(bObj$lsp0), Sl=bObj$Sl, weights=bObj$w, 
                             offset=bObj$offset, family=bObj$family, control=bObj$control, 
                             Mp=bObj$Mp, start=.init, needVb=(ctrl$loss=="cal" && ctrl$vtype=="b"))
          
          # In gamlss, we want to calibrate only the location and we need to reparametrize the coefficients
          .init <- .betas <- Sl.initial.repara(bObj$Sl, .fit$coef, inverse=TRUE, both.sides=FALSE)
          .betas <- .betas[.lpi[[1]]]
        } else {
          bObj$null.coef <- bObj$family$get.null.coef(bObj)$null.coef
          .fit <- .egamFit(x=bObj$X, y=bObj$y, sp=as.matrix(bObj$lsp0), Eb=bObj$Eb, UrS=bObj$UrS,
                           offset=bObj$offset, U1=bObj$U1, Mp=bObj$Mp, family = bObj$family, weights=bObj$w,
                           control=bObj$control, null.coef=bObj$null.coef, 
                           start=.init, needVb=(ctrl$loss == "cal" && ctrl$vtype == "b"))
          .init <- .betas <- .fit$coef
        }
        
        .mu <- pMat %*% .betas
        
        if( ctrl$loss == "cal" ){ # (1) Return standardized deviations from full data fit OR ... 
          if( ctrl$vtype == "b" ){ # (2) Use variance of bootstrap fit OR ...
            .Vp <- .getVp(.fit, bObj, bObj$lsp0, .lpi)
            .sdev <- sqrt(rowSums((pMat %*% .Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
          } else { # (2)  ... variance of the main fit
            .sdev <- sdev
          }
          z[[kk]] <- (.mu - mMU) / .sdev
        } else { # (1) ... out of sample observations minus their fitted values 
          z[[kk]] <- (bObj$y - .mu)[ !.wb ]
        }
        
        init[[kk]] <- .init
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