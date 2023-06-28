###############
#### Internal function that does the bootstrapping or cross-validation
###############
.tuneLearnBootstrapping <- function(lsig, form, fam, qu, ctrl, data, store, pMat, argGam, 
                                    multicore, cluster, ncores, paropts){
  
  n <- nrow(data)
  nt <- length(lsig)
  
  if( ctrl$sam == "boot" ){ # Create weights for K boostrap dataset OR...
    wb <- lapply(1:ctrl[["K"]], function(nouse) tabulate(sample(1:n, n, replace = TRUE), n))
  } else { # ... OR for K training sets for CV 
    tmp <- sample(rep(1:ctrl[["K"]], length.out = n), n, replace = FALSE)
    wb <- lapply(1:ctrl[["K"]], function(ii) tabulate(which(tmp != ii), n)) 
  }
  
  # Create gam object for bootstrap fits
  bObj <- do.call("gam", c(list("formula" = form, "family" = quote(elf(qu = qu, co = NA, theta = NA, link = ctrl$link)), "data" = quote(data), 
                                "sp" = if(length(store[[1]]$sp)){store[[1]]$sp}else{NULL}, fit = FALSE), argGam))
  
  # Preparing bootstrap object for gam.fit3
  bObj <- .prepBootObj(obj = bObj, eps = ctrl$epsB, control = argGam$control)
  
  # Internal function that fits the bootstrapped datasets and returns standardized deviations from full data fit. To be run in parallel.
  # GLOBALS: lsig, ctrl, store, pMat, bObj, argGam
  .getBootDev <- function(.wb)
  {   # # # # # # # # # .getBootDev START # # # # # # # # #
    y <- bObj$y 
    ns <- length(lsig); n  <- length(y)
    
    # Number of test observations
    nt <- ifelse(ctrl$loss == "cal", n, sum(!.wb)) 
    
    # Creating boot weights from boot indexes 
    bObj$w <- .wb
    
    init <- NULL
    .z <- vector("list", ns)
    for( ii in ns:1 )  # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harded)
    {   
      # In the gamlss case 'co' is a vector, and we have to take only those values that are in the boostrapped dataset.
      co <- store[[ii]]$co
      
      bObj$lsp0 <- log( store[[ii]]$sp )
      bObj$family$putCo( co )
      bObj$family$putTheta( lsig[ii] )
      
      init <- if(is.null(init)){ list(store[[ii]]$init) } else { list(init, store[[ii]]$init) }
      
      .offset <- bObj$offset
      
      bObj$null.coef <- bObj$family$get.null.coef(bObj)$null.coef
        fit <- .egamFit(x=bObj$X, y=bObj$y, sp=as.matrix(bObj$lsp0), Eb=bObj$Eb, UrS=bObj$UrS,
                        offset=bObj$offset, U1=bObj$U1, Mp=bObj$Mp, family = bObj$family, weights=bObj$w,
                        control=bObj$control, null.coef=bObj$null.coef, 
                        start=init, needVb=(ctrl$loss == "cal" && ctrl$vtype == "b"))
        init <- betas <- fit$coef
        
        if( is.null(.offset) ){ .offset <- numeric( nrow(pMat) )  }
        mu <- bObj$family$linkinv( pMat %*% betas + .offset )
    
      if( ctrl$loss == "cal" ){ # (1) Return standardized deviations from full data fit OR ... 
        if( ctrl$vtype == "b" ){ # (2) Use variance of bootstrap fit OR ...
          Vp <- fit$Vb
          sdev <- sqrt(rowSums((pMat %*% Vp) * pMat)) # same as sqrt(diag(pMat%*%Vp%*%t(pMat))) but (WAY) faster
        } else { # (2)  ... variance of the main fit
          sdev <- store[[ii]]$sdev
        }
        .z[[ii]] <- drop((mu - as.matrix(store[[ii]]$fit)[ , 1]) / sdev)
      } else { # (1) ... out of sample observations minus their fitted values 
        .z[[ii]] <- drop(y - mu)[ !.wb ]
      }
    }
    return( .z )
  }  # # # # # # # # # .getBootDev END # # # # # # # # #
  
  if( multicore ){ 
    # Making sure "qgam" is loaded on cluser
    paropts[[".packages"]] <- unique( c("qgam", paropts[[".packages"]]) )
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores) #, exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoParallel(cluster)
    
    # Exporting stuff. To about all environment being exported all the time, use .GlobalEnv  
    clusterExport(cluster, c("pMat", "bObj", "lsig", "ctrl", "store", "argGam", ".egamFit"), 
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
  .bindFun <- if( ctrl$loss == "cal" ) { "rbind" } else { "c" }
  z <- lapply(1:nt, function(.ii) do.call(.bindFun, lapply(z, function(.x) .x[[.ii]])))
  
  if( ctrl$loss == "cal"){ # ... calculate KL distance OR ...
    # KL distance for explanations see [*] below
    outLoss <- sapply(z, function(.x){ 
      .v <- .colVars(.x)
      return( mean( sqrt(.v + colMeans(.x)^2 - log(.v)) ) )
    })
    # E(z^2) = var(z) + E(z)^2 (var + bias)
    #outLoss <- sapply(z, function(.x) mean( (.colVars(.x) - 1)^2  + colMeans(.x)^2 ) )
    #outLoss <- sapply(z, function(.x) mean( colMeans(.x)^2 ) )  
    #outLoss <- sapply(z, function(.x) mean( (colMeans(.x^2) - 1)^2 ) ) 
    #outLoss <- sapply(z, function(.x) mean( apply(.x, 2, .adTest) ) ) # .adTest(as.vector(.x)))
    #outLoss <- sapply(z, function(.x) .adTest(as.vector(.x)) ) # .adTest(as.vector(.x)))
  } else { # ... pinball loss
    outLoss <- sapply(z, function(.x) .checkloss(.x, 0, qu = qu))
  }

  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return( outLoss )
  
}


##### [*] Code showing that KL distance is invariant to standardization
# mu1 <- rnorm(1)
# mu2 <- rnorm(1)
# v1 <- runif(1, 1, 2)
# v2 <- runif(1, 1, 2)
# 
# x <- rnorm(10000, mu1, sqrt(v1))
# 
# # KL distance between x ~ N(mu1, V1) and z ~ N(mu2, V2)
# v1/v2 + (mu1 - mu2)^2 / v2 + log(v2/v1)
# 
# # Empirical estimate of KL distance
# var(x)/v2 + (mean(x) - mu2)^2 / v2 + log(v2/var(x))
# 
# # Normalizing x using mu2 and V2, assume y is now N(0, 1)
# # and recalculate KL distance: the result must be the same
# y <- (x - mu2) / sqrt(v2)
# var(y) + (mean(y))^2 + log(1/var(y))
