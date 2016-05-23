#######
### Tune learning rate for quantile regression
#######
#' Learning rate tuning by calibration.
#' 
#' @param form a \code{gam} model formula.
#' @return output of \code{optimize}.
#'
#' @details XXX
#'         
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.   
#' @references Fasiolo and Wood XXX.           
#' @export
#' 

tuneLearnFast <- function(form, data, qu, err = 0.01, boot = NULL, 
                          multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                          control = list(), controlGam = list())
{ 
  n <- nrow(data)
  nq <- length(qu)
  
  # Setting up control parameter
  ctrl <- list( "init" = NULL,
                "brac" = log( c(1/5, 5) ), 
                "K" = 20,
                "tol" = .Machine$double.eps^0.25,
                "b" = 0,
                "verbose" = FALSE )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  tol <- ctrl[["tol"]]
  brac <- ctrl[["brac"]]
  
  # Sanity check
  if( tol > 0.1 * abs(diff(brac)) ) stop("tol > bracket_widths/10, choose smaller tolerance or larger bracket")
  
  # (Optional) create K boostrap dataset
  if( is.null(boot) ){
    tmp <- lapply(1:ctrl[["K"]], function(nouse) sample(1:n, n, replace = TRUE))
    boot <- lapply(tmp, function(ff) data[ff, ] )
  }
  
  # Gaussian fit, used for initialization
  if( is.formula(form) ) {
    fam <- "logF"
    gausFit <- gam(form, data = data, control = controlGam)
    varHat <- gausFit$sig2
  } else {
    fam <- "logFlss"
    gausFit <- gam(form, data = data, family = gaulss(b=ctrl[["b"]]), control = controlGam)
    varHat <- 1/gausFit$fit[ , 2]^2
  }  # Start = NULL in gamlss because it's not to clear how to deal with model for sigma 
  
  # Order quantiles so that those close to the median are dealt with first
  oQu <- order( abs(qu-0.5) )
  
  # (Optional) Initializing the search range for sigma
  if( is.null(ctrl[["init"]]) ){
    # We assume lam~0 and we match (5 times) the variance of a symmetric (median) Laplace density with that of the Gaussian fit.
    # This is an over-estimate for extreme quantiles, but experience suggests that it's better erring on the upper side.
    tmp <- 0.5 #qu[ oQu[1] ]
    if( !is.list(formula) ){
      isig <- log(sqrt( 5 * gausFit$sig2 * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
    } else {
      isig <- log(sqrt( 5 * (ctrl[["b"]]+exp(coef(gausFit)["(Intercept).1"]))^2 * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
    }
  } else {
    isig <- ctrl[["init"]]
  }
  
  # Create gam object for full data fits
  mObj <- gam(form, family = get(fam)(qu = NA, lam = NA, theta = NA), data = data, control = controlGam, fit = FALSE)
  
  # Create a gam object for each bootstrap sample
  bObj <- lapply(boot, function(bdat){
    out <- gam(form, family = get(fam)(qu = NA, lam = NA, theta = NA), data = bdat, 
               sp = gausFit$sp, control = controlGam, fit = FALSE)
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
                             control = ctrl, controlGam = controlGam)  
      
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
          # (unless the size of the bracket is < 10*tol)
          if( (abs(isig - mean(rans[kk, ])) < 0.25*wd) && (wd > 10*tol) ) brac <- brac / 2
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
      for(zz in 1:ii) segments(qu[oQu[zz]], rowMeans(rans)[oQu[zz]] - abs(diff(tmp[zz, ]))/4, 
                               qu[oQu[zz]], rowMeans(rans)[oQu[zz]] + abs(diff(tmp[zz, ]))/4, col = 1)
      plot(qu, efacts, xlab = "qu", "ylab" = "Bracket expansions")  
      plot(qu, errors)
    }
  }
  
  names(sigs) <- qu
  
  out <- list("lsigma" = sigs, "err" = errors, "ranges" = rans, "store" = store)
  
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
                           control, controlGam)
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
    mFit <- gam(G = mObj, in.out = initM[["in.out"]], start = initM[["start"]])
    mMU <- as.matrix(mFit$fit)[ , 1]
    mSP <- mFit$sp
    
    initM <- list("start" = coef(mFit), "in.out" = list("sp" = mSP, "scale" = 1))
    
    ## Function to be run in parallel (over boostrapped datasets)
    # GLOBAL VARS: bObj, pMat, initB, mSP, mMU, lam, lsig, qu
    .funToApply <- function(ind)
    {
      tmp <- lapply(ind, 
                    function(kk){
                      
                      .obj <- bObj[[kk]]; .pMat <- pMat[[kk]]; .init <- initB[[kk]]; 
                      
                      .obj$lsp0 <- log( mSP )
                      .obj$family$putQu( qu )
                      .obj$family$putLam( lam )
                      .obj$family$putTheta( lsig )
                      
                      fit <- gam(G = .obj, start = .init)
                      
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
                      sdev <- sqrt( diag( .pMat%*%Vp%*%t(.pMat) ) )
                      .z <- (mu - mMU) / sdev
                      
                      return( list("z" = .z, "init" = .init) )
                    })
      
      z <- as.vector(sapply(tmp, "[[", "z"))
      init <- lapply(tmp, "[[", "init")
      
      return( list("z" = z, "init" = init) )
    }
    
    if( !is.null(cluster) ){
      nc <- length(cluster)
      environment(.funToApply) <- .GlobalEnv
      clusterExport(cluster, c("initB", "mSP", "mMU", "lam", "lsig", "qu"), envir = environment())
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
    
    loss <- .adTest( as.vector(sapply(out, "[[", "z")) )
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





