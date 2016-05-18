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
  
  # Setting up control parameter
  ctrl <- list( "init" = NULL,
                "brac" = log( c(1/5, 5) ), 
                "K" = 20,
                "gausFit" = NULL,
                "tol" = .Machine$double.eps^0.25,
                "b" = 0,
                "verbose" = FALSE )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  if( multicore ){ 
    # Making sure "qgam" is loaded on clutser and collate user-specified packages
    paropts[[".packages"]] <- unique( c("qgam", paropts[[".packages"]]) )
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores) #, exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
  }
  
  tol <- ctrl[["tol"]]
  brac <- ctrl[["brac"]]
  
  # Sanity check
  if( tol > 0.1 * abs(diff(brac)) ) stop("tol > bracket_widths/10, choose smaller tolerance or larger bracket")
  
  # (Optional) create K boostrap dataset
  if( is.null(boot) ){
    tmp <- lapply(1:ctrl[["K"]], function(nouse) sample(1:n, n, replace = TRUE))
    boot <- lapply(tmp, function(ff) data[ff, ] )
  }
  
  # (Optional) Main Gaussian fit, used for initializations
  if( is.null(ctrl[["gausFit"]]) ) {
    if( plyr:::is.formula(form) ) {
      ctrl[["gausFit"]] <- gam(form, data = data, control = controlGam)
    } else {
      ctrl[["gausFit"]] <- gam(form, data = data, family = gaulss(b=ctrl[["b"]]), control = controlGam)
    } }
  
  # Order quantiles so that those close to the median are dealt with first
  oQu <- order( abs(qu-0.5) )
  
  # (Optional) Initializing the search range for sigma
  if( is.null(ctrl[["init"]]) ){
    # We assume lam~0 and we match (5 times) the variance of a symmetric (median) Laplace density with that of the Gaussian fit.
    # This is an over-estimate for extreme quantiles, but experience suggests that it's better erring on the upper side.
    tmp <- 0.5 #qu[ oQu[1] ]
    if( !is.list(formula) ){
      isig <- log(sqrt( 5 * ctrl$gausFit$sig2 * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
    } else {
      isig <- log(sqrt( 5 * (ctrl[["b"]]+exp(coef(ctrl$gausFit)["(Intercept).1"]))^2 * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
    }
  } else {
    isig <- ctrl[["init"]]
  }
  
  # Initialize
  nq <- length(qu)
  
  # Estimated learning rates, # of bracket expansions, error rates and bracket ranges used in bisection
  sigs <- efacts <- errors <- numeric(nq)
  rans <- matrix(NA, nq, 2)
  store <- vector("list", nq)
  names(sigs) <- names(errors) <- rownames(rans) <- qu
  
  # Here we need bTol > aTol otherwise we the new bracket will be too close to probable solution
  aTol <- 0.05
  bTol <- 0.2
  
  for(ii in 1:nq)
  {
    oi <- oQu[ii]
    
    ef <- 1
    
    repeat{
      
      # Compute bracket
      srange <- isig + ef * brac
      
      # Estimate log(sigma)
      res  <- .tuneLearnFast(form = form, data = data, qu = qu[oi], err = err, boot = boot, srange = srange, 
                             multicore = multicore, cluster = cluster, ncores = ncores, paropts = paropts,  
                             control = ctrl, controlGam = controlGam)  
      
      store[[oi]] <- cbind(store[[oi]], res[["store"]])
      lsig <- res$minimum
      
      # If solution not too close to boundary store results and determine bracket for next iter
      if( all(abs(lsig-srange) > aTol * abs(diff(srange))) ){ 
        
        sigs[oi] <- lsig
        rans[oi, ] <- srange
        efacts[oi] <- ef
        errors[oi] <- res$err
        
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
### Internal version, which works with scalar qu
########################################################################## 
.tuneLearnFast <- function(form, data, boot, qu, err, srange, 
                           multicore, cluster, ncores, paropts, 
                           control, controlGam)
{
  n <- nrow(data)
  nbo <- length(boot)
  
  gausFit <- control$gausFit
  
  # Gaussian fit, used for initializations 
  if( is.formula(form) ) {
    fam <- "logF"
    varHat <- gausFit$sig2
    initM <- list("start" = coef(gausFit) + c(qnorm(qu, 0, sqrt(gausFit$sig2)), rep(0, length(coef(gausFit))-1)), 
                  "in.out" = list("sp" = gausFit$sp, "scale" = 1)) 
  } else {
    fam <- "logFlss"
    varHat <- 1/gausFit$fit[ , 2]^2
    initM <- list("start" = NULL, "in.out" = list("sp" = gausFit$sp, "scale" = 1)) 
  }  # Start = NULL in gamlss because it's not to clear how to deal with model for sigma 
  
  # Create gam object for full data fits
  mObj <- gam(form, family = get(fam)(qu = qu, lam = NA, theta = NA), data = data, control = controlGam, fit = FALSE)
  
  # Create an object for each bootstrap sample
  bObj <- lapply(boot, function(bdat){
    out <- gam(form, family = get(fam)(qu = qu, lam = NA, theta = NA), data = bdat, 
               sp = gausFit$sp, control = controlGam, fit = FALSE)
    return( out )
  })
  
  # Objective function to be minimized
  # GLOBALS: varHat 
  objFun <- function(lsig, mObj, bObj, initM, initB, pMat)
  {
    lam <- err * sqrt(2*pi*varHat) / (2*log(2)*exp(lsig))
    mObj$family$putLam( lam )
    mObj$family$putTheta( lsig )
    
    mFit <- gam(G = mObj, in.out = initM[["in.out"]], start = initM[["start"]])
    
    initM <- list("start" = coef(mFit), "in.out" = list("sp" = mFit$sp, "scale" = 1))
    
    tmpDat <- as.data.frame( cbind(bObj, initB, pMat) )
    names(tmpDat) <- c(".obj", ".init", ".pMat")
    
    # Loop over bootstrap datasets to get standardized deviations from full data fit
    withCallingHandlers({
      out <- mlply(.data = tmpDat,
                   .fun = .funToApplyTuneLearnFast,
                   .parallel = multicore,
                   .inform = control[["verbose"]],
                   .paropts = paropts,
                   ### ... arguments start here
                   .mFit = mFit, .lam = lam, .lsig = lsig, .fData = data
      ) 
    }, warning = function(w) {
      # There is a bug in plyr concerning a useless warning about "..."
      if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
        invokeRestart("muffleWarning")
    })
    
    loss <- .adTest( as.vector(sapply(out, "[[", "z")) )
    initB <- lapply(out, "[[", "init")
    pMat <- lapply(out, "[[", "pMat")
    
    return( list("loss" = loss, "initM" = initM, "initB" = initB, "pMat" = pMat) )
    
  }
  
  init <- list("initM" = initM, "initB" = vector("list", nbo), "pMat" = vector("list", nbo))
  
  # If we get convergence error, we increase "err" up to 0.2. I the error persists (or if the 
  # error is of another nature) we throw an error
  repeat{
    res <- tryCatch(.brent(brac=srange, f=objFun, mObj = mObj, bObj = bObj, init = init, t = control$tol), error = function(e) e)
    #res <- tryCatch(.brent()optimize(obj, srange, tol = ctrl$tol), error = function(e) e)
    
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


###########
# Function that fits a single bootstrapped dataset (.obj) using a certain initialization (.init)
# and prediction matrix (.pMat). It returns standardized deviations from full data fit. 
# To be run in parallel (over boostrapped datasets) by .tuneLearnFast()
###########
.funToApplyTuneLearnFast <- function(.obj, .init, .pMat, .mFit, .lam, .lsig, .fData)
{
  # For some reason I have to do this
  .obj <- .obj[[1]]; .init <- .init[[1]]; .pMat <- .pMat[[1]]
  
  .obj$lsp0 <- log( .mFit$sp )
  .obj$family$putLam( .lam )
  .obj$family$putTheta( .lsig )
  
  fit <- gam(G = .obj, start = .init)
  
  .init <- betas <- coef(fit)
  Vp <- fit$Vp
  
  # Create prediction design matrix (only in first iteration)
  if( is.null(.pMat) ) { 
    .pMat <- predict.gam(fit, newdata = .fData, type = "lpmatrix") 
    lpi <- attr(.pMat, "lpi")
    if( !is.null(lpi) ){ 
      .pMat <- .pMat[ , lpi[[1]]] #"lpi" attribute lost here, re-inserted in next line 
      attr(.pMat, "lpi") <- lpi 
    }
  }
  
  # In the gamlss case, we are interested only in the calibrating the location model
  lpi <- attr(.pMat, "lpi")
  if( !is.null(lpi) ){
    betas <- betas[lpi[[1]]]
    Vp <- fit$Vp[lpi[[1]], lpi[[1]]]
  }
  
  mu <- .pMat %*% betas
  sdev <- sqrt( diag( .pMat%*%Vp%*%t(.pMat) ) )
  
  .z <- (mu - as.matrix(.mFit$fit)[ , 1]) / sdev
  
  return( list("z" = .z, "init" = .init, "pMat" = .pMat) )
}

