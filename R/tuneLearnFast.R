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
                          ncores = 1,  control = list())
{ 
  # Setting up control parameter
  ctrl <- list( "init" = NULL,
                "brac" = log(c(0.5, 2)), 
                "K" = 50,
                "gausFit" = NULL,
                "tol" = 1e-2,
                "verbose" = TRUE )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  tol <- ctrl[["tol"]]
  brac <- ctrl[["brac"]]
  
  # Sanity check
  if( tol > 0.1 * abs(diff(brac)) ) 
    stop("tol > bracket_widths/10, choose smaller tolerance or larger bracket")
  
  # (Optional) create K boostrap dataset
  if( is.null(boot) ){
    n <- nrow(data)
    tmp <- lapply(1:ctrl[["K"]], function(nouse) sample(1:n, n, replace = TRUE))
    boot <- lapply(tmp, function(ff) data[ff, ] )
  }
  
  # (Optional) Main Gaussian fit, used for initializations
  if( is.null(ctrl[["gausFit"]]) ) {
    if( plyr:::is.formula(form) ) {
      ctrl[["gausFit"]] <- gam(form, data = data)
    } else {
      ctrl[["gausFit"]] <- gam(form, data = data, family = gaulss)
    } }
  
  # Order quantiles so that those close to the median are dealt with first
  oQu <- order( abs(qu-0.5) )
  
  # (Optional) Initializing the search range for sigma
  if( is.null(ctrl[["init"]]) ){
    # We assume lam~0 and we match the variance of a symmetric Laplace density with that of the Gaussian fit.
    # We use the value of qu that is the closest to 0.5
    tmp <- qu[ oQu[1] ]
    isig <- log(sqrt(  ctrl$gausFit$sig2 * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
  } else {
    isig <- ctrl[["init"]]
  }
  
  # Initial search range
  brac <- ctrl$brac
  srange <- isig + brac

  # Initialize
  nt <- length(qu)

  # Estimated learning rates, # of bracket expansions, error rates and bracket ranges used in bisection
  sigs <- efacts <- errors <- numeric(nt)
  rans <- matrix(NA, nt, 2)
  names(sigs) <- names(errors) <- rownames(rans) <- qu
  
  # Here we need bTol > aTol otherwise we the new bracket will be too close to probable solution
  aTol <- 1.5 * tol
  bTol <- 2 * tol
  
  for(ii in 1:nt)
  {
    oi <- oQu[ii]
    
    ef <- 1
    
    repeat{
      
      # Compute bracket
      srange <- isig + ef * brac
      
      # Estimate log(sigma)
      res  <- .tuneLearnFast(form = form, data = data, qu = qu[oi], err = err,
                             boot = boot, srange = srange, ncores = ncores, ctrl = ctrl)  
      
      lsig <- res$minimum
      
      # If solution not too close to boundary store results and determine bracket for next iter
      if( all(abs(lsig-srange) > aTol) ){ 
        
        sigs[oi] <- lsig
        rans[oi, ] <- srange
        efacts[oi] <- ef
        errors[oi] <- res$err
        
        if(ii < nt)
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
      wd <- abs( diff(brac) )
      isig <- lsig + ifelse(lsig-srange[1] < aTol, - wd + bTol, wd - bTol)
      ef <- 2*ef
    }
    
    if( ctrl$verbose && (nt>1) )
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
  
  out <- list("lsigma" = sigs, "err" = errors, "ranges" = rans)
  
  return( out )
}

##########################################################################
### Internal version, which works with scalar qu
########################################################################## 
.tuneLearnFast <- function(form, data, boot, qu, err, srange, ncores, ctrl)
{
  n <- nrow(data)
  
  # Objective function to be minimized
  obj <- function(lsig)
  {
    lam <- err * sqrt(2*pi*ctrl$gausFit$sig2) / (2*log(2)*exp(lsig)) 
    
    # Extended Gam or ...
    if( plyr:::is.formula(form) ){
      
      mainFit <- gam(form, family = logF(qu = qu, lam = lam, theta = lsig), data = data)
      
      z <- sapply(boot, 
                  function(input)
                  {
                    fit <- gam(form, family = logF(qu = qu, lam = lam, theta = lsig), data = input, 
                               sp = mainFit$sp, start = coef(mainFit))
                    
                    pred <- predict(fit, newdata = data, se = TRUE)
                    
                    .z <- (pred$fit - mainFit$fit) / pred$se.fit
                    
                    return( .z )
                  })
      
    } else { #... Gamlss
      
      mainFit <- gam(form, family = logFlss2(qu = qu, lam = lam, offset = lsig), data = data)
      
      z <- sapply(boot, 
                  function(input)
                  {
                    fit <- gam(form, family = logFlss2(qu = qu, lam = lam, offset = lsig), 
                               data = input, sp = mainFit$sp)
                    
                    pred <- predict(fit, newdata = data, se = TRUE)
                    
                    .z <-  (pred$fit[ , 1] - mainFit$fit[ , 1]) / pred$se.fit[ , 1] 
                    
                    return( .z )
                  })
      
    }
    
    loss <- qgam:::.adTest( as.vector(z) )
    
    return( loss )
    
  }
  
  # If we get convergence error, we increase "err" up to 0.2. I the error persists (or if the 
  # error is of another nature) we throw an error
  repeat{
    res <- tryCatch(optimize(obj, srange, tol = ctrl$tol), error = function(e) e)
    
    if("error" %in% class(res)){
      if( grepl("can't correct step size", res) ) {
        if(err < 0.2){
          err <- min(2*err, 0.2)
          if(ctrl$verbose) message( paste("Increase \"err\" to ", err, " to get convergence") )
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

