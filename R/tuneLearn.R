####### Tuning the learning rate for Gibbs posterior

tuneLearn <- function(form, data, lsig, qu, err = 0.01, 
                      ncores = 1, control = list(), controlGam = list())
{ 
  if( length(qu) > 1 ) stop("length(qu) > 1, but this method works only for scalar qu")
  
  lsig <- sort( lsig )
  
  # Setting up control parameter
  ctrl <- list( "K" = 50, "b" = 0 )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  n <- nrow(data)
  nt <- length(lsig)
  
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
    
    fit <- gam(G = mainObj, in.out = initM[["in.out"]], start = initM[["start"]])
    
    initM <- list("start" = coef(fit), "in.out" = list("sp" = fit$sp, "scale" = 1))
    
    mainFit[[ii]] <- list("sp" = fit$sp, "fit" = fit$fitted, "lam" = lam)
  }
  
  # Fitting bootstrapped datasets
  z <- lapply( 1:ctrl[["K"]], # Loop over bootstrap datasets 
                 function(kk)
                 { 
                   init <- NULL
                   
                   # Create gam object
                   bObj <- gam(form, family = get(fam)(qu = qu, lam = NA, theta = NA), data = boot[[kk]], 
                               sp = mainFit[[1]]$sp, control = controlGam, fit = FALSE)
                   
                   .z <- matrix(NA, nt, n)
                   
                   for( ii in nt:1 )  # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harded)
                   {   
                     bObj$lsp0 <- log( mainFit[[ii]]$sp )
                     bObj$family$putLam( mainFit[[ii]]$lam )
                     bObj$family$putTheta( lsig[ii] )
                     
                     fit <- gam(G = bObj, start = init)
                     init <- coef(fit)
                     
                     pred <- predict(fit, newdata = data, se = TRUE)
                     
                     .z[ii, ] <- (as.matrix(pred$fit)[ , 1] - as.matrix(mainFit[[ii]]$fit)[ , 1]) / as.matrix(pred$se.fit)[ , 1]
                   }
                   
                   return( .z )
                 })

  z <- do.call("cbind", z)

  loss <- apply(z, 1, function(.x) .adTest(as.vector(.x)))
  names(loss) <- lsig
  
  return( list("loss" = loss) )
}














#   for( ii in nt:1 ) # START lsigma loop, from largest to smallest (because when lsig is small the estimation is harded)
#   {
#     if( is.formula(form) ) # EXTENDED GAM OR ....
#     {         
#       out <- lapply(1:ctrl[["K"]], 
#                     function(kk)
#                     {
#                       fit <- gam(form, family = logF(qu = qu, lam = mainFit[[ii]]$lam, theta = lsig[ii]), data = boot[[kk]], 
#                                  sp = mainFit[[ii]]$sp, 
#                                  mustart = if(is.null(initB[[kk]])) mainFit[[ii]]$fit[ind[[kk]]] else initB[[kk]], 
#                                  control = controlGam)
#                       
#                       pred <- predict(fit, newdata = data, se = TRUE)
#                       
#                       .z <- (pred$fit - mainFit[[ii]]$fit) / pred$se.fit
#                       
#                       return( list("z" = .z, "init" = fit$fitted) )
#                     })
#       
#     } else { # ... GAMLSS
#       
#       out <- lapply(1:ctrl[["K"]], 
#                     function(kk)
#                     {
#                       fit <- gam(form, family = logFlss(qu = qu, lam = mainFit[[ii]]$lam, offset = lsig[ii]), 
#                                  data = boot[[kk]], sp = mainFit[[ii]]$sp, control = controlGam, 
#                                  start = initB[[kk]])
#                       
#                       pred <- predict(fit, newdata = data, se = TRUE)
#                       
#                       .z <-  (pred$fit[ , 1] - mainFit$fit[ , 1]) / pred$se.fit[ , 1] 
#                       
#                       return( list("z" = .z, "init" = coef(fit)) )
#                     })
#       
#     }
#     
#     # Save stardardized residuals and initialization
#     initB <- lapply(out, "[[", "init")
#     z[[ii]] <- sapply(out, "[[", "z")
#     
#   }  # STOP lsigma loop, from largest to smallest