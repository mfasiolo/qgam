####### Tuning the learning rate for Gibbs posterior

tuneLearn <- function(form, data, nrep, lsig, qu, err = 0.01, ncores = 1)
{ 
  if( length(qu) > 1 ) stop("length(qu) > 1, but this method works only for scalar qu")
  
  n <- nrow(data)
  
  # Create bootstrapped datasets
  index <- lapply(1:nrep, function(nouse) sample(1:n, n, replace = TRUE))
  bootSets <- lapply(index, function(ff) data[ff, ] )
  
  # Main Gaussian fit, used for initializations
  if( is.formula(form) )
  {
   gauFit <- gam(form, data = data)
  } else {
   gauFit <- gam(form, data = data, family = gaulss)
  }
  
  lam <- err * sqrt(2*pi*gauFit$sig2) / (2*log(2)*exp(lsig)) 
  
  # Estimating coverage
  z <- lapply(1:length(lsig), function(ii){
        
    if( is.formula(form) )
    {
      
      mainFit <- gam(form, family = logF(qu = qu, lam = lam[ii], theta = lsig[ii]), data = data)
      
      out <- sapply(bootSets, 
                         function(input)
                         {
                           fit <- gam(form, family = logF(qu = qu, lam = lam[ii], theta = lsig[ii]), data = input, 
                                      sp = mainFit$sp, start = coef(mainFit))
                           
                           pred <- predict(fit, newdata = data, se = TRUE)
                           
                           .z <- (pred$fit - mainFit$fit) / pred$se.fit
                           
                           return( .z )
                         })
      
    } else {
      
    #stop("This does not work, yet")
    
    mainFit <- gam(form, family = logFlss2(qu = qu, lam = lam[ii], offset = lsig[ii]), data = data)
    
    out <- sapply(bootSets, 
                       function(input)
                       {
                         fit <- gam(form, family = logFlss2(qu = qu, lam = lam[ii], offset = lsig[ii]), 
                                    data = input, sp = mainFit$sp)
                         
                         pred <- predict(fit, newdata = data, se = TRUE)
                         
                         .z <-  (pred$fit[ , 1] - mainFit$fit[ , 1]) / pred$se.fit[ , 1] 
                         
                         return( .z )
                       })
    
    }
    
    return( out )
  })#, mc.cores = ncores)
  
  loss <- sapply(z, function(.x) .adTest(as.vector(.x)))
  
  return( list("loss" = loss, "lambda" = lam) )
}