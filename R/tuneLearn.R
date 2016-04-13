####### Tuning the learning rate for Gibbs posterior

tuneLearn <- function(form, data, nrep, sig, tau, err = 0.01, ncores = 1)
{ 
  n <- nrow(data)
  
  # Create bootstrapped datasets
  index <- lapply(1:nrep, function(nouse) sample(1:n, n, replace = TRUE))
  bootSets <- lapply(index, function(ff) data[ff, ] )
  
  gauFit <- gam(form, data = data)
  
  # Estimating coverage
  simul <- lapply(sig, function(inSig){
    
    if( is.formula(form) )
    {
      
      lam <- err * sqrt(2*pi*gauFit$sig2) / (2*log(2)*exp(inSig))  
      
      mainFit <- gam(form, family = logF(tau = tau, lam = lam, theta = inSig), data = data)
      
      coverage <- sapply(bootSets, 
                         function(input)
                         {
                           fit <- gam(form, family = logF(tau = tau, lam = lam, theta = inSig), data = input, 
                                      sp = mainFit$sp, start = coef(mainFit))
                           
                           pred <- predict(fit, newdata = data, se = TRUE)
                           
                           cover <- mean( (pred$fit+2*pred$se.fit > mainFit$fit) & (mainFit$fit > pred$fit-2*pred$se.fit) )
                           
                           return( cover )
                         })
      
    } else {
      
    stop("This does not work, yet")
    
    mainFit <- gam(form, family = logFlss2(tau = tau, lam = lam, offset = inSig), data = data)
    
    coverage <- sapply(bootSets, 
                       function(input)
                       {
                         fit <- gam(form, family = logFlss2(tau = tau, lam = lam, offset = inSig), 
                                    data = input, sp = mainFit$sp)
                         
                         pred <- predict(fit, newdata = data, se = TRUE)
                         
                         cover <- mean( (pred$fit[ , 1]+2*pred$se.fit[ , 1] > mainFit$fit[ , 1]) & 
                                        (mainFit$fit[ , 1] > pred$fit[ , 1]-2*pred$se.fit[ , 1]) )
                         
                         return( cover )
                       })
    
    }
    
    return( coverage )
  })#, mc.cores = ncores)
  
  cover <- sapply(simul, mean)
  
  return( cover )
}