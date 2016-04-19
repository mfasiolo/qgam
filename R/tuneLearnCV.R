####### Tuning the learning rate for Gibbs posterior

tuneLearnCV <- function(form, data, nrep, sig, qu, lamba = 0.1, ncores = 1, K = 10)
{ 
  n <- nrow(data)
  
  index <- sample(1:K, n, replace = TRUE)
  cvSets <- lapply(1:K, 
                 function(ii) 
                 { 
                   test <- data[index == ii, ]
                   train <- data[index != ii, ]
                   return( list("test" = test, "train" = train) )
                 })

  # Estimating coverage
  simul <- mclapply(sig, function(inSig){
    
    if( is.formula(form) )
    {
      
      mainFit <- gam(form, family = logF(qu = qu, lam = lam, theta = inSig), data = data)
      
      coverage <- sapply(cvSets, 
                         function(input)
                         {
                           fit <- gam(form, 
                                      family = logF(qu = qu, lam = lam, theta = inSig), 
                                      data = input[["train"]], sp = mainFit$sp)
                           
                           pred <- predict(fit, newdata = input[["test"]], se = TRUE)
                           
                           return( checkloss(input[["test"]], pred$fit, qu) )
                           
                         })
      
    } else {
      
      mainFit <- gam(form, family = logFlss2(qu = qu, lam = lam, offset = inSig), data = data)
      
      coverage <- sapply(sets, 
                         function(input)
                         {
                           fit <- gam(form, family = logFlss2(qu = qu, lam = lam, offset = inSig), 
                                      data = input[["train"]], sp = mainFit$sp)
                           
                           pred <- predict(fit, newdata = data, se = TRUE)
                           
                           return( checkloss(input[["test"]], pred$fit[ , 1], qu) )
                           
                         })
      
    }
    
    return( coverage )
  }, mc.cores = 4)
  
  cover <- sapply(simul, mean)
  
  return( cover )
}