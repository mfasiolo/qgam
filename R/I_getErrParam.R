##############
# Get "err" parameter for automatic loss smoothness selection 
##############
#
.getErrParam <- function(qu, gFit, varHat){
  
  # Estimated conditional mean and variance (latter could be constant)
  muHat <- as.matrix(gFit$fitted.values)[ , 1]

  # Raw residuals from Gaussian regression, normalized using estimated conditional SD
  r <- ( gFit$y - muHat ) / sqrt( varHat )

  n <- length( r )
  
  # Fixing dimension d to EDF of Gaussian fit. 
  # First use anova to find degrees of freedom of parametric terms in equation for location
  # Then find EDF of smooth terms in equation for location. unique() needed for "adaptive" smooths
  anv <- anova( gFit )
  d <- sum(anv$pTerms.df[!grepl("\\.1", rownames(anv$pTerms.table))]) + ("(Intercept)" %in% names(anv$p.coeff))
  d <- d + sum( unique(pen.edf(gFit)[!grepl("s\\.1|te\\.1|ti\\.1|t2\\.1", names(pen.edf(gFit)))]) )
  
  # Estimate parameters of shash density on standardized residuals
  parSH <- .fitShash( r )$par
  
  # Find probability p corresponding to the mode of shash density
  pmode  <- .shashCDF(.shashMode( parSH ), parSH)
  
  err <- qu * 0
  for(ii in 1:length(qu)){
    
    quX <- qu[ii]
    
    # If quantile qu is too close to mode, lower it or increase it by 0.05 to avoid f' = 0
    if( abs(quX - pmode) < 0.05 ){
      quX <- pmode + sign(quX - pmode) * 0.05
      quX <- max(min(quX, 0.99), 0.01) # To avoid going outside (0, 1)
    }
    
    # Quantile of shash at which derivatives should be estimated
    qhat <- .shashQf(quX, parSH)
    
    # Compure log(density) and log( abs(derivative of density) ) at quantile rqu
    # |Df / Dx| = |(Dlog(f) / Dx) * f|
    # log( |Df / Dx| ) = log( |Dlog(f) / Dx| ) + log( f )
    lf0 <- .llkShash(qhat, mu = parSH[1], tau = parSH[2], eps = parSH[3], phi = parSH[4])$l0
    lf1 <- - .llkShash(qhat, mu = parSH[1], tau = parSH[2], eps = parSH[3], phi = parSH[4], deriv = 1)$l1[1] # NB df/dx = -df/dmu
    lf1 <- log( abs(lf1) ) + lf0
    
    # f / f'^2 = exp( log(f) - 2 * log(|f'|) ) but we avoid dividing by almost zero
    h <- (d*9 / (n*pi^4))^(1/3) * exp(lf0/3 - 2*lf1/3)
    
    # Setting err too high might be problematic, so here it can be at most 1
    err[ii] <- min(h * 2 * log(2) / sqrt(2*pi), 1)
    
  }

  return( err )
  
}
