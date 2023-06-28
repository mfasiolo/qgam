##########
# Internal function which calculated calibration loss function based on the sandwich estimator
# wrt the regression coefficients. 
# INPUT
# - mFit: a gamObject fitted using the elflss family
# - X: model matrix
# - XFull: the full "lpmatrix" corresponding to both linear predictors
# - sdev: the posterior standard deviation of the first linear predicto 
# OUTPUT
# - a scalar indicating the loss
#
.sandwichLoss <- function(mFit, X, sdev, repar, alpha = NULL, VSim = NULL){
  
  # Posterior variance of fitted quantile (mu) 
  varOFI <- sdev ^ 2
  
  # Observed Fisher information and penalty matrix
  woW <- mFit$working.weights # NB: these can be negative
  OFI <- crossprod(sign(woW)*sqrt(abs(woW))*X, sqrt(abs(woW))*X)  
  P <- .getPenMatrix(q = ncol(X), UrS = repar$UrS, sp = log(mFit$sp), Mp = repar$Mp, U1 = repar$U1)
  
  # Calculate variance of the score
  grad <- .llkGrads(gObj = mFit, X = X)
  
  # Covariance matrix of the score
  V <- cov( grad ) 
  if( !is.null(alpha) ){
    V <- alpha * V + (1-alpha) * VSim
  }
  V <- V * nrow(X) 
  
  # Compute eigen-decomposition of V, get its rank and produce pseudo-inverse
  eV <- eigen( V )
  rv <- sum( eV$values > eV$values[1] * .Machine$double.eps )
  Q <- t( t(eV$vectors[ , 1:rv]) / sqrt(eV$values[1:rv]) ) 
  
  # Inverse 'Sandwich' posterior covariance
  iS <- (OFI %*% Q) %*% crossprod(Q, OFI) + P  
  
  # Computed the Cholesky factor
  C <- chol( iS )

  # Posterior variance using sandwich: var(mu) = diag( X %*% iS^-1 %*% t(X) )
  varSand <- rowSums(t(backsolve(C, forwardsolve(t(C), t(X)))) * X)
  
  bias <- numeric( nrow(X) ) * 0
  
  # Excluded observ where the variance is 0 (probably because all corresponding covariate values are 0)
  not0 <- which( !(varOFI == 0 & varSand == 0) )
  
  # Average distance KL(sand_i, post_i) on non-zeros
  outLoss <- mean( sqrt(varSand/varOFI + bias^2/varOFI + log(varOFI/varSand))[not0] )
  
  return( outLoss )
  
}