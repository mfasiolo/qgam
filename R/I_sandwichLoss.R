##########
# Internal function which calculated calibration loss function based on the sandwich estimator
# wrt the regression coefficients. 
# INPUT
# - mFit: a gamObject fitted using the elflss family
# - X: model matrix
# - sdev: the posterior standard deviation of the first linear predicto 
# OUTPUT
# - a scalar indicating the loss
#
.sandwichLoss <- function(mFit, X, sdev, repar, mObj, alpha = NULL, VSim = NULL){
  
  discrete <- !is.null(mObj$Xd)
  
  n <- length(mFit$y)
  
  # Posterior variance of fitted quantile (mu) 
  varOFI <- sdev ^ 2
  
  # Observed Fisher information and penalty matrix
  if(discrete){
    woW <- mFit$wt
    OFI <- XWXd(X=mObj$Xd,w=woW,k=mObj$kd,ks=mObj$ks,ts=mObj$ts,
                dt=mObj$dt,v=mObj$v,qc=mObj$qc,drop=mObj$drop)
  } else {
    woW <- mFit$working.weights # NB: these can be negative
    OFI <- crossprod(sign(woW)*sqrt(abs(woW))*X, sqrt(abs(woW))*X)  
  }
  
  P <- .getPenMatrix(q = length(mFit$coefficients), UrS = repar$UrS,
                     sp = log(mFit$sp), Mp = repar$Mp, U1 = repar$U1)
  
  # Calculate variance of the score
  if(discrete){
    sqrtw <- .llkGrads(gObj = mFit, X = NULL, mObj = mObj, type = "DllkDeta")$l1
    mXsqrtw <- XWyd(X=mObj$Xd,w=sqrtw,y=rep(1,n),k=mObj$kd,ks=mObj$ks,
                    ts=mObj$ts,dt=mObj$dt,v=mObj$v,qc=mObj$qc,drop=mObj$drop)/n
    V <- XWXd(X=mObj$Xd,w=sqrtw^2,k=mObj$kd,ks=mObj$ks,ts=mObj$ts,dt=mObj$dt,
              v=mObj$v,qc=mObj$qc,drop=mObj$drop)/(n-1) - mXsqrtw %*% t(mXsqrtw)
  } else {
    grad <- .llkGrads(gObj = mFit, X = X, mObj = mObj)
    V <- cov( grad )
  }
 
  if( !is.null(alpha) ){
    V <- alpha * V + (1-alpha) * VSim
  }
  V <- V * n 

  # Compute eigen-decomposition of V, get its rank and produce pseudo-inverse
  # Discard eigen-vectors corresponding to very LOW variance (low eigen-values) of the score,
  # the logic being that if the variance of the score of a parameter is LOW then
  # its variance will be HIGH --> we discard non-identifiable parameters.
  iVscore <- .invert_psdef_matrix(V, r = mFit$rank)
  
  # Inverse 'Sandwich' posterior covariance
  iS <- OFI %*% iVscore %*% OFI + P 
  
  # Compute eigen-decomposition of V_sand^{-1}, get its rank and produce pseudo-inverse
  # Discard LOW eigen-values of precision sandwich matrix --> they correspond to directions
  # of HIGH variance --> non-identifiable parameters
  VS <- .invert_psdef_matrix(iS, r = mFit$rank)

  # Posterior variance using sandwich: var(mu) = diag( X %*% iS^-1 %*% t(X) )
  if(discrete){
    varSand <- diagXVXd(X=mObj$Xd,V=VS,k=mObj$kd,ks=mObj$ks,
                        ts=mObj$ts,dt=mObj$dt,v=mObj$v,qc=mObj$qc,drop=mObj$drop)
  } else {
    varSand <- rowSums((X %*% VS) * X)
  }
  
  bias <- rep(0, n)
  
  # Excluded observ where the variance is 0 (probably because all corresponding covariate values are 0)
  not0 <- which( !(varOFI == 0 & varSand == 0) )
  
  # Average distance KL(sand_i, post_i) on non-zeros
  outLoss <- mean( sqrt(varSand/varOFI + bias^2/varOFI + log(varOFI/varSand))[not0] )
  
  return( outLoss )
  
}