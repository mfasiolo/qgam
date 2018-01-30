##########
# Internal function which calculated calibration loss function based on the sandwich estimator
# wrt the regression coefficients. 
# INPUT
# - mFit: a gamObject fitted using the elflss family
# - X: the submatrix of XFull corresponding to the first linear predictor (X = XFull[ , lpi[[2]]]).
#      Notice that X = XFull in the extended GAM case.
# - XFull: the full "lpmatrix" corresponding to both linear predictors
# - sdev: the posterior standard deviation of the first linear predicto 
# OUTPUT
# - a scalar indicating the loss
#
.sandwichLoss <- function(mFit, X, XFull, sdev, repar){
  
  lpi <- attr(X, "lpi")

  if( !is.null(lpi) ){ # GAMLSS version OR ...
    
    if( is.null(mFit$rp) ) { stop("mFit$rp is NULL, but a re-parametrization list is needed")  }
    
    mFit$Sl <- repar # Add reparametrization list

    # Extract observed Fisher information and invert the transformations
    OFI <- - mFit$lbb
    OFI <- Sl.repara(mFit$rp, OFI, inverse = TRUE)
    OFI <- Sl.initial.repara(mFit$Sl, OFI, inverse = TRUE, cov = FALSE)
    
    # Extract penalty matrix and invert transformations
    P <- mFit$St
    P <- Sl.repara(mFit$rp, P, inverse = TRUE)
    P <- Sl.initial.repara(mFit$Sl, P, inverse = TRUE, cov = FALSE)
    
    # Posterior variance of fitted quantile (mu) 
    varOFI <- sdev ^ 2
    
    # Calculate variance of the score
    grad <- .llkGrads(gObj = mFit, X = XFull, jj = lpi)
    
    V <- cov( grad ) * (nrow(X) - 1)         # Variance of the score
    S <- solve( OFI %*% solve(V, OFI) + P )[lpi[[1]], lpi[[1]]]  # 'Sandwich' posterior covariance

    # Posterior variance of mu using sandwich
    varSand <- rowSums((X %*% S) * X)
    
    bias <- numeric( nrow(X) ) * 0
    
    # Average distance KL(sand_i, post_i)
    outLoss <- mean( (varSand/varOFI + bias^2/varOFI + log(varOFI/varSand)) )
    
  } else { # Extended GAM version 
    
    # Working weights, observed Fisher information and penalty matrix
    w <- mFit$working.weights
    OFI <- t(X) %*% (w * X)
    P <- .getPenMatrix(q = ncol(X), UrS = repar$UrS, sp = log(mFit$sp), Mp = repar$Mp, U1 = repar$U1)
    
    # Posterior variance of fitted quantile (mu) 
    varOFI <- sdev ^ 2
    
    # Calculate variance of the score
    r <- mFit$residuals
    grad <- X * drop(w * r)
    V <- cov( grad ) * (length(w) - 1)
    S <- solve( OFI %*% solve(V, OFI) + P )  # 'Sandwich' posterior covariance
    
    # Posterior variance of mu using sandwich
    varSand <- rowSums((X %*% S) * X)
    
    bias <- numeric( length(w) ) * 0
    
    # Average distance KL(sand_i, post_i)
    outLoss <- mean( (varSand/varOFI + bias^2/varOFI + log(varOFI/varSand)) )
    
    initB <- NULL
    
  }
  
  return( outLoss )
  
}