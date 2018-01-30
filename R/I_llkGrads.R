##########
# Internal function which calculated the gradient of each element of the (unpenalized) log-likelihood
# wrt the regression coefficients. 
# INPUT
# - gObj: a gamObject fitted using the elflss family
# - X: the full "lpmatrix" corresponding to both linear predictors
# - jj: a list of indexes indicating the position of the coefficients corresponding to each linear predictor
# 
# OUTPUT
# - an n by p matrix where the i-th columns is the gradient of the i-th log-likelihood component
#
.llkGrads <- function(gObj, X, jj) {
  
  y <- gObj$y
  beta <- coef( gObj )
  fam <- gObj$family
  wt <- gObj$prior.weights
  offset <- gObj$offset
  
  n <- length( y )
  p <- length( beta )
  
  tau <- fam$getQu( )
  theta <- fam$getTheta( )
  lam <- fam$getLam( )
  
  # If offset is not null or a vector of zeros, give an error
  if( !is.null(offset[[1]]) && sum(abs(offset)) )  stop("offset not still available for this family")
  
  eta <- X[ , jj[[1]], drop=FALSE] %*% beta[jj[[1]]]
  mu <- fam$linfo[[1]]$linkinv( eta )
  eta1 <- X[ , jj[[2]], drop=FALSE] %*% beta[jj[[2]]] + theta
  sig <-  fam$linfo[[2]]$linkinv( eta1 ) 
  
  z <- (y - mu) / sig
  
  dl <- dlogis(y, mu, lam * sig)
  pl <- plogis(y, mu, lam * sig)
  
  # [1] Derivatives of llk wrt parameters
  l1 <- matrix(0, n, 2)
  l1[ , 1] <- wt * (pl - 1 + tau) / sig
  l1[ , 2] <- wt * (z * (pl - 1 + tau) - 1) / sig 
  
  # Derivative of link function 
  ig1 <- cbind(fam$linfo[[1]]$mu.eta(eta), fam$linfo[[2]]$mu.eta(eta1))
  
  # [2] Transform llk derivatives wrt mu to derivatives wrt linear predictor (eta)
  l1 <- l1 * ig1
  
  # [3] Transform into derivatives wrt regression coefficients
  # The i-th column of 'grads' is the score of the i-th likelihood component 
  grads <- matrix(0, n, p)
  for (i in 1:length(jj)) {
    grads[ , jj[[i]]] <- grads[ , jj[[i]]] + l1[ , i] * X[ , jj[[i]],  drop = FALSE]
  }
  
  return(grads)
  
}