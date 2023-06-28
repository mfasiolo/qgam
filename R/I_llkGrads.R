##########
# Internal function which calculates the gradient of each element of the (unpenalized) log-likelihood
# wrt the regression coefficients. 
# INPUT
# - gObj: a gamObject fitted using the elf or elflss family
# - X: the full "lpmatrix" corresponding to both linear predictors
# - type: set it to "DllkDb" if you want the derivative of the log-lik wrt the regression coeff, or 
#         to "DllkDeta" if you want those wrt the linear predictor.
# 
# OUTPUT
# - an n by p matrix where the i-th columns is the gradient of the i-th log-likelihood component
#
.llkGrads <- function(gObj, X, type = "DllkDb") {
  
  type <- match.arg(type, c("DllkDb", "DllkDeta"))
  
  y <- gObj$y
  beta <- coef( gObj )
  fam <- gObj$family
  wt <- gObj$prior.weights
  offset <- gObj$offset
  
  n <- length( y )
  p <- length( beta )
  
  tau <- fam$getQu( )
  theta <- fam$getTheta( )
  co <- fam$getCo( )
  
  sig <- exp( theta )
  lam <- co / sig
  
  if( is.null(offset) ){
    offset <- numeric( nrow(X) )  
  }
  
  eta <- X %*% beta + offset
  mu <- fam$linkinv( eta )
  
  z <- (y - mu) / sig
  
  pl <- plogis(y, mu, lam * sig)
  
  # [1] Derivatives of llk wrt parameters
  l1 <- numeric( n )
  l1 <- wt * (pl - 1 + tau) / sig
  
  # Derivative of link function 
  ig1 <- fam$mu.eta(eta)
  
  # [2] Transform llk derivatives wrt mu to derivatives wrt linear predictor (eta)
  l1 <- l1 * ig1
  if( type == "DllkDeta" ){ return(list("l1" = as.matrix(l1), "sig" = sig)) }
  
  # [3] Transform into derivatives wrt regression coefficients
  # The i-th column of 'grads' is the score of the i-th likelihood component 
  grads <- drop(l1) * X
  
  return(grads)
  
}
