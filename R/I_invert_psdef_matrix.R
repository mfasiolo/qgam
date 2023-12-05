# Function to invert positive definite matrix of rank <= r
.invert_psdef_matrix <- function(S, r){
  
  d <- diag(S)^(-0.5)
  D <- diag(d, ncol = length(d))
  S <- D %*% S %*% D
  eS <- eigen(S, symmetric = TRUE)
  r <- min(sum(eS$values > max(eS$values) * 1e-6), r)
  QiT <- (t(eS$vectors[ , 1:r, drop = FALSE]) / sqrt(eS$values[1:r])) %*% D
  return( crossprod(QiT) )
  
}