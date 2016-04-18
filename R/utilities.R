###########
# Utilities
##########

#### Vettorize check Loss function
checklossVett <- function(y, mu, p){
  
  n <- length( p )
  
  out <- sapply(1:n,
                function(ii){
                  return( checkloss(y, mu[ , ii], p[ii]) )
                })
  
  return( out )
}


#### Vettorized empirical cdf
qqVett <- function(y, mu){
  
  nq <- ncol( mu )
  nobs <- length( y )
  
  out <- sapply(1:nq,
                function(ii){
                  return( sum( (y - mu[ , ii]) < 0 ) / nobs )
                })
  
  return( out )
}


#### Andreson-Darling test for STANDARD normality
.adTest <- function(.x){
  
  n <- length(.x)
  .x <- sort(.x)
  
  # Cramer-von Mises statistic 
  # out <- 1/(12*n) + sum( ((2*1:n - 1)/(2*n) - pnorm(.x))^2 )
  
  logp1 <- pnorm(.x, log.p = TRUE)
  logp2 <- pnorm(.x, lower.tail = F, log.p = TRUE)
  
  h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
  out <- -n - mean(h)
  
  return( out )
}