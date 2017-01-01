####### Check loss function
.checkloss <- function(y, mu, qu, add = TRUE){
  
  tau <- 1 - qu
  
  d <- y - mu
  
  l <- d * 0
  
  l[d < 0] <- - tau*d[d<0]
  l[d > 0] <- - (tau-1)*d[d>0]
  
  if( add ) l <- sum(l)
  
  return( l )
  
} 

#### Vettorize check Loss function
.checklossVett <- function(y, mu, p){
  
  n <- length( p )
  
  out <- sapply(1:n,
                function(ii){
                  return( .checkloss(y, mu[ , ii], p[ii]) )
                })
  
  return( out )
}