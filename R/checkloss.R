

checkloss <- function(y, mu, qu){
  
  tau <- 1 - qu
  
  d <- y - mu

  l <- - sum( tau*d[d<0] ) - sum( (tau-1)*d[d>0] )
  
  return( l )
  
} 