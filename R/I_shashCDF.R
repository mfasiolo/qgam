####
# CDF of shash density
.shashCDF <- function(q, param) {

  mu <- param[1]
  sig <- exp(param[2])
  eps <- param[3]
  del <- exp(param[4])
  
  return( pnorm(sinh((asinh((q - mu)/(del * sig)) - eps/del) * del)) )
  
}