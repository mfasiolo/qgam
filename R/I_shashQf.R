###
# Quantile function of shash distribution
.shashQf <- function(p, param) {
  
  mu <- param[1]
  sig <- exp(param[2])
  eps <- param[3]
  del <- exp(param[4])
  
  return( mu + (del * sig) * sinh((1/del) * asinh(qnorm(p)) + (eps/del)) )
  
}