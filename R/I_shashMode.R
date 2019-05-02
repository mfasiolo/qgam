###
# Finding the mode of the shash density
#
.shashMode <- function(param)
{
  .objFun <- function(mu){
    - .llkShash(x = mu, 
                mu = param[1], tau = param[2], 
                eps = param[3], phi = param[4], deriv = 0)$l0
  }
  
  range <- c(.shashQf(0.001, param), .shashQf(0.999, param))
  
  mode <- optimize(f = .objFun, interval = range)$minimum
  
  return( mode )
}