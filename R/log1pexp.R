## Computes log(1+exp(x)) accurately
log1pexp <- function(x)
{
  big <- which(x > 18)
  
  out <- log1p( exp(x) )
  
  if( length(big) )
  {
    xb <- x[ big ]
    out[big] <- xb + exp(-xb)
  }
  
  return( out )
}