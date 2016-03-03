# Computes log-F density and its derivatives
dlf <- function(x, tau, mu, sig, lam, log = FALSE, deriv = 0)
{
  y <- (x - mu) / sig
  
  out <- tau * y - lam * log1pexp( y / lam ) - log( sig * lam * beta(lam*tau, (1-tau)*lam) ) 
  
  if( !log ) out <- exp(out)
  
  if( deriv > 0 )
  {
    out <- list("d" = out)
    
    dl <- dlogis(x, mu, lam*sig)
    pl <- plogis(x, mu, lam*sig)
    
    out$D <- c(  sum((pl - tau) / sig),  sum((y * (pl - tau) - 1) / sig) ) 
    
    if(deriv > 1)
    {
      H <- matrix(NA, 2, 2)
      
      H[1, 1] <- sum( - dl / sig )
      H[2, 2] <- sum( ( 2*y*(tau - pl - 0.5 * (x-mu)*dl) + 1 )/sig^2 )
      H[1, 2] <- H[2, 1] <- sum( - ((x-mu)*dl + pl - tau) / sig^2 )
      
      out$D2 <- H
      
      if(deriv > 2)
      {
        z <- y / lam
        der <- sigmoid(z, deriv = TRUE)
        
        out$D3m <- sum( der$D2 / (lam^2 * sig^3) )
        out$D4m <- sum( - der$D3 / (lam^3 * sig^4) )
        out$D2mDs <- sum( (z*der$D2 + 2*der$D1) / (lam * sig^3) )
        out$D3mDs <- sum( -(z*der$D3 + 3*der$D2) / (lam^2 * sig^4) )
        out$DmD2s <- sum( (2*(der$D0-tau) + 4*z*der$D1 + z^2*der$D2) / (sig^3) )
        out$D2mD2s <- sum( - (z^2*der$D3 + 6*z*der$D2 + 6*der$D1) / (lam*sig^4) )
      }
      
    }
    
  }
  
  return( out )
  
}


#############
# Checking if the density is normalized to 1
#############

# n <- 1e7
# x <- rt(n, df = 3)
# 
# dstud <- dt(x, df = 3)
# 
# mean( lf(x, tau = 0.5, mu = 0.5, sig = 1, lam = 0.1) / dstud )
# mean( lf(x, tau = 0.5, mu = 0.5, sig = 1, lam = 0.5) / dstud )
# mean( lf(x, tau = 0.2, mu = 0, sig = 1, lam = 0.1) / dstud )
# mean( lf(x, tau = 0.5, mu = 0, sig = 4, lam = 0.1) / dstud )