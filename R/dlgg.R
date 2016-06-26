# Computes log-F density and its derivatives
.dlgg <- function(y, mu, tau, phi, log = FALSE, deriv = 0)
{
  k <- exp( 2*phi )
  rk <- sqrt( k )
  sig <- exp( tau )
  
  z <- (y - mu) / sig
    
  out <- 2*(k - 0.5)*phi - tau - lgamma(k) + z*rk - k*exp(z/rk)
  
  if( !log ) out <- exp(out)
  
  if( deriv > 0 )
  {
    out <- list("d" = out)
    
    out$D <- cbind( rk/sig * (exp(z/rk) - 1), 
                    z*rk * ( exp(z/rk) - 1 ) - 1,
                    rk*z + exp(z/rk+phi) * (z-2*rk) + 4*k*phi + 2*k - 2*k*psigamma(k) - 1
                    ) 
        
    if(deriv > 1)
    {
      H <- matrix(NA, 3, 3)
      
      H[1, 1] <- sum( - exp(z/rk - 2*tau) )
      H[2, 2] <- sum( z*rk - z * exp(z/rk) * (z+rk) )
      H[3, 3] <- sum( z * (exp(z/rk)*(3*rk-z)+rk) - 4*k*(exp(z/rk) - 2*phi + psigamma(k) + k*psigamma(k, 1) - 2) )
      H[1, 2] <- H[2, 1] <- sum( -exp(z/rk - tau)*(z + rk) + rk/sig )
      H[1, 3] <- H[3, 1] <- sum( -exp(z/rk - tau)*(z - rk) - rk/sig )
      H[2, 3] <- H[3, 2] <- sum( -z * rk * ( exp(z/rk - phi) * (z - rk) + 1 ) )
      
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