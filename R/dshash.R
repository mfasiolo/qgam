sech <- function(x) 1/cosh(x)

# Computes log-F density and its derivatives
dshash <- function(y, mu, tau, eps, phi, log = FALSE, deriv = 0)
{
  sig <- exp( tau )
  del <- exp( phi )
  
  z <- (y - mu) / (sig*del)
  
  dTasMe <- del*asinh(z) - eps
  CC <- cosh( dTasMe )
  SS <- sinh( dTasMe )
  
  out <- -tau - 0.5*log(2*pi) + log(CC) - 0.5*log1p(z^2) - 0.5*SS^2
  
  if( !log ) out <- exp(out)
  
  if( deriv > 0 )
  {
    out <- list("d" = out)
    
    DpDeps <- tanh(-dTasMe) - 0.5*sinh(-2*dTasMe)
    DpDmu <- 1/(del*sig*sqrt(z^2+1))*(del*(DpDeps)+z/sqrt(z^2+1))
    DpDtau <- z*del*sig*DpDmu - 1
    DpDphi <- DpDtau + 1 - del*asinh(z)*DpDeps #-del*(asinh(z)-z/sqrt(z^2+1))*DpDeps + z^2/(z^2+1)
    
    out$D <- cbind( DpDmu, DpDtau, DpDeps, DpDphi )
                    
    if(deriv > 1)
    {
      
      g <- -dTasMe
      
      H <- matrix(NA, 4, 4)
      
      D2DmuDeps <- (sech(g)^2 - cosh(2*g)) / (sig*sqrt(z^2+1))
      D2DtauDeps <- z*del*sig*D2DmuDeps
      D2DmuDmu <- D2DmuDeps/(sig*sqrt(z^2+1)) + z*DpDeps/(sig^2*del*(z^2+1)^(3/2)) + (z^2-1)/(del*sig*del*sig*(z^2+1)^2)
      D2DmuDtau <- z*sig*del*D2DmuDmu - DpDmu
      D2Deps <- -2*cosh(g)^2 + sech(g)^2 + 1 
      D2Dtau <-  z*del*sig*D2DmuDtau
      D2DepsDphi <- D2DtauDeps - del*asinh(z)*D2Deps
      D2DphiDmu <- D2DmuDtau + DpDeps/(sig*sqrt(z^2+1)) - del*asinh(z)*D2DmuDeps
      D2DtauDphi <- z*del*sig*D2DphiDmu
      D2Dphi <- D2DtauDphi - del*asinh(z)*D2DepsDphi + del*(z/sqrt(z^2+1)-asinh(z))*DpDeps
      
      
      H[1, 1] <- sum( D2DmuDmu )
      H[2, 2] <- sum( D2Dtau  )
      H[3, 3] <- sum( D2Deps )
      H[4, 4] <- sum( D2Dphi )
      H[1, 2] <- H[2, 1] <- sum( D2DmuDtau )
      H[1, 3] <- H[3, 1] <- sum( D2DmuDeps )
      H[1, 4] <- H[4, 1] <- sum( D2DphiDmu )
      H[2, 3] <- H[3, 2] <- sum( D2DtauDeps ) 
      H[2, 4] <- H[4, 2] <- sum( D2DtauDphi ) 
      H[3, 4] <- H[4, 3] <- sum( D2DepsDphi ) 
  
      out$D2 <- H
      
      if(deriv > 2)
      {
        #D2 <- matrix
        
        Deee <-  -2*(sinh(2*g)+sech(g)^2*tanh(g))
        Deem <- Deee/(sig*sqrt(z^2+1))
        Dmme <- Deem/(sig*sqrt(z^2+1)) + z*D2Deps/(sig*sig*del*(z^2+1)^(3/2))
        Dmmm <- 2*z*D2DmuDeps/(sig*sig*del*sqrt(z^2+1)^3) + Dmme/(sig*sqrt(z^2+1)) + (2*z^2-1)*DpDeps/(sig^3*del^2*sqrt(z^2+1)^5) +
                2*z*(z^2-3)/((sig*del)^3*(z^2+1)^3)
        Dmmt <- z*sig*del*Dmmm - 2*D2DmuDmu
        Deet <- z*sig*del*Deem
        Dmte <- z*sig*del*Dmme - D2DmuDeps
        Dtte <- z*sig*del*Dmte
        Dttm <- z*sig*del*Dmmt - D2DmuDtau
        Dttt <- z*sig*del*Dttm
        Dmep <- Dmte + D2Deps/(sig*sqrt(z^2+1)) - del*asinh(z)*Deem
        Dpet <- z*sig*del*Dmep
        Deep <- Deet - del*asinh(z)*Deee
        Dppe <- Dpet - del*asinh(z)*Deep + del*( z/sqrt(z^2+1)-asinh(z) )*D2Deps
        Dmmp <- Dmmt + 2*D2DmuDeps/(sig*sqrt(z^2+1)) + 
          z*DpDeps/(del*sig*sig*sqrt(z^2+1)^3) - del*asinh(z)*Dmme
        Dmtp <- z*sig*del*Dmmp - D2DphiDmu
        Dttp <- z*sig*del*Dmtp
        Dppm <- Dmtp + D2DepsDphi/(sig*sqrt(z^2+1)) + 
                z^2*DpDeps/(sig*sqrt(z^2+1)^3) - 
                del*asinh(z)*Dmep + del*D2DmuDeps*(z/sqrt(z^2+1) - asinh(z))
        Dppt <- z*sig*del*Dppm
        Dppp <- Dppt - del*asinh(z)*Dppe + del*(z/sqrt(z^2+1)-asinh(z))*(2*D2DepsDphi + DpDeps) +
                del*(z/sqrt(z^2+1))^3 * DpDeps

      
        # mmm
        out$D3 <- c( sum( Dmmm ),
                     sum( Dmmt ),
                     sum( Dmme ),
                     sum( Dmmp ),

                     sum( Dttm ),
                     sum( Dttt ),
                     sum( Dtte ),
                     sum( Dttp ),
                     
                     sum( Deem ),
                     sum( Deet ),
                     sum( Deee ),
                     sum( Deep ),
                     
                     sum( Dppm ),
                     sum( Dppt ),
                     sum( Dppe ),
                     sum( Dppp ),
                     
                     sum( Dmte ),
                     sum( Dmep ), 
                     sum( Dpet ), 
                     sum( Dmtp ) )             
      }
      
    }
    
  }
  
  return( out )
  
}