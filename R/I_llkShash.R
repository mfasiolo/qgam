#######
# Log-likelihood of shash density and its derivatives
#######
.llkShash <- function(x, mu, tau, eps, phi, deriv = 0){
  
  sech <- function(.x){ 1 / cosh(.x) }
  
  sig <- exp( tau )
  del <- exp( phi )
  
  # 1) Calculate derivative of likelihood to appropriate order
  z <- (x - mu) / (sig*del)
  
  dTasMe <- del*asinh(z) - eps
  g <- -dTasMe
  CC <- cosh( dTasMe )
  SS <- sinh( dTasMe )
  
  l <- sum( -tau - 0.5*log(2*pi) + log(CC) - 0.5*log1pexp(2*log(abs(z))) - 0.5*SS^2 )
  
  out <- list("l0" = l)
  
  # Compute sqrt(x^2 + m) when |x| >> 0 and m is reasonably small (e.g. + 1 or - 1)
  sqrtX2pm <- function(x, m){ 
    x <- abs(x)
    kk <- which( x < 1e8 )
    if( length(kk) ){
      x[kk] <- sqrt(x[kk]^2 + m)
    }
    return(x)
  }
  
  # Compute (x^2 + m1) / (x^2 + m2)^2 when |x| >> 0 and m1, m2 are reasonably small (e.g. + 1 or - 1)
  x2m1DivX2m2SQ <- function(x, m1, m2){
    
    x <- abs(x)
    kk <- (x^2 + m1) < 0
    o <- x * 0
    if( any(kk) ){
      o[kk] <- (x[kk]^2 + m1) / (x[kk]^2 + m2)^2
    }
    if( sum(kk) < length(x) ){
      o[!kk] <-  ((sqrtX2pm(x[!kk], m1) / sqrtX2pm(x[!kk], m2)) / sqrtX2pm(x[!kk], m2))^2
    }
    
    return(o)
  }
  
  if( deriv > 0 )
  {
    zsd <- z*sig*del
    sSp1 <- sqrtX2pm(z, 1) # sqrt(z^2+1)
    asinhZ <- asinh(z)
    
    ## First derivatives 
    De <- tanh(g) - 0.5*sinh(2*g)
    Dm <- 1/(del*sig*sSp1)*(del*(De)+z/sSp1)
    Dt <- zsd*Dm - 1
    Dp <- Dt + 1 - del*asinhZ*De
    
    out$l1 <- c(sum(Dm), sum(Dt), sum(De), sum(Dp))
    
    if( deriv > 1 ){
      Dme <- (sech(g)^2 - cosh(2*g)) / (sig*sSp1)
      Dte <- zsd*Dme
      Dmm <- Dme/(sig*sSp1) + z*De/(sig^2*del*sSp1^3) + x2m1DivX2m2SQ(z, -1, 1)/(del*sig*del*sig)
      Dmt <- zsd*Dmm - Dm
      Dee <- -2*cosh(g)^2 + sech(g)^2 + 1 
      Dtt <-  zsd*Dmt
      Dep <- Dte - del*asinhZ*Dee
      Dmp <- Dmt + De/(sig*sSp1) - del*asinhZ*Dme
      Dtp <- zsd*Dmp
      Dpp <- Dtp - del*asinhZ*Dep + del*(z/sSp1-asinhZ)*De
      
      out$l2 <- matrix(c(sum(Dmm), sum(Dmt), sum(Dme), sum(Dmp), 
                         sum(Dmt), sum(Dtt), sum(Dte), sum(Dtp), 
                         sum(Dme), sum(Dte), sum(Dee), sum(Dep), 
                         sum(Dmp), sum(Dtp), sum(Dep), sum(Dpp)), 4, 4, byrow = TRUE)
      
    }
  }
    
    return( out )
    
}

