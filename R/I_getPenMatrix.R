#######################
# Function that extracts the penalty matrix from an extended gam object.
# It returns St expressed at the "user-level" parametrization (i.e. not the 
# parametrization used internally for fitting)
#
.getPenMatrix <- function(q, UrS, sp, Mp, U1)
{
  
  if (length(UrS)) { # Stable re-parameterization if needed....
    rp <- gam.reparam(UrS, sp, 0)
    T <- diag( q )
    T[1:ncol(rp$Qs), 1:ncol(rp$Qs)] <- rp$Qs
    T <- U1%*%T ## new params b'=T'b old params
    St <- rbind(cbind(rp$S,matrix(0,nrow(rp$S),Mp)),matrix(0,Mp,q))
    
    # Invert re-parametrization
    St <- T %*% St %*% t(T)
  } else { 
    T <- diag(q); 
    St <- matrix(0,q,q) 
  }
  
  return( St ) 

}
