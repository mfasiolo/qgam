###############################################
# Function that prepares the bootstrapped gamObject for use 
# with gam.fit4 or gam.fit5. It manly works out the re-parametrization
# to be used, which does not seem to depend on the smoothing parameters
# (hence it only needs to be calculated once).
#
.prepBootObj <- function(obj, eps, control)
{
  if( is.null(control) ){ control <- list() }
  ctrl <- do.call("gam.control", control)
  
  # Overwriting default tolerance, useful for using sloppy convergence test on 
  # bootstrapped fits
  if( !is.null(eps) ) { ctrl$epsilon <- eps }  
  obj$control <- ctrl
  obj$rS <- mini.roots(obj$S, obj$off, ncol(obj$X), obj$rank)
  Ssp <- totalPenaltySpace(obj$S,obj$H,obj$off,ncol(obj$X))
  obj$Eb <- Ssp$E       ## balanced penalty square root for rank determination purposes 
  obj$U1 <- cbind(Ssp$Y,Ssp$Z) ## eigen space basis
  obj$Mp <- ncol(Ssp$Z) ## null space dimension
  obj$UrS <- list()     ## need penalty matrices in overall penalty range space...
  if (length(obj$S)>0) for (i in 1:length(obj$S)) obj$UrS[[i]] <- t(Ssp$Y)%*%obj$rS[[i]] else i <- 0
  if (!is.null(obj$H)) { ## then the sqrt fixed penalty matrix H is needed for (RE)ML 
    obj$UrS[[i+1]] <- t(Ssp$Y)%*%mroot(obj$H)
  }
  obj$family <- fix.family.link(obj$family)
  
  return(obj)
}