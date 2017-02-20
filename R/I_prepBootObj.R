.prepBootObj <- function(obj, eps, control)
{
  if( is.null(control) ){ control <- list() }
  ctrl <- do.call("gam.control", control)
  ctrl$epsilon <- eps
  obj$control <- ctrl
  
  if (inherits(obj$family,"general.family")) {
    obj$Sl <- mgcv:::Sl.setup(obj) ## prepare penalty sequence
    obj$X <- mgcv:::Sl.initial.repara(obj$Sl,obj$X,both.sides=FALSE) ## re-parameterize accordingly
  }
  
  obj$rS <- mgcv:::mini.roots(obj$S, obj$off, ncol(obj$X), obj$rank)
  Ssp <- mgcv:::totalPenaltySpace(obj$S,obj$H,obj$off,ncol(obj$X))
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