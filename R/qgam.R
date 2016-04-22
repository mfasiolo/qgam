#######
### Fitting quantile gam model
#######
#' Fitting quantile gam model
#' 
#' @param \code{XXX} .
#' @return XXX.
#'
#' @details XXX
#'         
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.   
#' @references Fasiolo and Wood XXX.           
#' @export
#' 

qgam <- function(form, data, qu, lsigma = NULL, err = 0.01, ncores = 1, control = list(), controlGam = list())
{
  if( length(qu) > 1 ) stop("length(qu) > 1: you should use mqgam")
  
  # Initial Gaussian fit
  if( is.null(control[["gausFit"]]) )
  {
    if( is.formula(form) ){
      gausFit <- gam(form, data = data)
    } else {
      gausFit <- gam(form, data = data, family = gaulss(b=0))
    }
    control[["gausFit"]] <- gausFit
  }

  # Selecting the learning rate sigma
  if( is.null(lsigma) ) {  
    learn <- tuneLearnFast(form = form, data = data, err = err, qu = qu,
                           ncores = ncores, control = control, controlGam = controlGam)
    lsigma <- learn$lsigma
  }
  
  # Fit model
  if( is.formula(form) ){ # Extended Gam OR .....
    
    lam <- err * sqrt(2*pi*gausFit$sig2) / (2*log(2)*exp(lsigma))
    
    fit <- gam(form, family = logF(qu = qu, lam = lam, theta = lsigma), data = data, control = controlGam)
    
  } else { # .... Gamlss
    
    lam <- err * sqrt(2*pi/(gausFit$fit[ , 2]^2)) / (2*log(2)*exp(lsigma))
    
    fit <- gam(form, family = logFlss(qu = qu, lam = lam, offset = lsigma), data = data, control = controlGam)
    
  }
  
  return( fit )
}


