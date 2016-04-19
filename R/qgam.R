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
  nt <- length(qu)
  
  # Initial Gaussian fit
  gausFit <- gam(form, data = data, control = controlGam)
  control[["gausFit"]] <- gausFit
  
  # Output list
  out <- list()
  
  if( is.null(lsigma) ) { # Selecting the learning rate sigma OR ....
    learn <- tuneLearnFast(form = form, data = data, err = err, qu = qu, ncores = ncores, control = control)
    lsigma <- learn$lsigma
    out[["calibr"]] <- learn
  } else { # ... use the one provided by the user
    if( length(lsigma) == 1 ) {
      lsigma <- rep(lsigma, nt)
    } else {
      if( length(lsigma) != nt ) stop("lsigma should either be scalar of a vector of length(qu) ")
    } }
  

  # Selection lambda
  lam <- err * sqrt(2*pi*gausFit$sig2) / (2*log(2)*exp(lsigma))
  
  # Fitting a quantile model for each qu
  out[["fit"]] <- lapply(1:nt, function(ii){
    
    .out <- gam(form, family = logF(qu = qu[ii], lam = lam[ii], theta = lsigma[ii]), data = data)
    
    # Removing data and smooth matrix to reduce memory requirements. There quantities
    # are kept only inside the first fit ( qfit[[1]] )
    if(ii > 1){
      .out$model  <- NULL
      .out$smooth <- NULL 
    } 
        
    return( .out )
  })
  
  # Storing output list
  names(out[["fit"]]) <- qu
  out[["model"]] <- out[["fit"]][[1]][["model"]]
  out[["smooth"]] <- out[["fit"]][[1]][["smooth"]]
  out[["fit"]][[1]][["model"]] <- NULL
  out[["fit"]][[1]][["smooth"]] <- NULL
  
  out[["qu"]] <- qu
  out[["lambda"]] <- lam
  out[["lsigma"]] <- lsigma
  
  return( out )
}


