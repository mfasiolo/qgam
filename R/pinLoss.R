##########################
#' Pinball loss function
#' 
#' @description Evaluates the pinball loss.
#'  
#' @param y points at which the loss is evaluated. 
#' @param mu location parameter of the pinball loss.
#' @param qu quantile level of the loss.
#' @param add if TRUE the losses at which quantile level will be added up.
#' @return A numeric vector or matrix of evaluate losses.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' n <- 1000
#' x <- seq(0, 4, length.out = n)
#' plot(x, pinLoss(x, rep(2, n), qu = 0.9, add = FALSE), type = 'l', ylab = "loss")
#' 
pinLoss <- function(y, mu, qu, add = TRUE){
  
  # Recursive call for multiple quantiles
  if( length(qu) > 1 ){  
    n <- length( qu )
    l <- sapply(1:n,
                function(ii){
                  return( pinLoss(y, mu[ , ii], qu[ii], add = add) )
                })
    
    if( is.matrix(l) ){ colnames(l) <- qu } else { names(l) <- qu }
    
    return( l )
  }
  
  tau <- 1 - qu
  d <- y - mu
  l <- d * 0
  
  l[d < 0] <- - tau * d[ d < 0 ]
  l[d > 0] <- - (tau-1) * d[ d > 0 ]
  
  if( add ){ l <- sum(l) }
  
  return( l )
}