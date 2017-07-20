##########################
#' Calculating log(1+exp(x)) accurately
#' 
#' @description Calculates \code{log(1+exp(x))} in a numerically stable fashion.
#'  
#' @param x a numeric vector. 
#' @return A numeric vector where the i-th entry is equal to \code{log(1+exp(x[i]))}, but computed more stably.
#' @details We follow the recipe of Machler (2012), that is formula (10) page 7. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Machler, M. (2012). Accurately computing log(1-exp(-|a|)). 
#'             URL: \url{https://cran.r-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}.
#' @examples
#' set.seed(141)
#' library(qgam); 
#' x <- rnorm(100, 0, 100)
#' log1pexp(x) - log1p(exp(x))
log1pexp <- function(x)
{
  indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), T)
  
  kk <- which(indx==1)
  if( length(kk) ){  x[kk] <- exp(x[kk])  }
  
  kk <- which(indx==2)
  if( length(kk) ){  x[kk] <- log1p( exp(x[kk]) ) }
  
  kk <- which(indx==3)
  if( length(kk) ){  x[kk] <- x[kk] + exp(-x[kk]) }
  
  return(x)
}

