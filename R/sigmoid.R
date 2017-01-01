##########################
#' Sigmoid function and its derivatives
#' 
#' @description Calculates the sigmoid function and its derivatives.
#' @param y a numeric vector. 
#' @param deriv if \code{TRUE} alse the first three derivatives of the sigmoid function will be computed.
#' @return If \code{deriv==FALSE}, it returns a numeric vector equal to \code{1/(1+exp(-x))}. If
#'         \code{deriv==TRUE} it returns a list where the slot \code{$D0} contains \code{1/(1+exp(-x))}, 
#'         while \code{$D1}, \code{$D2} and \code{$D3} contain its first three derivatives.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' library(qgam)
#' set.seed(90)
#' h <- 1e-6
#' p <- rnorm(1e4, 0, 1e6)
#' sigmoid(p[1:50]) - 1/(1+exp(-p[1:50]))
#' 
#' ##### Testing sigmoid derivatives
#' e1 <- abs((sigmoid(p+h) - sigmoid(p-h)) / (2*h) - sigmoid(p, TRUE)[["D1"]]) / (2*h)
#' e2 <- abs((sigmoid(p+h, TRUE)$D1 - sigmoid(p-h, TRUE)$D1) / 
#'       (2*h) - sigmoid(p, TRUE)[["D2"]]) / (2*h)
#' e3 <- abs((sigmoid(p+h, TRUE)$D2 - sigmoid(p-h, TRUE)$D2) / 
#'       (2*h) - sigmoid(p, TRUE)[["D3"]]) / (2*h)
#' 
#' if( any(c(e1, e2, e3) > 1) ) stop("Sigmoid derivatives are not estimated accurately")
#'
#' @export sigmoid
#'
sigmoid <- function(y, deriv = FALSE)
{
  l0 <- plogis(y, 0, 1)
  if( deriv ){
    l1 <- l0 * (1-l0)
    l2 <- l1 - 2*l1*l0
    l3 <- l2 - 2*l2*l0 - 2*l1*l1
    out <- list("D0" = l0, "D1" = l1, "D2" = l2, "D3" = l3)
    return( out )
  } else {
    return( l0 )
  }
}