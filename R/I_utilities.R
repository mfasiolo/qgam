
.ctrlSetup <- function(innerCtrl, outerCtrl, verbose = TRUE)
{
  if(length(outerCtrl))
  {
    namOut <- names(outerCtrl)
    namIn <-  names(innerCtrl)
    
    if (verbose && length(noNms <- namOut[! namOut %in% namIn])) {
      
      warning("unknown names in control list: ", paste(noNms, collapse = ", "), ". They will not be used")
      
    }
    
    if(length(outerCtrl)) innerCtrl[namOut] <- outerCtrl
  }
  
  return(innerCtrl)
}




#### Andreson-Darling test for STANDARD normality
.adTest <- function(.x){
  
  n <- length(.x)
  .x <- sort(.x)
  
  # Cramer-von Mises statistic 
  # out <- 1/(12*n) + sum( ((2*1:n - 1)/(2*n) - pnorm(.x))^2 )
  
  logp1 <- pnorm(.x, log.p = TRUE)
  logp2 <- pnorm(.x, lower.tail = F, log.p = TRUE)
  
  h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
  out <- -n - mean(h)
  
  return( out )
}