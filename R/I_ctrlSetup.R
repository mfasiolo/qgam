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