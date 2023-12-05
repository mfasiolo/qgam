##
#
.should_we_use_discrete <- function(form, discrete){
  
  if( !discrete ){
    return( discrete )
  }

  if( !is.list(form) ){
    form <- list(form)
  }
  nlp <- length(form)
  
  form <- mgcv::interpret.gam(form)
  n_smooth <- sapply(form[1:nlp], function(.x) length(.x$smooth.spec))
  
  if( n_smooth[1] == 0 || sum(n_smooth) == 0 ){ # NOTE: this is specific for QGAMs!
    warning("The quantile model does not contain any smooth, ignoring `discrete = TRUE'")
    discrete <- FALSE
  }
  
  return(discrete)
}