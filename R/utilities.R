###########
# Utilities
##########

####### Check loss function
checkloss <- function(y, mu, qu){
  
  tau <- 1 - qu
  
  d <- y - mu
  
  l <- - sum( tau*d[d<0] ) - sum( (tau-1)*d[d>0] )
  
  return( l )
  
} 

#### Vettorize check Loss function
checklossVett <- function(y, mu, p){
  
  n <- length( p )
  
  out <- sapply(1:n,
                function(ii){
                  return( checkloss(y, mu[ , ii], p[ii]) )
                })
  
  return( out )
}


#### Vettorized empirical cdf
qqVett <- function(y, mu){
  
  nq <- ncol( mu )
  nobs <- length( y )
  
  out <- sapply(1:nq,
                function(ii){
                  return( sum( (y - mu[ , ii]) < 0 ) / nobs )
                })
  
  return( out )
}

#### Does stuff like "predict", "aic", etcetera
qdo <- function(obj, qu, fun, ...){
  
  if( !(qu %in% obj[["qu"]]) ) stop("qu is not in obj[[\"qu\"]].")
    
  tmpObj <- obj[["fit"]][[ which(obj[["qu"]] == qu) ]]
  
  tmpObj[["model"]] <- obj[["model"]]
  tmpObj[["smooth"]] <- obj[["smooth"]]
  
  out <- fun(tmpObj, ...)
  
  return( out )
}