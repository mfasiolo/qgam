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
  
  if( !(qu %in% names(obj[["fit"]])) ) stop("qu is not in obj[[\"qu\"]].")
    
  tmpObj <- obj[["fit"]][[ which(names(obj[["fit"]]) == qu) ]]
  
  tmpObj[["model"]] <- obj[["model"]]
  tmpObj[["smooth"]] <- obj[["smooth"]]
  
  out <- fun(tmpObj, ...)
  
  return( out )
}



####### Visual checks for mqgam()

checkMQGam <- function(obj)
{  
  cal <- obj$calibr
  est <- cal$store
  brac <- cal$ranges
  lsigma <- cal$lsigma
  errors <- cal$err
  
  qu <- as.numeric(names(cal$lsigma))
  nq <- length(qu)
  
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), 
         heights=c(2, 1))
  oldPar <- par(mai = c(1, 1, 0.1, 0.1))
  plot(qu, lsigma, ylim = range(as.vector(brac)), xlim = range(qu), col = 2, 
       ylab = expression("Log(" * sigma * ")"), xlab = "Quantile")
  points(qu, brac[ , 1], pch = 3)
  points(qu, brac[ , 2], pch = 3)
  points(qu, rowMeans(brac), pch = 3)
  for(zz in 1:nq) segments(qu[zz], mean(brac[zz, ]) - abs(diff(brac[zz, ]))/4, 
                           qu[zz], mean(brac[zz, ]) + abs(diff(brac[zz, ]))/4, col = 1)
  plot(qu, errors, xlab = "Quantile")
  
  readline(prompt = "Press <Enter> to see the next plot...")
  
  par(oldPar)
  
  pDim <- min( ceiling(sqrt(nq)), 2 )
  par(mfrow = c(pDim, pDim))
  for( ii in 1:nq )
  {
    plot(sort(est[[ii]][1, ]), est[[ii]][2, order(est[[ii]][1, ])], 
         main = substitute(Quantile == x, list(x = round(qu[ii], 3))), 
         ylab = "loss", xlab = expression(log(sigma)), type = 'b')
    abline(v = est[[ii]][1, which.min(est[[ii]][2, ])])
    
    if(ii %% (pDim^2) == 0) readline(prompt = "Press <Enter> to see the next plot...")
  }
  
  par(oldPar)
  
  return( invisible(NULL) )
  
}