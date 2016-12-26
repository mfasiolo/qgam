# Checking degrees of freedom: poor man's version of mgcv:::k.check: here there is no randomized
# test, because I am not sure it applies to quantile regression.
# This function only gives 1) maximum EDF 2) effective EDF
# Some of the code does useless stuff, 
.kcheck <- function(b) {
  ## function to check k in a gam fit... 
  ## does a randomization test looking for evidence of residual 
  ## pattern attributable to covariates of each smooth. 
  m <- length(b$smooth)
  if (m==0) return(NULL)

  kc <- edf<- rep(0,m)
  snames <- rep("",m)
  n <- nrow(b$model)
  
  for (k in 1:m) { ## work through smooths
    ok <- TRUE
    b$smooth[[k]]$by <- "NA" ## can't deal with by variables
    snames[k] <- b$smooth[[k]]$label
    ind <- b$smooth[[k]]$first.para:b$smooth[[k]]$last.para
    kc[k] <- length(ind)
    edf[k] <- sum(b$edf[ind]) 
  }
  k.table <- cbind(kc,edf)
  dimnames(k.table) <- list(snames, c("k\'","edf"))
  k.table
}