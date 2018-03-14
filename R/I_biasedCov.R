##########
# Internal function estimates the covariance matrix of the score function, under the assumption that
# Dllk/Deta is independent of the covariate vector x.
# INPUT
# - fit: a gamObject fitted using the elf or elflss family
# - EXXT: E(xx^T)
# - EXEXT: E(x)E(x)^T
# - X: the full "lpmatrix" corresponding to both linear predictors
# - lpi: a list of indexes indicating the position of the coefficients corresponding to each linear predictor
# 
# OUTPUT
# - the estimate of cov(DlDeta * X) and the mixture parameter alpha in (0, 1) to be used.
#
.biasedCov <- function(fit, X, EXXT, EXEXT, lpi)
{
  npar <- length( coef(fit) )
  DllDeta <- .llkGrads(gObj = fit, X = X, jj = lpi, type = "DllkDeta")
  ESS <- sum(abs(DllDeta$l1[ , 1]))^2 / sum(DllDeta$l1[ , 1]^2)
  
  # If DlDeta and X are independent: cov(DlDeta * X) = E(X*X^T) * E(DlDeta^2) - E(DlDeta)^2 * E(X)*E(X)^T 
  if( is.null(lpi) ){ # ELF OR...
    V <- EXXT * mean(DllDeta$l1^2) - mean(DllDeta$l1)^2 * EXEXT
  } else { # ... ELFLSS
    sig <- drop( DllDeta$sig )
    Deta <- DllDeta$l1 
    Deta[ , 1] <- Deta[ , 1] * sig
    
    j1 <- lpi[[1]]; j2 <- lpi[[2]]
    X[ , j1] <- X[ , j1] / sig
    EXXT <- crossprod(X, X) / nrow(X)
    EXEXT <- tcrossprod(colMeans(X), colMeans(X))
    
    Vmm <- mean(Deta[ , 1]^2) * EXXT[j1, j1] - mean(Deta[ , 1])^2 * EXEXT[j1, j1]
    Vss <- mean(Deta[ , 2]^2) * EXXT[j2, j2] - mean(Deta[ , 2])^2 * EXEXT[j2, j2]
    Vms <- mean(Deta[ , 1]*Deta[ , 2]) * EXXT[j1, j2] -  mean(Deta[ , 1]) * mean(Deta[ , 2]) * EXEXT[j1, j2]
    V <- cbind(rbind(Vmm, t(Vms)), rbind(Vms, Vss)) 
  }

  alpha <- min(ESS / npar^2, 1)

  return( list("V" = V, "alpha" = alpha ) )
}