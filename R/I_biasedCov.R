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
  DllDeta <- as.matrix( .llkGrads(gObj = fit, X = X, jj = lpi, type = "DllkDeta") )
  ESS <- sum(abs(DllDeta[ , 1]))^2 / sum(DllDeta[ , 1]^2)
  
  # If DlDeta and X are independent: cov(DlDeta * X) = E(X*X^T) * E(DlDeta^2) - E(DlDeta)^2 * E(X)*E(X)^T 
  if( is.null(lpi) ){ # ELF OR...
    V <- EXXT * mean(DllDeta^2) - mean(DllDeta)^2 * EXEXT
  } else { # ... ELFLSS
    j1 <- lpi[[1]]; j2 <- lpi[[2]]
    Vmm <- mean(DllDeta[ , 1]^2)*EXXT[j1, j1] - mean(DllDeta[ , 1])^2*EXEXT[j1, j1]  
    Vss <- mean(DllDeta[ , 2]^2)*EXXT[j2, j2] - mean(DllDeta[ , 2])^2*EXEXT[j2, j2]
    Vms <- mean(DllDeta[ , 1]*DllDeta[ , 2])*EXXT[j1, j2] -  mean(DllDeta[ , 1])*mean(DllDeta[ , 2])*EXEXT[j1, j2]
    V <- cbind(rbind(Vmm, t(Vms)), rbind(Vms, Vss)) 
  }

  alpha <- min(ESS / npar^2, 1)
  
  return( list("V" = V, "alpha" = alpha ) )
}