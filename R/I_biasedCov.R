##########
# Internal function estimates the covariance matrix of the score function, under the assumption that
# Dllk/Deta is independent of the covariate vector x.
# INPUT
# - fit: a gamObject fitted using the elf or elflss family
# - EXXT: E(xx^T)
# - EXEXT: E(x)E(x)^T
# - X: the model matrix
# 
# OUTPUT
# - the estimate of cov(DlDeta * X) and the mixture parameter alpha in (0, 1) to be used.
#
.biasedCov <- function(fit, X, EXXT, EXEXT, mObj)
{
  npar <- length( coef(fit) )
  DllDeta <- .llkGrads(gObj = fit, X = X, mObj = mObj, type = "DllkDeta")
  ESS <- sum(abs(DllDeta$l1[ , 1]))^2 / sum(DllDeta$l1[ , 1]^2)
  
  # If DlDeta and X are independent: cov(DlDeta * X) = E(X*X^T) * E(DlDeta^2) - E(DlDeta)^2 * E(X)*E(X)^T 
  V <- EXXT * mean(DllDeta$l1^2) - mean(DllDeta$l1)^2 * EXEXT
  
  alpha <- min(ESS / npar^2, 1)
  
  return( list("V" = V, "alpha" = alpha ) )
}