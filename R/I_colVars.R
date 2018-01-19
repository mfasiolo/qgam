
# Variance of columns of a matrix
.colVars <- function(.x){     
  .m <- colMeans(.x)
  .vr <- rowSums( (t(.x) - .m)^2 ) / (nrow(.x) - 1)  
  return(.vr)
}  