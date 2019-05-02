#######
# Fitting the SHASH density using optim.
#######
#
.fitShash <- function(x){
  
  shashF0 <- function(param, x){
    
    - .llkShash(x = x, mu = param[1], tau = param[2], eps = param[3], phi = param[4], deriv = 0)$l0
    
  } 
  
  shashF1 <- function(param, x){
    
    - .llkShash(x = x, mu = param[1], tau = param[2], eps = param[3], phi = param[4], deriv = 1)$l1
  } 
  
  shashF2 <- function(param, x){
    
    - .llkShash(x = x, mu = param[1], tau = param[2], eps = param[3], phi = param[4], deriv = 2)$l2

  } 
  
  # Seems to work OK, we could use Newton method in optimx but the Hessian was giving problems
  fitSH <- optim(par = c(mean(x), log(sd(x)), 0, 0), fn = shashF0, gr = shashF1, method = "BFGS", x = x)
  
  return( fitSH )
  
}

# 
# # Testing by finite differences
# library(numDeriv)
# y <- rchisq(1000, df = 3)
# parSH <- qgam:::.fitShash( y )$par
# 
# # Checking gradient
# jacobian(func = function(x){
#   - qgam:::.llkShash(x = y, mu = x[1], tau = x[2], eps = x[3], phi = x[4], deriv = 0)$l0
# }, x = parSH)
# 
# - qgam:::.llkShash(x = y, mu = parSH[1], tau = parSH[2], eps = parSH[3], phi = parSH[4], deriv = 1)$l1
# 
# par(mfrow = c(2, 2))
# # [1] Test of t-distributed
# x <- rt(1000, 5)
# parSH <- qgam:::.fitShash( x )$par
# 
# xseq <- seq(-4, 4, length.out = 1000)
# den <- unlist(sapply(xseq, qgam:::.llkShash, mu = parSH[1], tau = parSH[2], eps = parSH[3], phi = parSH[4]))
# 
# plot(xseq, exp(den), col = 1, ylim = range(exp(den), dt(xseq, 5)))
# lines(xseq, dt(xseq, 5), col = 2)
# abline(v = qgam:::.shashMode(parSH), lty = 2)
# 
# # [2] Test of chi-squ
# x <- rchisq(1000, df = 3)
# parSH <- qgam:::.fitShash( x )$par
# 
# xseq <- seq(0, 6, length.out = 1000)
# den <- unlist(sapply(xseq, qgam:::.llkShash, mu = parSH[1], tau = parSH[2], eps = parSH[3], phi = parSH[4]))
# 
# plot(xseq, exp(den), col = 1, ylim = range(exp(den), dt(xseq, 5)))
# lines(xseq, dchisq(xseq, df = 4), col = 2)
# abline(v = qgam:::.shashMode(parSH), lty = 2)
# 
# # [3] Test of gamma(3, 1)
# x <- rgamma(1000, 3, 1)
# parSH <- qgam:::.fitShash( x )$par
# 
# xseq <- seq(0, 6, length.out = 1000)
# den <- unlist(sapply(xseq, qgam:::.llkShash, mu = parSH[1], tau = parSH[2], eps = parSH[3], phi = parSH[4]))
# 
# plot(xseq, exp(den), col = 1, ylim = range(exp(den), dt(xseq, 5)))
# lines(xseq, dgamma(xseq, 3, 1), col = 2)
# abline(v = qgam:::.shashMode(parSH), lty = 2)
