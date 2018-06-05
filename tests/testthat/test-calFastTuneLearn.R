
context("calFastTuneLearn")

test_that("calFastTuneLearn", {
  
  set.seed(414)
  #par(mfrow = c(2, 2))
  #par(mar = c(5.1, 4.1, 0.1, 0.1))
  for(ii in 1:1){ #### !!!!!!!!!!   set to 1:4 to test also elfss
    if(ii == 1){
      ### 1) 4D Gaussian example
      dat <- gamSim(1, n=1000, dist="normal", scale=2, verbose=FALSE)
      form <- y ~ s(x0)+s(x1)+s(x2)+s(x3)
      lsig <- seq(-5.5, 4, length.out = 15)
      qus <- c(0.01, 0.5, 0.99)
      err <- 0.1
    }
    
    if(ii == 2){
      ### 2) 1D Gamma esample
      n <- 1000
      x <- seq(-4, 4, length.out = n)
      X <- cbind(1, x, x^2)
      beta <- c(0, 1, 1)
      sigma <- 1
      # sigma =  .1+(x+4)*.5 ## sigma grows with x
      f <- drop(X %*% beta)
      tauSim <- 0.9
      y <- f + rgamma(n, 3, 1)# rlf(n, 0, tau = tauSim, sig = sigma, lam)#  #  # rnorm(n, 0, sigma)
      form <- y ~ s(x, k = 30)
      dat <- data.frame(cbind(y, x))
      names(dat) <- c("y", "x")
      lsig <- seq(-5, 3, length.out = 15)
      qus <- c(0.01, 0.5, 0.95)
      err <- 0.1
    }
    
    if( ii == 3 ){
      ### 3) 3D Gamma esample
      n <- 1000
      x <- runif(n, -4, 4); z <- runif(n, -8, 8); w <- runif(n, -4, 4)
      X <- cbind(1, x, x^2, z, sin(z), w^3, cos(w))
      beta <- c(0, 1, 1, -1, 2, 0.1, 3)
      sigma <- 0.5
      f <- drop(X %*% beta)
      dat <- f + rgamma(n, 3, 1)
      dat <- data.frame(cbind(dat, x, z, w))
      names(dat) <- c("y", "x", "z", "w")
      bs <- "cr"
      formF <- y~s(x, k = 30, bs = bs) + s(z, k = 30, bs = bs) + s(w, k = 30, bs = bs)
      lsig <- seq(-3, 4, length.out = 15)
      qus <- c(0.01, 0.5, 0.95)
      err <- 0.1
    }
    
    if(ii == 4){
      ### 1) 4D Gaussian example BUT gamlss version
      dat <- gamSim(1, n=1000, dist="normal", scale=2, verbose=FALSE)
      form <- list(y ~ s(x0)+s(x1)+s(x2)+s(x3), ~ s(x0))
      lsig <- seq(-5.5, 4, length.out = 15)
      qus <- c(0.01, 0.5, 0.99)
      err <- 0.1
    }
    
    expect_error({
      calibr <- list("fast" = list(), "slow" = list())
      for(met in c("calFast", "cal")){
        calibr[[met]] <- lapply(qus, 
                                function(.q){
                                  tuneLearn(form,
                                            data = dat,
                                            lsig = lsig,
                                            qu = .q,
                                            err = err, 
                                            control = list("loss" = met, "progress" = "none"))})
      }
    } , NA)
    
    tmp <- cbind(sapply(calibr[["calFast"]], "[[", "loss"), sapply(calibr[["cal"]], "[[", "loss")) 
    matplot(lsig, log(tmp),  type = 'l', lty = c(1:3, 1:3), col = c(1, 1, 1, 2, 2, 2), ylab = "log-loss", 
            xlab = expression(log(sigma)))
    legend("topright", col = 1:2, legend = c("calFast", "cal"), lty = 1)
    legend("top", lty = 1:3, legend = qus)
  
  }
})
