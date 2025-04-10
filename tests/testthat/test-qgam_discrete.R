
context("qgam_discrete")

test_that("qgam_discrete", {
  
  set.seed(414)
  par(mfrow = c(2, 2))
  par(mar = c(5.1, 4.1, 1.1, 0.1))
  for(ii in 1:1){ #### !!!!!!!!!!   set to 1:4 to test also elfss
    if(ii == 1){
      ### 1) 4D Gaussian example
      dat <- gamSim(1, n=2000, dist="normal", scale=2, verbose=FALSE)
      form <- y ~ s(x0)+s(x1)+s(x2)+s(x3)
      lsig <- seq(-5.5, 4, length.out = 15)
      qus <- c(0.01, 0.5, 0.99)
    }
    
    if(ii == 2){
      ### 2) 1D Gamma esample
      n <- 2000
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
      qus <- c(0.01, 0.5, 0.95)
      
    }
    
    if( ii == 3 ){
      ### 3) 3D Gamma esample
      n <- 2000
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
      qus <- c(0.01, 0.5, 0.95)
      
    }
    
    if(ii == 4){
      ### 1) 4D Gaussian example BUT gamlss version
      dat <- gamSim(1, n=2000, dist="normal", scale=2, verbose=FALSE)
      form <- list(y ~ s(x0)+s(x1)+s(x2)+s(x3), ~ s(x0))
      qus <- c(0.01, 0.5, 0.99)
      
    }
    
    expect_error({
      calibr <- list("standard" = list(), "discrete" = list())
      for(met in c("standard", "discrete")){
        calibr[[met]] <- lapply(qus, 
                                function(.q){
                                  qgam(form,
                                       data = dat,
                                       discrete = (met == "discrete"),
                                       qu = .q,
                                       control = list("progress" = "none"))})
      }
    } , NA)
    
    tmp <- cbind(sapply(calibr[["standard"]], "[[", "fitted.values"), sapply(calibr[["discrete"]], "[[", "fitted.values")) 
    ylim <- range(tmp)

    plot(tmp[ , 1], tmp[ , 3+1], main = paste0("[",ii,"] qgam vs qgam discrete"), xlab = "standard fit", ylab = "discrete fit")
    abline(0, 1, col = 2)
    for(kk in 2:3){
     lines(tmp[ , 1], tmp[ , 3+1], col = kk)
    }
    legend("topleft", col = 1:3, lty = 1, legend = qus)
    
    expect_error({
      withCallingHandlers(
      {
      calibr <- list("standard" = list(), "discrete" = list())
      for(met in c("standard", "discrete")){
        calibr[[met]] <- mqgam(form,
                               data = dat,
                               discrete = (met == "discrete"),
                               qu = qus,
                               control = list("progress" = "none"))
        }
      },
        warning = function(w) {
          if (endsWith(conditionMessage(w), "algorithm did not converge") || endsWith(conditionMessage(w), "check results carefully"))
            invokeRestart("muffleWarning")
        })
      }, NA)
    
    tmp <- cbind(sapply(calibr[["standard"]]$fit, "[[", "fitted.values"), sapply(calibr[["discrete"]]$fit, "[[", "fitted.values")) 
    ylim <- range(tmp)
    
    plot(tmp[ , 1], tmp[ , 3+1], main = paste0("[",ii,"] mqgam vs mqgam discrete"), xlab = "standard fit", ylab = "discrete fit")
    abline(0, 1, col = 2)
    for(kk in 2:3){
      lines(tmp[ , 1], tmp[ , 3+1], col = kk)
    }
    legend("topleft", col = 1:3, lty = 1, legend = qus)

  }
})
