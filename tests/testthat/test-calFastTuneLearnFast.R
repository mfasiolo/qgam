
context("calFastTuneLearnFast")

test_that("calFastTuneLearnFast", {
  
  set.seed(41334)
  par(mfrow = c(2, 2))
  par(mar = c(5.1, 4.1, 1.1, 0.1))
  for(ii in 1:1){ #### !!!!!!!!!!   set to 1:4 to test also elfss
    if(ii == 1){
      ### 1) 4D Gaussian example
      dat <- gamSim(1, n=2000, dist="normal", scale=2, verbose=FALSE)
      form <- y ~ s(x0)+s(x1)+s(x2)+s(x3)
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
      qus <- c(0.1, 0.95, 0.99)
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
    
    # Checking that the loss evaluated by tuneLearn is close to that evaluated
    # by tuneLearnFast. They can't be exactly the same, because the order with which
    # the losses are evaluated are different (hence different initializations)
    expect_error({
      calibr <- list("fast" = list(), "slow" = list())
      calibr[["fast"]] <-  tuneLearnFast(form,
                                         data = dat,
                                         qu = qus,
                                         control = list("progress" = FALSE))
      
      calibr[["fast_discrete"]] <-  tuneLearnFast(form,
                                                  data = dat,
                                                  qu = qus,
                                                  discrete = TRUE,
                                                  control = list("progress" = FALSE))
      
      calibr[["slow"]] <- lapply(1:length(qus), 
                                 function(.kk){
                                   tuneLearn(form,
                                             data = dat,
                                             qu = qus[.kk],
                                             lsig = calibr[["fast"]]$store[[.kk]][1, ],
                                             control = list("progress" = FALSE))})
    }, NA)
    
    x <- lapply(calibr[["fast"]]$store, function(.inp) .inp[1, ])
    y1 <- lapply(calibr[["fast"]]$store, function(.inp) log(.inp[2, ]))
    y2 <- sapply(calibr[["slow"]], function(.inp) log(.inp$loss))
    
    plot(x[[1]], y1[[1]], col = 1, xlim = range(do.call("c", x)), ylim = range(c(do.call("c", y1), do.call("c", y2))), 
         ylab = "log-loss", xlab = expression(log(sigma)), main = "tuneLearnFast vs tuneLearn")
    points(x[[2]], y1[[2]], col = 2)
    points(x[[3]], y1[[3]], col = 3)
    lines(sort(x[[1]]), y2[[1]], col = 1)
    lines(sort(x[[2]]), y2[[2]], col = 2)
    lines(sort(x[[3]]), y2[[3]], col = 3)
    
    x1 <- lapply(calibr[["fast_discrete"]]$store, function(.inp) .inp[1, ])
    y1 <- lapply(calibr[["fast_discrete"]]$store, function(.inp) log(.inp[2, ]))
    y2 <- sapply(calibr[["slow"]], function(.inp) log(.inp$loss))
    
    plot(x1[[1]], y1[[1]], col = 1, xlim = range(do.call("c", x)), ylim = range(c(do.call("c", y1), do.call("c", y2))), 
         ylab = "log-loss", xlab = expression(log(sigma)), main = "tuneLearnFast_discrete vs tuneLearn")
    points(x1[[2]], y1[[2]], col = 2)
    points(x1[[3]], y1[[3]], col = 3)
    lines(sort(x[[1]]), y2[[1]], col = 1)
    lines(sort(x[[2]]), y2[[2]], col = 2)
    lines(sort(x[[3]]), y2[[3]], col = 3)
    
  }
  
})
