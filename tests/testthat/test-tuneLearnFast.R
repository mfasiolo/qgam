context("tuneLearnFast")

# test_that("tuneLearnFast_gamlss", {
# 
#   #set.seed(651)
#   n <- 1000
#   x <- seq(-4, 3, length.out = n)
#   X <- cbind(1, x, x^2)
#   beta <- c(0, 1, 1)
#   sigma =  1.2 + sin(2*x)
#   f <- drop(X %*% beta)
#   dat <- f + rnorm(n, 0, sigma)
#   dataf <- data.frame(cbind(dat, x))
#   names(dataf) <- c("y", "x")
#   form <- list(y~s(x, k = 30, bs = "cr"), ~ s(x, k = 30, bs = "cr"))
# 
#   QU <- 0.9
#   lossType <- rep(c("calFast", "cal", "pin"), each = 2)
# 
#   par(mfrow = c(3, 2))
#   par(mar = c(5.1, 4.1, 0.1, 2.1))
#   for(ii in c(1, 3, 5)){ # Set to 1:6 if you want to test all calibration methods
# 
#     expect_error({ # Actually we expect NO error!!
#       tun <- tuneLearnFast(form, data = dataf, qu = QU,
#                            control = list("loss" = lossType[ii], "progress" = FALSE, "K" = 20),
#                            multicore = ((ii %% 2) == 0), ncores = 2)
# 
#       fit <- qgam(form, qu = QU, lsig = tun$lsig, data = dataf)
# 
#       ylb <- if((ii %% 2) == 0) { paste(lossType[ii], "multicore") } else { lossType[ii] }
#       plot(x, dat, col = "grey", ylab = ylb)
#       tmp <- predict(fit, se = TRUE)
#       lines(x, tmp$fit[ , 1])
#       lines(x, tmp$fit[ , 1] + 3 * tmp$se.fit[ , 1], col = 2)
#       lines(x, tmp$fit[ , 1] - 3 * tmp$se.fit[ , 1], col = 2)
#       
#       plot(sort(tun$store[[1]][1, ]), tun$store[[1]][2, ][order(tun$store[[1]][1, ])], type = "b", 
#            ylab = "Calibration loss", xlab = expression(log(sigma)))
#     }
#     , NA)
# 
#   }
# 
# })




test_that("tuneLearnFast_egam", {
  
  set.seed(211)
  dataf <- gamSim(1,n=400,dist="normal",scale=2,verbose=FALSE)
  form <- y~s(x0)+s(x1)+s(x2)+s(x3)
  
  QU <- 0.9
  lossType <- rep(c("calFast", "cal", "pin"), each = 2)
  
  #par(mfrow = c(3, 2))
  for(ii in 1:2){
    
    expect_error({ # Actually we expect NO error!!
      tun <- tuneLearnFast(form, data = dataf, qu = QU,
                           control = list("loss" = lossType[ii], "K" = 20, "progress" = FALSE), 
                           multicore = ((ii %% 2) == 0), ncores = 2)
      
      fit <- qgam(form, qu = QU, lsig = tun$lsig, data = dataf)
      
      ylb <- if((ii %% 2) == 0) { paste(lossType[ii], "multicore") } else { lossType[ii] }
      plot(fit, select = 3, ylab = ylb)
      
      plot(sort(tun$store[[1]][1, ]), tun$store[[1]][2, ][order(tun$store[[1]][1, ])], type = "b", 
           ylab = "Calibration loss", xlab = expression(log(sigma)))
    }
    , NA)
    
  }
  
})

