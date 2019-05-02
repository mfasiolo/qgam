context("check_qgam")

# par(mfrow = c(1, 1))
# test_that("check_qgam_gamlss", {
#   
#   #set.seed(857758)
#   n <- 1000
#   x <- seq(-4, 3, length.out = n)
#   X <- cbind(1, x, x^2)
#   beta <- c(0, 1, 1)
#   sigma =  1.2 + sin(2*x)
#   f <- drop(X %*% beta)
#   dat <- f + rnorm(n, 0, sigma)
#   dataf <- data.frame(cbind(dat, x))
#   names(dataf) <- c("y", "x")
#   
#   expect_error({
#     
#     fit <- qgam(list(y~s(x, k = 15, bs = "cr"), ~ s(x, k = 15, bs = "cr")), data=dataf, err = 0.1, qu = 0.8, 
#                 control = list("progress" = FALSE))
#     invisible(capture.output( check(fit) ))
#     
#   } , NA)
#   
# })


test_that("check_qgam_egam", {
  
  set.seed(57576)
  dat <- gamSim(1,n=1000,dist="normal",scale=2, verbose = FALSE)
  
  expect_error({
    
    fit <- qgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, qu = 0.9, control = list("progress" = FALSE))
    invisible(capture.output( check(fit) ))
    
  } , NA)
  
})