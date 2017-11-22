
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/qgam)](https://cran.r-project.org/package=qgam)
[![Build Status](https://travis-ci.org/mfasiolo/qgam.svg?branch=master)](https://travis-ci.org/mfasiolo/qgam)

This R package offers methods for fitting additive quantile regression models based on splines, using the methods described in [Fasiolo et al., 2017](https://arxiv.org/abs/1707.03307).

The main functions are:
- `qgam` fits an additive quantile regression model to a single quantile. Very similar to `mgcv::gam`. It returns an `mgcv::gamObject`.
- `mqgam` fits the same additive quantile regression model to several quantiles. It is more efficient that calling `qgam` several time, 
          especially in terms of memory.
- `tuneLearn` useful for tuning the learning rate of the Gibbs posterior. It evaluates a calibration loss function on a grid of values 
              provided by the user. 
- `tuneLearnFast` similar to `tuneLearn`, but here the learning rate is selected by minimizing the calibration loss using Brent method.

## A first example: smoothing the motorcycle dataset

Let's start with a simple example. Here we are fitting a regression model with an adaptive spline basis to quantile 0.8 of the motorcycle dataset.
```R
# Maybe you need to install first: library(devtools); install_github("mfasiolo/qgam")
library(qgam); library(MASS)
require(RhpcBLASctl); blas_set_num_threads(1) #optional

set.seed(6436)
fit <- qgam(accel~s(times, k=20, bs="ad"), 
            data = mcycle, 
            err = 0.1, 
            qu = 0.8, 
            control = list("tol" = 0.01)) # <- sloppy tolerance to speed-up calibration 

# Plot the fit
xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
pred <- predict(fit, newdata = xSeq, se=TRUE)
plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
lines(xSeq$times, pred$fit, lwd = 1)
lines(xSeq$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
lines(xSeq$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)   
```

`qgam` automatically calls `tuneLearnFast` to select the learning rate. We can check whether the optimization 
succeded.
```R
check(fit$calibr)
```

The second plot suggests so. Alternatively, we could have selected the learning rate by evaluating the loss function on a grid.
```R
set.seed(6436)
cal <- tuneLearn(accel~s(times, k=20, bs="ad"), 
                 data = mcycle, 
                 err = 0.1, 
                 qu = 0.8,
                 lsig = seq(1, 3, length.out = 20)) #<- sequence of values for learning rate
                 
check(cal)
```


We might want to fit several quantiles at once. This can be done with `mqgam`.
```R
quSeq <- c(0.2, 0.4, 0.6, 0.8)
set.seed(6436)
fit <- mqgam(accel~s(times, k=20, bs="ad"), 
             data = mcycle, 
             err = 0.1, qu = quSeq, 
             control = list("tol" = 0.01)) # <- sloppy tolerance to speed-up calibration 
```

To save memory `mqgam` does not return one `mgcv::gamObject` for each quantile, but it avoids storing some redundant data (such as several copies of the design matrix). The output of `mqgam` can be manipulated using the `qdo` function.

```R
# Plot the data
xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))

# Predict each quantile curve and plot
for(iq in quSeq){
  pred <- qdo(fit, iq, predict, newdata = xSeq)
  lines(xSeq$times, pred, col = 2)
}

# Summary for quantile 0.4
qdo(fit, qu = 0.4, summary)

# Check learning rate optimization
check(fit$calibr)
```

## Dealing with heteroscedasticity

Let us simulate some data from an heteroscedastic model.
```R
set.seed(651)
n <- 5000
x <- seq(-4, 3, length.out = n)
X <- cbind(1, x, x^2)
beta <- c(0, 1, 1)
sigma =  1.2 + sin(2*x)
f <- drop(X %*% beta)
dat <- f + rnorm(n, 0, sigma)
dataf <- data.frame(cbind(dat, x))
names(dataf) <- c("y", "x")
   
qus <- seq(0.05, 0.95, length.out = 10)
plot(x, dat, col = "grey", ylab = "y")
for(iq in qus){ lines(x, qnorm(iq, f, sigma)) }
```

Let's fit some quantiles using an extended quantile GAM model with constant learning rate. To speed up things I've
pre-computed the learning rate. Just comment out the line `lsig = lsig,` if you want to re-computed it.
```R
lsig <- c(-0.96, -0.83, -0.69, -0.63, -0.76, -0.76, -0.89, -0.85, -0.99, -1.06)
fit <- mqgam(y~s(x, k = 30, bs = "cr"), 
             data = dataf,
             lsig = lsig,
             qu = qus, err = 0.05)
             
qus <- seq(0.05, 0.95, length.out = 10)
plot(x, dat, col = "grey", ylab = "y")
for(iq in qus){ 
 lines(x, qnorm(iq, f, sigma), col = 2)
 lines(x, qdo(fit, iq, predict))
}
legend("top", c("truth", "fitted"), col = 2:1, lty = rep(1, 2))
```

The fitted quantiles are close to the true ones, but their credible intervals don't vary much with x. Indeed, let's look at intervals for quantile 0.95.
```R
plot(x, dat, col = "grey", ylab = "y")
tmp <- qdo(fit, 0.95, predict, se = TRUE)
lines(x, tmp$fit)
lines(x, tmp$fit + 3 * tmp$se.fit, col = 2)
lines(x, tmp$fit - 3 * tmp$se.fit, col = 2)
```

We can do better by letting the learning rate vary with the covariate. Here I am fixing the intercept of 
the additive model for the learning rate, in order to avoid calibrating it.
```R
fit <- qgam(list(y~s(x, k = 30, bs = "cr"), ~ s(x, k = 30, bs = "cr")), 
            data = dataf, qu = 0.95, err = 0.05, lsig = -1.16)

plot(x, dat, col = "grey", ylab = "y")
tmp <- predict(fit, se = TRUE)
lines(x, tmp$fit[ , 1])
lines(x, tmp$fit[ , 1] + 3 * tmp$se.fit[ , 1], col = 2)
lines(x, tmp$fit[ , 1] - 3 * tmp$se.fit[ , 1], col = 2)
```

This looks much better. 


## Application to UK electricity load forecasting

Here we consider a UK electricity demand dataset, taken from the national grid [website](https://www.nationalgrid.com). The dataset covers the period January 2011 to June 2016 and it contains the following variables:

   - `NetDemand` net electricity demand between 11:30am and 12am.
   - `wM` instantaneous temperature, averaged over several English cities.
   - `wM_s95` exponential smooth of `wM`, that is `wM_s95[i] = a*wM[i] + (1-a)*wM_s95[i]` with `a=0.95`.
   - `Posan` periodic index in `[0, 1]` indicating the position along the year.
   - `Dow` factor variable indicating the day of the week.
   - `Trend` progressive counter, useful for defining the long term trend.
   - `NetDemand.48` lagged version of `NetDemand`, that is `NetDemand.48[i] = NetDemand[i-2]`.
   - `Holy` binary variable indicating holidays.
   - `Year` and `Date` should obvious, and partially redundant.
   
See [Fasiolo et al., 2017](https://arxiv.org/abs/1707.03307) for more details. This is how the demand over the period looks like:

```R
data("UKload")
EdfColors <- c("#FE5815", "#FFA02F", "#C4D600", "#509E2F", "#005BBB", "#001A70")

tmpx <- seq(UKload$Year[1], tail(UKload$Year, 1), length.out = nrow(UKload)) 
plot(tmpx, UKload$NetDemand, type = 'l', xlab = 'Year', ylab = 'Load')
```

Now we tune the learning rate on a grid, on 4 cores. We are interested the median.
```R
qu <- 0.5
form <- NetDemand~s(wM,k=20,bs='cr') + s(wM_s95,k=20,bs='cr') + 
        s(Posan,bs='ad',k=30,xt=list("bs"="cc")) + Dow + s(Trend,k=4) + NetDemand.48 + Holy

library(parallel)
cl <- makeCluster(4)
clusterEvalQ(cl, {library(RhpcBLASctl); blas_set_num_threads(1)})

tic <- proc.time()
set.seed(41241)
sigSeq <- seq(4, 8, length.out = 16)
closs <- tuneLearn(form = form, data = UKload, err = 0.1,
                   lsig = sigSeq, qu = qu, control = list("K" = 20), 
                   multicore = TRUE, cluster = cl)
proc.time() - tic
stopCluster(cl)

check(closs)
```

Let's fit the model with the learning rate corresponding to the lowest loss and let's look at the resulting smooth effects.
```R
lsig <- closs$lsig
fit <- qgam(form = form, data = UKload, err = 0.1, lsig = lsig, qu = qu)
plot(fit, scale = F, page = 1)
```

We can fit multiple quantiles using `mqgam`. Here we use pre-computed a learning-rate for each quantile, in order to save time.
```R
qus <- seq(0.05, 0.95, length.out = 20)
lsig <- rep(5, 20)

fitM <- mqgam(form = form, data = UKload, err = 0.1, lsig = lsig, qu = qus)

# Predict each quantile curve and plot
pmat <- list()
for(ii in 1:length(qus)){
  pmat[[ii]] <- qdo(fitM, qus[ii], predict)[1:100]
  #lines(pmat[[ii]], col = 2)
}
pmat <- do.call("rbind", pmat)

matplot(t(pmat),
        type = 'l', lty = 1, 
        col = colorRampPalette(EdfColors[1:4])(20),
        ylab = 'Load', xlab = "Day",)
points(UKload$NetDemand[1:100])
```

Let's look at how the conditional distribution of the response looks like on different days. 
```R
pmat <- t(t(pmat) - pmat[10, ])
pmat <- apply(pmat, 2, sort)

subDat <- pmat[ , 51:100]
par(mfrow = c(1, 2))
for(ii in 1:50)
{
qpos <- subDat[ , ii]

plot(qpos, qus, type = 'b', ylim = c(0, 1), "xlab" = "(centered) LOAD", 
     xaxt="n", xlim = range(subDat), ylab = "F(LOAD)")
axis(1, round(seq(min(subDat), max(subDat), length.out = 5), 1))
abline(h = 1, lty = 2)
abline(h = 0, lty = 2)

n <- length(qpos)
tmp <- diff(qus) / diff(qpos)
plot( rep(qpos, each = 2), c(0, rep(tmp, each = 2), 0), "xlab" = "(centered) LOAD", ylab = "f(LOAD)", 
       xaxt="n", xlim = range(subDat), type = 'l')
axis(1, round(seq(min(subDat), max(subDat), length.out = 5), 1))
polygon(x = rep(qpos, each = 2),
        y = c(0, rep(tmp, each = 2), 0),
        col = "gray",border="black")
abline(h = 0)
Sys.sleep(0.5)
}
```

References
----------------------------
  
  * Fasiolo, M., Goude Y., Nedellec R. and Wood, S. N. (2017). Fast calibrated additive quantile regression. URL: https://arxiv.org/abs/1707.03307