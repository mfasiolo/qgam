
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/qgam)](https://cran.r-project.org/package=qgam)
[![Build Status](https://travis-ci.org/mfasiolo/qgam.svg?branch=master)](https://travis-ci.org/mfasiolo/qgam)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mfasiolo/qgam?branch=master&svg=true)](https://ci.appveyor.com/project/mfasiolo/qgam)

# **qgam**: quantile GAMs

This R package offers methods for fitting additive quantile regression models based on splines, using the methods described in [Fasiolo et al., 2017](https://arxiv.org/abs/1707.03307).

See the [vignette](https://mfasiolo.github.io/qgam/articles/qgam.html) for an introduction to the most important functions:

- `qgam` fits an additive quantile regression model to a single quantile. Very similar to `mgcv::gam`. It returns an `mgcv::gamObject`.
- `mqgam` fits the same additive quantile regression model to several quantiles. It is more efficient that calling `qgam` several time, 
          especially in terms of memory.
- `tuneLearn` useful for tuning the learning rate of the Gibbs posterior. It evaluates a calibration loss function on a grid of values 
              provided by the user. 
- `tuneLearnFast` similar to `tuneLearn`, but here the learning rate is selected by minimizing the calibration loss using Brent method.

References
----------------------------
  
  * Fasiolo, M., Goude Y., Nedellec R. and Wood, S. N. (2017). Fast calibrated additive quantile regression. URL: https://arxiv.org/abs/1707.03307