#### Fit initial Gaussian models needed by QGAMs later on
.init_gauss_fit <- function(form, data, ctrl, argGam, qu, discrete){
  
  # We do not want initialisation or "sp" for qgam to be used to fit the gaussian gam
  argGam[c("coef", "start", "mustart", "etastart", "sp", "in.out")] <- NULL
  
  gam_name <- ifelse(discrete, "bam", "gam")
  
  if( is.formula(form) ) { # [A] Mean GAM 
    gFit <- do.call(gam_name, c(list("formula" = form, "data" = quote(data),
                                     "family" = gaussian(link=ctrl[["link"]]), "discrete" = discrete), argGam))
    varHat <- gFit$sig2
    formL <- form
  } else { # [B] Mean and variance GAM(s)
    if(discrete){ # B.a Discrete case: fit mean and variance separately
      gFit <- do.call(gam_name, c(list("formula" = form[[1]], "data" = quote(data),
                                          "family" = gaussian(link=ctrl[["link"]][[1]])), argGam, discrete = discrete))
      R <- residuals(gFit, type = "response")
      
      # Add variable to data with long name to avoid over-writing existing variables!
      data$my_abs_residuals_xyjz <- abs(R)
      
      RFit <- do.call(gam_name, c(list("formula" = as.formula(paste0(c("my_abs_residuals_xyjz", paste0(form[[2]])), collapse = " ")),
                                       "data" = quote(data)), argGam, discrete = discrete))
      
      varHat <- RFit$fitted.values^2 * (pi/2) # Convert MAD into variance
    } else { # B.b NOT Discrete case: fit mean and variance jointly with gaulss
      gFit <- do.call(gam_name, c(list("formula" = form, "data" = quote(data),
                                       "family" = gaulss(link=list(ctrl[["link"]], "logb"), b=ctrl[["b"]])), argGam))
      varHat <- 1/gFit$fit[ , 2]^2
    }
    formL <- form[[1]]
  }  
  
  # Provide initialisation for fitted quantiles
  initM <- lapply(qu, function(q){
    mustart <- qnorm(q, as.matrix(gFit$fitted.values)[, 1], sqrt(varHat))
    # Get initial regression coefficients, several cases to cover:
    # [A]: we do not model the variance and we have an intercept in the mean model
    # [B] We do model the variance and/or there is no intercept in the mean model
    # Under B we create synthetic responses, whose conditional mean is equal to the initial quantile mustart
    if(!is.list(form) & "(Intercept)" %in% names(coef(gFit))){ # [A]
      interc <- which(names(coef(gFit)) == "(Intercept)")
      coefstart <- coef(gFit)
      coefstart[interc] <- coefstart[interc] + mean(mustart - gFit$fitted.values)
    } else { # [B]
      data$my_new_mustart_for_coef_xyz <- mustart + residuals(gFit, type = "response") / sqrt(varHat)
      muFit <- do.call(gam_name, c(list("formula" = update(formL, my_new_mustart_for_coef_xyz ~ .), "data" = quote(data),
                                        "family" = gaussian(link=as.list(ctrl[["link"]])[[1]]), "discrete" = discrete), argGam))
      # Potentially pointless update to mustart, but this way should be coherent with coefstart
      mustart <- muFit$fitted.values
      coefstart <- coef(muFit)
    } 
    
    return( list("mustart" =  mustart, "coefstart" = coefstart,
                 "in.out" = NULL) ) # let gam() initialize sp via initial.spg())
  })
  
  return( list("gausFit" = gFit, "varHat" = varHat, "formL" = formL, "initM" = initM) )
}

