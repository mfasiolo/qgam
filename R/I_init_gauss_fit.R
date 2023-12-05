#### Fit initial Gaussian models needed by QGAMs later on
.init_gauss_fit <- function(form, data, ctrl, argGam, qu, discrete){
  
  gam_name <- ifelse(discrete, "bam", "gam")
  
  if( is.formula(form) ) { # [A] Mean GAM 
    gausFit <- if( is.null(ctrl[["gausFit"]]) ) {
      do.call(gam_name, c(list("formula" = form, "data" = quote(data),
                               "family" = gaussian(link=ctrl[["link"]]), "discrete" = discrete), argGam))
    } else {
      ctrl$gausFit
    }
    varHat <- gausFit$sig2
    formL <- form
  } else { # [B] Mean and variance GAM(s)
    if(is.null(ctrl[["gausFit"]])){ # B.1) Gaussian fit NOT available
      if(discrete){ # B.1.a Discrete case: fit mean and variance separately
        gausFit <- do.call(gam_name, c(list("formula" = form[[1]], "data" = quote(data),
                                            "family" = gaussian(link=ctrl[["link"]][[1]])), argGam, discrete = discrete))
        R <- residuals(gausFit, type = "response")
        
        # Add variable to data with long name to avoid over-writing existing variables!
        data$my_abs_residuals_xyjz <- abs(R)
        
        RFit <- do.call(gam_name, c(list("formula" = as.formula(paste0(c("my_abs_residuals_xyjz", paste0(form[[2]])), collapse = " ")),
                                         "data" = quote(data)), argGam, discrete = discrete))
        
        gausFit$my_var_fit_xyz <- RFit
        
        varHat <- RFit$fitted.values^2 * (pi/2)
      } else { # B.1.b NOT Discrete case: fit mean and variance jointly with gaulss
        gausFit <- do.call(gam_name, c(list("formula" = form, "data" = quote(data),
                                            "family" = gaulss(link=list(ctrl[["link"]], "logb"), b=ctrl[["b"]])), argGam))
        varHat <- 1/gausFit$fit[ , 2]^2
      }
    } else { # B.2 Gaussian fit available
      gausFit <- ctrl$gausFit
      if(discrete){
        varHat <- gausFit$my_var_fit_xyz$fitted.values^2 * (pi/2)
      }else{
        varHat <- 1/gausFit$fit[ , 2]^2
      }
    }
    formL <- form[[1]]
  }
  
  # Provide initialisation for fitted quantiles
  initM <- lapply(qu, function(q){
    list("mustart" =  qnorm(q, as.matrix(gausFit$fitted.values)[, 1], sqrt(varHat)), 
         "in.out" = NULL) # let gam() initialize sp via initial.spg()
  })
  
  return( list("gausFit" = gausFit, "varHat" = varHat, "formL" = formL, "initM" = initM) )
}

