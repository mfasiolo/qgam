## ******************************************************************** Fonctions
## permettant d'apprendre les modèles et de prévoir
## ********************************************************************

# Les modèles considérés sont : * Modèles par instants : modInstant pour fitter
# et predictModInstant pour prévoir - lm - rf - gbm - gam * rqInstant : rq par
# instant *



.modInstant <- function(data0, data1, eq, model = "gam", distribution = "gaussian", 
                       nodesize = 100) {
  
  firstInstant <- data1$Instant[1]
  minInstant <- min(data1$Instant)
  maxInstant <- max(data1$Instant)
  if (firstInstant == minInstant) {
    perm <- minInstant:maxInstant - minInstant + 1
  } else {
    perm <- c(firstInstant:maxInstant, minInstant:(firstInstant - 1)) - minInstant + 
      1
  }
  
  forecast <- numeric(nrow(data1))
  data0 <- split(data0, f = data0$Instant)
  mod <- lapply(data0, function(x) {
    if (model == "gam") {
      gam(formula(eq), data = x)
    } else if (model == "lm") {
      lm(formula(eq), data = x)
    } else if (model == "rf") {
      m <- randomForest(formula(eq), data = x, na.action = na.exclude, nodesize = nodesize, 
                        ntree = 1000)
      m$fitted <- m$predicted
      m
    } else if (model == "gbm") {
      m <- list()
      m$model <- gbm(formula(eq), data = x, n.trees = 7000, interaction.depth = 2, 
                     n.minobsinnode = 5, shrinkage = 0.005, distribution = distribution, 
                     verbose = FALSE)
      m$ntrees <- gbm.perf(m$model, method = "OOB", plot.it = FALSE)
      m$fitted <- predict(m$model, x, n.trees = m$ntrees, type = "response")
      m
    }
  })
  
  fitted <- lapply(1:24, function(x) {
    mod[[x]]$fitted
  })
  
  data1 <- split(data1, f = data1$Instant)
  prev <- lapply(1:24, function(x) {
    if (model == "gbm") {
      predict(mod[[x]]$model, data1[[x]], n.trees = mod[[x]]$ntrees, type = "response")
    } else {
      predict(mod[[x]], newdata = data1[[x]])
    }
  })
  
  prev0 <- do.call("rbind", prev[perm])
  fitted0 <- do.call("rbind", fitted[perm])
  list(forecast = as.numeric(prev0), fitted = as.numeric(fitted0), model = mod)
}

### prevision des modèles par instant

.predictModInstant <- function(g, data1, model = "gam") {
  
  mod <- g$model
  forecast <- numeric(nrow(data1))
  
  data1List <- split(data1, f = data1$Instant)
  prev <- lapply(1:24, function(x) {
    if (model == "gbm") {
      predict(mod[[x]]$model, data1List[[x]], n.trees = mod[[x]]$ntrees, type = "response")
    } else {
      predict(mod[[x]], newdata = data1List[[x]])
    }
  })
  
  firstInstant <- data1$Instant[1]
  minInstant <- min(data1$Instant)
  maxInstant <- max(data1$Instant)
  if (firstInstant == minInstant) {
    perm <- minInstant:maxInstant - minInstant + 1
  } else {
    perm <- c(firstInstant:maxInstant, minInstant:(firstInstant - 1)) - minInstant + 
      1
  }
  prev0 <- do.call("rbind", prev[perm])
  as.numeric(prev0)
}

## ************************************************* RQ par Instant
## *************************************************

.rqInstant <- function(eq, data0, data1, probs) {
  forecast_q <- matrix(0, nrow = nrow(data1), ncol = length(probs))
  fitted_q <- matrix(0, nrow = nrow(data0), ncol = length(probs))
  mod <- list()
  for (i in 0:23) {
    sel0 <- which(data0$Instant == i)
    sel1 <- which(data1$Instant == i)
    
    mod_i <- rq(eq, data = Data0[sel0, ], tau = probs)
    
    forecast_q[sel1, ] <- as.matrix(predict(mod_i, newdata = Data1[sel1, ]))
    fitted_q[sel0, ] <- as.matrix(mod_i$fitted.values)
    mod[[i + 1]] <- mod_i
  }
  return(list(forecast_q = forecast_q, fitted_q = fitted_q, mod = mod))
}



## ********************************************** Gam quantile
## **********************************************

## Gam quantile de base eq : equation des effets sur la moyenne eq.var : equation
## des effets sur la variance
.gam_quant2 <- function(eq, data0, data1, probs, estim = TRUE, eq.var = NULL) {
  
  # Définition des objets
  y0 <- as.character(terms(as.formula(eq))[[2]])  # variable à prévoir
  
  # Modèle GAM
  m_gam <- gam(as.formula(eq), data = data0)
  fitted_gam <- m_gam$fitted.values
  forecast_gam <- predict.gam(m_gam, newdata = data1)
  
  # Base de données pour quantiles
  data0_q <- data.frame(y = data0[, y0], predict(m_gam, type = "terms"))  # effets estimés
  
  data1_q <- data.frame(y = data1[, y0], predict(m_gam, data1, type = "terms"))  # effets validation
  
  # Modèle GAM sur la variance
  m_gam_var <- NULL
  if (!is.null(eq.var)) {
    var.gam <- as.character(terms(formula(eq.var))[[2]])  # 'res2'
    data0[, var.gam] <- abs(m_gam$residuals)^2  # on rajoute res2
    m_gam_var <- gam(as.formula(eq.var), data = data0)
    data0_q <- cbind(data0_q, predict(m_gam_var, type = "terms"))
    data1_q <- cbind(data1_q, predict(m_gam_var, data1, type = "terms"))
  }
  names(data0_q) <- c("y", paste("X", 2:ncol(data0_q), sep = "."))
  names(data1_q) <- names(data0_q)
  
  m_rq <- vector("list", length(probs))  # vecteur contenant pour chaque quantile le modèle de regression quantile
  
  fitted_q <- data.frame(matrix(NA, ncol = length(probs), nrow = length(data0_q$y)))
  forecast_q <- data.frame(matrix(NA, ncol = length(probs), nrow = length(data1_q$y)))
  
  for (j in 1:length(probs)) {
    m_rq[[j]] <- rq(y ~ ., tau = probs[j], data = data0_q, na.action = na.exclude)
    fitted_q[, j] <- m_rq[[j]]$fitted.values
    # c'est juste du ménage...
    m_rq[[j]]$fitted.values <- NULL
    m_rq[[j]]$x <- NULL
    m_rq[[j]]$y <- NULL
    m_rq[[j]]$residuals <- NULL
    forecast_q[, j] <- predict.rq(m_rq[[j]], newdata = data1_q)  # prévision quantile
  }
  
  fitted_q <- data.frame(fitted_q)
  names(fitted_q) <- probs
  forecast_q <- data.frame(forecast_q)
  names(forecast_q) <- probs
  
  return(list(forecast = forecast_gam, fitted = fitted_gam, forecast_q = forecast_q, 
              fitted_q = fitted_q, mod = list(eq = eq, gam = m_gam, rq = m_rq, gam.var = m_gam_var, 
                                              eq.var = eq.var)))
}

## Fonction prévision de gam quantile
.predict.gam_quant2 <- function(mod, newdata, quantile = "all") {
  nprobs <- ifelse(quantile == "all", 99, length(quantile))
  y0 <- as.character(terms(as.formula(mod[[1]]$eq))[[2]])
  q <- matrix(0, ncol = nprobs, nrow = nrow(newdata))
  for (i in 0:23) {
    sel_i <- which(newdata$Heure == i)
    
    data1_q_i <- data.frame(y = newdata[sel_i, y0], predict(mod[[i + 1]]$gam, 
                                                            newdata[sel_i, ], type = "terms"))
    if (!is.null(mod[[i + 1]]$eq.var)) {
      data1_q_i <- cbind(data1_q_i, predict(mod[[i + 1]]$gam.var, newdata[sel_i, 
                                                                          ], type = "terms"))
    }
    
    names(data1_q_i) <- c("y", paste("X", 2:ncol(data1_q_i), sep = "."))
    
    if (quantile == "all") {
      quantiles <- lapply(mod[[i + 1]]$rq, function(x) predict.rq(x, newdata = data1_q_i))
      q[sel_i, ] <- do.call("cbind", quantiles)
    } else {
      q[sel_i, ] <- predict.rq(mod$rq[[quantile]], newdata = data1_q_i)
    }
  }
  return(q)
}

## Fonction de prévision de gam quantile utilisée pour la création de scénarios
## météo Comprendre la différence avec predict.gam_quant2 !!!

.predict.gam_quant <- function(mod, newdata, quantile = "all") {
  y0 <- as.character(terms(as.formula(mod$eq))[[2]])
  data1_q <- data.frame(y = newdata[, y0], predict(mod$gam, newdata, type = "terms"))
  
  if (!is.null(mod$eq.var)) {
    data1_q <- cbind(data1_q, predict(mod$gam.var, newdata, type = "terms"))
  }
  
  names(data1_q) <- c("y", paste("X", 2:ncol(data1_q), sep = "."))
  
  if (quantile == "all") {
    quantiles <- lapply(mod$rq, function(x) predict.rq(x, newdata = data1_q))
    q <- do.call("cbind", quantiles)
  } else {
    q <- predict.rq(mod$rq[[quantile]], newdata = data1_q)
  }
  return(q)
}



## Gam quantile horaires (h = 0,...,23)
.gam_quant_instant <- function(eq, data0, data1, probs, eq.var = NULL) {
  forecast <- numeric(nrow(data1))
  fitted <- numeric(nrow(data0))
  forecast_q <- matrix(0, nrow = nrow(data1), ncol = length(probs))
  fitted_q <- matrix(0, nrow = nrow(data0), ncol = length(probs))
  mod <- list()
  for (i in 0:23) {
    cat(i, " ")
    sel0 <- which(data0$Instant == i)
    sel1 <- which(data1$Instant == i)
    
    mod_i <- gam_quant2(eq, data0[sel0, ], data1[sel1, ], probs, eq.var = eq.var)  ## test !! 
    forecast_q[sel1, ] <- as.matrix(mod_i$forecast_q)
    fitted_q[sel0, ] <- as.matrix(mod_i$fitted_q)
    forecast[sel1] <- mod_i$forecast
    fitted[sel0] <- mod_i$fitted
    mod[[i + 1]] <- mod_i$mod
  }
  return(list(forecast = forecast, fitted = fitted, forecast_q = forecast_q, fitted_q = fitted_q, 
              mod = mod))
}

## Version parallèle de gam_quant_instant
.gam_quant_instant.parallel <- function(eq, data0, data1, probs, eq.var = NULL) {
  forecast <- numeric(nrow(data1))  # résultats de prévision de la moyenne
  fitted <- numeric(nrow(data0))  # résultats d'estimation de la moyenne
  forecast_q <- matrix(0, nrow = nrow(data1), ncol = length(probs))  # résultats de prévision des quantiles (probs = 0.01,...)
  fitted_q <- matrix(0, nrow = nrow(data0), ncol = length(probs))  #... estimation ...
  mod.h <- list(NULL)  # variable temporaire pour mettre les modeles horaires
  H <- as.numeric(levels(as.factor(data0$Instant)))
  
  # modèles horaires (parallélisation)
  models <- foreach(h = H, .combine = c) %dopar% {
    sel0 <- which(data0$Instant == h)
    sel1 <- which(data1$Instant == h)
    
    mod_i <- gam_quant2(eq, data0[sel0, ], data1[sel1, ], probs, eq.var = eq.var)
    mod.h[[1]]$forecast_q <- as.matrix(mod_i$forecast_q)
    mod.h[[1]]$fitted_q <- as.matrix(mod_i$fitted_q)
    mod.h[[1]]$forecast <- mod_i$forecast
    mod.h[[1]]$fitted <- mod_i$fitted
    mod.h[[1]]$mod <- mod_i$mod
    mod.h
  }
  
  # on réagrège les données
  mod <- list()
  for (h in H) {
    sel0 <- which(data0$Instant == h)
    sel1 <- which(data1$Instant == h)
    forecast_q[sel1, ] <- models[[h + 1 - H[[1]]]]$forecast_q
    fitted_q[sel0, ] <- models[[h + 1 - H[[1]]]]$fitted_q
    forecast[sel1] <- models[[h + 1 - H[[1]]]]$forecast
    fitted[sel0] <- models[[h + 1 - H[[1]]]]$fitted
    mod[[h + 1 - H[[1]]]] <- models[[h + 1 - H[[1]]]]$mod
  }
  return(list(forecast = forecast, fitted = fitted, forecast_q = forecast_q, fitted_q = fitted_q, 
              mod = mod))
}



## ***************************************************** Fonctions GLM quantile
## *****************************************************


.glm_quant <- function(prob, X, Y, Y.fit, Y.forecast, Xnew, alpha = 1, nfold = 10, 
                      h) {
  Y.quant <- Y - Y.fit
  forecast_q <- matrix(0, ncol = length(prob), nrow = nrow(Xnew))
  for (i in c(1:length(prob))) {
    Kern <- exp(-(Y.quant - quantile(Y.quant, prob[i]))^2/(2 * h^2))
    Kern <- Kern/sum(Kern)
    cvfit.quant <- cv.glmnet(X, Y.quant, type.measure = "mse", nfolds = nfold, 
                             alpha = alpha, weights = Kern)
    YchapCV.quant <- predict(cvfit.quant, newx = Xnew)
    forecast_q[, i] <- exp(YchapCV.quant + Y.forecast)
    print(i)
  }
  return(forecast_q)
}


.glm_quant.par <- function(prob, X, Y, Y.fit, Y.forecast, Xnew, alpha = 1, nfold = 10, 
                          h) {
  Y.quant <- Y - Y.fit
  forecast_q <- foreach(i = c(1:length(prob)), .combine = cbind) %dopar% {
    Kern <- exp(-(Y.quant - quantile(Y.quant, prob[i]))^2/(2 * h^2))
    Kern <- Kern/sum(Kern)
    cvfit.quant <- cv.glmnet(X, Y.quant, type.measure = "mse", nfolds = nfold, 
                             alpha = alpha, weights = Kern, parallel = T)
    YchapCV.quant <- predict(cvfit.quant, newx = Xnew)
    exp(YchapCV.quant + Y.forecast)
  }
  return(forecast_q)
}



## ***************************************************** Autres modèles : à
## vérifier si ils sont utilisés
## ******************************************************

## Une version parallèle optimisée de Gam par instant
.gamH.parallel <- function(Data0, Data1, eq) {
  forecast <- array(0, dim = nrow(Data1))
  fitted <- array(0, dim = nrow(Data0))
  mod <- list()
  H <- as.numeric(levels(as.factor(Data0$Instant)))
  
  # On parallélise le calcul des modèles gam par tranche demi-horaire
  mod.h <- list()
  mod <- foreach(h = H, .combine = c) %dopar% {
    sel0 <- which(Data0$Instant == h)
    sel1 <- which(Data1$Instant == h)
    data_0h <- Data0[sel0, ]
    data_1h <- Data1[sel1, ]
    eq1 <- eval(parse(text = eq))
    mod.h[[1]] <- gam(eq1, data = data_0h, na.action = na.exclude)
    mod.h
  }
  # Prévisions et fittedimations de nos modèles
  for (h in H) {
    sel0 <- which(Data0$Instant == h)
    sel1 <- which(Data1$Instant == h)
    Data1_h <- Data1[sel1, ]
    fitted[sel0] <- mod[[h + 1 - H[[1]]]]$fitted.values
    forecast[sel1] <- predict(mod[[h + 1 - H[[1]]]], Data1_h)
  }
  l <- list(forecast = forecast, fitted = fitted, model = mod)
  return(l)
}


## une version horaire de quantregForest
.quantregForest_instant <- function(y, x0, x1, instant0, instant1, probs, ntree = 200, 
                                   nodesize = nodesize) {
  forecast <- numeric(nrow(x1))
  fitted <- numeric(nrow(x0))
  forecast_q <- matrix(0, nrow = nrow(x1), ncol = length(probs))
  fitted_q <- matrix(0, nrow = nrow(x0), ncol = length(probs))
  for (i in 0:23) {
    # cat(i, ' ')
    sel0 <- which(instant0 == i)
    sel1 <- which(instant1 == i)
    
    mod_i <- quantregForest(x0[sel0, ], y[sel0], ntree = ntree, nodesize = nodesize)
    forecast_q[sel1, ] <- as.matrix(predict(mod_i, newdata = x1[sel1, ], quantiles = probs))
    fitted_q[sel0, ] <- as.matrix(predict(mod_i, quantiles = probs))
  }
  return(list(forecast_q = forecast_q, fitted_q = fitted_q))
} 
