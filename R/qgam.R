#######
### Fitting quantile gam model
#######
#' Fitting quantile gam model
#' 
#' @param \code{XXX} .
#' @return XXX.
#'
#' @details XXX
#'         
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.   
#' @references Fasiolo and Wood XXX.           
#' @export
#' 

qgam <- function(form, data, tau, lsigma = NULL, err = 0.01, ncores = 1, control = list(), controlGam = list())
{
  nt <- length(tau)
  
  # Initial Gaussian fit
  gausFit <- gam(form, data = data, control = controlGam)
  control[["gausFit"]] <- gausFit
  
  # Output list
  out <- list()
  
  if( is.null(lsigma) ) { # Selecting the learning rate sigma OR ....
    learn <- tuneLearnFast(form = form, data = data, err = err, tau = tau, ncores = ncores, control = control)
    lsigma <- learn$lsigma
    out[["calibr"]] <- learn
  } else { # ... use the one provided by the user
    if( length(lsigma) == 1 ) {
      lsigma <- rep(lsigma, nt)
    } else {
      if( length(lsigma) != nt ) stop("lsigma should either be scalar of a vector of length(tau) ")
    } }
  

  # Selection lambda
  lam <- err * sqrt(2*pi*gausFit$sig2) / (2*log(2)*exp(lsigma))
  
  # Fitting a quantile model for each tau
  out[["fit"]] <- lapply(1:nt, function(ii){
    
    .out <- gam(form, family = logF(tau = tau[ii], lam = lam[ii], theta = lsigma[ii]), data = data)
    
    # Removing data and smooth matrix to reduce memory requirements. There quantities
    # are kept only inside the first fit ( qfit[[1]] )
    if(ii > 1){
      .out$model  <- NULL
      .out$smooth <- NULL 
    } 
        
    return( .out )
  })
  
  # Storing output list
  names(out[["fit"]]) <- tau
  out[["model"]] <- out[["fit"]][[1]][["model"]]
  out[["smooth"]] <- out[["fit"]][[1]][["smooth"]]
  out[["fit"]][[1]][["model"]] <- NULL
  out[["fit"]][[1]][["smooth"]] <- NULL
  
  out[["tau"]] <- tau
  out[["lambda"]] <- lam
  out[["lsigma"]] <- lsigma
  
  return( out )
}

# fit$model
# fit$smooth[[1]]$UZ   # fit$smooth

# 
# qFit <- lapply(1:nt, function(ii){
#   
#   lam <- err * sqrt(2*pi*gfit$sig2) / (2*log(2)*exp(sigs[ii])) 
#   
#   out <-  gam(form, family = logF(tau = tauSeq[ii], lam = lam, theta = sigs[ii]), data = data)
#   
#   return(out)
#   
# })
# 
# tuneLearnFast(form = y~s(x, k = 20), data = dataf, err = 0.02, tau = seq(0.1, 0.9, length.out = 5), ncores = 4, 
#               control = list("K" = 20))
# 
# 
# 
# 
# 
# 
# tic <- proc.time()
# set.seed(41241)
# nrep <- 20
# tmp <- lapply(1:nrep, function(nouse) sample(1:n, n, replace = TRUE))
# boot <- lapply(tmp, function(ff) dataf[ff, ] )
# tuneLearnFast(form = y~s(x), data = dataf, boot = boot, err = 0.01,
#               srange = c(-2, 1), tau = tau, ncores = 4)
# proc.time() - tic
# 
# tic <- proc.time()
# set.seed(41241)
# sigSeq <- seq(-2, 1, length.out = 16)
# closs <- tuneLearn(form = y~s(x), data = dataf, nrep = 20, err = 0.02,
#                    lsig = sigSeq, tau = tau, ncores = 4)
# proc.time() - tic
# 
# plot(sigSeq, closs$loss, type = 'b', lwd = 2, ylab = "Empirical Coverage", 
#      xlab = expression(sigma))
# 
# abline(v = fitREML$family$getTheta(F), col = 2, lwd = 2)
# abline(v = sigSeq[ which.min(closs$loss) ], lwd = 2)
# 
# 
# tauSeq <- c(0.05, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95)
# nt <- length(tauSeq)
# 
# nrep = 20
# tmp <- lapply(1:nrep, function(nouse) sample(1:n, n, replace = TRUE))
# boot <- lapply(tmp, function(ff) dataf[ff, ] )
# 
# tic <- proc.time()
# res <- list()
# sigs <- numeric(nt)
# for(ii in 1:nt)
# {
#   lsig <- tuneLearnFast(form = y~s(x), data = dataf, boot = boot, err = 0.02,
#                         srange = c(-2, 1), tau = tauSeq[[ii]], ncores = 4)$minimum
#   
#   res[[ii]] <-  gam(y~s(x, k = 30), 
#                     family = logF(tau = tauSeq[[ii]], lam = 0.1, theta = lsig), data = dataf)
#   sigs[ii] <- lsig
# }
# proc.time() - tic
# 
# plot(x, dat, main = "Coverage matching", pch = '.', ylab = "y")
# 
# for(ii in 1:nt){
#   truth <- f + qnorm(1-tauSeq[ii], 0, sigma)
#   fCV <- predict(res[[ii]], data.frame(x=x), se=TRUE)
#   
#   lines(x, truth, col = 3, lwd = 2)
#   lines(x, fCV$fit, lwd = 2)
# }
# 
# 
# #############################################################################
# tauSeq <- seq(0.1, 0.9, length.out = 5)
# 
# form <- y~s(x, k = 30) 
# data <- dataf
# 
# err <- 0.01
# 
# # Setting up control parameter
# ctrl <- list( "K" = 20,
#               "gauFit" = NULL,
#               "brac" = log( c(0.5, 2) ), 
#               "tol" = 1e-2,
#               "verbose" = TRUE )
# 
# ###### Function starts here
# 
# # (Optional) create K boostrap datasetd
# n <- nrow(data)
# tmp <- lapply(1:ctrl$K, function(nouse) sample(1:n, n, replace = TRUE))
# boot <- lapply(tmp, function(ff) data[ff, ] )
# 
# ###### Function starts here
# if( tol > 0.1 * abs(diff(brac)) ) stop("tol > bracket_widhts/10, choose smaller 
#                                        tolerance or larger bracket")
# 
# # Order quantiles so that those close to the median are dealt with first
# oTau <- order( abs(tauSeq-0.5) )
# 
# # Initial gaussian fit to get the scale
# gfit <- gam(form, data = data)
# 
# # We assume lam~0 and we match the variance of a symmetric Laplace density with that of the Gaussian fit.
# # We use the value of tau that is the closest to 0.5
# tmp <- tauSeq[ oTau[1] ]
# iSig <- log(sqrt(  gfit$sig2 * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
# 
# # Initialize
# nt <- length(tauSeq)
# srange <- iSig + brac
# 
# sigs <- efacts <- errors <- numeric(nt)
# rans <- matrix(NA, nt, 2)
# 
# # Here we need bTol > aTol otherwise we the new bracket will be too close to probable solution
# aTol <- 1.5 * tol
# bTol <- 2 * tol
# 
# for(ii in 1:nt)
# {
#   oi <- oTau[ii]
#   
#   ef <- 1
#   
#   repeat{
#     
#     # Compute bracket
#     srange <- iSig + ef * brac
#     
#     # Estimate log(sigma)
#     res  <- tuneLearnFast(form = form, data = data, tau = tauSeq[oi], err = err, boot = boot, 
#                           srange = srange, ncores = ncores, control = ctrl)
#     lsig <- res$minimum
#         
#     # If solution not too close to boundary store results and determine bracket for next iter
#     if( all(abs(lsig-srange) > aTol) ){ 
#       
#       sigs[oi] <- lsig
#       rans[oi, ] <- srange
#       efacts[oi] <- ef
#       errors[oi] <- res$err
#       
#       if(ii < nt)
#       {
#         kk <- oTau[ which.min(abs(tauSeq[oTau[ii+1]] - tauSeq[oTau[1:ii]])) ]
#         iSig <- sigs[kk] 
#         wd <- abs(diff(rans[kk, ]))
#         brac <- c(-1, 1) * wd / 2
#         
#         # If kk solution close to center of kk bracket, halve the bracket size 
#         # (unless the size of the bracket is < 10*tol)
#         if( (abs(iSig - mean(rans[kk, ])) < 0.25*wd) && (wd > 10*tol) ) brac <- brac / 2
#       }
#       
#       break
#     }
#     
#     # If solution is close to bracket boundaries, we shift bracket and expand it
#     wd <- abs( diff(brac) )
#     iSig <- lsig + ifelse(lsig-srange[1] < aTol, - wd + bTol, wd - bTol)
#     ef <- 2*ef
#   }
#   
#   if( verbose )
#   {
#     tseq <- oTau[1:ii]
#     tmp <- rans[tseq, , drop = FALSE]
#     layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
#            heights=c(2, 1))
#     par(mai = c(1, 1, 0.1, 0.1))
#     plot(tauSeq[tseq], sigs[oTau[1:ii]], ylim = range(as.vector(tmp)), xlim = range(tauSeq), col = 2, 
#          ylab = expression("Log(" * sigma * ")"), xlab = "tau")
#     points(tauSeq[tseq], tmp[ , 1], pch = 3)
#     points(tauSeq[tseq], tmp[ , 2], pch = 3)
#     points(tauSeq[tseq], rowMeans(tmp), pch = 3)
#     for(zz in 1:ii) segments(tauSeq[oTau[zz]], rowMeans(rans)[oTau[zz]] - abs(diff(tmp[zz, ]))/4, 
#                              tauSeq[oTau[zz]], rowMeans(rans)[oTau[zz]] + abs(diff(tmp[zz, ]))/4, col = 1)
#     plot(tauSeq, efacts, xlab = "tau", "ylab" = "Bracket expansions")  
#     plot(tauSeq, errors)
#   }
# }
# 
# list(  )
# 
# 
# 
# 
# plot(x, dat, main = "Coverage matching", pch = '.', ylab = "y")
# 
# for(ii in 1:nt){
#   truth <- f + qnorm(1-tauSeq[ii], 0, sigma)
#   fCV <- predict(qFit[[ii]], data.frame(x=x), se=TRUE)
#   
#   # lines(x, truth, col = 3, lwd = 2)
#   lines(x, fCV$fit, lwd = 2)
# }
# 
# ii = 1
# plot(x, dat, main = "Coverage matching", pch = '.', ylab = "y")
# fCV <- predict(qFit[[ii]], data.frame(x=x), se=TRUE)
# lines(x, fCV$fit, lwd = 2)
# lines(x, fCV$fit + 2*fCV$se.fit, lwd = 2, col = 2)
# lines(x, fCV$fit - 2*fCV$se.fit, lwd = 2, col = 2)
# ii <- ii + 1
# 
# 
# 
# 
# # iSig / ( lam^2 * ( trigamma(lam*tau) + trigamma(lam*(1-tau))) )
# # 
# # mu <- 0
# # tau <- 0.5
# # sig <- 0.36
# # lam <- 0.1
# # s <- rlf(1e6, mu, tau, sig, lam)
# # 
# # mean(s)
# # sig * lam * ( digamma(lam*tau) - digamma(lam*(1-tau)) ) + mu
# # 
# # var(s)
# # sig^2 * lam^2 * ( trigamma(lam*tau) + trigamma(lam*(1-tau)) )
# # 
# # 
# pr <- 0.5
# sig <- 0.36
# a <- rALD(1e4, mu = 0, sigma = sig, p = pr)
# var(a)
# varALD(mu = 0, sigma = sig, p = pr)
# sig^2 * (2*pr^2 - 2*pr + 1) / (pr^2*(1-pr)^2)
# 
# 
