## (c) Simon N. Wood & Matteo Fasiolo
## 2013-2015. Released under GPL2.

## extended families for mgcv, standard components. 
## family - name of family character string
## link - name of link character string
## linkfun - the link function
## linkinv - the inverse link function
## mu.eta - d mu/d eta function (derivative of inverse link wrt eta)
## note: for standard links this information is supplemented using 
##       function fix.family.link.extended.family with functions 
##       gkg where k is 2,3 or 4 giving the kth derivative of the 
##       link over the first derivative of the link to the power k.
##       for non standard links these functions muct be supplied.
## dev.resids - function computing deviance residuals.
## Dd - function returning derivatives of deviance residuals w.r.t. mu and theta. 
## aic - function computing twice - log likelihood for 2df to be added to.
## initialize - expression to be evaluated in gam.fit4 and initial.spg 
##              to initialize mu or eta.
## preinitialize - optional expression evaluated in estimate.gam to 
##                 e.g. initialize theta parameters (see e.g. ocat)
## postproc - expression to evaluate in estimate.gam after fitting (see e.g. betar)
## ls - function to evaluated log saturated likelihood and derivatives w.r.t.
##      phi and theta for use in RE/ML optimization. If deviance used is just -2 log 
##      lik. can njust return zeroes. 
## validmu, valideta - functions used to test whether mu/eta are valid.      
## n.theta - number of theta parameters.
## no.r.sq - optional TRUE/FALSE indicating whether r^2 can be computed for family
## ini.theta - function for initializing theta.
## putTheta, getTheta - functions for storing and retriving theta values in function 
##                      environment.
## rd - optional function for simulating response data from fitted model.
## residuals - optional function for computing residuals.
## predict - optional function for predicting from model, called by predict.gam.
## family$data - optional list storing any family specific data for use, e.g. in predict
##               function.

#######################
## Sinh-cosh density
#######################

.shashlss <- function (link = list("identity", "identity", "identity", "identity")) 
{ 
  ## Extended family object for modified log-F, to allow direct estimation of theta
  ## as part of REML optimization. Currently the template for extended family objects.
  ## length(theta)=2; log theta supplied. 
  ## Written by Matteo Fasiolo.
  npar <- 4
  if (length(link) != npar) stop("logFlss requires 4 links specified as character strings")
  okLinks <- lapply(1:npar, function(inp) c("identity"))
  stats <- list()
  param.names <- c("mu", "tau", "eps", "phi")
  for (i in 1:npar) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
      stop(link[[i]]," link not available for ", param.names[i]," parameter of shashlss")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
                          mu.eta=stats[[i]]$mu.eta),
                     class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  } 
  
  # variance <- function(mu) exp(get(".Theta"))  ##### XXX ##### Necessary?
  
  # validmu <- function(mu) all( is.finite(mu) )
  
  residuals <- function(object, type = c("deviance", "response")) {
    
    mu <- object$fitted[ , 1]
    tau <- object$fitted[ , 2]
    eps <- object$fitted[ , 3]
    phi <- object$fitted[ , 4]
    
    sig <- exp( tau )
    del <- exp( phi )
    
    type <- match.arg(type)
    
    # raw residuals  
    rsd <- object$y - mu - sig*del*exp(0.25)*(besselK(0.25, nu = (1/del+1)/2)+besselK(0.25, nu = (1/del-1)/2))/sqrt(8*pi) ####### XXX #######
      
    if (type=="response"){ 
      return(rsd)
    }
    else { ## compute deviance residuals
      sgn <- sign(rsd)
      
      z <- (object$y - mu) / (sig*del)
      
      dTasMe <- del*asinh(z) - eps
      CC <- cosh( dTasMe )
      SS <- sinh( dTasMe )
            
      l <- -tau - 0.5*log(2*pi) + log(CC) - 0.5*log1p(z^2) - 0.5*SS^2

#       stop("Matteo: I can't calculate the saturated likelihood")
      
      ls <- 0 ############~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!~#################
      
      rsd <- pmax(0, 2*(ls - l))
      
      rsd <- sqrt(rsd)*sgn
    }
    rsd
  } ## residuals
  
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv=0, d1b=0, d2b=0, Hp=NULL, rank=0, fh=NULL, D=NULL) {
    ## function defining the shash model log lik. 
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    ##        4 - everything.
    npar <- 4
    
    # If offset is not null or a vector of zeros, give an error
    if( !is.null(offset[[1]]) && sum(abs(offset)) )  stop("offset not still available for this family")
    
    jj <- attr(X,"lpi") ## extract linear predictor index
    
    eta <-  X[,jj[[1]],drop=FALSE] %*% coef[jj[[1]]]
    eta1 <- X[,jj[[2]],drop=FALSE] %*% coef[jj[[2]]]
    eta2 <- X[,jj[[3]],drop=FALSE] %*% coef[jj[[3]]]
    eta3 <- X[,jj[[4]],drop=FALSE] %*% coef[jj[[4]]]
    
    mu  <-  family$linfo[[1]]$linkinv(eta)
    tau <-  family$linfo[[2]]$linkinv(eta1)
    eps <-  family$linfo[[3]]$linkinv(eta2)
    phi <-  family$linfo[[4]]$linkinv(eta3)
    
    sig <- exp( tau )
    del <- exp( phi )
    
    n <- length(y)
    l1 <- matrix(0, n, npar)
    
    z <- (y - mu) / (sig*del)
    
    dTasMe <- del*asinh(z) - eps
    g <- -dTasMe
    CC <- cosh( dTasMe )
    SS <- sinh( dTasMe )
  
    l <- sum( -tau - 0.5*log(2*pi) + log(CC) - 0.5*log1p(z^2) - 0.5*SS^2 )
    
    if (deriv>0) {
      
      zsd <- z*sig*del
      sSp1 <- sqrt(z^2+1)
      asinhZ <- asinh(z)
      
      ## First derivatives
      De <- tanh(g) - 0.5*sinh(2*g)
      Dm <- 1/(del*sig*sSp1)*(del*(De)+z/sSp1)
      Dt <- zsd*Dm - 1
      Dp <- Dt + 1 - del*asinhZ*De
            
      l1[ , 1] <- Dm
      l1[ , 2] <- Dt
      l1[ , 3] <- De
      l1[ , 4] <- Dp
      
      ## the second derivatives  
      Dme <- (sech(g)^2 - cosh(2*g)) / (sig*sSp1)
      Dte <- zsd*Dme
      Dmm <- Dme/(sig*sSp1) + z*De/(sig^2*del*(z^2+1)^(3/2)) + (z^2-1)/(del*sig*del*sig*(z^2+1)^2)
      Dmt <- zsd*Dmm - Dm
      Dee <- -2*cosh(g)^2 + sech(g)^2 + 1 
      Dtt <-  zsd*Dmt
      Dep <- Dte - del*asinhZ*Dee
      Dmp <- Dmt + De/(sig*sSp1) - del*asinhZ*Dme
      Dtp <- zsd*Dmp
      Dpp <- Dtp - del*asinhZ*Dep + del*(z/sSp1-asinhZ)*De
      
      # Put them in matrix form
      l2 <- matrix(0, n, npar*(npar+1)/2)
      
      l2[ , 1] <- Dmm
      l2[ , 2] <- Dmt
      l2[ , 3] <- Dme
      l2[ , 4] <- Dmp
      l2[ , 5] <- Dtt
      l2[ , 6] <- Dte  
      l2[ , 7] <- Dtp 
      l2[ , 8] <- Dee 
      l2[ , 9] <- Dep 
      l2[ , 10] <- Dpp  
        
      ## need some link derivatives for derivative transform
      #ig1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(eta1))
      #g2 <- cbind(family$linfo[[1]]$d2link(mu), family$linfo[[2]]$d2link(sig))
    }
    
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults
    
    if (deriv>1) {
      
      ## the third derivatives
      Deee <-  -2*(sinh(2*g)+sech(g)^2*tanh(g))
      Dmee <- Deee/(sig*sSp1)
      Dmme <- Dmee/(sig*sSp1) + z*Dee/(sig*sig*del*(z^2+1)^(3/2))
      Dmmm <- 2*z*Dme/(sig*sig*del*sSp1^3) + Dmme/(sig*sSp1) + 
              (2*z^2-1)*De/(sig^3*del^2*sSp1^5) + 2*z*(z^2-3)/((sig*del)^3*(z^2+1)^3)
      Dmmt <- zsd*Dmmm - 2*Dmm
      Dtee <- zsd*Dmee
      Dmte <- zsd*Dmme - Dme
      Dtte <- zsd*Dmte
      Dmtt <- zsd*Dmmt - Dmt
      Dttt <- zsd*Dmtt
      Dmep <- Dmte + Dee/(sig*sSp1) - del*asinhZ*Dmee
      Dtep <- zsd*Dmep
      Deep <- Dtee - del*asinhZ*Deee
      Depp <- Dtep - del*asinhZ*Deep + del*( z/sSp1-asinhZ )*Dee
      Dmmp <- Dmmt + 2*Dme/(sig*sSp1) + z*De/(del*sig*sig*sSp1^3) - del*asinhZ*Dmme
      Dmtp <- zsd*Dmmp - Dmp
      Dttp <- zsd*Dmtp
      Dmpp <- Dmtp + Dep/(sig*sSp1) + z^2*De/(sig*sSp1^3) - 
              del*asinhZ*Dmep + del*Dme*(z/sSp1 - asinhZ)
      Dtpp <- zsd*Dmpp
      Dppp <- Dtpp - del*asinhZ*Depp + del*(z/sSp1-asinhZ)*(2*Dep + De) + del*(z/sSp1)^3 * De
    
      ## Put them in matrix form
      l3 <- matrix(0, n, (npar*(npar+3)+2)*npar/6) 
            
      l3[ , 1] <- Dmmm
      l3[ , 2] <- Dmmt
      l3[ , 3] <- Dmme
      l3[ , 4] <- Dmmp
      l3[ , 5] <- Dmtt
      l3[ , 6] <- Dmte
      l3[ , 7] <- Dmtp
      l3[ , 8] <- Dmee
      l3[ , 9] <- Dmep
      l3[ , 10] <- Dmpp
      l3[ , 11] <- Dttt
      l3[ , 12] <- Dtte
      l3[ , 13] <- Dttp
      l3[ , 14] <- Dtee
      l3[ , 15] <- Dtep
      l3[ , 16] <- Dtpp
      l3[ , 17] <- Deee
      l3[ , 18] <- Deep
      l3[ , 19] <- Depp
      l3[ , 20] <- Dppp
      
      #g3 <- cbind(family$linfo[[1]]$d3link(mu), family$linfo[[2]]$d3link(sig))
    }
    
    if (deriv>3) {
      ## the fourth derivatives
      
      #l4 <- matrix(0, n, 5) 
      
      stop("Matteo: don't have 4th order derivatives")
      
      #g4 <- cbind(family$linfo[[1]]$d4link(mu), family$linfo[[2]]$d4link(sig))
    }
        
    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4
      
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      #de <- mgcv:::gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)
      
      ## get the gradient and Hessian...
      ret <- mgcv:::gamlss.gH(X,jj,l1,l2,i2,l3=l3,i3=i3,l4=l4,i4=i4,
                              d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
      
    } else ret <- list()
    ret$l <- l; ret
  } ## end ll
  
  initialize <- expression({
    ## idea is to regress g(y) on model matrix for mean, and then 
    ## to regress the corresponding log absolute residuals on 
    ## the model matrix for log(sigma) - may be called in both
    ## gam.fit5 and initial.spg... note that appropriate E scaling
    ## for full calculation may be inappropriate for initialization 
    ## which is basically penalizing something different here.
    ## best we can do here is to use E only as a regularizer.
    n <- rep(1, nobs)
    ## should E be used unscaled or not?..
    use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      jj <- attr(x,"lpi")
      start <- rep(0,ncol(x))
      yt1 <- if (family$link[[1]]=="identity") y else 
        family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
      x1 <- x[,jj[[1]],drop=FALSE]
      e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
      #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
      if (use.unscaled) {
        qrx <- qr(rbind(x1,e1))
        x1 <- rbind(x1,e1)
        startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
        startji[!is.finite(startji)] <- 0       
      } else startji <- pen.reg(x1,e1,yt1)
      start[jj[[1]]] <- startji
      lres1 <- log(abs(y-family$linfo[[1]]$linkinv(x[,jj[[1]],drop=FALSE]%*%start[jj[[1]]])))
      x1 <-  x[,jj[[2]],drop=FALSE];e1 <- E[,jj[[2]],drop=FALSE]
      #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
      if (use.unscaled) {
        x1 <- rbind(x1,e1)
        startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
        startji[!is.finite(startji)] <- 0
      } else startji <- pen.reg(x1,e1,lres1)
      start[jj[[2]]] <- startji
      start[jj[[3]]] <- 0
      start[jj[[4]]] <- 0}
  }) ## initialize gaulss
  
  #   postproc <- expression({  ####### XXX ??? #######
  #    object$family$family <- 
  #      paste("logF(",round(object$family$getTheta(TRUE),3),")",sep="")
  #   })
  
  #   rd <- function(mu,wt,scale) {  ####### XXX TODO #######
  #     Theta <- exp(get(".Theta"))
  #     rnbinom(mu,size=Theta,mu=mu)
  #   }
  #   
  #   qf <- function(p,mu,wt,scale) {  ####### XXX TODO #######
  #     Theta <- exp(get(".Theta"))
  #     qnbinom(p,size=Theta,mu=mu)
  #   }
  
  
  structure(list(family="shashlss",ll=ll, link=paste(link), nlp=npar,
                 tri = mgcv:::trind.generator(npar), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 #postproc=postproc,
                 residuals=residuals,
                 #rd=rd,
                 #predict=predict,
                 linfo = stats, ## link information list
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
                 ls=1, ## signals that ls not needed here
                 available.derivs = 1 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## shashlss
