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
## Modified log-F density
#######################

logFlss <- function(link = list("identity", "log"), qu, lam, theta, remInter = TRUE) 
{ 
  # Some checks
  if( !remInter ){
    if( theta != 0 ){ stop("remInter == FALSE, but theta != 0") }
    theta <- 0 }
  
  if( !is.na(qu) && (findInterval(qu, c(0, 1) )!=1) ) stop("qu should be in (0, 1)")
  
  ## Extended family object for modified log-F, to allow direct estimation of theta
  ## as part of REML optimization. Currently the template for extended family objects.
  ## length(theta)=2; log theta supplied. 
  ## Written by Matteo Fasiolo.
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("logFlss requires 2 links specified as character strings")
  okLinks <- list(c("inverse", "log", "identity", "sqrt"), "log")
  stats <- list()
  param.names <- c("mu", "sigma")
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
      stop(link[[i]]," link not available for ", param.names[i]," parameter of logFlss")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
                          mu.eta=stats[[i]]$mu.eta),
                     class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  } 
  
  env <- new.env(parent = .GlobalEnv)
  
  assign(".qu", qu, envir = env)
  getQu <- function( ) get(".qu")
  putQu <- function(qu) assign(".qu", qu, envir=environment(sys.function()))
  
  assign(".lam", lam, envir = env)
  getLam <- function( ) get(".lam")
  putLam <- function(lam) assign(".lam", lam, envir=environment(sys.function()))
  
  assign(".theta", theta, envir = env)
  getTheta <- function( ) get(".theta")
  putTheta <- function(theta) assign(".theta", theta, envir=environment(sys.function()))
  
  # variance <- function(mu) exp(get(".Theta"))  ##### XXX ##### Necessary?
  
  # validmu <- function(mu) all( is.finite(mu) )
  
  residuals <- function(object, type = c("deviance", "response")) {
    
    tau <- 1 - get(".qu")
    theta <- get(".theta")
    lam <- get(".lam")
    
    mu <- object$fitted[ , 1]
    sig <- object$fitted[ , 2] * exp(theta)
    
    type <- match.arg(type)
    
    # raw residuals  
    rsd <- object$y - sig * lam * ( digamma(lam*tau) - digamma(lam*(1-tau)) ) - mu ####### XXX #######
    
    if (type=="response"){ 
      return(rsd)
    }
    else { ## compute deviance residuals
      sgn <- sign(rsd)
      
      z <- (object$y - mu) / sig
      
      dl <- dlogis(z-mu, 0, lam*sig)
      pl <- plogis(z-mu, 0, lam*sig)
      
      l <- tau * z - lam * .log1pexp( z / lam ) - log( sig * lam * beta(lam*tau, (1-tau)*lam) )
      
      ls <- tau*lam*log(tau) + lam*(1-tau)*log1p(-tau) - log(lam * sig * beta(lam*tau, lam*(1-tau)))
      
      rsd <- pmax(0, 2*(ls - l))
      
      rsd <- sqrt(rsd)*sgn
    }
    rsd
  } ## residuals
  
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv=0, d1b=0, d2b=0, Hp=NULL, rank=0, fh=NULL, D=NULL) {
    ## function defining the gamlss Gaussian model log lik. 
    ## N(mu,sigma^2) parameterized in terms of mu and log(sigma)
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    ##        4 - everything.
    tau <- 1 - get(".qu")
    theta <- get(".theta")
    lam <- get(".lam")
    
    if(!is.null(offset[[1]])) stop("offset not still available for this family")
    
    jj <- attr(X,"lpi") ## extract linear predictor index
    eta <- X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
    mu <- family$linfo[[1]]$linkinv(eta)
    eta1 <- X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]] + theta
    sig <-  family$linfo[[2]]$linkinv(eta1) 
    
    n <- length(y)
    l1 <- matrix(0,n,2)
    
    z <- (y - mu) / sig
    
    dl <- dlogis(z-mu, 0, lam*sig)
    pl <- plogis(z-mu, 0, lam*sig)
    
    l <- sum( tau * z - lam * .log1pexp( z / lam ) - log( sig * lam * beta(lam*tau, (1-tau)*lam) ) )
    
    if (deriv>0) {
      
      dl <- dlogis(y, mu, lam*sig)
      pl <- plogis(y, mu, lam*sig)
      
      l1[ , 1] <- (pl - tau) / sig
      l1[ , 2] <- (z * (pl - tau) - 1) / sig  
      
      ## the second derivatives
      
      l2 <- matrix(0, n, 3)
      
      ## order mm,ms,ss
      l2[ , 1] <- - dl / sig
      l2[ , 2] <- - ((y-mu)*dl + pl - tau) / sig^2
      l2[ , 3] <- (2*z*(tau - pl - 0.5 * (y-mu)*dl) + 1)/sig^2
      
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu), family$linfo[[2]]$d2link(sig))
    }
    
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults
    
    if (deriv>1) {
      
      zl <- z / lam
      der <- .sigmoid(zl, deriv = TRUE)
      
      ## the third derivatives
      ## order mmm,mms,mss,sss
      l3 <- matrix(0,n,4) 
      l3[ , 1] <- der$D2 / (lam^2 * sig^3)
      l3[ , 2] <- (zl*der$D2 + 2*der$D1) / (lam * sig^3)
      l3[ , 3] <- (2*(der$D0-tau) + 4*zl*der$D1 + zl^2*der$D2) / (sig^3)
      l3[ , 4] <- - 3*l2[ , 3]/sig + lam*zl^2/sig^3 * (3*der$D1 + zl*der$D2 + 1/(lam*zl^2)) 
      
      g3 <- cbind(family$linfo[[1]]$d3link(mu), family$linfo[[2]]$d3link(sig))
    }
    
    if (deriv>3) {
      ## the fourth derivatives
      ## order mmmm,mmms,mmss,msss,ssss
      
      l4 <- matrix(0, n, 5) 
      l4[,1] <- - der$D3 / (lam^3 * sig^4)
      l4[,2] <- -(zl*der$D3 + 3*der$D2) / (lam^2 * sig^4)
      l4[,3] <- -(zl^2*der$D3 + 6*zl*der$D2 + 6*der$D1) / (lam*sig^4)
      l4[,4] <- -3*l3[ , 3]/sig - zl/sig^4*(6*der$D1 + 6*zl*der$D2 + zl^2*der$D3 )
      l4[,5] <- -4*(2*l3[ , 4] + 3*l2[ , 3]/sig)/sig - lam*zl^3/sig^4 * (4*der$D2 + zl*der$D3 - 2/(lam*zl^3))
      
      g4 <- cbind(family$linfo[[1]]$d4link(mu), family$linfo[[2]]$d4link(sig))
    }
    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
      
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- mgcv:::gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)
      
      ## get the gradient and Hessian...
      ret <- mgcv:::gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
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
      x1 <-  x[,jj[[2]],drop=FALSE];
      e1 <- E[,jj[[2]],drop=FALSE]
      #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
      if (use.unscaled) {
        x1 <- rbind(x1,e1)
        startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
        startji[!is.finite(startji)] <- 0
      } else startji <- pen.reg(x1,e1,lres1)
      start[jj[[2]]] <- startji
    }
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
  
  environment(putTheta) <- environment(getTheta) <- environment(putLam) <- environment(getLam) <- 
    environment(ll) <- environment(residuals) <- environment(putQu) <- environment(getQu) <- env
  
  structure(list(family="logFlss",ll=ll,link=paste(link),nlp=2,
                 tri = mgcv:::trind.generator(2), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 drop.intercept = c(FALSE, remInter),
                 #postproc=postproc,
                 residuals=residuals,
                 getLam = getLam, putLam = putLam, getTheta = getTheta, putTheta = putTheta, putQu=putQu, getQu=getQu, 
                 #rd=rd,
                 #predict=predict,
                 linfo = stats, ## link information list
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
                 ls=1, ## signals that ls not needed here
                 available.derivs = 2 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## logF

