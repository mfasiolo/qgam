##########################
#' Extended log-F model with variable scale
#' 
#' @description The \code{elflss} family implements the Extended log-F (ELF) density of Fasiolo et al. (2017) and it is supposed
#'              to work in conjuction with the general GAM fitting methods of Wood et al. (2017), implemented by
#'              \code{mgcv}. It differs from the \code{elf} family, because here the scale of the density 
#'              (sigma, aka the learning rate) can depend of the covariates, while in 
#'              while in \code{elf} it is a single scalar. NB this function was use within the \code{qgam} function, but
#'              since \code{qgam} version 1.3 quantile models with varying learning rate are fitted using different methods
#'              (a parametric location-scale model, see Fasiolo et al. (2017) for details.).
#'              
#' @param link vector of two characters indicating the link function for the quantile location and for the log-scale.
#' @param qu parameter in (0, 1) representing the chosen quantile. For instance, to fit the median choose \code{qu=0.5}.
#' @param co positive vector of constants used to determine parameter lambda of the ELF density (lambda = co / sigma).
#' @param theta a scalar representing the intercept of the model for the log-scale log(sigma). 
#' @param remInter if TRUE the intercept of the log-scale model is removed. 
#' @return An object inheriting from mgcv's class \code{general.family}.
#' @details This function is meant for internal use only.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon N. Wood. 
#' @references Fasiolo, M., Wood, S.N., Zaffran, M., Nedellec, R. and Goude, Y., 2020. 
#'             Fast calibrated additive quantile regression. 
#'             Journal of the American Statistical Association (to appear).
#'             \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1725521}.
#'             
#'             Wood, Simon N., Pya, N. and Safken, B. (2017). Smoothing parameter and model selection for 
#'             general smooth models. Journal of the American Statistical Association.
#' @examples
#' \dontrun{
#' set.seed(651)
#' n <- 1000
#' x <- seq(-4, 3, length.out = n)
#' X <- cbind(1, x, x^2)
#' beta <- c(0, 1, 1)
#' sigma =  1.2 + sin(2*x)
#' f <- drop(X %*% beta)
#' dat <- f + rnorm(n, 0, sigma)
#' dataf <- data.frame(cbind(dat, x))
#' names(dataf) <- c("y", "x")
#' 
#' # Fit median using elflss directly: NOT RECOMMENDED
#' fit <- gam(list(y~s(x, bs = "cr"), ~ s(x, bs = "cr")), 
#'            family = elflss(theta = 0, co = rep(0.2, n), qu = 0.5), 
#'            data = dataf)
#' 
#' plot(x, dat, col = "grey", ylab = "y")
#' tmp <- predict(fit, se = TRUE)
#' lines(x, tmp$fit[ , 1])
#' lines(x, tmp$fit[ , 1] + 3 * tmp$se.fit[ , 1], col = 2)
#' lines(x, tmp$fit[ , 1] - 3 * tmp$se.fit[ , 1], col = 2) 
#' }     
#' 
#'
## (c) Simon N. Wood & Matteo Fasiolo
## 2013-2017. Released under GPL2.
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
elflss <- function(link = list("identity", "log"), qu, co, theta, remInter = TRUE) 
{ 
  # Some checks
  if( !remInter ){
    if( theta != 0 ){ stop("remInter == FALSE, but theta != 0") }
    theta <- 0 
  }
  
  if( !is.na(qu) && (findInterval(qu, c(0, 1) )!=1) ) stop("qu should be in (0, 1)")
  
  ## Extended family object for modified log-F, to allow direct estimation of theta
  ## as part of REML optimization. Currently the template for extended family objects.
  ## length(theta)=2; log theta supplied. 
  ## Written by Matteo Fasiolo.
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("elflss requires 2 links specified as character strings")
  okLinks <- list(c("inverse", "log", "identity", "sqrt"), "log")
  stats <- list()
  param.names <- c("mu", "sigma")
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
      stop(link[[i]]," link not available for ", param.names[i]," parameter of elflss")
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
  
  assign(".co", co, envir = env)
  getCo <- function( ) get(".co")
  putCo <- function(co) assign(".co", co, envir=environment(sys.function()))
  
  assign(".theta", theta, envir = env)
  getTheta <- function( ) get(".theta")
  putTheta <- function(theta) assign(".theta", theta, envir=environment(sys.function()))
  
  # variance <- function(mu) exp(get(".Theta"))  ##### XXX ##### Necessary?
  
  # validmu <- function(mu) all( is.finite(mu) )
  
  residuals <- function(object, type = c("deviance", "response")) {
    
    tau <- get(".qu")
    theta <- get(".theta")
    co <- get(".co")
    
    mu <- object$fitted[ , 1]
    sig <- object$fitted[ , 2] * exp(theta) # This will break if the link is not log!!
    lam <- co / sig
    
    type <- match.arg(type)
    
    # Raw residuals: y - E(y)
    rsd <- object$y - sig * lam * ( digamma(lam*(1-tau)) - digamma(lam*tau) ) - mu
    
    if (type=="response"){ 
      return(rsd)
    }
    else { ## compute deviance residuals
      sgn <- sign(rsd)
      
      z <- (object$y - mu) / sig
      
      dl <- dlogis(z-mu, 0, lam*sig)
      pl <- plogis(z-mu, 0, lam*sig)
      
      l <- (1-tau) * z - lam * log1pexp( z / lam ) - log( sig * lam * beta(lam*(1-tau), lam*tau) )
      
      ls <- (1-tau)*lam*log1p(-tau) + lam*tau*log(tau) - log(lam * sig * beta(lam*(1-tau), lam*tau))
      
      rsd <- pmax(0, 2*(ls - l))
      
      rsd <- sqrt(rsd)*sgn
    }
    rsd
  } ## residuals
  
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv=0, d1b=0, d2b=0, Hp=NULL, rank=0, fh=NULL, D=NULL) {
    ## function defining the gamlss Gaussian model log lik. 
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    ##        4 - everything.
    tau <- get(".qu")
    theta <- get(".theta")
    co <- get(".co")
    
    if( !is.null(offset) ){ offset[[3]] <- 0 } # Not sure whether this is needed
    
    discrete <- is.list(X)
    
    jj <- attr(X,"lpi") ## extract linear predictor index
    eta <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[1]]) else X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
    if( !is.null(offset[[1]]) ){ eta <- eta + offset[[1]] }
    mu <- family$linfo[[1]]$linkinv(eta)
    
    eta1 <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[2]]) else X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]]
    if( !is.null(offset[[2]]) ){ eta1 <- eta1 + offset[[2]] }
    sig <-  family$linfo[[2]]$linkinv(eta1) 
    lam <- co / sig 
    
    n <- length(y)
    if( is.null(wt) ) { wt <- numeric(n) + 1 }
    l1 <- matrix(0, n, 2)
    
    z <- (y - mu) / sig
    zc <- (y - mu) / co
    lpxp <- log1pexp( zc ) 
    
    l <- drop(crossprod(wt, (1-tau) * z - co*lpxp/sig - log( co*beta(co*(1-tau)/sig, co*tau/sig) )))

    if (deriv>0) {
      
      # D logBeta / D lam;  D^2 logBeta / D lam^2 
      dBdL <- (1-tau) * digamma(lam*(1-tau)) + tau * digamma(lam*tau) - digamma(lam)  
      d2BdL2 <- (1-tau)^2 * trigamma(lam*(1-tau)) + tau^2 * trigamma(lam*tau) - trigamma(lam)  

      # D lam / D sig;  D^2 lam / D sig^2 
      dLdS <- - co / sig^2
      d2LdS2 <- 2 * co / sig^3

      # D logBeta / D sig;  D^2 logBeta / D sig^2
      dBdS <- dLdS * dBdL
      d2BdS2 <- d2LdS2*dBdL + (dLdS)^2*d2BdL2
      
      gLog <- (tau-1)*(y-mu) + co * lpxp
      
      dl <- dlogis(y, mu, co)
      pl <- plogis(y, mu, co)
      
      dLLKdmu <- (pl - 1 + tau) / sig
      l1[ , 1] <- wt * dLLKdmu
      l1[ , 2] <- wt * ( gLog/sig^2 - dBdS )
      
      ## the second derivatives
      D2LLKdmu2 <- - dl / sig
      l2 <- matrix(0, n, 3)
      ## order mm,ms,ss
      l2[ , 1] <- wt * D2LLKdmu2
      l2[ , 2] <- wt * ( - dLLKdmu / sig )
      l2[ , 3] <- wt * ( - 2*gLog/sig^3 - d2BdS2  )
      
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu), family$linfo[[2]]$d2link(sig))
      
    }
    
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults
    
    if (deriv>1) {
      # D^3 logBeta / D lam^3 ;  D^3 lam / D sig^3;  D^3 logBeta / D sig^3
      d3BdL3 <- (1-tau)^3 * psigamma(lam*(1-tau), 2) + tau^3 * psigamma(lam*tau, 2) - psigamma(lam, 2)
      d3LdS3 <- - 6 * co / sig^4
      d3BdS3 <- d3LdS3*dBdL + 3*dLdS*d2LdS2*d2BdL2 + dLdS^3*d3BdL3
      
      der <- sigmoid(zc, deriv = TRUE)
      
      D3LLKdmu3 <- der$D2 / (sig*co^2)
      
      ## the third derivatives
      ## order mmm,mms,mss,sss
      l3 <- matrix(0,n,4) 
      l3[ , 1] <- wt * D3LLKdmu3
      l3[ , 2] <- wt * ( - D2LLKdmu2 / sig )
      l3[ , 3] <- wt * ( 2 * dLLKdmu / sig^2 )
      l3[ , 4] <- wt * ( 6*gLog/sig^4 - d3BdS3 )
      
      g3 <- cbind(family$linfo[[1]]$d3link(mu), family$linfo[[2]]$d3link(sig))
    }
    
    if (deriv>3) {
      # D^4 logBeta / D lam^4 ;  D^4 lam / D sig^4;  D^4 logBeta / D sig^4
      d4BdL4 <- (1-tau)^4 * psigamma(lam*(1-tau), 3) + tau^4 * psigamma(lam*tau, 3) - psigamma(lam, 3)
      d4LdS4 <- 24 * co / sig^5
      d4BdS4 <- d4LdS4*dBdL + 3*d2LdS2^2*d2BdL2 + 4*dLdS*d3LdS3*d2BdL2 + 6*(dLdS)^2*d2LdS2*d3BdL3 + dLdS^4*d4BdL4 
      
      ## the fourth derivatives, order: mmmm,mmms,mmss,msss,ssss
      l4 <- matrix(0, n, 5) 
      l4[ , 1] <- wt * ( - der$D3 / (sig*co^3) )
      l4[ , 2] <- wt * ( - D3LLKdmu3 / sig )
      l4[ , 3] <- wt * ( 2 * D2LLKdmu2 / sig^2 )
      l4[ , 4] <- wt * ( - 6 * dLLKdmu / sig^3 )
      l4[ , 5] <- wt * ( -24*gLog/sig^5 - d4BdS4 )
      
      g4 <- cbind(family$linfo[[1]]$d4link(mu), family$linfo[[2]]$d4link(sig))
    }
    if (deriv) {
      i2 <- family$tri$i2;    i3 <- family$tri$i3;    i4 <- family$tri$i4
      
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)
      
      ## get the gradient and Hessian...
      ret <- gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                              d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
      
    } else ret <- list()
    ret$l <- l; ret
  } ## end ll
  
  initialize <- expression({ # COPIED EXACTLY FROM gaulss()
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
      if (!is.null(offset)) offset[[3]] <- 0
      yt1 <- if (family$link[[1]]=="identity") y else 
        family$linfo[[1]]$linkfun(abs(y)+max(y)*1e-7)
      if (!is.null(offset[[1]])) yt1 <- yt1 - offset[[1]]
      if (is.list(x)) { ## discrete case
        start <- rep(0,max(unlist(jj)))
        R <- suppressWarnings(chol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[1]])+crossprod(E[,jj[[1]]]),pivot=TRUE))
        Xty <- XWyd(x$Xd,rep(1,length(y)),yt1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[1]])
        piv <- attr(R,"pivot")
        rrank <- attr(R,"rank") 
        startji <- rep(0,ncol(R))
        if (rrank<ncol(R)) {
          R <- R[1:rrank,1:rrank]
          piv <- piv[1:rrank]
        }
        startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
        startji[!is.finite(startji)] <- 0
        start[jj[[1]]] <- startji
        eta1 <- Xbd(x$Xd,start,k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,drop=x$drop,lt=x$lpid[[1]])
        lres1 <- log(abs(y-family$linfo[[1]]$linkinv(eta1)))
        if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
        R <- suppressWarnings(chol(XWXd(x$Xd,w=rep(1,length(y)),k=x$kd,ks=x$ks,ts=x$ts,dt=x$dt,v=x$v,qc=x$qc,nthreads=1,drop=x$drop,lt=x$lpid[[2]])+crossprod(E[,jj[[2]]]),pivot=TRUE))
        Xty <- XWyd(x$Xd,rep(1,length(y)),lres1,x$kd,x$ks,x$ts,x$dt,x$v,x$qc,x$drop,lt=x$lpid[[2]])
        piv <- attr(R,"pivot")
        rrank <- attr(R,"rank")
        startji <- rep(0,ncol(R))
        if (rrank<ncol(R)) {
          R <- R[1:rrank,1:rrank]
          piv <- piv[1:rrank]
        }
        startji[piv] <- backsolve(R,forwardsolve(t(R),Xty[piv]))
        start[jj[[2]]] <- startji
      } else { ## regular case
        start <- rep(0,ncol(x))
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
        if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
        x1 <-  x[,jj[[2]],drop=FALSE];e1 <- E[,jj[[2]],drop=FALSE]
        #ne1 <- norm(e1); if (ne1==0) ne1 <- 1
        if (use.unscaled) {
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(lres1,rep(0,nrow(E))))   
          startji[!is.finite(startji)] <- 0
        } else startji <- pen.reg(x1,e1,lres1)
        start[jj[[2]]] <- startji
      }  
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
  
  environment(putTheta) <- environment(getTheta) <- environment(putCo) <- environment(getCo) <- 
    environment(ll) <- environment(residuals) <- environment(putQu) <- environment(getQu) <- env
  
  structure(list(family="elflss",ll=ll,link=paste(link),nlp=2,
                 tri = trind.generator(2), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 drop.intercept = c(FALSE, remInter),
                 #postproc=postproc,
                 residuals=residuals,
                 getCo = getCo, putCo = putCo, getTheta = getTheta, putTheta = putTheta, putQu=putQu, getQu=getQu, 
                 #rd=rd,
                 #predict=predict,
                 linfo = stats, ## link information list
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
                 ls=1, ## signals that ls not needed here
                 discrete.ok = TRUE,
                 available.derivs = 2 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## elflss
