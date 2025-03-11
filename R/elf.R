##########################
#' Extended log-F model with fixed scale
#' 
#' @description The \code{elf} family implements the Extended log-F density of Fasiolo et al. (2017) and it is supposed
#'              to work in conjuction with the extended GAM methods of Wood et al. (2017), implemented by
#'              \code{mgcv}. It differs from the \code{elflss} family, because here the scale of the density (sigma, aka the learning rate) is a single scalar, 
#'              while in \code{elflss} it can depend on the covariates. At the moment the family is mainly intended for internal use, 
#'              use the \code{qgam} function to fit quantile GAMs based on ELF.
#'  
#' @param theta a scalar representing the log-scale log(sigma). 
#' @param link the link function between the linear predictor and the quantile location.
#' @param qu parameter in (0, 1) representing the chosen quantile. For instance, to fit the median choose \code{qu=0.5}.
#' @param co positive constant used to determine parameter lambda of the ELF density (lambda = co / sigma).
#'           Can be vector valued.
#' @return An object inheriting from mgcv's class \code{extended.family}.
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
#' library(qgam)
#' set.seed(2)
#' dat <- gamSim(1,n=400,dist="normal",scale=2)
#' 
#' # Fit median using elf directly: FAST BUT NOT RECOMMENDED
#' fit <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), 
#'            family = elf(co = 0.1, qu = 0.5), data = dat)
#' plot(fit, scale = FALSE, pages = 1)     
#' 
#' # Using qgam: RECOMMENDED
#' fit <- qgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, qu = 0.8)
#' plot(fit, scale = FALSE, pages = 1)      
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
elf <- function (theta = NULL, link = "identity", qu, co) { 
  
  # Some checks
  if( !is.na(qu) && (findInterval(qu, c(0, 1) )!=1) ) stop("qu should be in (0, 1)")
  
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt")) stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else stop(linktemp, " link not available for elf family; available links are \"identity\", \"log\" and \"sqrt\"")
  }
  ## Theta <-  NULL;
  n.theta <- 1
  if ( !is.null(theta) ) {
    
    iniTheta <- theta ## fixed log theta supplied
    
    n.theta <- 0 ## signal that there are no theta parameters to estimate
    
  } else iniTheta <- 0 ## inital log theta value
  
  env <- new.env(parent = environment(elf)) #.GlobalEnv) ##########!!!!!!!!!!!!!!!!~########################
  
  assign(".Theta", iniTheta, envir = env)
  getTheta <- function(trans=FALSE) if (trans) exp(get(".Theta")) else get(".Theta")
  putTheta <- function(theta) assign(".Theta", theta, envir=environment(sys.function()))
  
  assign(".qu", qu, envir = env)
  getQu <- function( ) get(".qu")
  putQu <- function(qu) assign(".qu", qu, envir=environment(sys.function()))
  
  assign(".co", co, envir = env)
  getCo <- function( ) get(".co")
  putCo <- function(co) assign(".co", co, envir=environment(sys.function()))
  
  # variance <- function(mu) exp(get(".Theta"))  ##### XXX ##### Necessary?
  
  validmu <- function(mu) all( is.finite(mu) )
  
  dev.resids <- function(y, mu, wt, theta=NULL) {        ##### XXX #####
    if( is.null(theta) ) theta <- get(".Theta")
    tau <- get(".qu")
    co <- get(".co")
    
    mu <- drop(mu)
    sig <- exp(theta)
    lam <- co
    sig <- sig * lam / mean(lam)

    term <-
      (1 - tau)*lam*log1p(-tau) +
      lam*tau*log(tau) -
      (1 - tau)*(y - mu) +
      lam*log1pexp((y - mu) / lam)

    2 * wt * term / sig
  }
  
  Dd <- function(y, mu, theta, wt, level=0) {
    
    tau <- get(".qu")
    co <- get(".co")
    mu <- drop(mu)
    
    ## derivatives of the deviance...
    sig <- exp(theta)
    
    lam <- co
    sig <- sig * lam / mean(lam)

    dl <- dlogis(y-mu, 0, lam)
    pl <- plogis(y-mu, 0, lam)
    r <- list()
    ## get the quantities needed for IRLS. 
    ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
    ## Dmu is deriv w.r.t. mu once, etc...
    r$Dmu <- - 2 * wt * ( (pl - 1 + tau) / sig )
    r$Dmu2 <- 2 * wt * ( dl / sig )
    # r$EDmu2 <- 2 * wt * ((1-tau)*tau / (lam + 1)) / sig^2 ## exact (or estimated) expected weight #### XXX ####
    r$EDmu2 <- r$Dmu2 # It make more sense using the observed information everywhere
    if (level>0) { ## quantities needed for first derivatives
      der <- sigmoid((y - mu) / lam, deriv = TRUE)

      term <-
        (1 - tau)*lam*log1p(-tau) +
        lam*tau*log(tau) -
        (1 - tau)*(y - mu) +
        lam*log1pexp((y - mu) / lam)

      r$Dth <- -2 * wt * term / sig
      r$Dmuth <- -r$Dmu
      r$Dmu3 <- -(2 * wt * der$D2) / (sig * lam^2)
      r$Dmu2th <- -r$Dmu2
    }
    if (level>1) { ## whole damn lot
      r$Dmu4 <- (2 * wt * der$D3) / (sig * lam^3)
      r$Dth2 <- -r$Dth
      r$Dmuth2 <- r$Dmu
      r$Dmu2th2 <- r$Dmu2
      r$Dmu3th <- -r$Dmu3
    }
    r
  }
  
  aic <- function(y, mu, theta=NULL, wt, dev) { 
    if (is.null(theta)) theta <- get(".Theta")
    sig <- exp(theta)
    tau <- get(".qu")
    co <- get(".co")
    
    lam <- co
    sig <- sig * lam / mean(lam)
    mu <- drop(mu)

    term <-
      -(1 - tau) * (y - mu) / sig +
      lam * log1pexp( (y - mu) / lam ) / sig +
      log(lam * beta(lam * (1 - tau) / sig, tau * lam / sig))
    2 * sum(term * wt)
  }
  
  ls <- function(y, w, theta, scale) {
    tau <- get(".qu")
    co <- get(".co")
    sig <- exp(theta)
    
    lam <- co
    sig <- sig * lam / mean(lam)
    ## the log saturated likelihood function.
    ls <- sum(w * (
      (1 - tau) * lam * log1p(-tau) / sig +
        lam * tau * log(tau) / sig -
        log(lam) - lbeta(lam * (1 - tau) / sig, lam * tau / sig)
    ))

    lsth <-
      sum(w * (
        -(1 - tau) * log1p(-tau) -
          tau * log(tau) +
          (1 - tau) * digamma(lam * (1 - tau) / sig) +
          tau * digamma(lam * tau / sig) -
          digamma(lam / sig)
      ) * lam / sig)
    
    lsth2 <-
      -lsth -
      sum(w * (
        (1 - tau)^2 * trigamma(lam * (1 - tau) / sig) +
          tau^2 * trigamma(lam * tau / sig) -
          trigamma(lam / sig)
      ) * lam^2 / sig^2)
    
    list(ls=ls, ## saturated log likelihood
         lsth1=lsth, ## first deriv vector w.r.t theta - last element relates to scale, if free
         lsth2=lsth2) ## Hessian w.r.t. theta, last row/col relates to scale, if free
  }
  
  initialize <- expression({
    
    mustart <- quantile(y, family$getQu()) + y * 0 # this ---> y <--- is very bad idea
    
  })
  
  #postproc <- expression({  ####### XXX ??? #######
  #  object$family$family <- 
  #    paste("elf(",round(object$family$getTheta(TRUE),3),")",sep="")
  #})
  
  #   rd <- function(mu,wt,scale) {  ####### XXX TODO #######
  #     Theta <- exp(get(".Theta"))
  #     rnbinom(mu,size=Theta,mu=mu)
  #   }
  #   
  #   qf <- function(p,mu,wt,scale) {  ####### XXX TODO #######
  #     Theta <- exp(get(".Theta"))
  #     qnbinom(p,size=Theta,mu=mu)
  #   }
  
  get.null.coef <- function(G,start=NULL,etastart=NULL,mustart=NULL,...) {
    ## Get an estimate of the coefs corresponding to maximum reasonable deviance...
    y <- G$y
    weights <- G$w
    nobs <- G$n ## ignore codetools warning!!
    ##start <- etastart <- mustart <- NULL
    family <- G$family
    eval(family$initialize) ## have to do this to ensure y numeric
    y <- as.numeric(y)
    mum <- quantile(y, get(".qu")) + 0*y
    etam <- family$linkfun(mum)
    null.coef <- qr.coef(qr(G$X), etam)
    null.coef[is.na(null.coef)] <- 0;
    ## get a suitable function scale for optimization routines
    null.scale <- sum(family$dev.resids(y, mum, weights))/nrow(G$X) 
    list(null.coef=null.coef,null.scale=null.scale)
  }
  
  
  #  environment(rd)<- environment(qf) <- environment(variance) <- 
  environment(dev.resids) <- environment(ls) <- environment(aic) <- environment(Dd) <- 
    environment(getTheta) <- 
    environment(putTheta) <- environment(putCo) <- environment(getCo) <-
    environment(putQu) <- environment(getQu) <- environment(get.null.coef) <- env
  
  structure(list(family = "elf", link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,
                 #variance=variance,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 #postproc=postproc,
                 ls=ls,
                 validmu = validmu, valideta = stats$valideta, n.theta=n.theta, 
                 ini.theta = iniTheta, putTheta=putTheta,getTheta=getTheta, 
                 putQu=putQu, getQu=getQu, 
                 putCo=putCo,getCo=getCo, get.null.coef=get.null.coef,
                 use.wz=TRUE
                 #, rd=rd,qf=qf
  ),
  class = c("extended.family","family"))
} ## elf
