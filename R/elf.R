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
#' @param lam parameter lambda of the ELF density, it must be positive. See Fasiolo et al. (2017) for details.
#' @return An object inheriting from mgcv's class \code{extended.family}.
#' @details This function is meant for internal use only.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon N. Wood. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2017). Fast calibrated additive quantile regression. 
#'             Available at \url{https://github.com/mfasiolo/qgam/blob/master/draft_qgam.pdf}.
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
#'            family = elf(theta = 0, lam = 0.5, qu = 0.5), data = dat)
#' plot(fit, scale = FALSE, pages = 1)     
#' 
#' # Using qgam: RECOMMENDED
#' fit <- qgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, err = 0.05, qu = 0.8)
#' plot(fit, scale = FALSE, pages = 1)      
#'
#' @export elf
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
elf <- function (theta = NULL, link = "identity", qu, lam) { 
  
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
    else stop(linktemp, " link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"")
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
  
  assign(".lam", lam, envir = env)
  getLam <- function( ) get(".lam")
  putLam <- function(lam) assign(".lam", lam, envir=environment(sys.function()))
  
  # variance <- function(mu) exp(get(".Theta"))  ##### XXX ##### Necessary?
  
  validmu <- function(mu) all( is.finite(mu) )
  
  dev.resids <- function(y, mu, wt, theta=NULL) {        ##### XXX #####
    if( is.null(theta) ) theta <- get(".Theta")
    tau <- 1 - get(".qu")
    lam <- get(".lam")
    
    sig <- exp(theta)
    
    
    z <- (y - drop(mu)) / sig
    
    term <- tau*lam*log(tau) + lam*(1-tau)*log1p(-tau) - tau*z + lam*log1pexp( z / lam )
    
    2 * wt * term
  }
  
  Dd <- function(y, mu, theta, wt, level=0) {
    
    tau <- 1 - get(".qu")
    lam <- get(".lam")
    mu <- drop(mu)
    
    ## derivatives of the deviance...
    sig <- exp(theta)
    
    z <- (y - mu) / sig
    
    dl <- dlogis(y-mu, 0, lam*sig)
    pl <- plogis(y-mu, 0, lam*sig)
    
    r <- list()
    ## get the quantities needed for IRLS. 
    ## Dmu2eta2 is deriv of D w.r.t mu twice and eta twice,
    ## Dmu is deriv w.r.t. mu once, etc...
    r$Dmu <- - 2 * wt * ( (pl - tau) / sig )
    r$Dmu2 <- 2 * wt * ( dl / sig )
    r$EDmu2 <- 2 * wt * (tau*(1-tau) / (lam + 1)) / sig^2 ## exact (or estimated) expected weight #### XXX ####
    if (level>0) { ## quantities needed for first derivatives
      zl <- z / lam
      der <- sigmoid(zl, deriv = TRUE)
      
      r$Dth <- - 2 * wt * sig * ( z * (pl - tau) / sig ) 
      #r$Dmuth <- 2 * wt * sig * ( ((y-mu)*dl + pl - tau) / sig^2 )
      r$Dmuth <- 2 * wt * ( ((y-mu)*dl + pl - tau) / sig )
      #r$Dmu3 <- - 2 * wt * ( der$D2 / (lam^2 * sig^3) )  
      r$Dmu3 <- - (2 * wt * ( der$D2 / (lam * sig) )) / (lam*sig^2) 
      
      #D2mDt <- (zl*der$D2 + 2*der$D1) / (lam*sig^3)
      D2mDt <- ((zl*der$D2 + 2*der$D1) / (lam*sig)) / (sig^2)
      r$Dmu2th <- - 2 * wt * sig * D2mDt
    } 
    if (level>1) { ## whole damn lot
      # r$Dmu4 <- 2 * wt * ( der$D3 / (lam^3 * sig^4) )
      r$Dmu4 <- (2 * wt * ( der$D3 / (lam * sig^2) )) / (lam^2 * sig^2)
      #       r$Dth2 <- - 2 * wt * sig * ( z * (pl - tau) / sig + 
      #                                    sig * (2*z*(tau - pl - 0.5 * (y-mu)*dl)/sig^2) )
      r$Dth2 <- - 2 * wt *  ( z * (pl - tau)  + (2*z*(tau - pl - 0.5 * (y-mu)*dl)) )
      #       r$Dmuth2 <- 2 * wt * sig * ( ((y-mu)*dl + pl - tau) / sig^2 - 
      #                                   sig * (2*(der$D0-tau) + 4*zl*der$D1 + zl^2*der$D2) / (sig^3) )
      r$Dmuth2 <- 2 * wt * ( ((y-mu)*dl + pl - tau) / sig - 
                               (2*(der$D0-tau) + 4*zl*der$D1 + zl^2*der$D2) / sig )
      #       r$Dmu2th2 <- - 2 * wt * sig * (D2mDt - sig * (zl^2*der$D3 + 6*zl*der$D2 + 6*der$D1) / (lam*sig^4))
      r$Dmu2th2 <- - 2 * wt * (sig * D2mDt - (zl^2*der$D3 + 6*zl*der$D2 + 6*der$D1) / (lam*sig^2))
      #       r$Dmu3th <- 2 * wt * sig * ( (zl*der$D3 + 3*der$D2) / (lam^2 * sig^4) )
      r$Dmu3th <- 2 * wt * ( (zl*der$D3 + 3*der$D2) / (lam * sig) ) / (lam * sig^2)
    }
    r
  }
  
  aic <- function(y, mu, theta=NULL, wt, dev) { ###### XXX only likelihood, no df? 
    if (is.null(theta)) theta <- get(".Theta")
    sig <- exp(theta)
    
    tau <- 1 - get(".qu")
    lam <- get(".lam")
    
    z <- (y - drop(mu)) / sig
    
    term <- - tau * z + lam * log1pexp( z / lam ) + log( sig * lam * beta(lam*tau, (1-tau)*lam) )
    
    2 * sum(term * wt)
  }
  
  ls <- function(y, w, n, theta, scale) { ##### XXX n is number of observations?
    tau <- 1 - get(".qu")
    lam <- get(".lam")
    ## the log saturated likelihood function.
    sig <- exp(theta)
    
    ls <- sum( w * (tau*lam*log(tau) + lam*(1-tau) * log1p(-tau) - log(lam * sig * beta(lam*tau, lam*(1-tau)))) )
    
    #lsth <- - sig * sum(w / sig)
    lsth <- - sum(w)
    
    #     lsth2 <- sig * sum(w / sig^2)
    lsth2 <- sum(w / sig)
    
    list(ls=ls, ## saturated log likelihood
         lsth1=lsth, ## first deriv vector w.r.t theta - last element relates to scale, if free
         lsth2=lsth2) ## Hessian w.r.t. theta, last row/col relates to scale, if free
  }
  
  initialize <- expression({
    
    mustart <- y
    
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
  
  
  environment(dev.resids) <- environment(ls) <- environment(aic) <- environment(Dd) <- environment(getTheta) <-
    #  environment(rd)<- environment(qf) <- environment(variance) <- 
    environment(putTheta) <- environment(putLam) <- environment(getLam) <-
    environment(putQu) <- environment(getQu) <- env
  
  structure(list(family = "elf", link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, dev.resids = dev.resids,Dd=Dd,
                 #variance=variance,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 #postproc=postproc,
                 ls=ls,
                 validmu = validmu, valideta = stats$valideta, n.theta=n.theta, 
                 ini.theta = iniTheta, putTheta=putTheta,getTheta=getTheta, 
                 putQu=putQu, getQu=getQu, 
                 putLam=putLam,getLam=getLam, 
                 use.wz=TRUE
                 #, rd=rd,qf=qf
  ),
  class = c("extended.family","family"))
} ## elf