.egamFit <- function(x, y, sp, Eb, UrS,
                    weights, start, offset, U1, Mp, family, 
                    control, null.coef, needVb) {
  ## Routine for fitting GAMs beyond exponential family.
  ## Inputs as gam.fit3 except that family is of class "extended.family", while
  ## sp contains the vector of extended family parameters, followed by the log smoothing parameters,
  
  scoreType <- "REML"
  deriv <- 0
  theta <- family$getTheta()
  ## penalized <- if (length(UrS)>0) TRUE else FALSE
  
  x <- as.matrix(x)  
  nSp <- length(sp) 
  rank.tol <- .Machine$double.eps*100 ## tolerance to use for rank deficiency
  q <- ncol(x)
  n <- nobs <- nrow(x)  
  
  xnames <- dimnames(x)[[2]]
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  
  ## Now a stable re-parameterization is needed....
  if (length(UrS)) {
    rp <- mgcv:::gam.reparam(UrS, sp, deriv)
    T <- diag(q)
    T[1:ncol(rp$Qs),1:ncol(rp$Qs)] <- rp$Qs
    T <- U1%*%T ## new params b'=T'b old params
    
    null.coef <- t(T)%*%null.coef  
    
    # Start is a list of vectors, each is a different possible initialization
    for(jj in 1:length(start)){ start[[jj]] <- t(T)%*%start[[jj]] }
    
    ## form x%*%T in parallel 
    x <- .Call(mgcv:::C_mgcv_pmmult2,x,T,0,0,control$nthreads)
    rS <- list()
    for (i in 1:length(UrS)) {
      rS[[i]] <- rbind(rp$rS[[i]],matrix(0,Mp,ncol(rp$rS[[i]])))
    } ## square roots of penalty matrices in current parameterization
    Eb <- Eb%*%T ## balanced penalty matrix
    rows.E <- q-Mp
    Sr <- cbind(rp$E,matrix(0,nrow(rp$E),Mp))
    St <- rbind(cbind(rp$S,matrix(0,nrow(rp$S),Mp)),matrix(0,Mp,q))
  } else { 
    T <- diag(q); 
    St <- matrix(0,q,q) 
    rSncol <- rows.E <- Eb <- Sr <- 0   
    rS <- list(0)
    rp <- list(det=0,det1 = 0,det2 = 0,fixed.penalty=FALSE)
  }
  
  ## re-parameterization complete. Initialization....
  
  nvars <- ncol(x)
  if (nvars==0) stop("emtpy models not available")
  if (is.null(weights)) weights <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  
  linkinv <- family$linkinv
  valideta <- family$valideta
  validmu <- family$validmu
  dev.resids <- family$dev.resids
  
  ## need an initial `null deviance' to test for initial divergence...
  ## if (!is.null(start)) null.coef <- start - can be on edge of feasible - not good
  null.eta <- as.numeric(x%*%null.coef + as.numeric(offset))
  
  ## call the families initialization code...
  mustart <- NULL
  eval(family$initialize)
  
  coefold <- mu <- eta <- NULL
  old.pdev <- null.pdev <- sum(dev.resids(y, linkinv(null.eta), weights, theta)) + t(null.coef)%*%St%*%null.coef 
  
  # Calculating pen deviance using each initialization and choosing the best
  tmp <- lapply(start, function(.st){
    if (length(.st) != nvars){stop("Length of start should equal ", nvars, " and correspond to initial coefs for ", deparse(xnames))}
    .eta <- offset + as.vector(if (NCOL(x) == 1) x * .st else x %*% .st)
    .mu <- linkinv(.eta)
    .pdev <- sum(dev.resids(y, .mu, weights, theta)) + t(.st)%*%St%*%.st
    return( list("eta"=.eta, "mu"=.mu, "pdev"=.pdev, "start"=.st) )
  })
  tmp <- tmp[[ which.min(sapply(tmp, "[[", "pdev")) ]]
  coefold <- start <- tmp$start
  eta <- etaold <- tmp$eta
  mu <- tmp$mu
  pdev <- tmp$pdev
  
  # BAD start (it's worse than null.coef) reset everything
  if (pdev>old.pdev){ 
    start <- coefold <- eta <- etaold <- mu <- NULL
  } else { # GOOD start
    old.pdev <- pdev
  } 
  
  ## Initialization of mu and eta (and coefold) if "start" did not do it
  if(is.null(coefold)){ coefold <- null.coef }
  if (is.null(eta)){ 
    eta <- family$linkfun(mustart) 
    mu <- linkinv(eta) 
    etaold <- eta
  }
  
  conv <-  boundary <- FALSE
  
  for (iter in 1:control$maxit) { ## start of main fitting iteration 
    if (control$trace) cat(iter," ")
    dd <- mgcv:::dDeta(y,mu,weights,theta,family,0) ## derivatives of deviance w.r.t. eta
    
    # good <- is.finite(dd$Deta.Deta2)
    
    w <- dd$Deta2 * .5;
    wz <- w*(eta-offset) - .5*dd$Deta
    z <- (eta-offset) - dd$Deta.Deta2
    good <- is.finite(z)&is.finite(w)
    if (control$trace&sum(!good)>0) cat("\n",sum(!good)," not good\n")
    if (sum(!good)) {
      use.wy <- TRUE
      good <- is.finite(w)&is.finite(wz)
      z[!is.finite(z)] <- 0 ## avoid NaN in .C call - unused anyway
    } else use.wy <- family$use.wz
    
    oo <- .C(mgcv:::C_pls_fit1,   
             y=as.double(z[good]),X=as.double(x[good,]),w=as.double(w[good]),wy = as.double(wz[good]),
             E=as.double(Sr),Es=as.double(Eb),n=as.integer(sum(good)),
             q=as.integer(ncol(x)),rE=as.integer(rows.E),eta=as.double(z),
             penalty=as.double(1),rank.tol=as.double(rank.tol),
             nt=as.integer(control$nthreads),use.wy=as.integer(use.wy))
    if (oo$n<0) { ## then problem is indefinite - switch to +ve weights for this step
      if (control$trace) cat("**using positive weights\n")
      # problem is that Fisher can be very poor for zeroes  
      
      ## index weights that are finite and positive 
      good <- is.finite(dd$Deta2)
      good[good] <- dd$Deta2[good]>0 
      w[!good] <- 0
      wz <- w*(eta-offset) - .5*dd$Deta
      z <- (eta-offset) - dd$Deta.Deta2
      good <- is.finite(z)&is.finite(w) 
      if (sum(!good)) {
        use.wy <- TRUE
        good <- is.finite(w)&is.finite(wz)
        z[!is.finite(z)] <- 0 ## avoid NaN in .C call - unused anyway
      } else use.wy <- family$use.wz
      
      oo <- .C(mgcv:::C_pls_fit1, ##C_pls_fit1,
               y=as.double(z[good]),X=as.double(x[good,]),w=as.double(w[good]),wy = as.double(wz[good]),
               E=as.double(Sr),Es=as.double(Eb),n=as.integer(sum(good)),
               q=as.integer(ncol(x)),rE=as.integer(rows.E),eta=as.double(z),
               penalty=as.double(1),rank.tol=as.double(rank.tol),
               nt=as.integer(control$nthreads),use.wy=as.integer(use.wy))
    }
    # if(control$epsilon == Inf){ ### MATTEO ###################################
    #   startOld <- drop(start)
    #   lprOld <- drop(x%*%startOld)
    # }
    start <- oo$y[1:ncol(x)] ## current coefficient estimates
    penalty <- oo$penalty ## size of penalty
    
    eta <- drop(x%*%start) ## the linear predictor (less offset)
    
    if (any(!is.finite(start))) { ## test for breakdown
      conv <- FALSE
      warning("Non-finite coefficients at iteration ", 
              iter)
      return(list(REML=NA)) ## return immediately signalling failure
    }        
    
    mu <- linkinv(eta <- eta + offset)
    dev <- sum(dev.resids(y, mu, weights,theta)) 
    
    ########################## MATTEO ###################################
    # if(control$epsilon == Inf)
    # {
    #   .myObj <- function(.alpha){
    #     .param <- .alpha * start + (1-.alpha) * startOld
    #     .lpr <- .alpha * eta + (1-.alpha) * lprOld
    #     .mu <- linkinv(.lpr)
    #     .dev <- sum(dev.resids(y, .mu, weights, theta))
    #     return(.dev)
    #   }
    #   
    #   .opt <- optimize(.myObj, c(0, 2))$minimum
    #   cat(" ", .opt)
    #   
    #   start <- .opt * start + (1-.opt) * startOld
    #   eta <- .opt * eta + (1-.opt) * lprOld              
    #   mu <- linkinv(eta)
    #   dev <- sum(dev.resids(y, mu, weights,theta))
    #   
    #   boundary <- TRUE
    #   penalty <- t(start)%*%St%*%start ## reset penalty too
    # }
    ########################## MATTEO ###################################
    
    ## now step halve under non-finite deviance...
    if (!is.finite(dev)) {
      if (is.null(coefold)) {
        if (is.null(null.coef)) 
          stop("no valid set of coefficients has been found:please supply starting values", 
               call. = FALSE)
        ## Try to find feasible coefficients from the null.coef and null.eta
        coefold <- null.coef
        etaold <- null.eta
      }
      #warning("Step size truncated due to divergence", 
      #            call. = FALSE)
      ii <- 1
      while (!is.finite(dev)) {
        if (ii > control$maxit) 
          stop("inner loop 1; can't correct step size")
        ii <- ii + 1
        start <- (start + coefold)/2
        eta <- (eta + etaold)/2               
        mu <- linkinv(eta)
        dev <- sum(dev.resids(y, mu, weights,theta))
        
      }
      boundary <- TRUE
      penalty <- t(start)%*%St%*%start ## reset penalty too
      if (control$trace) 
        cat("Step halved: new deviance =", dev, "\n")
    } ## end of infinite deviance correction
    
    ## now step halve if mu or eta are out of bounds... 
    if (!(valideta(eta) && validmu(mu))) {
      #warning("Step size truncated: out of bounds", 
      #         call. = FALSE)
      ii <- 1
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > control$maxit) 
          stop("inner loop 2; can't correct step size")
        ii <- ii + 1
        start <- (start + coefold)/2
        eta <- (eta + etaold)/2 
        mu <- linkinv(eta)
      }
      boundary <- TRUE
      dev <- sum(dev.resids(y, mu, weights))
      penalty <- t(start)%*%St%*%start ## need to reset penalty too
      if (control$trace) 
        cat("Step halved: new deviance =", dev, "\n")
    } ## end of invalid mu/eta handling
    
    ## now check for divergence of penalized deviance....
    
    pdev <- dev + penalty  ## the penalized deviance 
    if (control$trace) cat("penalized deviance =", pdev, "\n")
    
    div.thresh <- 10*(.1+abs(old.pdev))*.Machine$double.eps^.5
    
    if (pdev-old.pdev>div.thresh) { ## solution diverging
      ii <- 1 ## step halving counter
      if (iter == 1 && (pdev-null.pdev>div.thresh)) { ## Doing worse than null.coef at 1st iterat -> shrink towards zero
        etaold <- null.eta; coefold <- null.coef; old.pdev <- null.pdev
      }
      while (pdev - old.pdev > div.thresh)  { ## step halve until pdev <= old.pdev
        if (ii > 100) 
          stop("inner loop 3; can't correct step size")
        ii <- ii + 1
        start <- (start + coefold)/2 
        eta <- (eta + etaold)/2               
        mu <- linkinv(eta)
        dev <- sum(dev.resids(y, mu, weights,theta))
        
        pdev <- dev + t(start)%*%St%*%start ## the penalized deviance
        if (control$trace) 
          cat("Step halved: new penalized deviance =", pdev, "\n")
      }
    } ## end of pdev divergence
    
    ## convergence testing...
    
    if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) {
      ## Need to check coefs converged adequately, to ensure implicit differentiation
      ## ok. Testing coefs unchanged is problematic under rank deficiency (not guaranteed to
      ## drop same parameter every iteration!)       
      grad <- 2 * t(x[good,])%*%((w[good]*(x%*%start)[good]-wz[good]))+ 2*St%*%start 
      if (max(abs(grad)) > control$epsilon*max(abs(start+coefold))/2) {
        old.pdev <- pdev  ## not converged quite enough
        coef <- coefold <- start
        etaold <- eta 
        ##muold <- mu
      } else { ## converged
        conv <- TRUE
        coef <- start
        break 
      }
    } else { ## not converged
      old.pdev <- pdev
      coef <- coefold <- start
      etaold <- eta 
    }
  } ## end of main loop
  
  ## so at this stage the model has been fully estimated
  coef <- as.numeric(T %*% coef)
  
  ## now obtain derivatives, if these are needed...
  # check.derivs <- FALSE
  # while (check.derivs) { ## debugging code to check derivatives
  #   eps <- 1e-7
  #   fmud.test(y,mu,weights,theta,family,eps = eps)
  #   fetad.test(y,mu,weights,theta,family,eps = eps)
  # }   
  
  ##################### Calculate the Bayesian covariance matrix
  if( needVb )
  {
    dd <- mgcv:::dDeta(y,mu,weights,theta,family,deriv)
    w <- dd$Deta2 * .5
    z <- (eta-offset) - dd$Deta.Deta2 ## - .5 * dd$Deta[good] / w
    wf <- pmax(0,dd$EDeta2 * .5) ## Fisher type weights 
    wz <- w*(eta-offset) - 0.5*dd$Deta ## Wz finite when w==0
    
    gdi.type <- if (any(abs(w)<.Machine$double.xmin*1e20)||any(!is.finite(z))) 1 else 0   
    good <- is.finite(wz)&is.finite(w)   
    
    residuals <- z - (eta - offset)
    residuals[!is.finite(residuals)] <- NA 
    z[!is.finite(z)] <- 0 ## avoid passing NA etc to C code  
    
    ntot <- length(theta) + length(sp)
    rSncol <- unlist(lapply(UrS,ncol))
    ## Now drop any elements of dd that have been dropped in fitting...
    if (sum(!good)>0) { ## drop !good from fields of dd, weights and pseudodata
      z <- z[good]; w <- w[good]; wz <- wz[good]; wf <- wf[good]
      dd$Deta <- dd$Deta[good];dd$Deta2 <- dd$Deta2[good] 
      dd$EDeta2 <- dd$EDeta2[good]
      if (deriv>0) dd$Deta3 <- dd$Deta3[good]
      if (deriv>1) dd$Deta4 <- dd$Deta4[good]
      if (length(theta)>1) {
        if (deriv>0) {  
          dd$Dth <- dd$Dth[good,]; 
          dd$Detath <- dd$Detath[good,]; dd$Deta2th <- dd$Deta2th[good,]
          if (deriv>1) {  
            dd$Detath2 <- dd$Detath2[good,]; dd$Deta3th <- dd$Deta3th[good,]
            dd$Deta2th2 <- dd$Deta2th2[good,];dd$Dth2 <- dd$Dth2[good,]
          }
        }
      } else {
        if (deriv>0) { 
          dd$Dth <- dd$Dth[good]; 
          dd$Detath <- dd$Detath[good]; dd$Deta2th <- dd$Deta2th[good]
          if (deriv>1) {
            dd$Detath2 <- dd$Detath2[good]; dd$Deta3th <- dd$Deta3th[good]
            dd$Deta2th2 <- dd$Deta2th2[good]; dd$Dth2 <- dd$Dth2[good]
          }
        } 
      }
    }
    
    oo <- .C(mgcv:::C_gdi2,
             X=as.double(x[good,]),E=as.double(Sr),Es=as.double(Eb),rS=as.double(unlist(rS)),
             U1 = as.double(U1),sp=as.double(exp(sp)),theta=as.double(theta),
             z=as.double(z),w=as.double(w),wz=as.double(wz),wf=as.double(wf),Dth=as.double(dd$Dth),
             Det=as.double(dd$Deta),
             Det2=as.double(dd$Deta2),Dth2=as.double(dd$Dth2),Det.th=as.double(dd$Detath),
             Det2.th=as.double(dd$Deta2th),Det3=as.double(dd$Deta3),Det.th2 = as.double(dd$Detath2),
             Det4 = as.double(dd$Deta4),Det3.th=as.double(dd$Deta3th), Deta2.th2=as.double(dd$Deta2th2),
             beta=as.double(coef),b1=as.double(rep(0,ntot*ncol(x))),w1=as.double(rep(0,ntot*length(z))),
             D1=as.double(rep(0,ntot)),D2=as.double(rep(0,ntot^2)),
             P=as.double(0),P1=as.double(rep(0,ntot)),P2 = as.double(rep(0,ntot^2)),
             ldet=as.double(1-2*(scoreType=="ML")),ldet1 = as.double(rep(0,ntot)), 
             ldet2 = as.double(rep(0,ntot^2)),
             rV=as.double(rep(0,ncol(x)^2)),
             rank.tol=as.double(.Machine$double.eps^.75),rank.est=as.integer(0),
             n=as.integer(sum(good)),q=as.integer(ncol(x)),M=as.integer(nSp),
             n.theta=as.integer(length(theta)), Mp=as.integer(Mp),Enrow=as.integer(rows.E),
             rSncol=as.integer(rSncol),deriv=as.integer(deriv),
             fixed.penalty = as.integer(rp$fixed.penalty),nt=as.integer(control$nthreads),
             type=as.integer(gdi.type),dVkk=as.double(rep(0,nSp^2)))
    
    rV <- matrix(oo$rV,ncol(x),ncol(x)) ## rV%*%t(rV)*scale gives covariance matrix 
    rV <- T %*% rV
    Vb <- rV%*%t(rV)
  } else {
    Vb <- NULL
  }
  
  list("coefficients" = coef, "Vb" = Vb)
  
} ## gam.fitF
