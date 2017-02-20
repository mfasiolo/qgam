.gamlssFit <- function(x,y,lsp,Sl,weights,offset,family,control,Mp,start,needVb){
  ## NOTE: offset handling - needs to be passed to ll code
  ## fit models by general penalized likelihood method, 
  ## given doubly extended family in family. lsp is log smoothing parameters
  ## Stabilization strategy:
  ## 1. Sl.repara
  ## 2. Hessian diagonally pre-conditioned if +ve diagonal elements
  ##    (otherwise indefinite anyway)
  ## 3. Newton fit with perturbation of any indefinite hessian
  ## 4. At convergence test fundamental rank on balanced version of 
  ##    penalized Hessian. Drop unidentifiable parameters and 
  ##    continue iteration to adjust others.
  ## 5. All remaining computations in reduced space.
  ##    
  ## Idea is that rank detection takes care of structural co-linearity,
  ## while preconditioning and step 1 take care of extreme smoothing parameters
  ## related problems. 
  deriv <- 0
  
  penalized <- if (length(Sl)>0) TRUE else FALSE
  
  nSp <- length(lsp)
  q <- ncol(x)
  nobs <- length(y)
  
  if (penalized) {
    Eb <- attr(Sl,"E") ## balanced penalty sqrt
    
    ## the stability reparameterization + log|S|_+ and derivs... 
    rp <- mgcv:::ldetS(Sl,rho=lsp,fixed=rep(FALSE,length(lsp)),np=q,root=TRUE) 
    x <- mgcv:::Sl.repara(rp$rp,x) ## apply re-parameterization to x
    Eb <- mgcv:::Sl.repara(rp$rp,Eb) ## root balanced penalty 
    St <- crossprod(rp$E) ## total penalty matrix
    E <- rp$E ## root total penalty
    attr(E,"use.unscaled") <- TRUE ## signal initialization code that E not to be further scaled   
    for(jj in 1:length(start)){ start[[jj]] <- as.numeric(mgcv:::Sl.repara(rp$rp, start[[jj]])) } ## re-para start
    ## NOTE: it can be that other attributes need re-parameterization here
    ##       this should be done in 'family$initialize' - see mvn for an example. 
    
  } else { ## unpenalized so no derivatives required
    deriv <- 0 
    rp <- list(ldetS=0,rp=list())
    St <- matrix(0,q,q)
    E <- matrix(0,0,q) ## can be needed by initialization code
  }
  
  if (is.null(weights)) weights <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  
  # Choose best initialization
  llf <- family$ll
  tmp <- sapply(start, function(.str){
    .llp <- llf(y,x,.str,weights,family,offset=offset,deriv=0)$l - (t(.str)%*%St%*%.str)/2 
    return( .llp )
  })
  coef <- start[[ which.max(tmp) ]]
  
  ## get log likelihood, grad and Hessian (w.r.t. coefs - not s.p.s) ...
  ll <- llf(y,x,coef,weights,family,offset=offset,deriv=1) 
  ll0 <- ll$l - (t(coef)%*%St%*%coef)/2
  rank.checked <- FALSE ## not yet checked the intrinsic rank of problem 
  rank <- q;drop <- NULL
  eigen.fix <- FALSE
  converged <- FALSE
  check.deriv <- FALSE; eps <- 1e-5 
  drop <- NULL;bdrop <- rep(FALSE,q) ## by default nothing dropped
  perturbed <- 0 ## counter for number of times perturbation tried on possible saddle
  for (iter in 1:(2*control$maxit)) { ## main iteration
    ## get Newton step... 
    if (check.deriv) {
      fdg <- ll$lb*0; fdh <- ll$lbb*0
      for (k in 1:length(coef)) {
        coef1 <- coef;coef1[k] <- coef[k] + eps
        ll.fd <- llf(y,x,coef1,weights,family,offset=offset,deriv=1)
        fdg[k] <- (ll.fd$l-ll$l)/eps
        fdh[,k] <- (ll.fd$lb-ll$lb)/eps
      }
    }
    grad <- ll$lb - St%*%coef 
    Hp <- -ll$lbb+St
    D <- diag(Hp)
    indefinite <- FALSE
    if (sum(D <= 0)) { ## Hessian indefinite, for sure
      D <- rep(1,ncol(Hp))
      if (eigen.fix) {
        eh <- eigen(Hp,symmetric=TRUE);
        ev <- abs(eh$values)
        Hp <- eh$vectors%*%(ev*t(eh$vectors))
      } else {
        Ib <- diag(rank)*abs(min(D))
        Ip <- diag(rank)*abs(max(D)*.Machine$double.eps^.5)
        Hp <- Hp  + Ip + Ib
      }
      indefinite <- TRUE
    } else { ## Hessian could be +ve def in which case Choleski is cheap!
      D <- D^-.5 ## diagonal pre-conditioner
      Hp <- D*t(D*Hp) ## pre-condition Hp
      Ip <- diag(rank)*.Machine$double.eps^.5   
    }
    L <- suppressWarnings(chol(Hp,pivot=TRUE))
    mult <- 1
    while (attr(L,"rank") < rank) { ## rank deficient - add ridge penalty 
      if (eigen.fix) {
        eh <- eigen(Hp,symmetric=TRUE);ev <- eh$values
        thresh <- max(min(ev[ev>0]),max(ev)*1e-6)*mult
        mult <- mult*10
        ev[ev<thresh] <- thresh
        Hp <- eh$vectors%*%(ev*t(eh$vectors)) 
        L <- suppressWarnings(chol(Hp,pivot=TRUE))
      } else {
        L <- suppressWarnings(chol(Hp+Ip,pivot=TRUE))
        Ip <- Ip * 100 ## increase regularization penalty
      }
      indefinite <- TRUE
    }
    
    piv <- attr(L,"pivot")
    ipiv <- piv;ipiv[piv] <- 1:ncol(L)
    step <- D*(backsolve(L,forwardsolve(t(L),(D*grad)[piv]))[ipiv])
    
    c.norm <- sum(coef^2)
    if (c.norm>0) { ## limit step length to .1 of coef length
      s.norm <- sqrt(sum(step^2))
      c.norm <- sqrt(c.norm)
      if (s.norm > .1*c.norm) step <- step*0.1*c.norm/s.norm
    }
    ## try the Newton step...
    coef1 <- coef + step 
    ll <- llf(y,x,coef1,weights,family,offset=offset,deriv=1) 
    ll1 <- ll$l - (t(coef1)%*%St%*%coef1)/2
    khalf <- 0;fac <- 2
    while ((!is.finite(ll1)||ll1 < ll0) && khalf < 25) { ## step halve until it succeeds...
      step <- step/fac;coef1 <- coef + step
      ll <- llf(y,x,coef1,weights,family,offset=offset,deriv=0)
      ll1 <- ll$l - (t(coef1)%*%St%*%coef1)/2
      if (ll1>=ll0) {
        ll <- llf(y,x,coef1,weights,family,offset=offset,deriv=1)
      } else { ## abort if step has made no difference
        if (max(abs(coef1-coef))==0) khalf <- 100
      }
      khalf <- khalf + 1
      if (khalf>5) fac <- 5
    } ## end step halve
    
    if (!is.finite(ll1) || ll1 < ll0) { ## switch to steepest descent... 
      step <- -.5*drop(grad)*mean(abs(coef))/mean(abs(grad))
      khalf <- 0
    }
    
    while ((!is.finite(ll1)||ll1 < ll0) && khalf < 25) { ## step cut until it succeeds...
      step <- step/10;coef1 <- coef + step
      ll <- llf(y,x,coef1,weights,family,offset=offset,deriv=0)
      ll1 <- ll$l - (t(coef1)%*%St%*%coef1)/2
      if (ll1>=ll0) {
        ll <- llf(y,x,coef1,weights,family,offset=offset,deriv=1)
      } else { ## abort if step has made no difference
        if (max(abs(coef1-coef))==0) khalf <- 100
      }
      khalf <- khalf + 1
    }
    
    if ((is.finite(ll1)&&ll1 >= ll0)||iter==control$maxit) { ## step ok. Accept and test
      coef <- coef + step
      ## convergence test...
      ok <- (iter==control$maxit||(abs(ll1-ll0) < control$epsilon*abs(ll0) 
                                   && max(abs(grad)) < .Machine$double.eps^.5*abs(ll0))) 
      if (ok) { ## appears to have converged
        if (indefinite) { ## not a well defined maximum
          if (perturbed==5) stop("indefinite penalized likelihood in gam.fit5 ")
          if (iter<4||rank.checked) {
            perturbed <- perturbed + 1
            coef <- coef*(1+(runif(length(coef))*.02-.01)*perturbed) + 
              (runif(length(coef)) - 0.5 ) * mean(abs(coef))*1e-5*perturbed 
            ll <- llf(y,x,coef,weights,family,offset=offset,deriv=1) 
            ll0 <- ll$l - (t(coef)%*%St%*%coef)/2
          } else {        
            rank.checked <- TRUE
            if (penalized) {
              Sb <- crossprod(Eb) ## balanced penalty
              Hb <- -ll$lbb/norm(ll$lbb,"F")+Sb/norm(Sb,"F") ## balanced penalized hessian
            } else Hb <- -ll$lbb/norm(ll$lbb,"F")
            ## apply pre-conditioning, otherwise badly scaled problems can result in
            ## wrong coefs being dropped...
            D <- abs(diag(Hb))
            D[D<1e-50] <- 1;D <- D^-.5
            Hb <- t(D*Hb)*D
            qrh <- qr(Hb,LAPACK=TRUE)
            rank <- Rrank(qr.R(qrh))
            if (rank < q) { ## rank deficient. need to drop and continue to adjust other params
              drop <- sort(qrh$pivot[(rank+1):q]) ## set these params to zero 
              bdrop <- 1:q %in% drop ## TRUE FALSE version
              ## now drop the parameters and recompute ll0...
              lpi <- attr(x,"lpi")
              xat <- attributes(x)
              xat$dim <- xat$dimnames <- NULL
              coef <- coef[-drop]
              St <- St[-drop,-drop]
              x <- x[,-drop] ## dropping columns from model matrix
              if (!is.null(lpi)) { ## need to adjust column indexes as well
                ii <- (1:q)[!bdrop];ij <- rep(NA,q)
                ij[ii] <- 1:length(ii) ## col i of old model matrix is col ij[i] of new 
                
                for (i in 1:length(lpi)) {
                  lpi[[i]] <- ij[lpi[[i]][!(lpi[[i]]%in%drop)]] # drop and shuffle up
                }
              } ## lpi adjustment done
              for (i in 1:length(xat)) attr(x,names(xat)[i]) <- xat[[i]]
              attr(x,"lpi") <- lpi
              attr(x,"drop") <- drop ## useful if family has precomputed something from x
              ll <- llf(y,x,coef,weights,family,offset=offset,deriv=1) 
              ll0 <- ll$l - (t(coef)%*%St%*%coef)/2
            } 
          }
          
        } else { ## not indefinite really converged
          converged <- TRUE
          break
        }
      } else ll0 <- ll1 ## step ok but not converged yet
    } else { ## step failed.
      converged  <- FALSE
      if (is.null(drop)) bdrop <- rep(FALSE,q)
      warning(paste("step failed: max abs grad =",max(abs(grad))))
      break
    }
  } ## end of main fitting iteration
  
  ## at this stage the Hessian (of pen lik. w.r.t. coefs) should be +ve definite,
  ## so that the pivoted Choleski factor should exist...
  if (iter == 2*control$maxit&&converged==FALSE) 
    warning(gettextf("iteration limit reached: max abs grad = %g",max(abs(grad))))
  
  ldetHp <- 2*sum(log(diag(L))) - 2 * sum(log(D)) ## log |Hp|
  
  if (!is.null(drop)) { ## create full version of coef with zeros for unidentifiable 
    fcoef <- rep(0,length(bdrop));fcoef[!bdrop] <- coef
  } else fcoef <- coef
  
  if( !needVb ){ 
    
    return( list("coefficients" = mgcv:::Sl.repara(rp$rp,fcoef,inverse=TRUE)) ) ## undo re-parameterization of coef 
    
  } else {
    
    dVkk <- d1l <- d2l <- d1bSb <- d2bSb <- d1b <- d2b <- d1ldetH <- d2ldetH <- d1b <- d2b <- NULL
    
    ## get grad and Hessian of REML score...
    REML <- -as.numeric(ll$l - t(coef)%*%St%*%coef/2 + rp$ldetS/2  - ldetHp/2  + Mp*log(2*pi)/2)
    
    REML1 <- if (deriv<1) NULL else -as.numeric( # d1l # cancels
      - d1bSb/2 + rp$ldet1/2  - d1ldetH/2 ) 
    
    if (control$trace) {
      cat("\niter =",iter,"  ll =",ll$l,"  REML =",REML,"  bSb =",t(coef)%*%St%*%coef/2,"\n")
      cat("log|S| =",rp$ldetS,"  log|H+S| =",ldetHp,"  n.drop =",length(drop),"\n")
      if (!is.null(REML1)) cat("REML1 =",REML1,"\n")
    }
    REML2 <- if (deriv<2) NULL else -( d2l - d2bSb/2 + rp$ldet2/2  - d2ldetH/2 ) 
    ## bSb <- t(coef)%*%St%*%coef
    lpi <- attr(x,"lpi")
    if (is.null(lpi)) { 
      linear.predictors <- if (is.null(offset)) as.numeric(x%*%coef) else as.numeric(x%*%coef+offset)
      fitted.values <- family$linkinv(linear.predictors) 
    } else {
      fitted.values <- linear.predictors <- matrix(0,nrow(x),length(lpi))
      if (!is.null(offset)) offset[[length(lpi)+1]] <- 0
      for (j in 1:length(lpi)) {
        linear.predictors[,j] <- as.numeric(x[,lpi[[j]],drop=FALSE] %*% coef[lpi[[j]]])
        if (!is.null(offset[[j]])) linear.predictors[,j] <-  linear.predictors[,j] + offset[[j]]
        fitted.values[,j] <- family$linfo[[j]]$linkinv( linear.predictors[,j]) 
      }
    }
    coef <- mgcv:::Sl.repara(rp$rp,fcoef,inverse=TRUE) ## undo re-parameterization of coef 
    
    if (!is.null(drop)&&!is.null(d1b)) { ## create full version of d1b with zeros for unidentifiable 
      db.drho <- matrix(0,length(bdrop),ncol(d1b));db.drho[!bdrop,] <- d1b
    } else db.drho <- d1b
    ## and undo re-para...
    if (!is.null(d1b)) db.drho <- t(mgcv:::Sl.repara(rp$rp,t(db.drho),inverse=TRUE,both.sides=FALSE)) 
    
    ret <- list(coefficients=coef,family=family,y=y,prior.weights=weights,
                fitted.values=fitted.values, linear.predictors=linear.predictors,
                scale.est=1, ### NOTE: needed by newton, but what is sensible here? 
                REML= REML,REML1= REML1,REML2=REML2,
                rank=rank,aic = -2*ll$l, ## 2*edf needs to be added
                ##deviance = -2*ll$l,
                l= ll$l,## l1 =d1l,l2 =d2l,
                lbb = ll$lbb, ## Hessian of log likelihood
                L=L, ## chol factor of pre-conditioned penalized hessian
                bdrop=bdrop, ## logical index of dropped parameters
                D=D, ## diagonal preconditioning matrix
                St=St, ## total penalty matrix
                rp = rp$rp,
                db.drho = db.drho, ## derivative of penalty coefs w.r.t. log sps.
                #bSb = bSb, bSb1 =  d1bSb,bSb2 =  d2bSb,
                S1=rp$ldet1,
                #S=rp$ldetS,S1=rp$ldet1,S2=rp$ldet2,
                #Hp=ldetHp,Hp1=d1ldetH,Hp2=d2ldetH,
                #b2 = d2b)
                niter=iter,H = ll$lbb,dH = ll$d1H,dVkk=dVkk)#,d2H=llr$d2H)
    ## debugging code to allow components of 2nd deriv of hessian w.r.t. sp.s 
    ## to be passed to deriv.check.... 
    #if (!is.null(ll$ghost1)&&!is.null(ll$ghost2)) { 
    #  ret$ghost1 <- ll$ghost1; ret$ghost2 <- ret$ghost2
    #} 
    return( ret )
  }
} ## end of gam.fit5