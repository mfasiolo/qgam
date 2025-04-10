##########################
#' Fast learning rate calibration for the Gibbs posterior
#' 
#' @description The learning rate (sigma) of the Gibbs posterior is tuned either by calibrating the credible intervals for the fitted
#'              curve, or by minimizing the pinball loss on out-of-sample data. This is done by bootrapping or by k-fold cross-validation. 
#'              Here the loss function is minimized, for each quantile, using a Brent search.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
#' @param qu The quantile of interest. Should be in (0, 1).
#' @param discrete If TRUE then covariate discretisation is used for faster model fitting. See \code{mgcv::}\link[mgcv]{bam} for details.
#' @param err An upper bound on the error of the estimated quantile curve. Should be in (0, 1). 
#'            Since qgam v1.3 it is selected automatically, using the methods of Fasiolo et al. (2017).
#'            The old default was \code{err=0.05}.
#' @param multicore If TRUE the calibration will happen in parallel.
#' @param ncores Number of cores used. Relevant if \code{multicore == TRUE}.
#' @param cluster An object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster.
#' @param paropts a list of additional options passed into the foreach function when parallel computation is enabled. 
#'                This is important if (for example) your code relies on external data or packages: 
#'                use the .export and .packages arguments to supply them so that all cluster nodes 
#'                have the correct environment set up for computing. 
#' @param control A list of control parameters for \code{tuneLearn} with entries: \itemize{
#'                   \item{\code{loss} = loss function use to tune log(sigma). If \code{loss=="cal"} is chosen, then log(sigma) is chosen so that
#'                                       credible intervals for the fitted curve are calibrated. See Fasiolo et al. (2017) for details.
#'                                       If \code{loss=="pin"} then log(sigma) approximately minimizes the pinball loss on the out-of-sample
#'                                       data.}
#'                   \item{\code{sam} = sampling scheme use: \code{sam=="boot"} corresponds to bootstrapping and \code{sam=="kfold"} to k-fold
#'                                      cross-validation. The second option can be used only if \code{ctrl$loss=="pin"}.}
#'                   \item{\code{vtype} = type of variance estimator used to standardize the deviation from the main fit in the calibration.
#'                                        If set to \code{"m"} the variance estimate obtained by the full data fit is used, if set to \code{"b"}
#'                                        than the variance estimated produced by the bootstrap fits are used. By default \code{vtype="m"}.}
#'                   \item{\code{epsB} = positive tolerance used to assess convergence when fitting the regression coefficients on bootstrap data.  
#'                                       In particular, if \code{|dev-dev_old|/(|dev|+0.1)<epsB} then convergence is achieved. 
#'                                       Default is \code{epsB=1e-5}.}
#'                   \item{\code{K} = if \code{sam=="boot"} this is the number of boostrap datasets, while if \code{sam=="kfold"} this is the 
#'                                    number of folds. By default \code{K=50}.}
#'                   \item{\code{init} = an initial value for the log learning rate (log(sigma)). 
#'                                       By default \code{init=NULL} and the optimization is initialized by other means.}
#'                   \item{\code{brac} = initial bracket for Brent method. By default \code{brac=log(c(0.5, 2))}, so the initial 
#'                                       search range is \code{(init + log(0.5), init + log(2))}.}
#'                   \item{\code{tol} = tolerance used in the Brent search. By default \code{tol=.Machine$double.eps^0.25}.
#'                                      See \code{?optimize} for details.}
#'                   \item{\code{aTol} = Brent search parameter. If the solution to a Brent get closer than 
#'                                       \code{aTol * abs(diff(brac))} to one of the extremes of the bracket, the optimization is
#'                                       stop and restarted with an enlarged and shifted bracket. \code{aTol=0.05} should be > 0 and values > 0.1
#'                                       don't quite make sense. By default \code{aTol=0.05}.}
#'                   \item{\code{redWd} = parameter which determines when the bracket will be reduced.
#'                                        If \code{redWd==10} then the bracket is halved if the nearest solution
#'                                        falls within the central 10\% of the bracket's width. By default \code{redWd = 10}.}
#'                   \item{\code{b} = offset parameter used by the mgcv::gauslss, which we estimate to initialize the quantile
#'                                    fit (when a variance model is used). By default \code{b=0}.}
#'                   \item{\code{link} = Link function to be used. See \code{?elf} and \code{?elflss} for defaults.}
#'                   \item{\code{verbose} = if TRUE some more details are given. By default \code{verbose=FALSE}.}
#'                   \item{\code{progress} = if TRUE progress in learning rate estimation is reported via printed text.
#'                                           \code{TRUE} by default.}
#' }
#' @param argGam A list of parameters to be passed to \code{mgcv::gam}. This list can potentially include all the arguments listed
#'               in \code{?gam}, with the exception of \code{formula}, \code{family} and \code{data}.
#' @return A list with entries: \itemize{
#'                   \item{\code{lsig} = a vector containing the values of log(sigma) that minimize the loss function, 
#'                                       for each quantile.}
#'                   \item{\code{err} = the error bound used for each quantile. Generally each entry is identical to the
#'                                      argument \code{err}, but in some cases the function increases it to enhance stability.}
#'                   \item{\code{ranges} = the search ranges by the Brent algorithm to find log-sigma, for each quantile. }
#'                   \item{\code{store} = a list, where the i-th entry is a matrix containing all the locations (1st row) at which
#'                                        the loss function has been evaluated and its value (2nd row), for the i-th quantile.}
#'                   \item{\code{final_fit} = a list, where the i-th entry is a list with the estimated conditional quantile 
#'                                             (\code{mustart}), smoothing parameters (\code{in.out$sp}) and scale parameter 
#'                                             (\code{in.out$scale}, should always be equal to 1) at the optimal value of log(sigma).}
#' }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Wood, S.N., Zaffran, M., Nedellec, R. and Goude, Y., 2020. 
#'             Fast calibrated additive quantile regression. 
#'             Journal of the American Statistical Association (to appear).
#'             \doi{10.1080/01621459.2020.1725521}.
#' @examples
#' library(qgam); library(MASS)
#' 
#' ###
#' # Single quantile fit
#' ###
#' # Calibrate learning rate on a grid
#' set.seed(5235)
#' tun <- tuneLearnFast(form = accel~s(times,k=20,bs="ad"), 
#'                      data = mcycle, 
#'                      qu = 0.2)
#' 
#' # Fit for quantile 0.2 using the best sigma
#' fit <- qgam(accel~s(times, k=20, bs="ad"), data = mcycle, qu = 0.2, lsig = tun$lsig)
#' 
#' pred <- predict(fit, se=TRUE)
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", 
#'      ylim = c(-150, 80))
#' lines(mcycle$times, pred$fit, lwd = 1)
#' lines(mcycle$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
#' lines(mcycle$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2) 
#' 
#' ###
#' # Multiple quantile fits
#' ###
#' # Calibrate learning rate on a grid
#' quSeq <- c(0.25, 0.5, 0.75)
#' set.seed(5235)
#' tun <- tuneLearnFast(form = accel~s(times, k=20, bs="ad"), 
#'                      data = mcycle, 
#'                      qu = quSeq)
#' 
#' # Fit using estimated sigmas
#' fit <- mqgam(accel~s(times, k=20, bs="ad"), data = mcycle, qu = quSeq, lsig = tun$lsig)
#' 
#' # Plot fitted quantiles
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", 
#'      ylim = c(-150, 80))
#' for(iq in quSeq){
#'   pred <- qdo(fit, iq, predict)
#'   lines(mcycle$times, pred, col = 2)
#' }  
#' 
#' \dontrun{
#' # You can get a better fit by letting the learning rate change with "accel"
#' # For instance
#' tun <- tuneLearnFast(form = list(accel ~ s(times, k=20, bs="ad"), ~ s(times)), 
#'                       data = mcycle, 
#'                       qu = quSeq)
#' 
#' fit <- mqgam(list(accel ~ s(times, k=20, bs="ad"), ~ s(times)),
#'              data = mcycle, qu = quSeq, lsig = tun$lsig)
#' 
#' # Plot fitted quantiles
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", 
#'      ylim = c(-150, 80))
#' for(iq in quSeq){
#'   pred <- qdo(fit, iq, predict)
#'   lines(mcycle$times, pred, col = 2)
#' }
#' } 
#'
tuneLearnFast <- function(form, data, qu, discrete = FALSE, err = NULL,
                          multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                          control = list(), argGam = NULL)
{ 
  discrete <- .should_we_use_discrete(form = form, discrete = discrete)
  
  gam_name <- ifelse(discrete, "bam", "gam")
  
  # Removing all NAs, unused variables and factor levels from data
  data <- .cleanData(.dat = data, .form = form, .drop = argGam$drop.unused.levels)
  
  n <- nrow(data)
  nq <- length(qu)
  
  if( !is.null(err) && length(err) != nq ){
    if(length(err) == 1) {
      err <- rep(err, nq)
    } else {
      stop("\"err\" should either be a scalar or a vector of the same length as \"qu\".")
    }
  }
  
  # Setting up control parameter
  ctrl <- list( "loss" = "calFast", "sam" = "boot", "vtype" = "m", "epsB" = 1e-5,
                "init" = NULL, "brac" = log( c(1/2, 2) ),  "K" = 50,
                "redWd" = 10, "tol" = .Machine$double.eps^0.25, "aTol" = 0.05, "b" = 0,
                "init_qgam" = NULL,
                "link" = "identity",
                "verbose" = FALSE, "progress" = TRUE )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  if( ctrl$progress == "none" ) { ctrl$progress <- FALSE }
  if( !(ctrl$vtype%in%c("m", "b")) ) stop("control$vtype should be either \"m\" or \"b\" ")
  if( !(ctrl$loss%in%c("calFast", "cal", "pin")) ) stop("control$loss should be either \"cal\", \"pin\" or \"calFast\" ")
  if( !(ctrl$sam%in%c("boot", "kfold")) ) stop("control$sam should be either \"boot\" or \"kfold\" ")
  if( (ctrl$loss=="cal") && (ctrl$sam=="kfold")  ) stop("You can't use control$sam == \"kfold\" when ctrl$loss==\"cal\" ")
  if(discrete && ctrl$loss != "calFast"){ stop("discrete = TRUE can be used only with control$loss != \"calFast\"") }
  
  tol <- ctrl[["tol"]]
  brac <- ctrl[["brac"]]
  
  if( length(argGam$sp) && ctrl$loss != c("calFast") ){ stop("Cannot fix smoothing parameters unless control$loss == \"calFast\".") }
  
  # Sanity check
  if( tol > 0.1 * abs(diff(brac)) ) stop("tol > bracket_widths/10, choose smaller tolerance or larger bracket")
  
  if( ctrl$sam == "boot" ){ # Create weights for K boostrap dataset OR...
    wb <- lapply(1:ctrl[["K"]], function(nouse) tabulate(sample(1:n, n, replace = TRUE), n))
  } else { # ... OR for K training sets for CV 
    tmp <- sample(rep(1:ctrl[["K"]], length.out = n), n, replace = FALSE)
    wb <- lapply(1:ctrl[["K"]], function(ii) tabulate(which(tmp != ii), n)) 
  }
  
  tmp <- ctrl$init_qgam
  if(is.null(tmp)){
    tmp <- .init_gauss_fit(form = form, data = data, ctrl = ctrl, argGam = argGam, qu = qu, discrete = discrete)
  }
  gausFit <- tmp$gausFit
  formL <- tmp$formL
  varHat <- tmp$varHat
  initM <- tmp$initM
  
  # Order quantiles so that those close to the median are dealt with first
  oQu <- order( abs(qu - 0.5) )
  
  # (Optional) Initializing the search range for sigma
  if( is.null(ctrl[["init"]]) ){
    # We assume lam~0 and we match the variance of a symmetric (median) Laplace density with that of the Gaussian fit.
    # This is an over-estimate for extreme quantiles, but experience suggests that it's better erring on the upper side.
    tmp <- 0.5 #qu[ oQu[1] ]
    isig <- log(sqrt( mean(varHat) * (tmp^2*(1-tmp)^2) / (2*tmp^2-2*tmp+1) ))
  } else {
    isig <- ctrl[["init"]]
  }
  
  # Create gam object for full data fits
  mObj <- do.call(gam_name, c(list("formula" = formL, "family" = quote(elf(qu = NA, co = NA, theta = NA, link = ctrl$link)),
                                   "data" = quote(data), "fit" = FALSE, "discrete" = discrete), argGam))
  
  # Remove "sp" as it is already been fixed
  argGam <- argGam[ names(argGam) != "sp" ]
  
  # Create gam object for bootstrap fits
  bObj <- do.call(gam_name, c(list("formula" = formL, "family" = quote(elf(qu = NA, co = NA, theta = NA, link = ctrl$link)), "data" = quote(data),
                                   "sp" = if(length(gausFit$sp)){gausFit$sp}else{NULL}, fit = FALSE, "discrete" = discrete), argGam))
  
  # Preparing bootstrap object for gam.fit3
  bObj <- .prepBootObj(obj = bObj, eps = ctrl$epsB, control = argGam$control)
  
  # Preparing reparametrization list and hide it within mObj. This will be needed by the sandwich calibration
  if( ctrl$loss == "calFast" ){
    mObj$hidRepara <- .prepBootObj(obj = mObj, eps = NULL, control = argGam$control)[ c("UrS", "Mp", "U1") ] 
  }
  
  # Create prediction design matrices for each bootstrap sample or CV fold
  class( mObj ) <- c("gam", "glm", "lm") 
  mObj$coefficients <- rep(0, ncol(mObj$X))  # Needed to fool predict.gam
  
  # Stuff needed for the sandwich estimator
  if(discrete){
    pMat <- NULL
    colmX <- XWyd(X=mObj$Xd,w=rep(1,n),y=rep(1,n),k=mObj$kd,ks=mObj$ks,
                  ts=mObj$ts,dt=mObj$dt,v=mObj$v,qc=mObj$qc,drop=mObj$drop)/n
    sandStuff <- list("XFull" = pMat,
                      "EXXT" = XWXd(X=mObj$Xd,w=rep(1,n),k=mObj$kd,ks=mObj$ks,
                                    ts=mObj$ts,dt=mObj$dt,v=mObj$v,qc=mObj$qc,
                                    drop=mObj$drop)/n,  # E(xx^T)
                      "EXEXT" = tcrossprod( colmX, colmX)) # E(x)E(x)^T
  }else{
    pMat <- predict.gam(mObj, newdata = data, type = "lpmatrix")
    sandStuff <- list("XFull" = pMat,
                      "EXXT" = crossprod(pMat, pMat) / n,                    # E(xx^T)
                      "EXEXT" = tcrossprod( colMeans(pMat), colMeans(pMat))) # E(x)E(x)^T
  }
  
  if( multicore ){ 
    # Create cluster
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores) #, exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoParallel(cluster)
    
    # Load "qgam" and user-specified packages
    tmp <- unique( c("qgam", paropts[[".packages"]]) )
    clusterExport(cluster, "tmp", envir = environment())
    clusterEvalQ(cluster, { lapply(tmp, library, character.only = TRUE) })
    paropts[[".packages"]] <- NULL
    
    # Export bootstrap objects, prediction matrix and user-defined stuff
    tmp <- unique( c("bObj", "pMat", "wb", "ctrl", "argGam", ".egamFit", paropts[[".export"]]) )
    clusterExport(cluster, tmp, envir = environment())
    paropts[[".export"]] <- NULL
  }
  
  # Estimated learning rates, num of bracket expansions, error rates and bracket ranges used in bisection
  sigs <- efacts <- errors <- numeric(nq)
  rans <- matrix(NA, nq, 2)
  store <- vector("list", nq)
  final_fit <- vector("list", nq)
  names(sigs) <- names(errors) <- rownames(rans) <- qu
  
  # Here we need bTol > aTol, otherwise the new bracket will be too close to the probable solution
  bTol <- 4*ctrl$aTol
  
  if( is.null(err) ){ err <- .getErrParam(qu = qu, gFit = gausFit, varHat = varHat) }
  
  if(ctrl$progress){ cat("Estimating learning rate. Each dot corresponds to a loss evaluation. \n") }
  for(ii in 1:nq)
  {
    oi <- oQu[ii]
    
    ef <- 1
    
    if(ctrl$progress){ cat("qu =", qu[oi]) }
    
    repeat{
      
      # Compute bracket
      srange <- isig + ef * brac
      
      # Estimate log(sigma) using brent methods with current bracket (srange)
      res  <- .tuneLearnFast(mObj = mObj, bObj = bObj, pMat = pMat, sandStuff = sandStuff, wb = wb, qu = qu[oi], err = err[oi],
                             srange = srange, gausFit = gausFit, varHat = varHat, initM = initM[[oi]],
                             multicore = multicore, cluster = cluster, ncores = ncores, paropts = paropts,  
                             control = ctrl, argGam = argGam)  
      
      # Store loss function evaluations
      store[[oi]] <- cbind(store[[oi]], res[["store"]])
      final_fit[[oi]] <- res[["final_fit"]]
      lsig <- res$minimum
      
      # If solution not too close to boundary store results and determine bracket for next iteration
      if( all(abs(lsig-srange) > ctrl$aTol * abs(diff(srange))) ){ 
        
        sigs[oi] <- lsig
        rans[oi, ] <- srange
        efacts[oi] <- ef
        errors[oi] <- res$err
        
        # Determine what quantile needs to be dealt with next, then choose bracket and initialization using old results
        if(ii < nq)
        {
          kk <- oQu[ which.min(abs(qu[oQu[ii+1]] - qu[oQu[1:ii]])) ]
          isig <- sigs[kk] 
          wd <- abs(diff(rans[kk, ]))
          brac <- c(-1, 1) * wd / 2
          
          # If kk solution close to center of kk bracket, halve the bracket size 
          # (unless the size of the bracket is < 10*tol or the bracket has been expanded in the old iteration)
          if( (abs(isig - mean(rans[kk, ])) < wd/ctrl$redWd) && (wd > 10*tol) && (efacts[kk] == 1)) brac <- brac / 2
        }
        
        break
      }
      
      # If solution is close to bracket boundaries, we shift bracket and expand it
      # This (- wd + bTol*wd)/2 is divided by 2 to make the algorithm more reluctant to reduce lsig
      wd <- abs( diff(brac) )
      isig <- lsig + ifelse(lsig-srange[1] < ctrl$aTol*wd, (- wd + bTol*wd)/2, wd - bTol*wd)
      ef <- 2*ef
    }
    
    if( ctrl$verbose && (nq>1) )
    {
      tseq <- oQu[1:ii]
      tmp <- rans[tseq, , drop = FALSE]
      layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
             heights=c(2, 1))
      par(mai = c(1, 1, 0.1, 0.1))
      plot(qu[tseq], sigs[oQu[1:ii]], ylim = range(as.vector(tmp)), xlim = range(qu), col = 2, 
           ylab = expression("Log(" * sigma * ")"), xlab = "qu")
      points(qu[tseq], tmp[ , 1], pch = 3)
      points(qu[tseq], tmp[ , 2], pch = 3)
      points(qu[tseq], rowMeans(tmp), pch = 3)
      for(zz in 1:ii) segments(qu[oQu[zz]], rowMeans(rans)[oQu[zz]] - abs(diff(tmp[zz, ]))/ctrl$redWd, 
                               qu[oQu[zz]], rowMeans(rans)[oQu[zz]] + abs(diff(tmp[zz, ]))/ctrl$redWd, col = 1)
      plot(qu, efacts, xlab = "qu", "ylab" = "Bracket expansions")  
      plot(qu, errors)
    }
    
    if(ctrl$progress){ cat("done \n") }
    
  }
  
  if( any(errors > err) ){ message("We had to increase \`err\` for some of the quantiles. See fit$calibr$err") }
  
  names(sigs) <- qu
  
  out <- list("lsig" = sigs, "err" = errors, "ranges" = rans, "store" = store, "final_fit" = final_fit)
  attr(out, "class") <- "learnFast"
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return( out )
}

##########################################################################
### Internal version, which works for a single quantile qu
########################################################################## 
.tuneLearnFast <- function(mObj, bObj, pMat, sandStuff, wb, qu, err,
                           srange, gausFit, varHat, initM,
                           multicore, cluster, ncores, paropts, 
                           control, argGam)
{

  init <- list("initM" = initM, "initB" = vector("list", control$K))
  
  # If we get convergence error, we increase "err" up to 0.2. If the error persists (or if the 
  # error is of another nature) we throw an error
  repeat{
    res <- tryCatch(.brent(brac=srange, f=.objFunLearnFast, mObj = mObj, bObj = bObj, wb = wb, init = init, 
                           pMat = pMat, SStuff = sandStuff, qu = qu, ctrl = control, varHat = varHat, err = err, argGam = argGam,
                           multicore = multicore, paropts = paropts, cluster = cluster, t = control$tol, aTol = control$aTol), 
                    error = function(e) e)
    
    if("error" %in% class(res)){
      if( grepl("can't correct step size", res) ) {
        if(err < 0.2){
          err <- min(2*err, 0.2)
          if(control$verbose) message( paste("Increase \"err\" to ", err, " to get convergence") )
        } else {
          stop("I increased \"err\" up to 0.2, but still didn't get convergence.")
        }
      } else {
        stop( res )
      }
    } else { break } }
  
  res[["err"]] <- err
  
  return( res )
}





