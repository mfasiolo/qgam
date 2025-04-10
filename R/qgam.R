##########################
#' Fit a smooth additive quantile regression model
#' 
#' @description This function fits a smooth additive regression model for a single quantile.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
#' @param qu The quantile of interest. Should be in (0, 1).
#' @param discrete If TRUE then covariate discretisation is used for faster model fitting. See \code{mgcv::}\link[mgcv]{bam} for details.
#' @param lsig The value of the log learning rate used to create the Gibbs posterior. By defauls \code{lsig=NULL} and this
#'             parameter is estimated by posterior calibration described in Fasiolo et al. (2017). Obviously, the function is much faster
#'             if the user provides a value. 
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
#' @param control A list of control parameters. The only one relevant here is \code{link}, which is the link function
#'                used (see \code{?elf} and \code{?elflss} for defaults). All other control parameters are used by 
#'                \code{tuneLearnFast}. See \code{?tuneLearnFast} for details.
#' @param argGam A list of parameters to be passed to \code{mgcv::gam}. This list can potentially include all the arguments listed
#'               in \code{?gam}, with the exception of \code{formula}, \code{family} and \code{data}.
#' @return A \code{gamObject}. See \code{?gamObject}.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Wood, S.N., Zaffran, M., Nedellec, R. and Goude, Y., 2020. 
#'             Fast calibrated additive quantile regression. 
#'             Journal of the American Statistical Association (to appear).
#'             \doi{10.1080/01621459.2020.1725521}.
#' @references Fasiolo, M., Wood, S.N., Zaffran, M., Nedellec, R. and Goude, Y., 2021. 
#'             qgam: Bayesian Nonparametric Quantile Regression Modeling in R. 
#'             Journal of Statistical Software, 100(9), 1-31, \doi{10.18637/jss.v100.i09}.
#' @examples
#
#' #####
#' # Univariate "car" example
#' ####
#' library(qgam); library(MASS)
#' 
#' # Fit for quantile 0.5 using the best sigma
#' set.seed(6436)
#' fit <- qgam(accel~s(times, k=20, bs="ad"), data = mcycle, qu = 0.5)
#' 
#' # Plot the fit
#' xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
#' pred <- predict(fit, newdata = xSeq, se=TRUE)
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
#' lines(xSeq$times, pred$fit, lwd = 1)
#' lines(xSeq$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
#' lines(xSeq$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)   
#' 
#' \dontrun{
#' # You can get a better fit by letting the learning rate change with "accel"
#' # For instance
#' fit <- qgam(list(accel ~ s(times, k=20, bs="ad"), ~ s(times)),
#'            data = mcycle, qu = 0.8)
#' 
#' pred <- predict(fit, newdata = xSeq, se=TRUE)
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
#' lines(xSeq$times, pred$fit, lwd = 1)
#' lines(xSeq$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
#' lines(xSeq$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)  
#' }
#' 
#' #####
#' # Multivariate Gaussian example
#' ####
#' library(qgam)
#' set.seed(2)
#' dat <- gamSim(1,n=400,dist="normal",scale=2)
#' 
#' fit <- qgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, qu = 0.5)
#' plot(fit, scale = FALSE, pages = 1)      
#' 
#' ######
#' # Heteroscedastic example
#' ######
#' \dontrun{
#' set.seed(651)
#' n <- 2000
#' x <- seq(-4, 3, length.out = n)
#' X <- cbind(1, x, x^2)
#' beta <- c(0, 1, 1)
#' sigma =  1.2 + sin(2*x)
#' f <- drop(X %*% beta)
#' dat <- f + rnorm(n, 0, sigma)
#' dataf <- data.frame(cbind(dat, x))
#' names(dataf) <- c("y", "x")
#' 
#' fit <- qgam(list(y~s(x, k = 30, bs = "cr"), ~ s(x, k = 30, bs = "cr")), 
#'             data = dataf, qu = 0.95)
#' 
#' plot(x, dat, col = "grey", ylab = "y")
#' tmp <- predict(fit, se = TRUE)
#' lines(x, tmp$fit)
#' lines(x, tmp$fit + 2 * tmp$se.fit, col = 2)
#' lines(x, tmp$fit - 2 * tmp$se.fit, col = 2)
#' }
#'
qgam <- function(form, data, qu, discrete = FALSE, lsig = NULL, err = NULL, 
                 multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                 control = list(), argGam = NULL)
{
  if( length(qu) > 1 ) stop("length(qu) > 1, so you should use mqgam()")
  
  discrete <- .should_we_use_discrete(form = form, discrete = discrete)
  
  gam_name <- ifelse(discrete, "bam", "gam")
  
  # Removing all NAs, unused variables and factor levels from data
  data <- .cleanData(.dat = data, .form = form, .drop = argGam$drop.unused.levels)
  
  # Setting up control parameter (mostly used by tuneLearnFast)
  ctrl <- list("verbose" = FALSE, "b" = 0, "link" = "identity")
  
  # Checking if the control list contains unknown names entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control, verbose = FALSE)
  
  tmp <- ctrl$init_qgam
  if(is.null(tmp)){
    tmp <- .init_gauss_fit(form = form, data = data, ctrl = ctrl, argGam = argGam, qu = qu, discrete = discrete)
  }
  ctrl$init_qgam <- tmp
  gausFit <- tmp$gausFit
  formL <- tmp$formL
  varHat <- tmp$varHat
  initM <- tmp$initM
  
  # Get loss smoothness
  if( is.null(err) ){ err <- .getErrParam(qu = qu, gFit = gausFit, varHat = varHat) }

  # Selecting the learning rate sigma
  learn <- NULL
  if( is.null(lsig) ) {  
    learn <- tuneLearnFast(form = form, data = data, qu = qu, discrete = discrete, err = err, multicore = multicore, cluster = cluster, 
                           ncores = ncores, paropts = paropts, control = ctrl, argGam = argGam)
    lsig <- learn$lsig
    err <- learn$err # Over-writing err parameter!
    argGam$mustart <- learn$final_fit[[1]]$mustart
    # Annoyingly, initial coeffs are supplied via "coef" argument in bam() and "start" in gam()
    argGam[[ ifelse(discrete, "coef", "start") ]] <- learn$final_fit[[1]]$coefstart 
    argGam$in.out <- learn$final_fit[[1]]$in.out
  } else {
    if( is.null(argGam$coef) && is.null(argGam$start) && is.null(argGam$mustart)  ) {
      argGam$mustart <- initM$mustart
      argGam[[ ifelse(discrete, "coef", "start") ]] <- initM$coefstart
    }
  }
  
  co <- err * sqrt(2*pi*varHat) / (2*log(2))
  
  # Fit model for fixed log-sigma
  fit <- do.call(gam_name, c(list("formula" = formL, "family" = quote(elf(qu = qu, co = co, theta = lsig, link = ctrl$link)), 
                               "data" = quote(data), "discrete" = discrete), argGam))
  
  fit$calibr <- learn
  
  class(fit) <- c("qgam", class(fit))
  
  return( fit )
  
}


