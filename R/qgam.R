##########################
#' Fit a smooth additive quantile regression model
#' 
#' @description This function fits a smooth additive regression model for a single quantile.
#' 
#' @param form A GAM formula, or a list of formulae. See ?mgcv::gam details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula.
#'             By default the variables are taken from environment(formula): typically the environment from which gam is called.
#' @param qu The quantile of interest. Should be in (0, 1).
#' @param lsig The value of the log learning rate used to create the Gibbs posterior. By defauls \code{lsig=NULL} and this
#'             parameter is estimated by posterior calibration described in Fasiolo et al. (2016). Obviously, the function is much faster
#'             if the user provides a value. 
#' @param err An upper bound on the error of the estimated quantile curve. Should be in (0, 1). See Fasiolo et al. (2016) for details.
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
#' @param ... additional arguments passed to \code{mgcv::gam}.
#' @return A \code{gamObject}. See \code{?gamObject}.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2017). Fast calibrated additive quantile regression. Available at
#'             \url{https://arxiv.org/abs/1707.03307}.
#' @examples
#
#' #####
#' # Univariate "car" example
#' ####
#' library(qgam); library(MASS)
#' 
#' # Fit for quantile 0.5 using the best sigma
#' set.seed(6436)
#' fit <- qgam(accel~s(times, k=20, bs="ad"), data = mcycle, err = 0.05, qu = 0.5)
#' 
#' # Plot the fit
#' xSeq <- data.frame(cbind("accel" = rep(0, 1e3), "times" = seq(2, 58, length.out = 1e3)))
#' pred <- predict(fit, newdata = xSeq, se=TRUE)
#' plot(mcycle$times, mcycle$accel, xlab = "Times", ylab = "Acceleration", ylim = c(-150, 80))
#' lines(xSeq$times, pred$fit, lwd = 1)
#' lines(xSeq$times, pred$fit + 2*pred$se.fit, lwd = 1, col = 2)
#' lines(xSeq$times, pred$fit - 2*pred$se.fit, lwd = 1, col = 2)   
#' 
#' #####
#' # Multivariate Gaussian example
#' ####
#' library(qgam)
#' set.seed(2)
#' dat <- gamSim(1,n=400,dist="normal",scale=2)
#' 
#' fit <- qgam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat, err = 0.05, qu = 0.5)
#' plot(fit, scale = FALSE, pages = 1)      
#' 
#' ######
#' # Heteroscedastic example
#' ######
#' \dontrun{
#' set.seed(651)
#' n <- 5000
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
#'             data = dataf, qu = 0.95, lsig = -1.16)
#' 
#' plot(x, dat, col = "grey", ylab = "y")
#' tmp <- predict(fit, se = TRUE)
#' lines(x, tmp$fit[ , 1])
#' lines(x, tmp$fit[ , 1] + 3 * tmp$se.fit[ , 1], col = 2)
#' lines(x, tmp$fit[ , 1] - 3 * tmp$se.fit[ , 1], col = 2)
#' }
#'
qgam <- function(form, data, qu, lsig = NULL, err = 0.05, 
                 multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1, paropts = list(),
                 control = list(), argGam = NULL)
{
  if( length(qu) > 1 ) stop("length(qu) > 1, so you should use mqgam()")
  
  # Removing all NAs and unused levels from data
  if( inherits(data, "groupedData") ) { data <- as.data.frame( data ) }
  data <- droplevels( na.omit( data ) )
  
  # Setting up control parameter (mostly used by tuneLearnFast)
  ctrl <- list("gausFit" = NULL, "verbose" = FALSE, "b" = 0, "link" = if(is.formula(form)){"identity"}else{list("identity", "log")})
  
  # Checking if the control list contains unknown names entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control, verbose = FALSE)
  
  # Gaussian fit, used for initializations 
  if( is.formula(form) ) {
    fam <- "elf"
    if( is.null(ctrl[["gausFit"]]) ) { ctrl$gausFit <- do.call("gam", c(list("formula" = form, "data" = data), argGam)) }
    varHat <- ctrl$gausFit$sig2
  } else {
    fam <- "elflss"
    if( is.null(ctrl[["gausFit"]]) ) { ctrl$gausFit <- do.call("gam", c(list("formula" = form, "data" = data, "family" = gaulss(b=ctrl[["b"]])), argGam)) }
    varHat <- 1/ctrl$gausFit$fit[ , 2]^2
  }   
  
  # Selecting the learning rate sigma
  learn <- NULL
  if( is.null(lsig) ) {  
    learn <- tuneLearnFast(form = form, data = data, qu = qu, err = err, multicore = multicore, cluster = cluster, 
                           ncores = ncores, paropts = paropts, control = ctrl, argGam = argGam)
    lsig <- learn$lsig
  }
  
  # Fit model for fixed log-sigma
  # Do not use 'start' gausFit in gamlss case because it's not to clear how to deal with model for sigma
  if( fam=="elf" && is.null(argGam$start) ) { argGam$start <- coef(ctrl$gausFit) + c(qnorm(qu, 0, sqrt(varHat)), rep(0, length(coef(ctrl$gausFit))-1))  }
  co <- err * sqrt(2*pi*varHat) / (2*log(2))
  fit <- do.call("gam", c(list("formula" = form, "family" = get(fam)(qu = qu, co = co, theta = lsig), "data" = data), argGam))
  
  fit$calibr <- learn
  
  class(fit) <- c("qgam", class(fit))
  
  return( fit )
}


