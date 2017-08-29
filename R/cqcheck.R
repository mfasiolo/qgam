##########################
#' Visually checking a fitted quantile model
#' 
#' @description Given an additive quantile model, fitted using \code{qgam}, \code{cqcheck} provides some plots
#'              that allow to check what proportion of responses, \code{y}, falls below the fitted quantile.
#'  
#' @param obj the output of a \code{qgam} call. 
#' @param v if a 1D plot is required, \code{v} should be either a single character or a numeric vector. In the first case
#'          \code{v} should be the names of one of the variables in the dataframe \code{X}. In the second case, the length
#'          of \code{v} should be equal to the number of rows of \code{X}. If a 2D plot is required, \code{v} should be 
#'          either a vector of two characters or a matrix with two columns.  
#' @param X a dataframe containing the data used to obtain the conditional quantiles. By default it is NULL, in which
#'          case predictions are made using the model matrix in \code{obj$model}.
#' @param y vector of responses. Its i-th entry corresponds to the i-th row of X.  By default it is NULL, in which
#'          case it is internally set to \code{obj$y}.
#' @param nbin a vector of integers of length one (1D case) or two (2D case) indicating the number of bins to be used
#'             in each direction. Used only if \code{bound==NULL}.
#' @param bound in the 1D case it is a numeric vector whose increasing entries represent the bounds of each bin.
#'              In the 2D case a list of two vectors should be provided. \code{NULL} by default. 
#' @param lev the significance levels used in the plots, this determines the width of the confidence 
#'            intervals. Default is 0.05.
#' @param scatter if TRUE a scatterplot is added (using the \code{points} function). FALSE by default.
#' @param ... extra graphical parameters to be passed to \code{plot()}.
#' @return Simply produces a plot.
#' @details Having fitted an additive model for, say, quantile \code{qu=0.4} one would expect that about 40% of the 
#'          responses fall below the fitted quantile. This function allows to visually compare the empirical number
#'          of responses (\code{qu_hat}) falling below the fit with its theoretical value (\code{qu}). In particular, 
#'          the responses are binned, which the bins being constructed along one or two variables (given be arguments
#'          \code{v}). Let (\code{qu_hat[i]}) be the proportion of responses below the fitted quantile in the ith bin.
#'          This should be approximately equal to \code{qu}, for every i. In the 1D case, when \code{v} is a single
#'          character or a numeric vector, \code{cqcheck} provides a plot where: the horizontal line is \code{qu}, 
#'          the dots correspond to \code{qu_hat[i]} and the grey lines are confidence intervals for \code{qu}. The
#'          confidence intervals are based on \code{qbinom(lev/2, siz, qu)}, if the dots fall outside them, then 
#'          \code{qu_hat[i]} might be deviating too much from \code{qu}. In the 2D case, when \code{v} is a vector of two
#'          characters or a matrix with two columns, we plot a grid of bins. The responses are divided between the bins
#'          as before, but now don't plot the confidence intervals. Instead we report the empirical proportions \code{qu_hat[i]}
#'          for the non-empty bin, and with colour the bins in red if \code{qu_hat[i]<qu} and in green otherwise. If       
#'          \code{qu_hat[i]} falls outside the confidence intervals we put an * next to the numeric \code{qu_hat[i]} and
#'          we use more intense colours.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' #######
#' # Bivariate additive model y~1+x+x^2+z+x*z/2+e, e~N(0, 1)
#' #######
#' \dontrun{
#' library(qgam)
#' set.seed(15560)
#' n <- 500
#' x <- rnorm(n, 0, 1); z <- rnorm(n)
#' X <- cbind(1, x, x^2, z, x*z)
#' beta <- c(0, 1, 1, 1, 0.5)
#' y <- drop(X %*% beta) + rnorm(n) 
#' dataf <- data.frame(cbind(y, x, z))
#' names(dataf) <- c("y", "x", "z")
#' 
#' #### Fit a constant model for median
#' qu <- 0.5
#' fit <- qgam(y~1, qu = qu, data = dataf)
#' 
#' # Look at what happens along x: clearly there is non linear pattern here
#' cqcheck(obj = fit, v = c("x"), X = dataf, y = y) 
#' 
#' #### Add a smooth for x
#' fit <- qgam(y~s(x), qu = qu, err = 0.05, data = dataf)
#' cqcheck(obj = fit, v = c("x"), X = dataf, y = y) # Better!
#' 
#' # Lets look across x and z. As we move along z (x2 in the plot) 
#' # the colour changes from green to red
#' cqcheck(obj = fit, v = c("x", "z"), X = dataf, y = y, nbin = c(5, 5))
#' 
#' # The effect look pretty linear
#' cqcheck(obj = fit, v = c("z"), X = dataf, y = y, nbin = c(10))
#' 
#' #### Lets add a linear effect for z 
#' fit <- qgam(y~s(x)+z, qu = qu, data = dataf)
#' 
#' # Looks better!
#' cqcheck(obj = fit, v = c("z"))
#' 
#' # Lets look across x and y again: green prevails on the top-left to bottom-right
#' # diagonal, while the other diagonal is mainly red.
#' cqcheck(obj = fit, v = c("x", "z"), nbin = c(5, 5))
#' 
#' ### Maybe adding an interaction would help?
#' fit <- qgam(y~s(x)+z+I(x*z), qu = qu, data = dataf)
#' 
#' # It does! The real model is: y ~ 1 + x + x^2 + z + x*z/2 + e, e ~ N(0, 1)
#' cqcheck(obj = fit, v = c("x", "z"), nbin = c(5, 5))
#' }
#'
cqcheck <- function(obj, v, X = NULL, y = NULL, nbin = c(10, 10), bound = NULL, lev = 0.05, scatter = FALSE, ...)
{
  #### Set up
  if( is.null(X) ){ 
    X <- obj$model
    if( is.null(y) ){ y <- obj$y }
  } else {
    if( is.null(y) ){ stop("If you provide X you must provide also the corresponding vector of responses y") }
  }
  
  if( length(y)!=nrow(X) ){ stop("length(y)!=nrow(X)") }
  
  ####### Setting up 1D and 2D cases
  if( is.character(v) ){ # Name(s) of variable(s) in X provided OR ...
    if(length(v) == 1){ ## 1D CASE ##
      if( !(v %in% names(X)) ) stop("(v %in% names(X)) == FALSE")
      x1 <- X[[v]]
      x2 <- NULL
    } else {
      if(length(v) == 2){ ## 2D CASE ##
        if( !(v[1] %in% names(X)) ) stop("(v[1] %in% names(X)) == FALSE")
        if( !(v[2] %in% names(X)) ) stop("(v[2] %in% names(X)) == FALSE")
        x1 <- X[[v[1]]]
        x2 <- X[[v[2]]]
      } else { 
        stop("If is.character(v)==TRUE, then length(v) should be either 1 or 2.") 
      }
    }
  } else { # ... actual numeric value of the variable(s) provided 
    if(is.vector(v)){ ## 1D CASE ##
      x1 <- v
      x2 <- NULL
      if(length(v) != nrow(X)){ stop("length(v) != ncol(X)") }
    } else {
      if(is.matrix(v)){ ## 2D CASE ##
        if(ncol(v)!=2){ stop("In the 2D case, v should be a matrix with 2 columns or a vector of 2 characters") } 
        x1 <- v[ , 1]
        x2 <- v[ , 2]
      } else {
        stop("In the 2D case, v should be a matrix with 2 columns or a vector of 2 characters")
      }
    }
  } 
  
  # Discard NAs from X, y, x1 and x2. We don't do this on v, hence if v is numeric it is dangerous to use it from here onwards
  good <- complete.cases(X, y, x1, x2)
  y <- y[ good ]
  X <- X[good, , drop = FALSE]
  x1 <- x1[ good ]
  x2 <- x2[ good ]
  
  # Calculating proportion of observation falling below estimate quantile curve
  n <- nrow(X)
  mu <- as.matrix(predict(obj, newdata = X))[ , 1]
  res <- (mu - y) > 0
  qu <- obj$family$getQu()
  
  # Now branching for main computation
  if( is.null(x2) ) ################ ONE VARIABLE
  {
    if(length(x1) != n) stop("length(x1) != ncol(X)")
    
    if( is.null(bound )){ # Create bounds OR ...
      nbin1 <- nbin[1]
      bound <- seq(min(x1), max(x1), length.out = nbin1 + 1)
    } else { # ... use those already given
      nbin1 <- length(bound-1)
    }
    
    # For each bin: count number of responses that are smaller than the fitted quantile
    indx <- as.factor( .bincode(x1, bound, T, T) )       # Attribute data to bins
    levels(indx) <- 1:nbin1
    bsize <- as.vector( table(indx) )                    # Count number of data in each bin
    indx <- as.integer(indx)
    bins <- numeric(nbin1)             
    for(ii in 1:nbin1){ bins[ii] <- sum(res[indx==ii]) } # Count number of 1s in each bin
    
    # Remove empty bins
    while( any(bsize == 0) )
    {
      bad <- which(bsize == 0)[1]
      bsize <- bsize[-bad]
      bins <- bins[-bad]
      bound <- bound[ -min(bad+1, nbin1) ]
      if(bad<nbin1){
        indx[indx>bad] <- indx[indx>bad]-1
      }
      nbin1 <- nbin1 - 1
    }
    
    ub <- qbinom(lev/2, bsize, qu, lower.tail = FALSE) / bsize
    lb <- qbinom(lev/2, bsize, qu) / bsize
    
    #svpar <- par(no.readonly = TRUE) 
    par(mar = c(5.1, 4.6, 4.1, 2.1))
    x <- sort(x1)
    tmp <- rep(bins/bsize, bsize)
    plot(x, tmp, ylim = range(ub, lb, tmp), type = 'l', col = "white", ylab = expression(hat(P)(y<hat(mu))), ...)
    abline(h = qu, lty = 2)
    for(ii in 1:nbin1){
      prop <- bins[ii]/bsize[ii]
      lines(c(bound[ii], bound[ii+1]), rep(ub[ii], 2), col = "grey")
      points(c(bound[ii], bound[ii+1]), rep(ub[ii], 2), pch = 3, lwd = 0.2)
      lines(c(bound[ii], bound[ii+1]), rep(lb[ii], 2), col = "grey")
      points(c(bound[ii], bound[ii+1]), rep(lb[ii], 2), pch = 3, lwd = 0.2)
      points((bound[ii] + bound[ii+1])/2, prop, 
             col = ifelse(prop < ub[ii] && prop > lb[ii], 1, 2))
      rug(x1[indx==ii], col = ifelse(ii%%2, 3, 4))
    }
    #par(svpar) # Warning: this causes problems with shiny
  } else { ################ ... TWO VARIABLES  
    
    if(length(x1) != n) stop("length(x1) != ncol(X)")
    if(length(x2) != n) stop("length(x2) != ncol(X)")
    
    if( is.null(bound) ){ # Bounds created OR ... 
      if( length(nbin) != 2 ){ stop("In the 2D case, nbin should be a vector of length 2") }
      bound1 <- seq(min(x1), max(x1), length.out = nbin[1] + 1)
      bound2 <- seq(min(x2), max(x2), length.out = nbin[2] + 1)
    } else { # ... already provided
      if( length(bound) != 2 ){ stop("In the 2D case, bound a list of two numeric vectors") }
      bound1 <- bound[[1]]
      bound2 <- bound[[2]]
      nbin <- c(length(bound1)-1, length(bound2)-1)
    }
    
    # For each bin: count number of responses that are smaller than the fitted quantile
    bins <- bsize <- matrix(0, nbin[2], nbin[1])
    for(ir in 1:nbin[2]){ # From bottom upward
      inRow <- which( x2 >= bound2[ir] & x2 <= bound2[ir+1] )
      if(length(inRow)){
        x1In <- x1[inRow]
        resIn <- res[inRow]
        for(ic in 1:nbin[1]){ # From left to right
          tmp <- which( x1In >= bound1[ic] & x1In <= bound1[ic+1] )
          bsize[ir, ic] <- length(tmp)   
          bins[ir, ic] <- sum( resIn[tmp] )
        }}}
    
    # Plot!
    plot(x1, x2, pch = ".", col = "white", ylim = range(bound2), xlim = range(bound1), 
         main = expression(hat(P)(y<hat(mu))), ...)
    for(tmp in bound2){ segments(x0 = min(bound1), x1 = max(bound1), y0 = tmp, y1 = tmp) }
    for(tmp in bound1){ segments(y0 = min(bound2), y1 = max(bound2), x0 = tmp, x1 = tmp) }
    
    prop <- bins / bsize 
    for(ir in 1:nbin[2]){
      for(ic in 1:nbin[1]){
        if( is.finite(prop[ir, ic]) )
        {
          pr <- prop[ir, ic]
          siz <- bsize[ir, ic]
          b <- qbinom(lev/2, siz, qu, lower.tail = (pr<=qu)) / siz
          
          sig <- if(pr>qu){ pr>b } else{ pr<b }
          polygon(x = c(bound1[ic], bound1[ic], bound1[ic+1], bound1[ic+1], bound1[ic]),
                  y = c(bound2[ir], bound2[ir+1], bound2[ir+1], bound2[ir], bound2[ir]),
                  col = rgb(pr<qu, pr>qu, 0, alpha = ifelse(sig, 0.8, 0.3)))
          
          text(x = (bound1[ic]+bound1[ic+1])/2, y = (bound2[ir]+bound2[ir+1])/2, 
               paste(round(pr, 3), ifelse(sig, "*", ''), sep = ""))
        }
      }}
    
    rug(x1, side = 1)
    rug(x2, side = 2)
    if(scatter){ points(x1, x2, pch = ".") }
  }
  
   return( invisible(NULL) )
}