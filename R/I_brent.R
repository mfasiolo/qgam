
#*****************************************************************************
# Shamelessly copied (and translated) from Burkardt's 
# website https://people.sc.fsu.edu/~jburkardt/m_src/brent/local_min.m
#
## LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
#
#  Discussion:
#
#    The method used is a combination of golden section search and
#    successive parabolic interpolation.  Convergence is never much slower
#    than that for a Fibonacci search.  If F has a continuous second
#    derivative which is positive at the minimum (which is not at A or
#    B), then convergence is superlinear, and usually of the order of
#    about 1.324....
#
#    The values EPSI and T define a tolerance TOL = EPSI * abs ( X ) + T.
#    F is never evaluated at two points closer than TOL.
#
#    If F is a unimodal function and the computed values of F are always
#    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
#    LOCAL_MIN approximates the abscissa of the global minimum of F on the
#    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
#
#    If F is not unimodal, then LOCAL_MIN may approximate a local, but
#    perhaps non-global, minimum to the same accuracy.
#
#    Thanks to Jonathan Eggleston for pointing out a correction to the 
#    golden section step, 01 July 2013.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    13 April 2008
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    MATLAB version by John Burkardt.
#    R vesion from Matteo Fasiolo
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization Without Derivatives,
#    Dover, 2002,
#    ISBN: 0-486-41998-3,
#    LC: QA402.5.B74.
#
#  Parameters:
#
#    Input, real A, B, the endpoints of the interval.
#
#    Input, real EPSI, a positive relative error tolerance.
#    EPSI should be no smaller than twice the relative machine precision,
#    and preferably not much less than the square root of the relative
#    machine precision.
#
#    Input, real T, a positive absolute error tolerance.
#
#    Input, function value = F ( x ), the name of a user-supplied
#    function whose local minimum is being sought.
#
#    Output, real X, the estimated value of an abscissa
#    for which F attains a local minimum value in [A,B].
#
#    Output, real FX, the value F(X).
#

.brent <- function(brac, f, mObj, bObj, init, pMat, qu, varHat, cluster, t = .Machine$double.eps^0.25, ...)
{
  brac <- sort(brac)
  a <- brac[1]
  b <- brac[2]
  
  # Relative tolerance, as in ?optimize. No need to touch it, I think.
  epsi = sqrt(.Machine$double.eps)
  
  # cc is the square of the inverse of the golden ratio.
  cc = 0.5 * ( 3.0 - sqrt(5.0) )
  
  sa = a
  sb = b
  x = sa + cc * ( b - a )
  w = x
  v = w
  e = 0.0
  feval = f(lsig = x, mObj = mObj, bObj = bObj, initM = init[["initM"]], initB = init[["initB"]], 
            pMat = pMat, qu = qu, varHat = varHat, cluster = cluster, ...)
  fx = feval$loss
  fw = fx
  fv = fw
  
  # Storing all evaluations points and function values
  jj <- 1
  store <- list()
  store[[jj]] <- list("x" = x, "f" = fx, "initM" = feval[["initM"]], "initB" = feval[["initB"]])
  jj <- jj + 1
  
  while( TRUE ) {
    
    m = 0.5 * ( sa + sb )
    tol = epsi * abs ( x ) + t
    t2 = 2.0 * tol
    
    #  Check the stopping criterion.
    if( abs(x-m) <= (t2 - 0.5 * (sb-sa)) ) { break }
    
    #  Fit a parabola.
    r = 0.0
    q = r
    p = q
    
    if ( tol < abs(e) ) {
      
      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0 * ( q - r )
      
      if( 0.0 < q ) { p = - p }
      
      q = abs ( q )
      r = e
      e = d
    }
    
    #  Take the parabolic interpolation step OR ...
    if ( (abs(p) < abs(0.5 * q * r)) && 
         (q * ( sa - x )) < p &&
         (p < q * ( sb - x )) ) {
      
      d = p / q
      u = x + d
      
      #  F must not be evaluated too close to A or B.
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 ) {
        
        if ( x < m ) { d = tol } else { d = - tol }
        
      }
      
    } else { # ...  a golden-section step. 
      
      if ( x < m ){ e = sb - x } else { e = sa - x }
      
      d = cc * e
    }
    
    #  F must not be evaluated too close to X.
    if ( tol <= abs( d ) ){
      u = x + d
    } else {
      if ( 0.0 < d ) { u = x + tol } else { u = x - tol }
    }
    
    init <- store[[ which.min(abs(u - sapply(store, "[[", "x"))) ]]
    feval = f(lsig = u, mObj = mObj, bObj = bObj, initM = init[["initM"]], initB = init[["initB"]], 
              pMat = pMat, qu = qu, varHat = varHat, cluster = cluster, ...)
    fu = feval$loss
    store[[jj]] <- list("x" = u, "f" = fu, "initM" = feval[["initM"]], "initB" = feval[["initB"]])
    jj <- jj + 1
    
    #  Update A, B, V, W, and X.
    if ( fu <= fx ){
      
      if ( u < x ) { sb = x } else { sa = x }
      
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      
    } else { 
      
      if ( u < x ) { sa = u } else { sb = u }
      
      if ( (fu <= fw) || (w == x) ) {
        v = w
        fv = fw
        w = u
        fw = fu
      } else {
        if ( (fu <= fv) || (v == x) || (v == w) ){
          v = u
          fv = fu
        } 
      }
    }
  }
  
  store <- rbind( sapply(store, "[[", "x"), sapply(store, "[[", "f") )
  
  return( list("minimum" = x, "objective" = fx, "store" = store) )
}

###################
######### TEST
###################

# # Test 1
# f <- function (x, k) (x - k)^2
# xmin <- optimize(f, c(0, 1), tol = 0.0001, k = 1/3)
# xmin
# 
# qgam:::.brent(brac = c(0, 1), f = f, t = 1e-4, k = 1/3)
# 
# # Test 2
# f  <- function(x) ifelse(x > -1, ifelse(x < 4, exp(-1/abs(x - 1)), 10), 10)
# fp <- function(x) { print(x); f(x) }
# 
# plot(f, -2,5, ylim = 0:1, col = 2)
# optimize(fp, c(-7, 20))   # ok
# qgam:::.brent(c(-7, 20), f = f)

