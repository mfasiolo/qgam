# Simulated random variables from log-F density
rlf <- function(n, mu, tau, sig, lam)
{
  x <- pmax( rgamma(n, lam*tau, 1), 1e-300)
  y <- pmax( rgamma(n, lam*(1-tau), 1), 1e-300)
  
  out <- log(x) - log(y)
  
  out <- lam * sig * out + mu
  
  return( out )
}

# Test
# mu <- 3
# tau <- 0.2
# sig <- 1
# lam <- 0.1
# s <- rlf(1e4, mu, tau, sig, lam)
# truth <- dlf(sort(s), tau, mu, sig, lam)
# plot(sort(s), truth, type = 'l')
# lines(density(s, bw = 0.1), col = 2)
# 
# mu <- 3
# tau <- 0.1
# sig <- 0.3
# lam <- 0.4
# s <- rlf(1e6, mu, tau, sig, lam)
# 
# mean(s)
# sig * lam * ( digamma(lam*tau) - digamma(lam*(1-tau)) ) + mu
# 
# var(s)
# sig^2 * lam^2 * ( trigamma(lam*tau) + trigamma(lam*(1-tau)) )