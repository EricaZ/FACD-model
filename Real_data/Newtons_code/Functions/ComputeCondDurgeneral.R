ComputeCondDurgeneral <- function(x, p = 1, q = 1, 
                           w = 0.1, a.vec = 0.1, b.vec = 0.7) {
  # Computes the conditional durations
  #
  # Args:
  #   x:
  #   p:
  #   q:
  #   a.vec:
  #   b.vec:
  #
  # Returns:
  #   a vector of conditonal durations
  
  # Initialization
  num.n <- length(x)
  xbar <- mean(x)
  x.vec <- c(rep(xbar, p), x)
  dur.vec <- c(rep(xbar, q), rep(0, num.n))
  
  # Compute conditional durations
  for (i in 1: num.n) {
    apart <- sum(a.vec * x.vec[(p+i-1):i])
    bpart <- sum(b.vec * dur.vec[(q+i-1):i])
    dur.vec[q+i] <- w + apart + bpart
  }
  
  dur <- dur.vec[(q+1): (q+num.n)]  # an n.num vector
  
  return(dur)
}