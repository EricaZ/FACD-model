ComputeCondDur <- function(x, coeff) {
  # Computes the conditional durations
  #
  # Args:
  #   x:
  #   coeff: vector (w, a, b)
  #
  # Returns:
  #   a vector of conditonal durations
  
  # Initialization
  w <- coeff[1]
  a <- coeff[2]
  b <- coeff[3]
  num.n <- length(x)
  xbar <- mean(x)
  x.vec <- c(xbar, x)
  dur.vec <- c(xbar, rep(0, num.n))
  
  # Compute conditional durations
  for (i in 1: num.n) {
    dur.vec[1+i] <- w + a * x.vec[i] + b * dur.vec[i]
  }
  
  return(dur.vec[2: (1+num.n)])
}