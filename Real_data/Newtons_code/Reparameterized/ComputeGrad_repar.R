ComputeGrad_repar <- function(b, x, dur) {
  # Computes the gradiants of conditional durations.
  #
  # Args:
  #   b:
  #   x:
  #   dur:
  #
  # Returns:
  #   a matrix (p+q rows, num.n columns) of the gradiant
  
  xbar <- mean(x)
  num.n <- length(x)
  x.vec <- c(xbar, x)
  dur.vec <- c(xbar, dur)
  grad.mat <- matrix(0, nrow = 2, ncol = 1 + num.n)
  
  for(i in 1: num.n) {
    grad.mat[, 1+i] = c(x.vec[i], dur.vec[i]) + grad.mat[, i] * b   
  }

  return(grad.mat[, 2:(1+num.n)])
}