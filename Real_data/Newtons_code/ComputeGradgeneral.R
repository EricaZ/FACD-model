ComputeGradgeneral <- function(b.vec, x, p, q, dur) {
  # Computes the gradiants of conditional durations.
  #
  # Args:
  #   b.vec:
  #   x:
  #   p:
  #   q:
  #   dur:
  #
  # Returns:
  #   a matrix (1+p+q rows, num.n columns) of the gradiant
  
  xbar <- mean(x)
  num.n <- length(x)
  x.vec <- c(rep(xbar, p), x)
  dur.vec <- c(rep(xbar, q), dur)
  grad.mat <- matrix(0, nrow = 1 + p + q, ncol = q + num.n)
  
  for(i in 1: num.n) {
    term1 <- as.matrix(c(1, x.vec[(p+i-1):i], dur.vec[(q+i-1):i]))
    term2 <- grad.mat[, (q+i-1):i] %*% as.matrix(b.vec)
    grad.mat[, q+i] = term1 + term2   
  }
  
  gradient <- grad.mat[, (q+1):(q+num.n)]
  
#   write.csv(gradient, "grad.csv")
  return(gradient)
}