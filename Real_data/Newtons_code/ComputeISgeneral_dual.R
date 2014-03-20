ComputeISgeneral_dual <- function(c1, dur, gradient) {
  # Computes the information matrix and the score function, when gamma is fixed 
  # (i.e. only estimate 3 parameters)
  #
  # Args:
  #   c1:
  #   dur: 
  #   gradient:
  #
  # Returns:
  #   list(i.mat, s.vec)
  #     infomat: information matrix
  #     score: a vector of the score function
  
  k1 <- mean(c1^2)
  num.coeff <- nrow(gradient)
  num.n <- length(dur)
  scaled.grad <- matrix(nrow = num.coeff, ncol = num.n)
  
  for (i in 1: num.n) {
    scaled.grad[, i] <- gradient[, i] / dur[i]
  }

  B <- (scaled.grad %*% t(scaled.grad)) / num.n
  
  infomat <- k1 * B
  score <- scaled.grad %*% c1 / num.n
  
  return(list(infomat = infomat, score = score))
}