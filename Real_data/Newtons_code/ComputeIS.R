ComputeIS <- function(IS2 = FALSE, c, dur, gradient) {
  # Computes the information matrix and the score function.
  #
  # Args:
  #   c:
  #   dur: 
  #   gradient:
  #   IS2: TRUE if output the information matrix for (w,a,b) and the scaled gradiant matrix
  #
  # Returns:
  #   list(i.mat, s.vec)
  #     infomat: information matrix
  #     score: a vector of the score function
  
  c1 <- c$c1
  c2 <- c$c2
  
  k1 <- mean(c1^2)
  k2 <- mean(c2^2)
  k3 <- mean(c1 * c2)
  num.n <- length(dur)
  scaled.grad <- matrix(nrow = 3, ncol = num.n)
  
  A <- matrix(0, nrow = 3, ncol = 1)
  
  for (i in 1: num.n) {
    scaled.grad[, i] <- gradient[, i] / dur[i]
    A <- A + scaled.grad[, i]
  }
  
  A <- A / num.n
  B <- (scaled.grad %*% t(scaled.grad)) / num.n
  
  if (IS2 == FALSE) {
    infomat <- rbind(c(k2, k3*t(A)), cbind(k3*A, k1*B))
    score <- rbind(mean(c2), scaled.grad %*% c1 / num.n)
    return(list(infomat = infomat, score = score))
  } else {
    infomat1 <- k1 * B - k3^2 / k2 * A %*% t(A)
    return(list(infomat1 = infomat1, scaled.grad = scaled.grad))
  }
}