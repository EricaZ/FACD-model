ComputeIS_repar <- function(c, dur, gradient) {
  # Computes the information matrix and the score function.
  #
  # Args:
  #   c:
  #   dur: 
  #   gradient:
  #
  # Returns:
  #   list(i.mat, s.vec)
  #     infomat: information matrix
  #     score: a vector of the score function
  
  c1 <- c$c1
  c2 <- c$c2  # 2 by n matrix
  
  num.n <- length(dur)
  k1 <- mean(c1^2)
  k2 <- c2 %*% t(c2) / num.n  # 2 by 2 matrix
  k3 <- c2 %*% c1 / num.n  # 2 by 1 vector

  scaled.grad <- matrix(nrow = 2, ncol = num.n)
  
  A <- matrix(0, nrow = 2, ncol = 1)
  
  for (i in 1: num.n) {
    scaled.grad[, i] <- gradient[, i] / dur[i]
    A <- A + scaled.grad[, i]
  }
  
  A <- A / num.n
  B <- (scaled.grad %*% t(scaled.grad)) / num.n
  
  infomat <- rbind(cbind(k2, k3%*%t(A)), cbind(t(k3%*%t(A)), k1*B))
  score <- rbind(mean(c2[1,]), mean(c2[2,]), scaled.grad %*% c1 / num.n)
  
  return(list(infomat = infomat, score = score))
}