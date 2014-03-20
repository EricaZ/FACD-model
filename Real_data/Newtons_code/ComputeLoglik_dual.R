ComputeLoglik_dual <- function(r, x, err, dur, distrib) {
  # Computes log likelihood.
  #
  # Args:
  #   r: shape parameter of the modified frechet distribution
  #   x: 
  #   err: 
  #   dur
  #   distrib: 
  #
  # Return:
  #   value of the log likelihood function
  
  if (distrib == "frechet") {
    a <- gamma(1 - 1 / r)
    cga <- a^(-r)
    loglik <- r * log(dur) - cga * err^(-r) - (1 + r) * log(x) + log(r * cga)
  }
  
  if (distrib == "modifrechet") {
    p <- 1 - 1 / r
    g <- gamma(p)
    q <- 9.21^(-1 / r)  # 9.21 = -log(delta) if delta = 0.0001
    bga <- 1 / (g - q)
    aga <- -q * bga
    loglik <- log(r) - log(bga) - (1 + r) * (log(err - aga) - log(bga)) - ((err - aga) / bga)^(-r) - log(dur)
  }
  
  return(sum(loglik))
}