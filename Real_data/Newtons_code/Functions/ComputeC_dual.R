ComputeC_dual <- function(r, err, distrib) {
  # Computes c1 for estimation of 3 by 3 variance-covariance matrix.
  #
  # Args:
  #   r: shape parameter of frechet distribution
  #   err: 
  #   distrib: 
  #
  # Returns:
  #   c1: numeric
  
  log.err <- log(err)
  
  if (distrib == "frechet") {
    a <- gamma(1 - 1 / r)
    cga <- a^(-r)
    dr.cga <- -digamma(1 - 1 / r) * a^(-r + 1) * r^(-2) * log(a)
    power.err <- err^(-r)    
    c1 <- r * (1 - cga * power.err) 
  }
  
  if (distrib == "modifrechet") {
    p <- 1 - 1 / r
    g <- gamma(p)
    q <- 9.21^(-1 / r)  # 9.21 = -log(delta) if delta = 0.0001
    bga <- 1 / (g - q)
    aga <- -q * bga
    power.err <- ((err - aga) / bga)^(-r)  
    c1 <- (err + aga / r - err * power.err) * r / (err - aga)
  }

  return(c1)
}