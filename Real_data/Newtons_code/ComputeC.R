ComputeC <- function(r, err, distrib) {
  # Computes c1 and c2 for estimation of variance-covariance matrix.
  #
  # Args:
  #   r: shape parameter of frechet distribution
  #   err: 
  #   distrib: 
  #
  # Returns:
  #   list(c1, c2)
  #     c1 and c2 are vectors
  
  log.err <- log(err)
  
  if (distrib == "frechet") {
    a <- gamma(1 - 1 / r)
    cga <- a^(-r)
    dr.cga <- -digamma(1 - 1 / r) * a^(-r + 1) * r^(-2) * log(a)
    power.err <- err^(-r)
    
    c1 <- r * (1 - cga * power.err)
    c2 <- cga * power.err * log.err - log.err - dr.cga * power.err + 1 / r + dr.cga / cga  
  }
  
  if (distrib == "weibull") {
    a <- gamma(1 + 1 / r)
    cga <- a^r
    dr.cga <- -digamma(1 + 1 / r) * a^(r + 1) * r^(-2) * log(a)
    power.err <- err^r
    
    c1 <- -r * (1 - cga * power.err)
    c2 <- -cga * power.err * log.err + log.err - dr.cga * power.err + 1 / r + dr.cga / cga 
  }

  if (distrib == "modifrechet") {
    p <- 1 - 1 / r
    g <- gamma(p)
    q <- 9.21^(-1 / r)  # 9.21 = -log(delta) = m if delta = 0.0001
    bga <- 1 / (g - q)
    aga <- -q * bga
    power.err <- ((err - aga) / bga)^(-r)  
    c1 <- (err + aga / r - err * power.err) * r / (err - aga)
    
    dr.aga <- g*q/(r*q-r*g)^2*(digamma(p)-log(9.21))
    dr.bga <- (r*q-r*g)^(-2)*(q*log(9.21)-g*digamma(p))
    c2 <- 1 / r + dr.aga / (err - aga) + (r * dr.bga / bga - log((err - aga) / bga) + r * dr.aga / (err - aga)) / (1 - power.err)
  }
  
  return(list(c1 = c1, c2 = c2))
}