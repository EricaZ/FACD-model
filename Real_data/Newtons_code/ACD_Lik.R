ACD_Lik <- function(param, x, p, q, distrib="frechet", lik = TRUE){
#   browser()
  r <- param[1]
  w <- param[2]
  a.vec <- param[3:(2+p)]
  b.vec <- param[(3+p):(2+p+q)]
  
  # Initialization
  num.n <- length(x)
  xbar <- mean(x)
#   browser()
  x.vec <- c(rep(xbar, p), x)
  dur.vec <- c(rep(xbar, q), rep(0, num.n))
  
  # Compute conditional durations
  for (i in 1: num.n) {
    apart <- sum(a.vec * x.vec[(p+i-1):i])
    bpart <- sum(b.vec * dur.vec[(q+i-1):i])
    dur.vec[q+i] <- w + apart + bpart
  }
  
  dur <- dur.vec[(q+1): (q+num.n)]  # an n.num vector
  err <- x / dur
    
  if (lik == TRUE) {
    if (distrib == "frechet") {
      a <- gamma(1 - 1 / r)
      cga <- a^(-r)
      loglik <- r * log(dur) - cga * err^(-r) - (1 + r) * log(x) + log(r * cga)
    }
    print(param)
    return(sum(loglik))  
  } else {
    return(list(res = err, param = param))
  }
}
