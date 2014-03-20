library(evd)

SimulateACDgeneral <- function(param = list(w = 0.2, a.vec = 0.2, b.vec = 0.6, r = 2), 
                        offset = 200, num.n = 1000, num.rep = 1000, distrib = "frechet") {
  # Simulates num.n repliacations of ACD time series using ACD model.
  # 
  # Args: 
  #   param: a list of true parameters
  #          list(w, avec, bvec, r)
  #            w: constant term
  #            a.vec: a vector containing coefficients of lagged x's
  #            b.vec: a vector containing coefficients of lagged conditional durations
  #            r: shape parameter of error distribution, if weibull and frechet distribution
  #   offset:
  #   num.n:
  #   num.rep:
  #   distrib: 
  #
  # Returns:
  #   a matrix (num.n rows, num.rep columns) containing the simulated time series.   
  
  p <- length(param$a.vec)
  q <- length(param$b.vec)
  num.rn <- (offset + num.n) * num.rep
  
  # Generate random numbers 
  if (distrib == "frechet") {
    rand.vec <- rfrechet(num.rn, shape = param$r, scale = 1 / gamma(1 - 1 / param$r)) 
  }
  
  if (distrib == "exp") {
    rand.vec <- rexp(num.rn, rate = 1) 
  }

  if (distrib == "weibull") {
    rand.vec <- rweibull(num.rn, param$r, scale = 1 / gamma(1 + 1 / param$r))
  }
 
  rand.mat <- matrix(rand.vec, nrow = offset + num.n, ncol = num.rep)
  
  # Initialize variables
  x.vec <- rep(0, p + offset + num.n)
  dur.vec <- rep(0, q + offset + num.n)  
  x.mat <- matrix(nrow = num.n, ncol = num.rep)
  
  # Compute durations
  for (j in 1: num.rep) {
    for (i in 1: (offset + num.n)) {
      apart <- sum(param$a.vec * x.vec[(p+i-1):i])
      bpart <- sum(param$b.vec * dur.vec[(q+i-1):i])
      dur.vec[q+i] <- param$w + apart + bpart
      x.vec[p+i] <- dur.vec[q+i] *  rand.mat[i, j] 
    }
    x.mat[, j] <- x.vec[(p+offset+1): (p+offset+num.n)]
#   write.csv(rand.mat[, j], "random.csv")
  }
  
  return(x.mat)  
}

