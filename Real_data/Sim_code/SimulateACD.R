library(evd)

SimulateACD <- function(param, distrib, offset = 200, num.n = 1000, num.rep = 1000) {
  # Simulates num.n repliacations of ACD time series using ACD model.
  # 
  # Args: 
  #   param: a list of true parameters (r, w, a, b) or (w, a, b)
  #   distrib:
  #   offset:
  #   num.n:
  #   num.rep: 
  #
  # Returns:
  #   a matrix (num.n rows, num.rep columns) containing the simulated time series.   
  
  # Initialization 
  if (distrib == "exp") {
    w <- param[1]
    a <- param[2]
    b <- param[3]
  }
  
  if ((distrib == "weibull")|(distrib == "frechet")) {
    r <- param[1]
    w <- param[2]
    a <- param[3]
    b <- param[4]
  }
  
  num.rn <- (offset + num.n) * num.rep
  
  # Generate random numbers 
  if (distrib == "frechet") {
    rand.vec <- rfrechet(num.rn, shape = r, scale = 1 / gamma(1 - 1 / r)) 
  }
  
  if (distrib == "exp") {
    rand.vec <- rexp(num.rn, rate = 1) 
  }

  if (distrib == "weibull") {
    rand.vec <- rweibull(num.rn, r, scale = 1 / gamma(1 + 1 / r))
  }
  
#   for (i in 1: num.rn)
#     cat(rand.vec[i], "\n", file = "rand.csv",sep = ",", append=TRUE)
  
  rand.mat <- matrix(rand.vec, nrow = offset + num.n, ncol = num.rep)
  
  # Initialize variables
  x.vec <- rep(0, 1 + offset + num.n)
  dur.vec <- rep(0, 1 + offset + num.n)  
  x.mat <- matrix(nrow = num.n, ncol = num.rep)
  
  # Compute durations
  for (j in 1: num.rep) {
    for (i in 1: (offset + num.n)) {
      dur.vec[1+i] <- w + a * x.vec[i] + b * dur.vec[i]
      x.vec[1+i] <- dur.vec[1+i] *  rand.mat[i, j] 
    }
    x.mat[, j] <- x.vec[(offset+2): (1+offset+num.n)]
  }
  
  return(x.mat)  
}

