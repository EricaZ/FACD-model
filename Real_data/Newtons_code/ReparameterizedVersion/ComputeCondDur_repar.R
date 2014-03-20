ComputeCondDur_repar <- function(x, a, b) {

  # Initialization
  num.n <- length(x)
  xbar <- mean(x)
  x.vec <- c(xbar, x)
  dur.vec <- c(xbar, rep(0, num.n))
  
  # Compute conditional durations
  for (i in 1: num.n) {
    dur.vec[1+i] <- 1 + a * x.vec[i] + b * dur.vec[i]
  }
  
  return(dur.vec[2: (1+num.n)])
}