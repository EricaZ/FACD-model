FitWACD <- function(x, id.str = "test_run", maxit = 500, portmanteau = FALSE,
                   init.param = c(2, 0.1, 0.2, 0.6)) {
  # Fits the ACD model
  #
  # Args:
  #   x: sample
  #   id.str: character string, identifier of the run, for naming of the csv output files
  #   maxit: maximum iteration number
  #   init.param: a vector of initial parameters (r, w, a, b)
  
#   traceparamfile <- paste(id.str, "_trace.csv", sep = "") 
  finalparamfile <- paste(id.str, "_record.csv", sep = "") 
  
  num.n <- length(x)
  itr = 0
  param <- init.param
  
  # Optimization
  repeat {
    
    r <- param[1]
    w <- param[2]
    a <- param[3]
    b <- param[4]
    
    itr = itr + 1
    
    # Compute conditional durations
    dur <- ComputeCondDur(x = x, coeff = c(w, a, b))
    
    # Compute errors
    err <- x / dur
    
    # Compute gradient of conditional durations
    gradient <- ComputeGrad(b = b, x = x, dur = dur)

    # Compute c1 and c2
    c <- ComputeC(r = r, err = err, distrib = "weibull")
    
    # Compute information matrix and score function
    infoscore <- ComputeIS(c = c, dur = dur, gradient = gradient)
    
    # Compute parameters for next iteration
    inv.infomat <- solve(infoscore$infomat)
    adj <- inv.infomat %*% infoscore$score 
    new.param <- param + adj
    
    cat("\n", "Iteration", itr, "\n")
    
    print(new.param)
    
#     for(i in 1: 4) {
#       cat(new.param[i], " ", sep = ",", file = traceparamfile, append = TRUE)
#     }
#     cat(" ", infoscore$score, "\n", sep = ",", file = traceparamfile, append = TRUE)
    
    # ==========================================================================
    # Impose constraints on new parameters
    ## Lower bound = 0 for w, a.vec, b.vec
    for (i in 2: 4) {
      if (new.param[i] < 0) {
        new.param[i] <- init.param[i]
      }
    }
    
    ## Lower bound = 0 for r, for weibull distrib
    if(new.param[1] < 0) {
      new.param[1] <- init.param[1]
    }

    ## Upper bound = 1 for sum of b.vec
    if (new.param[4] > 1) {
      new.param[4] <- max(new.param[4] / 1.3, init.param[4])
    }

    # ===========================================================================
    # Check stopping condition
    if (max(abs(new.param - param)) < 0.0001) {
      break
    } else {
      if (itr == maxit) {
        break
      }
    }
    
    # Renew parameters for next iteration
    param <- new.param
  }
  
  # Compute the estimated variance of the parameters
  param.var <- diag(inv.infomat) / num.n

  # Compute the residual ACFs
  res.acf <- vector(length = 6)
  for (k in 1: 6) {
    res.acf[k] <- sum((err[(k+1):num.n]-1)*(err[1:(num.n-k)]-1)) / sum((err-1)^2)  
  }

  # Compute the estimated variance of the residual ACFs
  Omega.input <- ComputeIS(c = c, dur = dur, gradient = gradient, IS2 = TRUE)
  Omega.mat <- ComputeOmega(input = Omega.input, err = err)
  acf.var <- diag(Omega.mat) / num.n

  # Compute the portmanteau test statistic
  if (portmanteau == TRUE) {
    Q <- num.n * t(res.acf) %*% solve(Omega.mat) %*% res.acf
  }

  cat("\n")

  cat("WACD", "\n", sep = ",", file = finalparamfile, append = TRUE)
  cat(param, " ", param.var, " ", res.acf, " ", acf.var, " ", sep = ",", file = finalparamfile, append = TRUE)

  if (portmanteau == TRUE) {
    cat(" ", Q, " ", sep = ",", file = finalparamfile, append = TRUE) 
  }

  cat(" ", itr, "\n", sep = ",", file = finalparamfile, append = TRUE) 

  return(list(res = err, param = param))
}