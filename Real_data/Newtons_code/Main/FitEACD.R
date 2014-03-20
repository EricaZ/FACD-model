FitEACD <- function(x, id.str, maxit = 500, init.param = c(0.1, 0.2, 0.6)) {
  # Fits the exponential ACD model
  #
  # Args:
  #   x: sample
  #   id.str: character string, identifier of the run, for naming of the csv output files
  #   maxit: maximum iteration number
  #   init.param: a vector of initial parameters (w, a, b)
  
#  traceparamfile <- paste(id.str, "_trace.csv", sep = "") 
  finalparamfile <- paste(id.str, "_record.csv", sep = "") 
  
  num.n <- length(x)
  itr = 0
  param <- init.param 
  
  # Optimization
  repeat {
    itr = itr + 1
    
    # Compute conditional durations
    dur <- ComputeCondDur(x = x, coeff = param)  
    
    # Compute errors
    err <- x / dur
    
    # Compute gradient of conditional durations
    gradient <- ComputeGrad(b = param[3], x = x, dur = dur)

    # Compute information matrix and score function
    k <- mean((err - 1)^2)
    scaled.grad <- matrix(nrow = 3, ncol = num.n)
    
    for (i in 1: num.n) {
      scaled.grad[, i] <- gradient[, i] / dur[i]
    }

    B <- (scaled.grad %*% t(scaled.grad)) / num.n

    infomat <- k * B
    score <- scaled.grad %*% (err - 1) / num.n

    # Compute parameters for next iteration
    inv.infomat <- solve(infomat)
    adj <- inv.infomat %*% score 
    new.param <- param + adj
    
    cat("Iteration", itr, "\n")  

    print(new.param)
    
#     for(i in 1: 3) {
#       cat(new.param[i], " ", sep = ",", file = traceparamfile, append = TRUE)
#     }   
#     cat(" ", score, "\n", sep = ",", file = traceparamfile, append = TRUE)
    
    # =====================================================================
    # Impose constraints on new parameters
    ## Lower bound = 0 for w, a, b
    for (i in 1: 3) {
      if (new.param[i] < 0) {
        new.param[i] <- init.param[i]
      }
    }

    ## Upper bound = 1 for b
    if (new.param[3] > 1) {
      new.param[3] <- max(new.param[3] / 1.3, init.param[3])  
    }

    # ======================================================================
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
  
  cat("\n")
  
  cat("EACD", "\n", sep = ",", file = finalparamfile, append = TRUE)
  cat(param, " ", score, " ", itr, "\n", sep = ",", file = finalparamfile, append = TRUE) 

  return(list(res = err, param = param))
}