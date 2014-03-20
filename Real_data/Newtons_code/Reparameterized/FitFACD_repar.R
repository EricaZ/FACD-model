FitFACD_repar <- function(x, id.str = "test_run", maxit = 500, 
                    init.param = c(2, 1, 0.2, 0.6)) {
  
  finalparamfile <- paste(id.str, "_record.csv", sep = "") 
  
  num.n <- length(x)
  itr = 0 
  param <- init.param
  
  # Optimization
  repeat {
    
    r <- param[1]
    s <- param[2]
    a <- param[3]
    b <- param[4]
    
    itr = itr + 1
    
    # Compute conditional durations
    dur <- ComputeCondDur_repar(x = x, a = a, b = b)
    
    # Compute errors
    err <- x / dur
    
    # Compute gradient of conditional durations
    gradient <- ComputeGrad_repar(b = b, x = x, dur = dur)
    
    # Compute c1, c2
    c <- ComputeC_repar(r = r, s = s, err = err)
    
    # Compute information matrix and score function
    infoscore <- ComputeIS_repar(c = c, dur = dur, gradient = gradient)
    
    # Compute parameters for next iteration
    inv.infomat <- solve(infoscore$infomat)
    adj <- inv.infomat %*% infoscore$score 
    new.param <- param + adj
    
    cat("\n", "Iteration", itr, "\n")
    
    print(new.param)
    
    for(i in 1: 4) {
      cat(new.param[i], " ", sep = ",", file = finalparamfile, append = TRUE)
    }  
    cat(" ", infoscore$score, "\n", sep = ",", file = finalparamfile, append = TRUE)
    
    # Impose constraints on new parameters
    ## Lower bound = 0 for s, a, b
    for (i in 2: 4) {
      if (new.param[i] < 0) {
        new.param[i] <- init.param[i]
      }
    }
    
    ## Lower bound = 1 for r, for frechet distrib
    if(new.param[1] < 1) {
      new.param[1] <- init.param[1]
    }
    
    ## Upper bound = 1 for b
    if (new.param[4] > 1) {
      new.param[4] <- max(new.param[4] / 1.3, init.param[4])
    }
    
    # =======================================================================
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
  
  cat("\n")
  
  cat("FACD_repar \n", sep = ",", file = finalparamfile, append = TRUE)
  cat(param, " ", param.var, " ", itr, "\n", sep = ",", file = finalparamfile, append = TRUE)
  
  return(list(res = err, param = param))
}