FitACD_dual <- function(x, id.str = "test_run", maxit = 500,
                        gamma.vec = seq(from = 1.1, to = 1.5, by = 0.1),
                        init.param = c(0.2, 0.2, 0.6),
                        distrib = "frechet") {
  # Fits the ACD model
  #
  # Args:
  #   x: sample
  #   id.str: character string, identifier of the run, for naming of the csv output files
  #   maxit: maximum iteration number
  #   init.param: a vector of initial parameters (w, a, b)
  #   distrib: "frechet" or "modifrechet"
  
#   traceparamfile <- paste(id.str, "_trace.csv", sep = "") 
#   finalparamfile <- paste(id.str, "_record.csv", sep = "")  
  optparamfile <- paste(id.str, "_record.csv", sep = "") 
  
  num.ga <- length(gamma.vec)
  max.loglik <- 0  # initialize the max(loglik) among estimations with different gamma's
  
  for (j in 1: num.ga) {
     
    itr = 0 
    r <- gamma.vec[j]
    param <- c(r, init.param)
    
    # Optimization
    repeat {
      
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
      
      # Compute c1
      c1 <- ComputeC_dual(r = r, err = err, distrib = distrib)

      # Compute information matrix and score function
      infoscore <- ComputeIS_dual(c1 = c1, dur = dur, gradient = gradient)
      
      # Compute parameters for next iteration
      inv.infomat <- solve(infoscore$infomat)
      adj <- inv.infomat %*% infoscore$score 
      new.param <- param + c(0, adj)
      
      cat("\n", "Iteration", itr, "\n")
      
      print(new.param)
      
#       for(i in 1: 4) {
#         cat(new.param[i], " ", sep = ",", file = traceparamfile, append = TRUE)
#       }
#       
#       cat(" ", infoscore$score, "\n", sep = ",", file = traceparamfile, append = TRUE)
      
      # ===========================================================================================
      # Impose constraints on new parameters
      ## Lower bound = 0 for w, a, b
      for (i in 2: 4) {
        if (new.param[i] < 0) {
          new.param[i] <- init.param[i-1]
        }
      }
      
      ## Upper bound = 1 for a and b
      for (i in 3: 4) {
        if (new.param[i] > 1) {
          new.param[i] <- init.param[i-1]
        }
      }
      
      # =============================================================================================
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
    
#     for(i in 1: 4) {
#       cat(param[i], " ", sep = ",", file = finalparamfile, append = TRUE)
#     }
    
    # Compute the log likelihood of the estimated model
    loglik <- ComputeLoglik_dual(r = r, x = x, err = err, dur = dur, distrib = distrib)

#     cat(" ", infoscore$score, " ", loglik, " ", itr, "\n", sep = ",", file = finalparamfile, append = TRUE) 
      
    # Update max(loglik)
    if ((j == 1)|(max.loglik < loglik)) {
      max.loglik <- loglik
      opt.param <- param
    } 
  }
  
  # ==================== end of outer for-loop =======================
  
  cat("mFACD", "\n", sep = ",", file = optparamfile, append = TRUE)

  for(i in 1: 4) {
    cat(opt.param[i], " ", sep = ",", file = optparamfile, append = TRUE)
  }
  
  cat(" ", max.loglik, "\n", sep = ",", file = optparamfile, append = TRUE)
  
  return(list(res = err, param = opt.param))
  
}




