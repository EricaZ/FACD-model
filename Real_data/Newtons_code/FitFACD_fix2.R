FitFACD_fix <- function(x, id.str = "test_run", maxit = 500, r = 2, portmanteau = FALSE,
                        init.param = c(0.1, 0.2, 0.6), distrib = "frechet") {
  # Fits the ACD model
  #
  # Args:
  #   x: sample
  #   id.str: character string, identifier of the run, for naming of the csv output files
  #   maxit: maximum iteration number
  #   init.param: a vector of initial parameters (w, a, b)
  
  finalparamfile <- paste(id.str, "_record.csv", sep = "") 
  
  num.n <- length(x)
  itr = 0
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
  
  # Compute log likelihood
  loglik <- ComputeLoglik_dual(r = r, x = x, err = err, dur = dur, distrib = distrib) 
  
  # Compute the estimated variance of the parameters
  param.var <- diag(inv.infomat) / num.n
  
  # Compute the residual ACFs
  res.acf <- vector(length = 6)
  for (k in 1: 6) {
    res.acf[k] <- sum((err[(k+1):num.n]-1)*(err[1:(num.n-k)]-1)) / sum((err-1)^2)  
  }
  
  res.acf2 <- vector(length = 12)
  for (k in 1: 12) {
    res.acf2[k] <- sum((err[(k+1):num.n]-1)*(err[1:(num.n-k)]-1)) / sum((err-1)^2)  
  }
  
  res.acf3 <- vector(length = 18)
  for (k in 1: 18) {
    res.acf3[k] <- sum((err[(k+1):num.n]-1)*(err[1:(num.n-k)]-1)) / sum((err-1)^2)  
  }
  
  
  # Compute the estimated variance of the residual ACFs
  c <- ComputeC(r = r, err = err, distrib = distrib)
  Omega.input <- ComputeIS(c = c, dur = dur, gradient = gradient, IS2 = TRUE)
  Omega.mat <- ComputeOmega(input = Omega.input, err = err)
  acf.var <- diag(Omega.mat) / num.n
  
  Omega.mat2 <- ComputeOmega2(input = Omega.input, err = err)
  acf.var2 <- diag(Omega.mat2) / num.n
  
  Omega.mat3 <- ComputeOmega3(input = Omega.input, err = err)
  acf.var3 <- diag(Omega.mat3) / num.n
  
  # Compute the estimated variance of the 3 parameters (except gamma)
  param.var2 <- diag(solve(Omega.input$infomat1)) / num.n
  
  # Compute the portmanteau test statistic
  if (portmanteau == TRUE) {
    Q <- num.n * t(res.acf) %*% solve(Omega.mat) %*% res.acf
    Q2 <- num.n * t(res.acf2) %*% solve(Omega.mat2) %*% res.acf2
    Q3 <- num.n * t(res.acf3) %*% solve(Omega.mat3) %*% res.acf3 
  }
  
  cat("\n")
  
  cat("fixFACD", distrib, "\n", sep = ",", file = finalparamfile, append = TRUE)
  cat(param, " ", param.var, " ", param.var2, " ", loglik, " ", res.acf3, " ", acf.var3, " ", sep = ",", file = finalparamfile, append = TRUE)
  
  if (portmanteau == TRUE) {
    cat(" ", Q, " ", Q2, " ", Q3, sep = ",", file = finalparamfile, append = TRUE) 
  }
  
  cat(" ", itr, "\n", sep = ",", file = finalparamfile, append = TRUE) 

  return(list(res = err, param = param))
  
}



