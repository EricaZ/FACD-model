FitACDgeneral_fix <- function(x, id.str, p = 1, q = 1, maxit = 100, r = 2, portmanteau = FALSE,
                   init.param = list(w = 0.2, a.vec = 0.2, b.vec = 0.6),
                   distrib = "frechet") {
  # Fits the general ACD model
  #
  # Args:
  #   x: sample
  #   id.str: character string, identifier of the run, for naming of the csv output files
  #   p: order of lagged x's
  #   q: order of lagged conditional durations
  #   maxit: maximum iteration number
  #   init.param: a list of initial parameters
  #               list(w, a.vec, b.vec, r)
  #                 w: constant term
  #                 a.vec: a vector containing coefficients of lagged x's
  #                 b.vec: a vector containing coefficients of lagged conditional durations
  #                 r: shape parameter of error distribution  
  #   distrib: "frechet" or "weibull"
  #
  
#   traceparamfile <- paste(id.str, "_trace.csv", sep = "") 
  finalparamfile <- paste(id.str, "_record.csv", sep = "") 
  
  # Check validity of inputs
  if (!any(length(init.param$a.vec) == p, length(init.param$b.vec) == q)) {
    stop("Mismatch of orders.")
  }
  
  num.n <- length(x)
  itr = 0
  
  # Convert parameters into a vector
  initparam.vec <- c(r, init.param$w, init.param$a.vec, init.param$b.vec)
  param <- init.param  # a list containing w, a.vec, b.vec
  param.vec <- initparam.vec
  
  # Optimization
  repeat {
    itr = itr + 1
    
    # Compute conditional durations
    dur <- ComputeCondDurgeneral(x = x, p = p, q = q, 
                          w = param$w, a.vec = param$a.vec, b.vec = param$b.vec)  

    # Compute errors
    err <- x / dur
    
    # Compute gradient of conditional durations
    gradient <- ComputeGradgeneral(b.vec = param$b.vec, x = x, p = p, q = q, dur = dur)

    # Compute c1
    c1 <- ComputeC_dual(r = r, err = err, distrib = distrib)
    
    # Compute information matrix and score function
    infoscore <- ComputeISgeneral_dual(c1 = c1, dur = dur, gradient = gradient)
    
    # Compute parameters for next iteration
    inv.infomat <- solve(infoscore$infomat)
    adj <- inv.infomat %*% infoscore$score 
    newparam.vec <- param.vec + c(0, adj)
    
    cat("Iteration", itr, "\n")  
    print(newparam.vec)
    
#     for(i in 1:(2+p+q)) {
#       cat(newparam.vec[i], " ", sep = ",", file = traceparamfile, append = TRUE)
#     }
    
    # ===========================================================================================
    # Impose constraints on new parameters
    ## Lower bound = 0 for w, a.vec, b.vec
    for (i in 2:(2+p+q)) {
      if (newparam.vec[i] < 0) {
        newparam.vec[i] <- initparam.vec[i]
      }
    }
    
    ## Lower bound = 1 for r, for frechet distrib
    if((distrib == "frechet") && (newparam.vec[1] < 1)) {
      newparam.vec[1] <- initparam.vec[1]
    }
    
    ## Lower bound = 0 for r, for weibull distrib
    if((distrib == "weibull") && (newparam.vec[1] < 0)) {
      newparam.vec[1] <- initparam.vec[1]
    }

    ## Upper bound = 1 for sum of b.vec
    b.sum <- sum(newparam.vec[(3+p):(2+p+q)])
    if (b.sum > 1) {
      newparam.vec[(3+p):(2+p+q)] <- pmax(newparam.vec[(3+p):(2+p+q)] / (b.sum * 1.3), initparam.vec[(3+p):(2+p+q)])
    }

    # =============================================================================================
    # Check stopping condition
    if (max(abs(newparam.vec - param.vec)) < 0.0001) {
      break
    } else {
      if (itr == maxit) {
        break
      }
    }
    
#     cat("\n", sep = ",", file = traceparamfile, append = TRUE)
    
    # Renew parameters for next iteration
    param.vec <- newparam.vec
    param <- list(w = newparam.vec[2], 
                  a.vec = newparam.vec[3:(p+2)], 
                  b.vec = newparam.vec[(p+3):(p+q+2)])
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
Omega.input <- ComputeISgeneral(c = c, dur = dur, gradient = gradient, IS2 = TRUE)

Omega.mat <- ComputeOmega(input = Omega.input, err = err)
acf.var <- diag(Omega.mat) / num.n

Omega.mat2 <- ComputeOmega2(input = Omega.input, err = err)
acf.var2 <- diag(Omega.mat2) / num.n

Omega.mat3 <- ComputeOmega3(input = Omega.input, err = err)
acf.var3 <- diag(Omega.mat3) / num.n

# Compute the portmanteau test statistic
if (portmanteau == TRUE) {
  Q <- num.n * t(res.acf) %*% solve(Omega.mat) %*% res.acf
  Q2 <- num.n * t(res.acf2) %*% solve(Omega.mat2) %*% res.acf2
  Q3 <- num.n * t(res.acf3) %*% solve(Omega.mat3) %*% res.acf3 
}

cat("\n")

cat(paste("fix_general ", distrib, " ACD with p = ", p, " and q = ", q, sep=""), "\n", sep = ",", file = finalparamfile, append = TRUE)
cat(param.vec, " ", param.var, " ", res.acf3, " ", acf.var3, " ", sep = ",", file = finalparamfile, append = TRUE)

if (portmanteau == TRUE) {
  cat(" ", Q, " ", sep = ",", file = finalparamfile, append = TRUE) 
  cat(" ", Q2, " ", sep = ",", file = finalparamfile, append = TRUE) 
  cat(" ", Q3, " ", sep = ",", file = finalparamfile, append = TRUE) 
}

cat(" ", itr, "\n", sep = ",", file = finalparamfile, append = TRUE) 

  return(list(res = err, param = param.vec))
}