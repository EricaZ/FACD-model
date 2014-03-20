FitEACDgeneral <- function(x, id.str, p = 1, q = 1, maxit = 100,
                   init.param = list(w = 0.2, a.vec = 0.2, b.vec = 0.6)) {
  # Fits the exponential ACD model
  #
  # Args:
  #   x: sample
  #   id.str: character string, identifier of the run, for naming of the csv output files
  #
  # 
  #   p: order of lagged x's
  #   q: order of lagged conditional durations
  #   maxit: maximum iteration number
  #   init.param: a list of initial parameters
  #               list(w, a.vec, b.vec)
  #                 w: constant term
  #                 a.vec: a vector containing coefficients of lagged x's
  #                 b.vec: a vector containing coefficients of lagged conditional durations
  #
  # Returns:
  #   an S4 object of the class "acdModel" containing all fitted results
  
  traceparamfile <- paste(id.str, "_trace.csv", sep = "") 
  finalparamfile <- paste(id.str, "_final.csv", sep = "") 
  itrnumfile <- paste(id.str, "_itrnum.csv", sep = "") 
  scorefile <- paste(id.str, "_score.csv", sep = "") 
  
  num.n <- length(x)
  
  # Check validity of inputs
  if (!any(length(init.param$a.vec) == p, length(init.param$b.vec) == q)) {
    stop("Mismatch of orders.")
  }
  
  itr = 0
  
  # Convert parameters into a vector
  initparam.vec <- c(init.param$w, init.param$a.vec, init.param$b.vec)
  param <- init.param  # a list containing w, a.vec, b.vec
  param.vec <- initparam.vec
  
  # Optimization
  repeat {
    itr = itr + 1
    
    # Compute conditional durations
    dur <- ComputeCondDur(x = x, p = p, q = q, 
                          w = param$w, a.vec = param$a.vec, b.vec = param$b.vec)  
    
#     write.csv(dur, "condDur.csv")
    
    # Compute errors
    err <- x / dur
    
    # Compute gradient of conditional durations
    gradient <- ComputeGrad(b.vec = param$b.vec, x = x, p = p, q = q, dur = dur)
    
#     write.csv(gradient, "gradient.csv")

    # Compute information matrix and score function
    k <- mean((err - 1)^2)
    scaled.grad <- matrix(nrow = 1 + p + q, ncol = num.n)
    
    for (i in 1: num.n) {
      scaled.grad[, i] <- gradient[, i] / dur[i]
    }

    B <- (scaled.grad %*% t(scaled.grad)) / num.n

    infomat <- k * B
    score <- scaled.grad %*% (err - 1) / num.n
    
    cat(score, "\n", sep = ",", file = scorefile, append = TRUE)

    # Compute parameters for next iteration
    inv.infomat <- solve(infomat)
    adj <- inv.infomat %*% score 
    newparam.vec <- param.vec + adj
    
    cat("Iteration", itr, "\n")  

#     print(infomat)
#     print(inv.infomat)
#     print(score)
       
#     cat("Iteration", itr, "\n", sep = ",", file = "infomat.csv", append = TRUE)
#     cat("Iteration", itr, "\n", sep = ",", file = "inverse.csv", append = TRUE)
#     for (i in 1:(1+p+q)) {
#       cat(infomat[i, ], "\n", sep = ",", file = "infomat.csv", append = TRUE)
#       cat(inv.infomat[i, ], "\n", sep = ",", file = "inverse.csv", append = TRUE)   
#     }
#     
#     cat(score, "\n", sep = ",", file = "score.csv", append = TRUE)
    
    print(newparam.vec)
    
    for(i in 1:(1+p+q)) {
      cat(newparam.vec[i], " ", sep = ",", file = traceparamfile, append = TRUE)
    }
    
    # ===========================================================================================
    # Impose constraints on new parameters
    ## Lower bound = 0 for w, a.vec, b.vec
    for (i in 1: (1+p+q)) {
      if (newparam.vec[i] < 0) {
        newparam.vec[i] <- initparam.vec[i]
        
        cat("Hit lower bound of parameter", i)
        cat("Hit lower bound of parameter", i, " ", sep = ',', file = traceparamfile, append = TRUE)
      }
    }

    ## Upper bound = 1 for sum of b.vec
    b.sum <- sum(newparam.vec[(2+p):(1+p+q)])
    if (b.sum > 1) {
      newparam.vec[(2+p):(1+p+q)] <- pmax(newparam.vec[(2+p):(1+p+q)] / (b.sum * 1.3), initparam.vec[(2+p):(1+p+q)])
      
      print("Hit upper bound of parameters b.vec")
      cat("Hit upper bound of parameters b.vec ", sep = ',', file = traceparamfile, append = TRUE)
    }

    # =============================================================================================
    # Check stopping condition
    if (max(abs(newparam.vec - param.vec)) < 0.0001) {
      break
    } else {
      if (itr == maxit) {
        print("Maximum iteration number is reached.")
        cat("Maximum iteration number is reached.", sep = ',', file = traceparamfile, append = TRUE)
        break
      }
    }
    
    cat("\n", sep = ",", file = traceparamfile, append = TRUE)
    
    # Renew parameters for next iteration
    param.vec <- newparam.vec
    param <- list(w = newparam.vec[1], 
                  a.vec = newparam.vec[2:(1+p)], 
                  b.vec = newparam.vec[(2+p):(1+p+q)])
  }
  
  cat("\n", sep = ',', file = traceparamfile, append = TRUE)

#   cat("\n", sep = ',', file = "infomat.csv", append = TRUE)
#   cat("\n", sep = ',', file = "inverse.csv", append = TRUE)

  for(i in 1:(1+p+q)) {
    cat(param.vec[i], " ", sep = ",", file = finalparamfile, append = TRUE)
  }
  cat("\n", sep = ',', file = finalparamfile, append = TRUE)
  
  cat(itr, "\n", sep = ",", file = itrnumfile, append = TRUE)

#   std.vec <- sqrt(diag(inv.infomat))
#   pval.vec <- 2*(1 - pt(abs(param.vec / std.vec), length(x) - length(std.vec)))
#   
#   stdParam <- list(w = std.vec[2], 
#                    a.vec = std.vec[3:(p+2)], 
#                    b.vec = std.vec[(p+3):(p+q+2)])
#   
#   pValues <- list(w = pval.vec[2], 
#                   a.vec = pval.vec[3:(p+2)], 
#                   b.vec = pval.vec[(p+3):(p+q+2)])
#   
#   modelOut <- new("ACDModel",
#                 x = x,
#                 p = p,
#                 q = q,                  
#                 param = param,
#                 stdParam = stdParam, 
#                 pValues = pValues,
#                 condDur= dur,
#                 distrib = distrib,
#                 timeRun = date()             
#                 )
# 
#   return(modelOut)  
}