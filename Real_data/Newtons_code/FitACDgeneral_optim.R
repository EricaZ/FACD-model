FitACDgeneral_optim <- function(x, id.str, p = 1, q = 1, maxit = 500, init.param = c(2,0.1,0.2,0.6), 
                          distrib = "frechet") {
  
  finalparamfile <- paste(id.str, "_record.csv", sep = "") 
  
  # Optimization
  options(warn=-1)
  param <- init.param
  fitted.param <- optim(par = param,
                      fn = ACD_Lik, x = x, p = p, q = q, distrib = distrib,
                      method ="BFGS" , hessian = TRUE,
                      control = list(maxit=maxit,fnscale=-1))
  options(warn=1)
  
#   stdParam<-sqrt(diag(solve(-fitted.param$hessian)))
#   pValues<-2*(1-pt(abs(fitted.param$par/stdParam),length(x)-length(stdParam)))
#   
  print(fitted.param$par)
   
  cat("\n")
  cat(paste("general FACD with p = ", p, " and q = ", q, sep=""),"\n", sep = ",", file = finalparamfile, append = TRUE)
  cat(fitted.param$par, "\n", sep = ",", file = finalparamfile, append = TRUE)
  
  return(ACD_Lik(fitted.param$par, x, p, q, distrib="frechet", lik = FALSE))

}