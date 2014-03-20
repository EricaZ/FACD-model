ComputeOmega2 <- function(input, err) {
  # Computes the Omega matrix for the residual ACFs (lags 1 ~ 6)
  
  num.n <- length(err)
  H <- matrix(nrow = nrow(input$infomat1), ncol = 12)
  for (k in 1: 12) {
    H[, k] <- -input$scaled.grad[, (k+1): num.n] %*% (err[1: (num.n-k)] - 1) / num.n
  }
  
  inv.infomat1 <- solve(input$infomat1)
  return(diag(12) - (mean((err - 1)^2))^(-2) * t(H) %*% inv.infomat1 %*% H)
}