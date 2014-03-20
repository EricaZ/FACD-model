ComputeC_repar <- function(r, s, err) {
  # Computes c1 and c2 for estimation of variance-covariance matrix.
  #
  # Args:
  #   r: shape parameter of frechet distribution
  #   s: scale parameter of frechet distribution
  #   err:
  #
  # Returns:
  #   list(c1, c2)
  #     c1 is vector; c2 is matrix
  
  scaled.err <- err / s
  common <- 1 - scaled.err ^ (-r)
  
  c1 <- r * common
  c2 <- rbind(1 / r - log(scaled.err) * common, r / s * common)
  
  return(list(c1 = c1, c2 = c2))
}