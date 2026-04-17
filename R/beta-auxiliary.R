## ML EVALUATOR ================================================================
#' @keywords internal
.ml_beta_ll <- function(b, y, x, z, w = NULL, ...)
{
  k1 <- ncol(x) # Number of coefficients for the mean.
  k <- k1 + ncol(z) # Total number of coefficients.
  
  if(is.null(w))
    w <- rep(1, nrow(x))
  
  beta <- b[1:k1]
  delta <- b[(k1+1):k]
  
  # Useful operations
  xb <- as.vector(x %*% cbind(beta))
  zd <- as.vector(z %*% cbind(delta))
  
  # Clamping for stability
  up_bound <- 700
  lo_bound <- -100
  
  zd_clamp <- pmax(pmin(zd, up_bound), lo_bound)
  
  mu_y <- plogis(xb)
  mu_n <- plogis(xb, lower.tail = TRUE)
  phi <- exp(zd_clamp)
  
  ll <- w * (lgamma(phi) - lgamma(mu_y * phi) - lgamma(mu_n * phi) +
               (mu_y * phi -1) * log(y) + (mu_n * phi - 1) * log(1 - y))
  
  # Gradient
  d_mu <- dlogis(xb)
  
  gb <- w * phi * d_mu * (digamma(mu_n * phi) - digamma(mu_y * phi) +
                             log(y) - log(1 - y)) * x
  gd <- w * phi * (digamma(phi) - digamma(phi * mu_y) * mu_y -
                      digamma(phi * mu_n) * mu_n + mu_y * log(y) +
                      mu_n * log(1 - y)) * z
  
  g <- c(gb, gd)
  
  # Hessian
  s_bb <- w * phi * d_mu *  ((1 - 2 * mu_y) * (digamma(mu_n * phi) -
                  digamma(mu_y * phi) + log(y) - log(1 - y)) - phi * d_mu *
                  (trigamma(mu_n * phi) + trigamma(mu_y * phi)))
  s_bd <- w * phi * d_mu * (digamma(mu_n * phi) - digamma(mu_y * phi) + log(y) -
                            log(1 - y) + phi * (trigamma(mu_n * phi) * mu_n -
                            trigamma(mu_y * phi) * mu_y))
  s_dd <- w * phi * (digamma(phi) - digamma(phi * mu_y) * mu_y -
                     digamma(phi * mu_n) * mu_n + mu_y * log(y) +
                     mu_n * log(1 - y) + phi * (trigamma(phi) -
                     trigamma(mu_y * phi) * mu_y^2 - trigamma(mu_n * phi) * mu_n ^2))
  
  H_bb <- crossprod(x * s_bb, x)
  H_bd <- crossprod(x * s_bd, z)
  H_dd <- crossprod(z * s_dd, z)
  
  H <- rbind(cbind(H_bb, H_bd),
             cbind(t(H_bd), H_dd))
  
  # Set the attribute in ll to pass it back to maxLik
  attr(ll, "gradient") <- cbind(gb, gd)
  
  # Set the attribute in ll to pass it back to maxLik
  attr(ll, "hessian") <- H
  
  return(ll)
}