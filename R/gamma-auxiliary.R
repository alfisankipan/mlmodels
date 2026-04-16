.ml_gamma_ll <- function(b, y, x, z, w = NULL, ...)
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
  zd_xb_clamp <- pmax(pmin(zd - xb, up_bound), lo_bound)
  nu <- exp(zd_clamp)
  e_zd_xb <- exp(zd_xb_clamp)
  
  ll <- w * (- lgamma(nu) + nu * (zd - xb) + (nu - 1) * log(y) - y * e_zd_xb)
  
  # Gradient
  gb <- w * (y * e_zd_xb - nu) * x
  gd <- w * (nu * (- digamma(nu) + zd - xb + log(y) + 1) - y * e_zd_xb) * z
  
  g <- c(gb, gd)
  
  # Hessian
  s_bb <- w * (- y * e_zd_xb)
  s_bd <- w * (y * e_zd_xb - nu)
  s_dd <- w * (nu * (- digamma(nu) - trigamma(nu) * nu + zd - xb + log(y) + 2) - y * e_zd_xb)
  
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