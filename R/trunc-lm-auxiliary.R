## ML EVALUATOR ================================================================
#' @keywords internal
.ml_trunc_lm_ll <- function(b, y, x, z, left, right, w = NULL, lognormal = FALSE)
{
  # The last coefficients in b are the coefficients for the natural log of sigma
  k1 <- ncol(x) # Number of coefficients for the mean.
  k <- k1 + ncol(z) # Total number of coefficients.
  
  # If we don't have weights we set them to ones to be able to use the same
  # formulas underneath.
  if(is.null(w))
    w <- rep(1, nrow(x))
  
  # Extract the coefficients
  beta <- b[1:k1]
  delta <- b[(k1+1):k]
  
  # Useful operations
  xb <- as.vector(x %*% cbind(beta))
  zd <- as.vector(z %*% cbind(delta))
  s <- exp(zd)
  z_y <- (y - xb) / s
  a <- (left - xb) / s
  g <- (right - xb) / s
  phi_a <- dnorm(a)
  phi_g <- dnorm(g)
  Phi_a <- pnorm(a)
  Phi_g <- pnorm(g)
  
  denom <- Phi_g - Phi_a
  
  # Calculate components
  g_phi_g <- g * phi_g
  a_phi_a <- a * phi_a
  g2_phi_g <- (g^2 - 1) * Phi_g
  a2_phi_a <- (a^2 - 1) * Phi_a
  g_g2_phi_g <- g * phi_g * (1 - g^2)
  a_a2_phi_a <- a * phi_a * (1 - a^2)
  
  # Clean Inf * 0 cases
  g_phi_g[is.infinite(g) & phi_g == 0] <- 0
  a_phi_a[is.infinite(a) & phi_a == 0] <- 0
  g2_phi_g[is.infinite(g) & phi_g == 0] <- 0
  a2_phi_a[is.infinite(a) & phi_a == 0] <- 0
  g_g2_phi_g[is.infinite(g) & phi_g == 0] <- 0
  a_a2_phi_a[is.infinite(a) & phi_a == 0] <- 0
  
  # 1. Define the components for the ratios
  num_lam  <- phi_g - phi_a
  num_wlam <- g_phi_g - a_phi_a
  num_bd   <- g2_phi_g - a2_phi_a
  num_dd   <- g_g2_phi_g - a_a2_phi_a
  
  # 2. Calculate the ratios
  lam      <- num_lam / denom
  wlam     <- num_wlam / denom
  ratio_bd <- num_bd / denom
  ratio_dd <- num_dd / denom
  
  xs <- x / s                                                                   # Standardized to use as matrix in gradient and Hessian
  
  ## LL
  ll <- dnorm(z_y, log = TRUE) - zd - log(denom)
  if(lognormal) ll <- ll - y                                                    # y because it's already log-transformed
  ll <- ll * w
  
  ## GRADIENT
  # Partial with respect to beta.
  gb <- w * as.vector(z_y + lam) * xs
  
  # Partial with respect to delta.
  gd <- w * as.vector(z_y^2 + wlam - 1) * z
  
  ## HESSIAN
  
  s_bb <- w * as.vector(wlam + lam^2 - 1)
  
  s_bd <- w * as.vector(ratio_bd + lam * wlam - 2 * z_y)
  
  s_dd <- w * as.vector(wlam^2 - 2 * z_y^2 - ratio_dd)
  
  
  H_bb <- crossprod(xs * s_bb, xs)
  H_bd <- crossprod(xs * s_bd, z)
  H_dd <- crossprod(z * s_dd, z)
  
  H <- rbind(cbind(H_bb, H_bd),
             cbind(t(H_bd), H_dd))
  
  # Set the attribute in ll to pass it back to maxLik
  attr(ll, "gradient") <- cbind(gb, gd)
  
  # Set the attribute in ll to pass it back to maxLik
  attr(ll, "hessian") <- H
  
  return(ll)
}