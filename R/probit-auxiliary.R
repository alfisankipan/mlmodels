## ML EVALUATOR ==============================================================
#' @keywords internal
.ml_probit_ll <- function(b, y, x, z = NULL, w = NULL)
{
  if(is.null(w))
    w <- rep(1, nrow(x))
  if (is.null(z)) {
    # -- Homoskedastic binary probit -------------------------------------------
    xb <- as.vector(x %*% cbind(b))
    py <- pnorm(xb)
    pn <- pnorm(-xb)
    
    # Weighted log-likelihood
    ll <- ((1-y) * pnorm(-xb, log.p = TRUE) + y * pnorm(xb, log.p = TRUE)) * w
    
    # Gradient
    g <- w * as.vector(y * dnorm(xb) / py - (1 - y) * dnorm(xb) / pn) * x
    
    # Hessian
    scalar <- w * ((1 - y) * dnorm(xb) / pn * (xb - dnorm(xb) / pn) -
            y * dnorm(xb) / py * (xb + dnorm(xb) / py))
    
    H   <- crossprod(x * scalar, x)
  } else {
    # -- Heteroskedastic binary probit ------------------------------------------
    k1   <- ncol(x)
    k    <- k1 + ncol(z)
    beta <- b[1:k1]
    delta <- b[(k1+1):k]
    
    zd <- as.vector(z %*% cbind(delta))
    xz <- x / exp(zd)
    xb <- as.vector(xz %*% cbind(beta)) # <- this is actually xb / sigma.
    py <- pnorm(xb)
    pn <- pnorm(-xb)
    den <- dnorm(xb)
    
    # Weighted log-likelihood
    ll <- ((1-y) * pnorm(-xb, log.p = TRUE) + y * pnorm(xb, log.p = TRUE)) * w
    
    # Gradient
    gb <- w * as.vector(y * den / py - (1 - y) * den / pn) * xz
    gd <- w * as.vector((1 - y) * den / pn - y * den / py) * xb * z
    g  <- cbind(gb, gd)
    
    # Hessian
    s_bb <- as.vector(w * ((1 - y) * den / pn * (xb - den / pn) -
                             y * den / py * (xb + den / py)))
    
    s_bd <- as.vector(w * (y * den / py * (xb * (xb + den / py) - 1) -
                             (1 - y) * den / pn * (xb * (xb - den / pn) - 1)))
    
    s_dd <- as.vector(w * xb^2 * ((1 - y) * den / pn * (xb - den / pn - 1) -
                                    y * den / py * (xb + den / py - 1)))
    
    # Remember, xz has the standardized x. :D
    H_bb <- crossprod(xz * s_bb, xz)
    H_bd <- crossprod(xz * s_bd, z)
    H_dd <- crossprod(z * s_dd, z)
    
    H <- rbind(cbind(H_bb, H_bd),
               cbind(t(H_bd), H_dd))
  }
  
  # Attach gradient and Hessian as attributes
  attr(ll, "gradient") <- g
  attr(ll, "hessian")  <- H
  
  return(ll)
}