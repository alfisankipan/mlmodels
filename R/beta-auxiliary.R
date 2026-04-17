## Gradients by Observation ====================================================
#' @keywords internal
.ml_beta_gradientObs <- function(object)
{
  if (!inherits(object, "ml_beta"))
    cli::cli_abort("`object` must be a model of class 'ml_beta'.",
                   call = NULL)
  
  b <- coef(object)
  y <- object$model$value$outcomes[[1]]
  x <- as.matrix(object$model$value$predictors)
  z <- as.matrix(object$model$scale$predictors)
  w <- if(is.null(object$model$weights))
    rep(1, nrow(x))
  else
    object$model$weights
  k1 <- ncol(x)
  k <- k1 + ncol(z)
  
  if (length(b) != k)
    cli::cli_abort("The length of the coefficients ({length(b)}) \\
                   does not match with the number of parameters ({k}).",
                   call = NULL)
  
  beta <- b[1:k1]
  delta <- b[(k1+1):k]
  
  xb <- as.vector(x %*% cbind(beta))
  zd <- as.vector(z %*% cbind(delta))
  
  mu_y <- plogis(xb)
  mu_n <- plogis(xb, lower.tail = FALSE)
  phi <- exp(zd)
  d_mu <- dlogis(xb)
  
  gb <- w * phi * d_mu * (digamma(mu_n * phi) - digamma(mu_y * phi) +
                            log(y) - log(1 - y)) * x
  gd <- w * phi * (digamma(phi) - digamma(phi * mu_y) * mu_y -
                     digamma(phi * mu_n) * mu_n + mu_y * log(y) +
                     mu_n * log(1 - y)) * z
  
  g <- cbind(gb, gd)
  return(g)
}

# Hessians by Observation ======================================================
.ml_beta_hessianObs <- function(object)
{
  if (!inherits(object, "ml_beta"))
    cli::cli_abort("`object` must be a model of class 'ml_beta'.",
                   call = NULL)
  
  b <- coef(object)
  y <- object$model$value$outcomes[[1]]
  x <- as.matrix(object$model$value$predictors)
  z <- as.matrix(object$model$scale$predictors)
  w <- if(is.null(object$model$weights))
    rep(1, nrow(x))
  else
    object$model$weights
  k1 <- ncol(x)
  k <- k1 + ncol(z)
  
  if (length(b) != k)
    cli::cli_abort("The length of the coefficients ({length(b)}) \\
                   does not match with the number of parameters ({k}).",
                   call = NULL)
  
  beta <- b[1:k1]
  delta <- b[(k1+1):k]
  
  xb <- as.vector(x %*% cbind(beta))
  zd <- as.vector(z %*% cbind(delta))
  
  mu_y <- plogis(xb)
  mu_n <- plogis(xb, lower.tail = FALSE)
  phi <- exp(zd)
  d_mu <- dlogis(xb)
  
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
  
  H_stacked <- matrix(0, nrow = nrow(x) * k, ncol = k)
  
  for(i in seq_len(nrow(x)))
  {
    # Extracting the elements for the observation we need.
    xi <- cbind(x[i, ])
    zi <- cbind(z[i, ])
    
    # Second partial with respect both times to beta.
    h_bb <- s_bb[i] * tcrossprod(xi)
    
    # Second partial first with respect to beta and then to delta
    h_bd <- s_bd[i] * tcrossprod(xi,zi)
    
    # Transpose that.
    h_db <- t(h_bd)
    
    # Second partial with respect both times to lnsigma.
    h_dd <- s_dd[i] * tcrossprod(zi)
    
    start_row <- (i - 1) * k + 1
    end_row <- i * k
    
    # Form the observation's Hessian
    H_stacked[start_row:end_row, ] <- rbind(cbind(h_bb, h_bd),
                                            cbind(h_db, h_dd))
  }
  # Stack all individual Hessians
  return(H_stacked)
}

# Log-likelihood by Observations ===============================================
.ml_beta_loglikeObs <- function(object)
{
  if (!inherits(object, "ml_beta"))
    cli::cli_abort("`object` must be a model of class 'ml_beta'.",
                   call = NULL)
  
  b <- coef(object)
  y <- object$model$value$outcomes[[1]]
  x <- as.matrix(object$model$value$predictors)
  z <- as.matrix(object$model$scale$predictors)
  w <- if(is.null(object$model$weights))
    rep(1, nrow(x))
  else
    object$model$weights
  k1 <- ncol(x)
  k <- k1 + ncol(z)
  
  if (length(b) != k)
    cli::cli_abort("The length of the coefficients ({length(b)}) \\
                   does not match with the number of parameters ({k}).",
                   call = NULL)
  
  beta <- b[1:k1]
  delta <- b[(k1+1):k]
  
  xb <- as.vector(x %*% cbind(beta))
  zd <- as.vector(z %*% cbind(delta))
  
  mu_y <- plogis(xb)
  mu_n <- plogis(xb, lower.tail = FALSE)
  phi <- exp(zd)
  
  ll <- w * (lgamma(phi) - lgamma(mu_y * phi) - lgamma(mu_n * phi) +
               (mu_y * phi - 1) * log(y) + (mu_n * phi - 1) * log(1 - y))
  
  return(ll)
}

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
  mu_n <- plogis(xb, lower.tail = FALSE)
  phi <- exp(zd_clamp)
  
  ll <- w * (lgamma(phi) - lgamma(mu_y * phi) - lgamma(mu_n * phi) +
               (mu_y * phi - 1) * log(y) + (mu_n * phi - 1) * log(1 - y))
  
  # Gradient
  d_mu <- dlogis(xb)
  
  gb <- w * phi * d_mu * (digamma(mu_n * phi) - digamma(mu_y * phi) +
                             log(y) - log(1 - y)) * x
  gd <- w * phi * (digamma(phi) - digamma(phi * mu_y) * mu_y -
                      digamma(phi * mu_n) * mu_n + mu_y * log(y) +
                      mu_n * log(1 - y)) * z
  
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