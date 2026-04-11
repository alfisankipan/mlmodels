## NB1 FUNCTIONS ===============================================================

# NB1 Hessian Observations
#' @keywords internal
.ml_negbin_nb1_hessianObs <- function(object)
{
  if (!inherits(object, "ml_negbin"))
    cli::cli_abort("`object` must be a model of class 'ml_negbin' (from ml_negin).",
                   call = NULL)
  
  # Add a check that the negbin type is NB1, when I add that to the estimated object
  
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
  mu_alpha = exp(xb - zd)
  
  p_disp <- plogis(zd)
  # p_ndisp <- plogis(zd, lower.tail = FALSE) # Check if this is needed. I think only for gradient.
  l_p_ndisp <- - plogis(zd, lower.tail = FALSE, log.p = TRUE)
  
  psi_mu_alpha_y <- digamma(mu_alpha + y)
  psi_mu_alpha <- digamma(mu_alpha)
  
  psi_1_mu_alpha <- trigamma(mu_alpha)
  psi_1_mu_alpha_y <- trigamma(mu_alpha + y)
  d_p_disp <- dlogis(zd)
  
  s_bb <- w * mu_alpha * (psi_mu_alpha_y - psi_mu_alpha - l_p_ndisp
                          + mu_alpha * (psi_1_mu_alpha_y - psi_1_mu_alpha))
  
  s_bd <- w * mu_alpha * (mu_alpha * ( psi_1_mu_alpha - psi_1_mu_alpha_y)
                          + psi_mu_alpha - psi_mu_alpha_y + l_p_ndisp - p_disp)
  
  s_dd <- w * (mu_alpha * (psi_mu_alpha_y - psi_mu_alpha + p_disp - l_p_ndisp
                           + mu_alpha * ( psi_1_mu_alpha_y - psi_1_mu_alpha)
                           + p_disp^2) - y * d_p_disp)
  
  # H_list <- vector("list", nrow(x))
  
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
  H_stacked
}




# NB1 Likelihood Evaluator -----------------------------------------------------
#' @keywords internal
.ml_negbin_nb1_ll <- function(b, y, x, z, w = NULL, ...)
{
  k1 <- ncol(x) # Number of coefficients for the mean.
  k <- k1 + ncol(z) # Total number of coefficients.
  
  if(is.null(w))
    w <- rep(1, nrow(x))
  
  # Extract the coefficients for the mean.
  beta <- b[1:k1]
  delta <- b[(k1+1):k]
  
  # Useful operations
  xb <- as.vector(x %*% cbind(beta))
  zd <- as.vector(z %*% cbind(delta))
  
  mu_alpha = exp(xb - zd)
  p_disp <- plogis(zd)
  p_ndisp <- plogis(zd, lower.tail = FALSE)
  
  # ln(1 + exp(zd)) = ln(1 / (1 / (1 + exp(zd)))) = ln(1 / p_no) = ln(1) - ln(p_no) = - ln(p_no)
  l_p_ndisp <- - plogis(zd, lower.tail = FALSE, log.p = TRUE)
  
  ## LL
  ll <- w * (lgamma(mu_alpha + y) - lgamma(mu_alpha) - lgamma(y + 1)
             - l_p_ndisp * (mu_alpha + y) + y * zd)
  
  ## GRADIENT
  
  psi_mu_alpha_y <- digamma(mu_alpha + y)
  psi_mu_alpha <- digamma(mu_alpha)


  # Partial with respect to beta.
  gb <- w * mu_alpha * (psi_mu_alpha_y - psi_mu_alpha - l_p_ndisp) * x

  # Partial with respect to delta.
  gd <- w * (mu_alpha * ( psi_mu_alpha - psi_mu_alpha_y - p_disp
                          + l_p_ndisp) + y * p_ndisp ) * z

  ## HESSIAN
  psi_1_mu_alpha <- trigamma(mu_alpha)
  psi_1_mu_alpha_y <- trigamma(mu_alpha + y)
  d_p_disp <- dlogis(zd)

  s_bb <- w * mu_alpha * (psi_mu_alpha_y - psi_mu_alpha - l_p_ndisp
                           + mu_alpha * (psi_1_mu_alpha_y - psi_1_mu_alpha))

  s_bd <- w * mu_alpha * (mu_alpha * ( psi_1_mu_alpha - psi_1_mu_alpha_y)
                           + psi_mu_alpha - psi_mu_alpha_y + l_p_ndisp - p_disp)

  s_dd <- w * (mu_alpha * (psi_mu_alpha_y - psi_mu_alpha + p_disp - l_p_ndisp
               + mu_alpha * ( psi_1_mu_alpha_y - psi_1_mu_alpha) + p_disp^2)
               - y * d_p_disp)


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

.ml_negbin_nb2_ll <- function(b, y, x, z, w = NULL)
{
  #
}
