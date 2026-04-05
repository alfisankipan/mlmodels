## HESSIANS ====================================================================
#' Observed Hessian by observation for ml_lm models
#'
#' Returns the Hessian matrix evaluated at each observation. Used internally
#' by the Information Matrix test (`IMtest`).
#'
#' @param object A fitted `ml_lm` object.
#'
#' @return Matrix with one row per observation and columns corresponding to the
#'   second derivatives with respect to all parameters.
#'
#' @keywords internal
.ml_lm_hessianObs <- function(object)
{
  if (!inherits(object, "ml_lm"))
    cli::cli_abort("`object` must be a model of class 'ml_lm' (from ml_lm).",
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
  
  xb <- x %*% cbind(beta)
  zd <- z %*% cbind(delta)
  s <- exp(zd)
  u <- y - xb
  
  # For each observation we form the hessian and we add it as a new row on an
  # overall matrix to return.
  H <- NULL
  for(i in 1:nrow(x))
  {
    # Extracting the elements for the observation we need.
    xi <- cbind(x[i, ])
    zi <- cbind(z[i, ])
    si <- s[i]
    ui <- u[i]
    wi <- w[i]
    
    # Second partial with respect both times to beta.
    hbb <- - wi * si^(-2) * tcrossprod(xi)
    
    # Second partial first with respect to beta and then to s
    hbs <- -2 * wi * (ui / si^2) * tcrossprod(xi,zi)
    
    # Transpose that.
    hsb <- t(hbs)
    
    # Second partial with respect both times to lnsigma.
    hss <- -2 * wi * (ui / si)^2 * tcrossprod(zi)
    
    # Form the observation's Hessian
    h <- rbind(cbind(hbb, hbs),
               cbind(hsb, hss))
    # Add it to the final Hessian.
    H <- rbind(H, h)
  }
  return(H)
}

## ML EVALUATOR ================================================================
#' Log-likelihood, gradient and Hessian calculation for ml_lm models
#'
#' This is the function passed to [maxLik::maxLik()] for Gaussian linear models.
#' It returns the log-likelihood and attaches the gradient and Hessian as
#' attributes (required by maxLik when analytical derivatives are available).
#'
#' @param b Numeric vector of all coefficients (value + scale).
#' @param y Numeric vector of outcomes.
#' @param x Design matrix for the value equation.
#' @param z Design matrix for the scale equation.
#' @param w Numeric vector of weights (can be `NULL`).
#'
#' @return Numeric vector of log-likelihood values (one per observation).
#'   The returned object has attributes `gradient` and `hessian`.
#'
#' @keywords internal
.ml_lm_ll <- function(b, y, x, z, w = NULL, lognormal = FALSE)
{
  # The last coefficient in b is the coefficient for the natural log of sigma
  k1 <- ncol(x) # Number of coefficients for the mean.
  k <- k1 + ncol(z) # Total number of coefficients.

  # If we don't have weights we set them to ones to be able to use the same
  # formulas underneath.
  if(is.null(w))
    w <- rep(1, nrow(x))

  # Extract the coefficients for the mean.
  beta <- b[1:k1]
  delta <- b[(k1+1):k]

  # Useful operations
  xb <- x %*% cbind(beta)
  zd <- z %*% cbind(delta)
  s <- exp(zd)
  u <- y - xb

  ## LL
  ll <- dnorm(u / s, log = TRUE) - zd
  if(lognormal) ll <- ll - y # y because it's already log-transformed
  ll <- ll * w

  ## GRADIENT
  # Partial with respect to beta.
  gb <- w * as.vector(u / s^2) * x

  # Partial with respect to delta.
  gd <- w * as.vector((u / s)^2 - 1) * z

  # Set the attribute in ll to pass it back to maxLik
  attr(ll, "gradient") <- cbind(gb, gd)

  ## HESSIAN

  # For the Hessian we have to calculate the matrix by observation and then
  # add it to an overall matrix that I set with zero values.
  H <- matrix(0, nrow = k, ncol = k)

  for(i in 1:nrow(x))
  {
    # Extracting the elements for the observation we need.
    xi <- cbind(x[i, ])
    zi <- cbind(z[i, ])
    si <- s[i]
    ui <- u[i]
    wi <- w[i]

    # Second partial with respect both times to beta.
    hbb <- - wi * si^(-2) * tcrossprod(xi)

    # Second partial first with respect to beta and then to s
    hbs <- -2 * wi * (ui / si^2) * tcrossprod(xi,zi)

    # Transpose that.
    hsb <- t(hbs)

    # Second partial with respect both times to lnsigma.
    hss <- -2 * wi * (ui / si)^2 * tcrossprod(zi)

    # Form the observation's Hessian
    h <- rbind(cbind(hbb, hbs),
               cbind(hsb, hss))
    # Add it to the final Hessian.
    H <- H + h
  }

  # Set the attribute in ll to pass it back to maxLik
  attr(ll, "hessian") <- H

  return(ll)
}

# VARIANCE BOOTSTRAP -----------------------------------------------------------
#' Bootstrap variance-covariance matrix for ml_lm models
#'
#' Specialized (faster) bootstrap method for `ml_lm` objects. It calls
#' `.ml_lm.fit()` directly instead of going through the generic `update()`
#' method.
#'
#' @param object A fitted `ml_lm` object.
#' @param repetitions Number of bootstrap replications.
#' @param seed Random seed for reproducibility.
#' @param cl_var Clustering variable (if doing clustered bootstrap).
#' @param progress Logical. Whether to show a progress bar.
#' @param ... Not currently used.
#'
#' @keywords internal
#' @keywords internal
.vcov_boot.ml_lm <- function(object,
                             repetitions = 999,
                             seed = NULL,
                             cl_var = NULL,
                             progress = TRUE,
                             ...)
{
  
  # --- 0. Validity Checks -----------------------------------------------------
  if(!inherits(object, "ml_lm"))
    cli::cli_abort("`object` needs to be of class 'ml_lm'.")

  if (is.null(seed)) seed <- sample.int(1e6, 1)
  set.seed(seed)

  if (is.null(object$model$value$outcomes) ||
      is.null(object$model$value$predictors) ||
      is.null(object$model$scale$predictors))
    cli::cli_abort("The sample data was not stored properly.")

  # --- 1. Sample data extraction. ---------------------------------------------
  y <- object$model$value$outcomes[[1]]
  x <- as.matrix(object$model$value$predictors)
  z <- as.matrix(object$model$scale$predictors)
  n <- nrow(x)
  if(is.null(object$model$weights))
    w <- rep(1,n)
  else
    w <- object$model$weights

  is_clustered <- !is.null(cl_var)
  if (is_clustered) {
    cluster_ids <- unique(cl_var[object$model$sample])
    n_cluster   <- length(cluster_ids)
  }
  
  # --- 2. Bootstrap area ------------------------------------------------------
  if (progress) {
    if (is_clustered) {
      cli::cli_alert_info("Clustered bootstrap with {.val {repetitions}} repetitions and {.val {n_cluster}} clusters.")
    } else {
      cli::cli_alert_info("Bootstrap with {.val {repetitions}} repetitions.")
    }
    cat(cli::col_blue(" 0"))
    for (i in seq(10, 50, by = 10)) cat(cli::col_blue(sprintf("%10d", i)))
    cat("\n")
    cat(cli::col_blue(strrep("=", 52), "\n"))
  }

  success     <- logical(repetitions)
  coef_matrix <- matrix(NA_real_, nrow = repetitions, ncol = length(coef(object)))

  if (!is.null(object$model$constraints$maxLik)) {
    cli::cli_warn(
      c("Bootstrap variance with constraints may be unreliable.",
        "i" = "Different bootstrap samples often produce infeasible log-likelihoods at the supplied starting values.",
        "i" = "Equality constraints in particular lead to very low convergence rates.",
        "i" = "Consider using `type = 'robust'` or `type = 'cluster'` (with `cl_var`) instead.")
    )
  }
  
  # --- 2.1 Bootstrap loop -----------------------------------------------------
  for (i in seq_len(repetitions)) {
    if (progress && i %% 50 == 1 && i > 1) cat("\n ")
    else if(progress && i == 1) cat(" ")

    tryCatch({
      if (is_clustered) {
        sampled_clusters <- sample(cluster_ids, size = n_cluster, replace = TRUE)
        boot_idx <- unlist(lapply(sampled_clusters, function(cid) {
          which(cl_var[object$model$sample] == cid)
        }))
      } else {
        boot_idx <- sample(n, n, replace = TRUE)
      }
      
      y_boot <- y[boot_idx]
      x_boot <- x[boot_idx, , drop = FALSE]
      z_boot <- z[boot_idx, , drop = FALSE]
      w_boot <- w[boot_idx]
      
      updated <- .ml_lm.fit(y = y_boot,
                            x = x_boot,
                            z = z_boot,
                            w = w_boot,
                            lognormal = object$model$log_info$value$is_log,
                            constraints = object$model$constraints$maxLik,
                            start       = object$model$start,
                            method      = object$model$method,
                            control     = object$model$control)
      
      if (updated$code %in% c(0L, 1L, 2L, 8L)) {
        if (progress) cat(cli::col_green("."))
        success[i] <- TRUE
        coef_matrix[i, ] <- coef(updated)
      } else {
        if (progress) cat(cli::col_red("x"))
        success[i] <- FALSE
        coef_matrix[i, ] <- NA_real_
      }
    }, error = function(e) {
      if (progress) cat(cli::col_red("x"))
      success[i] <- FALSE
      coef_matrix[i, ] <- NA_real_
    })
  }

  if (progress) {
    cat("\n")
    cat(cli::col_blue(strrep("=", 52), "\n"))
  }

  # --- 3. Final reporting -----------------------------------------------------
  if (progress) {
    cat("\n")
    cli::cli_text("Bootstrapping finished - {round(mean(success) * 100, 1)}% of replications converged.")
  }
  
  if (mean(success) < 0.7) {
    cli::cli_warn("Low convergence rate - bootstrap results may be unreliable.")
  }

  # Variance from successful replications only
  valid_rows <- complete.cases(coef_matrix)
  vcov_boot  <- var(coef_matrix[valid_rows, , drop = FALSE])

  dimnames(vcov_boot) <- list(names(coef(object)), names(coef(object)))
  vcov_boot
}
