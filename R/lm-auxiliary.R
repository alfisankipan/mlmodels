# ML EVALUATOR -----------------------------------------------------------------
# Function to evaluate the log-likelihood, gradient and Hessian.
# The first argument in these functions always has to be the vector with ALL
# the coefficients (parameters) that we're estimating, because that's what
# maxLik() requires. You can then add any quantity of other arguments you
# may need for the function, and you will pass those arguments to maxLik(), so
# that it can then pass them through to this function.
ml_lm_ll <- function(b, y, x, z, w = NULL)
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

  ## GRADIENT
  # Partial with respect to beta.
  gb <- as.vector(u / s^2) * x

  # Partial with respect to delta.
  gd <- as.vector((u / s)^2 - 1) * z

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

    # Second partial with respect both times to beta.
    hbb <- -si^(-2) * tcrossprod(xi)

    # Second partial first with respect to beta and then to s
    hbs <- -2 * (ui / si^2) * tcrossprod(xi,zi)

    # Transpose that.
    hsb <- t(hbs)

    # Second partial with respect both times to lnsigma.
    hss <- -2 * (ui / si)^2 * tcrossprod(zi)

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

# HESSIANS ---------------------------------------------------------------------
# Function to return the hessians of all observations of the linear model.
# The must be a maxLik object where we've added a term called model
# which is a list of two elements: value and scale. Each of them is the object
# returned by hardhat's mold() function, for its respective equation.
ml_lm_hessianObs <- function(object)
{
  if (!inherits(object, "ml_lm"))
    cli::cli_abort("`object` must be a model of class 'ml_lm' (from ml_lm).",
                   call = NULL)

  b <- coef(object)
  y <- object$model$value$outcomes[[1]]
  x <- as.matrix(object$model$value$predictors)
  z <- as.matrix(object$model$scale$predictors)
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

    # Second partial with respect both times to beta.
    hbb <- -si^(-2) * tcrossprod(xi)

    # Second partial first with respect to beta and then to s
    hbs <- -2 * (ui / si^2) * tcrossprod(xi,zi)

    # Transpose that.
    hsb <- t(hbs)

    # Second partial with respect both times to lnsigma.
    hss <- -2 * (ui / si)^2 * tcrossprod(zi)

    # Form the observation's Hessian
    h <- rbind(cbind(hbb, hbs),
               cbind(hsb, hss))
    # Add it to the final Hessian.
    H <- rbind(H, h)
  }
  return(H)
}

# VARIANCE BOOTSTRAP -----------------------------------------------------------
#' @keywords internal
vcov_boot.ml_lm <- function(object,
                            repetitions = 999,
                            seed = NULL,
                            cl_var = NULL,
                            progress = TRUE,
                            ...)
{
  if(!inherits(object, "ml_lm"))
    cli::cli_abort("`object` needs to be of class 'ml_lm'.")

  if (is.null(seed)) seed <- sample.int(1e6, 1)
  set.seed(seed)

  if (is.null(object$model$value$outcomes) ||
      is.null(object$model$value$predictors) ||
      is.null(object$model$scale$predictors))
    cli::cli_abort("The sample data was not stored properly.")

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

  # ── Bootstrap loop ───────────────────────────────────────────────────
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

  for (i in seq_len(repetitions)) {
    if (progress && i %% 50 == 1 && i > 1) cat("\n ")
    else if(progress && i == 1) cat(" ")

    res <- tryCatch({
      if (is_clustered) {
        # Clustered bootstrap - sample whole clusters
        sampled_clusters <- sample(cluster_ids, size = n_cluster, replace = TRUE)
        boot_idx <- unlist(lapply(sampled_clusters, function(cid) {
          which(cl_var[object$model$sample] == cid)
        }))
      } else {
        # Regular (individual) bootstrap
        boot_idx <- sample(nrow(used_data), nrow(used_data), replace = TRUE)
      }

      y_boot    <- y[boot_idx]
      x_boot    <- x[boot_idx, , drop = FALSE]
      z_boot    <- z[boot_idx, , drop = FALSE]
      w_boot    <- w[boot_idx]                     # subset weights to bootstrap sample

      # Pass weights to update()
      updated <- .ml_lm.fit(y = y_boot,
                            x = x_boot,
                            z = z_boot,
                            w = w_boot)

      if (updated$code %in% c(1L, 2L, 8L)) {
        if (progress) cat(cli::col_green("."))
        success[i] <- TRUE
        coef(updated)
      } else {
        if (progress) cat(cli::col_red("x"))
        success[i] <- FALSE
        rep(NA_real_, length(coef(object)))
      }
    }, error = function(e) {
      if (progress) cat(cli::col_red("x"))
      success[i] <- FALSE
      rep(NA_real_, length(coef(object)))
    })

    coef_matrix[i, ] <- res
  }

  if (progress) {
    cat("\n")
    cat(cli::col_blue(strrep("=", 52), "\n"))
  }

  # Final summary
  success_rate <- mean(success) * 100
  cli::cli_text("Bootstrapping finished — {round(success_rate, 1)}% of replications converged.")

  if (success_rate < 70) {
    cli::cli_warn("Low convergence rate — bootstrap results may be unreliable.")
  }

  # Variance from successful replications only
  valid_rows <- complete.cases(coef_matrix)
  vcov_boot  <- var(coef_matrix[valid_rows, , drop = FALSE])

  dimnames(vcov_boot) <- list(names(coef(object)), names(coef(object)))
  vcov_boot
}
