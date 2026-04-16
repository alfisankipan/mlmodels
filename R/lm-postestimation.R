## PREDICT =====================================================================
#' Predictions for ml_lm objects
#'
#' @details
#' ### ml_lm prediction types
#' The `type` argument controls what quantity is returned. Behavior differs
#' depending on whether the outcome was modeled in logs (\code{log(y)}).
#'
#' | Type | Normal (linear) case | Lognormal case (\code{log(y)}) | Notes |
#' |---------------------|---------------------------------------|------------------------------------------------------|-------|
#' | `link` | Linear predictor for scale (zd) | Linear predictor on log scale (μ_log) | Scale equation |
#' | `fitted` | xb (mean predictor) | xb (original log-scale predictor) | Mean equation |
#' | `response`, `mean` | xb (E\code{[y]}) | E\code{[y]} = exp(μ_log + σ²/2) - shift | Proper expected value on original scale |
#' | `median` | xb (same as mean) | exp(μ_log) - shift | Median of y |
#' | `sigma`, `sd` | σ (sd of ε) | σ (sd of \code{log(y)}) | On log scale |
#' | `sigma_y`, `sd_y` | same as `sigma` | sd(y) | Only meaningful in lognormal case |
#' | `variance` | σ² | σ² (variance of \code{log(y)}) | On log scale |
#' | `variance_y` | same as `variance` | Var(y) = exp(2μ_log + σ²)(exp(σ²) - 1) | Only meaningful in lognormal case |
#' | `zd` | Linear predictor for scale (zd) | Linear predictor for scale (zd) | Alias for `link` |
#'
#' When the outcome is log-transformed, `response` (or `mean`) returns the
#' correct lognormal expected value on the original scale of y. The `median`
#' is the simple exponential back-transform.
#'
#' Coefficient names in the object use the prefixes `value::` and `scale::`.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @rdname predict.mlmodel
#' @export
predict.ml_lm <- function(object,
                          newdata = NULL,
                          type = "response",
                          se.fit = FALSE,
                          vcov = NULL,
                          vcov.type = "oim",
                          cl_var = NULL,
                          repetitions = 999,
                          seed = NULL,
                          progress = FALSE,
                          ...)
{
  if (!inherits(object, "ml_lm"))
    cli::cli_abort("`object` must be of class 'ml_lm'.")
  type <- rlang::arg_match(type, c("response", "mean", "median", "fitted",
                                   "sigma", "sd", "variance",
                                   "sigma_y", "sd_y", "variance_y",
                                   "link", "zd"))
  # ── Prepare predictors (using hardhat) ─────────────────────────────
  is_heteroskedastic <- !is.null(object$model$scale_formula)
  log_info <- object$model$log_info$value
  if (is.null(newdata)) {
    val_predictors <- object$model$value$predictors
    scale_predictors <- if (is_heteroskedastic)
      object$model$scale$predictors
    else
      matrix(1, nrow = nrow(val_predictors), ncol = 1)
    sample_idx <- object$model$sample
    n_orig <- length(sample_idx)
  } else {
    val_bp <- object$model$value$blueprint
    forged <- hardhat::forge(newdata, blueprint = val_bp, outcomes = FALSE)
    val_predictors <- forged$predictors
    scale_predictors <- if (is_heteroskedastic)
    {
      scale_bp <- object$model$scale$blueprint
      forged <- hardhat::forge(newdata, blueprint = scale_bp, outcomes = FALSE)
      forged$predictors
    }
    else
      matrix(1, nrow = nrow(val_predictors), ncol = 1)
    sample_idx <- rep(TRUE, nrow(val_predictors))
    n_orig <- nrow(val_predictors)
  }
  # ── Extract coefficients and compute linear predictors ─────────────
  coefs <- coef(object)
  k_mean <- ncol(val_predictors)
  beta <- coefs[1:k_mean]
  delta <- coefs[(k_mean + 1):length(coefs)]
  X <- as.matrix(val_predictors)
  Z <- as.matrix(scale_predictors)
  xb <- as.vector(X %*% beta)
  zd <- as.vector(Z %*% delta)
  sigma <- exp(zd)
  # ── Point predictions ──────────────────────────────────────────────
  if (!log_info$is_log) {
    # Normal case
    out <- switch(type,
                  "link" = ,
                  "zd" = zd,
                  "fitted" = xb,
                  "response" = ,
                  "mean" = ,
                  "median" = xb,
                  "sigma" = ,
                  "sd" = ,
                  "sigma_y" = ,
                  "sd_y" = sigma,
                  "variance" = ,
                  "variance_y" = sigma^2)
  } else {
    # Lognormal case
    mu_log <- xb - log(log_info$multiplier)
    if (log_info$shift != 0) {
      cli::cli_warn(
        "Outcome was transformed as log(y + d). Back-transformation is approximate."
      )
    }
    Ey <- exp(mu_log + sigma^2 / 2) - log_info$shift
    My <- exp(mu_log) - log_info$shift
    Vy <- exp(2 * mu_log + sigma^2) * (exp(sigma^2) - 1)
    out <- switch(type,
                  "link" = mu_log,
                  "zd" = zd,
                  "fitted" = xb,
                  "response" = ,
                  "mean" = Ey,
                  "median" = My,
                  "sigma" = ,
                  "sd" = sigma,
                  "sigma_y" = ,
                  "sd_y" = sqrt(Vy),
                  "variance" = sigma^2,
                  "variance_y" = Vy)
  }
  
  # ── Align in-sample predictions to original data length ────────────
  if (is.null(newdata) && any(!sample_idx)) {
    full_out <- rep(NA_real_, n_orig)
    full_out[sample_idx] <- out
    out <- full_out
  }
  
  if (!se.fit) return(out)
  # ── Delta-method standard errors ───────────────────────────────────
  n_obs   <- length(xb)
  n_beta  <- length(beta)
  n_delta <- length(delta)
  g <- matrix(0, nrow = n_obs, ncol = n_beta + n_delta)
  
  if (!log_info$is_log) {
    # ── Normal case: pre-compute all possible gradients (cheap) ─────
    g_mean_beta   <- X
    g_mean_delta  <- matrix(0, n_obs, n_delta)
    g_link_beta   <- matrix(0, n_obs, n_beta)
    g_link_delta  <- if (is_heteroskedastic) Z else matrix(1, n_obs, 1)
    g_sigma_beta  <- matrix(0, n_obs, n_beta)
    g_sigma_delta <- if (is_heteroskedastic) sigma * Z else sigma
    g_var_beta    <- matrix(0, n_obs, n_beta)
    g_var_delta   <- if (is_heteroskedastic) 2 * sigma^2 * Z else 2 * sigma^2
    
    g[, 1:n_beta] <- switch(type,
                            "link" = ,
                            "zd"   = g_link_beta,
                            "fitted" = ,
                            "response" = ,
                            "mean" = ,
                            "median" = g_mean_beta,
                            "sigma" = ,
                            "sd" = ,
                            "sigma_y" = ,
                            "sd_y" = g_sigma_beta,
                            "variance" = ,
                            "variance_y" = g_var_beta,
                            g_mean_beta) # default = mean
    
    g[, (n_beta + 1):(n_beta + n_delta)] <- switch(type,
                                                   "link" = ,
                                                   "zd"   = g_link_delta,
                                                   "fitted" = ,
                                                   "response" = ,
                                                   "mean" = ,
                                                   "median" = g_mean_delta,
                                                   "sigma" = ,
                                                   "sd" = ,
                                                   "sigma_y" = ,
                                                   "sd_y" = g_sigma_delta,
                                                   "variance" = ,
                                                   "variance_y" = g_var_delta,
                                                   g_mean_delta) # default = mean
  } else {
    # Lognormal case — all moments computed unconditionally (independent of type)
    mu_log   <- xb - log(log_info$multiplier)
    sigma_i  <- sigma
    Ey       <- exp(mu_log + sigma_i^2 / 2) - log_info$shift # true mean
    My       <- exp(mu_log) - log_info$shift                 # median
    Vy       <- exp(2 * mu_log + sigma_i^2) * (exp(sigma_i^2) - 1)
    
    # Gradients (your original derivatives, now using the correct Ey/My/Vy)
    g_mean_beta   <- Ey * X
    g_mean_delta  <- if (is_heteroskedastic) Ey * sigma_i^2 * Z else Ey * sigma_i^2
    g_med_beta    <- My * X
    g_med_delta   <- matrix(0, n_obs, n_delta)
    g_var_beta    <- 2 * Vy * X
    g_var_delta   <- if (is_heteroskedastic)
      2 * sigma_i^2 * exp(2 * xb + sigma_i^2) * (2 * exp(sigma_i^2) - 1) * Z
    else
      2 * sigma_i^2 * exp(2 * xb + sigma_i^2) * (2 * exp(sigma_i^2) - 1)
    
    g[, 1:n_beta] <- switch(type,
                            "link" = ,
                            "fitted" = X,
                            "response" = ,
                            "mean" = g_mean_beta,
                            "median" = g_med_beta,
                            "variance_y" = g_var_beta,
                            "variance" = matrix(0, n_obs, n_beta),
                            matrix(0, n_obs, n_beta)) # default for sigma, sd, etc.
    
    g[, (n_beta + 1):(n_beta + n_delta)] <- switch(type,
                                                   "link" = ,
                                                   "fitted" = 0,
                                                   "response" = ,
                                                   "mean" = g_mean_delta,
                                                   "median" = 0,
                                                   "variance_y" = g_var_delta,
                                                   "variance" = if (is_heteroskedastic) Z else 1,
                                                   if (is_heteroskedastic) Z else 1) # default for sigma, sd, etc.
  }
  
  full_vcov <- .process_vcov(object,
                             vcov = vcov,
                             vcov.type = vcov.type,
                             cl_var = cl_var,
                             repetitions = repetitions,
                             seed = seed,
                             progress = progress)
  
  # ── Check for unusable variance matrix ─────────────────────────────
  if (any(!is.finite(full_vcov)) || any(is.na(full_vcov))) {
    cli::cli_warn(
      c("Variance matrix is unusable (contains NAs or non-finite values).",
        "i" = "This usually happens with bootstrap when constraints are present.",
        "i" = "Standard errors will be returned as NA.")
    )
    se_fit <- rep(NA_real_, length(out))
  } else {
    # ── Delta-method standard errors ─────────────────────────────────
    se_fit <- sqrt(rowSums(g * (g %*% full_vcov)))
  }
  
  # Align in-sample SEs
  if (is.null(newdata) && any(!sample_idx)) {
    full_se <- rep(NA_real_, n_orig)
    full_se[sample_idx] <- se_fit
    se_fit <- full_se
  }
  list(fit = out, se.fit = se_fit)
}

## RESIDUALS ===================================================================
#' Extract Residuals from ml_lm objects
#'
#' Returns the residuals (observed minus fitted values) from the value equation.
#'
#' @param object A fitted `ml_lm` object.
#' @param ... Not currently implemented.
#'
#' @return A numeric vector of residuals with length equal to the number of
#'   observations used in estimation.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @export
residuals.ml_lm <- function(object, ...)
{
  if (!inherits(object, "ml_lm"))
    cli::cli_abort("`object` must be a model of class 'ml_lm' (from ml_lm).",
                   call = NULL)
  if (!is.null(object$model$residuals))
    return(object$model$residuals)

  # Compute on the fly if not stored
  y_clean <- object$model$value$outcomes[[1]]          # already clean (length = n_used)

  # Get fitted values only for the used observations
  fit_full  <- predict(object, type = "response")
  fit_clean <- fit_full[object$model$sample]           # length = n_used

  return(y_clean - fit_clean)
}

## PRINT SUMMARY ===============================================================
#' Print the summary statistics from an `ml_lm` estimation.
#'
#' @param x An object of class `summary.ml_lm`, usually from [summary.ml_lm()][mlmodels::summary.ml_lm()].
#'
#' @param digits A numeric scalar with the number of decimal places to use in
#'    the statistics
#'
#' @param ... Currently not implemented.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @export
print.summary.ml_lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!inherits(x, "summary.ml_lm"))
    cli::cli_abort("`x` needs to be a `summary.ml_lm` object.")
  
  cat("\nMaximum Likelihood Model\n")
  cat(" Type:", x$model_type, "\n")
  cat("---------------------------------------\n")
  
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  
  if (!x$converged) {
    cat("WARNING: Model did NOT converge!\n")
    cat("Convergence code:", x$code %||% "???", "-", x$message %||% "", "\n\n")
  } else {
    # Log-Likelihood + Joint tests (only when converged)
    cat("Log-Likelihood:", format(x$logLik, nsmall = 2, digits = digits + 1), "\n\n")
    cat("Wald significance tests:\n")
    
    any_test_printed <- FALSE
    for (test in c("all", "mean", "scale")) {
      w <- x$significance[[test]]
      if (is.null(w) || isTRUE(w$singular) || !is.finite(w$pval)) {
        next   # skip silently (happens in homoskedastic case or useless variance)
      }
      any_test_printed <- TRUE
      p_str <- if (w$pval < 1e-8) "< 1e-8" else sprintf("%.4f", w$pval)
      cat(sprintf(" %s: Chisq(%d) = %.3f, Pr(>Chisq) = %s\n",
                  tools::toTitleCase(test), w$df, w$waldstat, p_str))
    }
    
    if (!any_test_printed) {
      cat(" Tests were not computable (singular or not finite variance).\n")
    }
  }
  
  if(!is.null(x$vcov.type))
    vcov_type <- switch (x$vcov.type,
                         "oim" = "Original Information Matrix",
                         "opg" = "Outer Product of Gradients (BHHH)",
                         "robust" = if(is.null(x$vcov.cluster)) "Robust" else "Cluster-Robust",
                         "boot" = if(is.null(x$vcov.cluster)) paste0("Bootstrap (",
                                                                     x$boot.reps, " repetitions)")
                             else paste0("Cluster Bootstrap (",
                                         x$boot.reps, " repetitions)"),
                         "jack" = if(is.null(x$vcov.cluster)) "Jackknife" else "Cluster Jackknife",
                         x$vcov.type
    )
  else
    vcov_type <- "User Supplied (Unknown)"
  
  cat("\nVariance type:", vcov_type)
  if (!is.null(x$vcov.cluster)) {
    cat(" | Clusters:", x$vcov.cluster$n_cluster)
    if (!is.null(x$vcov.cluster$var_name))
      cat(" (", x$vcov.cluster$var_name, ")", sep = "")
  }
  cat("\n---------------------------------------\n")
  
  old_pen <- getOption("scipen")
  options(scipen = .mlmodels_get_default("scipen"))
  
  # Capture the whole output of printCoefmat into a vector of strings.
  captured <- capture.output(printCoefmat(x$coefficients,
                                          digits = digits,
                                          signif.legend = TRUE))
  
  # Get the number of coefficients in each equation.
  k1 <- sum(grepl("^value::", rownames(x$coefficients)))
  k2 <- sum(grepl("^scale::", rownames(x$coefficients)))
  
  # The first row is the header of the table, have to indent it to align it with
  # the coefficients after.
  cat("  ", captured[1],"\n")
  # Value header. Depends if we have the dependent's variable name.
  val_head <- if (!is.null(x$response_name) && 
                  nzchar(trimws(x$response_name))) {
    paste0("Value (", trimws(x$response_name), "):")
  } else {
    "Value:"
  }
  cat(val_head)
  cat("  ", captured[2:(k1+1)],
      sep = "\n  ")
  cat("Scale (log(sigma)):")
  cat("  ", captured[(k1+2):(k1+k2+1)],
      sep = "\n  ")
  # It seems as if there is an empty line between the coefficients and the legend
  # so we need to add one more line and start at k1+k2+3, to avoid that empty one.
  cat("---------------------------------------",
      captured[(k1+k2+3):length(captured)],
      sep = "\n")
  
  options(scipen = old_pen)
  
  if (x$converged) {
    cat("---\n")
    cat("Number of observations:", x$nobs, "\n")
    if (!is.null(x$df.residual))
      cat("Residual degrees of freedom:", x$df.residual, "\n")
    if (!is.null(x$sigma))
      cat("Residual standard error (sigma):", format(x$sigma, digits = digits), "\n")
    cat("Multiple R-squared: ", format(x$r.squared, digits = digits),
        " Adjusted R-squared: ", format(x$adj.r.squared, digits = digits), "\n", sep = "")
    cat("AIC:", format(x$AIC, nsmall = 2, digits = digits + 1),
        " BIC:", format(x$BIC, nsmall = 2, digits = digits + 1), "\n")
  } else {
    cat("\nGoodness-of-fit statistics not available (model did not converge).\n")
  }
  
  invisible(x)
}

## SUMMARY =====================================================================
#' @rdname summary.mlmodel
#' @export
summary.ml_lm <- function(object,
                          correlation = FALSE,
                          vcov = NULL,           # User-supplied variance matrix
                          vcov.type = "oim",
                          cl_var = NULL,
                          repetitions = 999,
                          seed = NULL,
                          progress = FALSE,
                          ...)
{
  if (!inherits(object, "ml_lm"))
    cli::cli_abort("`object` must be a model of class 'ml_lm'.")

  converged <- object$code %in% c(0, 1, 2, 8)
  # Get variance-covariance matrix once
  vcov_mat <- .process_vcov(object,
                            vcov = vcov,
                            vcov.type   = vcov.type,
                            cl_var      = cl_var,
                            repetitions = repetitions,
                            seed        = seed,
                            progress    = progress)

  # Check if the variance matrix is usable
  usable_vcov <- TRUE
  if (any(!is.finite(vcov_mat)) || any(is.na(vcov_mat))) {
    usable_vcov <- FALSE
    warn_msg <- "Variance matrix is not usable (contains NAs or non-finite values)."
  } else if (!.is_invertible(vcov_mat)) {
    usable_vcov <- FALSE
    warn_msg <- "Variance matrix is not invertible (likely singular or nearly singular)."
  }
  if(!usable_vcov)
    cli::cli_warn(
      c(warn_msg,
        "i" = "This can happen with bootstrap under constraints or in models with high collinearity.",
        "i" = "Joint significance tests and correlation matrix will be skipped.",
        "i" = "Consider using `vcov.type = 'robust'` for more stable results.")
    )

  n <- object$model$n_used
  k_total <- length(coef(object))
  k_scale <- ncol(object$model$scale$predictors)
  is_heteroskedastic <- !is.null(object$model$scale_formula)
  k_mean <- k_total - k_scale

  # Start building the summary object
  s <- list()
  
  # Store the response's variable name
  s$response_name <- object$model$response_name
  # Read metadata from attributes (preferred source)
  s$vcov.type <- attr(vcov_mat, "vcov.type")

  # Clustered variance handling
  if (!is.null(attr(vcov_mat, "clustered")) && attr(vcov_mat, "clustered")) {
    s$vcov.cluster <- .vcov_cluster_info(object, attr(vcov_mat, "cluster.var"))
  } else if (vcov.type %in% c("cluster", "robust") && !is.null(cl_var)) {
    # Fallback when attributes are missing (should rarely happen)
    s$vcov.cluster <- .vcov_cluster_info(object, cl_var)
  } else {
    s$vcov.cluster <- NULL
  }
  
  # Boostrap variance handling.
  if(s$vcov.type == "boot")
    s$boot.reps <- attr(vcov_mat, "rep")

  # Basic information
  s$logLik <- as.numeric(object$maximum %||% NA_real_)
  s$call           <- object$call                    # ← Now using root-level call
  s$formula        <- object$model$formula
  s$scale_formula  <- object$model$scale_formula
  s$nobs           <- n
  s$df.residual    <- n - k_mean
  s$converged      <- converged
  s$is_heteroskedastic <- is_heteroskedastic

  # Coefficient table
  se <- sqrt(diag(vcov_mat))

  s$coefficients <- cbind(
    Estimate   = coef(object),
    `Std. Error` = se,
    `z value`  = coef(object) / se,
    `Pr(>|z|)` = 2 * pnorm(abs(coef(object) / se), lower.tail = FALSE)
  )

  # Stats if converged
  if (converged) {
    y <- object$model$value$outcomes[[1]]
    yhat <- object$model$fitted.values

    s$r.squared      <- cor(y, yhat)^2
    s$adj.r.squared  <- 1 - (1 - s$r.squared) * (n - 1) / (n - k_mean)

    ll <- s$logLik
    s$AIC            <- -2 * ll + 2 * k_total
    s$BIC            <- -2 * ll + log(n) * k_total

    if (!is_heteroskedastic) {
      s$sigma <- exp(coef(object)[k_mean + 1])
    }


    if(usable_vcov)
    {
      # Joint significance tests (reuse vcov_mat)
      idx_mean <- if (object$model$value$blueprint$intercept) 2:k_mean else 1:k_mean

      if (is_heteroskedastic) {
        idx_scale <- if (object$model$scale$blueprint$intercept) (k_mean + 2):k_total
        else (k_mean + 1):k_total

        s$significance <- list(
          all  = waldtest(object, indices = c(idx_mean, idx_scale), vcov = vcov_mat),
          mean = waldtest(object, indices = idx_mean, vcov = vcov_mat),
          scale = waldtest(object, indices = idx_scale, vcov = vcov_mat)
        )
      } else {
        s$significance <- list(
          all  = waldtest(object, indices = idx_mean, vcov = vcov_mat),
          mean = NULL,
          scale = NULL
        )
      }
    }
    else
    {
      s$significance <- list(
        all  = NULL,
        mean = NULL,
        scale = NULL
      )
    }

  } else {
    s$r.squared <- s$adj.r.squared <- s$AIC <- s$BIC <- s$sigma <- s$significance <- NULL
  }
  
  
  if(correlation && converged && usable_vcov)
    s$correlation <- cov2cor(vcov_mat)
  else
    s$correlation <- NULL
  
  s$model_type <- object$model$description
  class(s) <- c("summary.ml_lm", "summary.mlmodel", "summary")
  s
}

## UPDATE ======================================================================
#' @rdname update.mlmodel
#' @export
update.ml_lm <- function(object,
                         formula. = NULL,
                         scale. = NULL,
                         data = NULL,
                         weights = NULL,
                         ...,
                         evaluate = TRUE)
{
  if (is.null(call <- object$call))
    cli::cli_abort("`object` does not contain a `call` component.", call = NULL)
  
  # Update value formula if explicitly requested (rare for bootstrap)
  if (!is.null(formula.)) {
    call$value <- update.formula(formula(object), formula.)
  }
  
  # Update scale formula if explicitly requested
  if (!is.null(scale.)) {
    if (identical(scale., ~1) || identical(scale., ~0)) {
      call$scale <- NULL
    } else {
      call$scale <- scale.
    }
  }
  # If scale. is not supplied → keep the original scale formula (crucial for heteroskedastic models)
  
  # Update data (the main thing vcovBS passes)
  if (!is.null(data)) {
    call$data <- data
  }
  
  # Update weights — this is the important part you asked about
  if (!is.null(weights)) {
    call$weights <- weights          # Yes, this correctly adds/replaces weights
  }
  # If weights is NULL (most common case from vcovBS), we do NOTHING → original weights stay in the call
  
  # Forward any other arguments the user might pass explicitly
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0) {
    for (arg in names(extras)) {
      call[[arg]] <- extras[[arg]]
    }
  }
  
  # Evaluate or return the updated call
  if (evaluate) {
    eval(call, envir = parent.frame())
  } else {
    call
  }
}
