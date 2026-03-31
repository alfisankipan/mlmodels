## PREDICT ---------------------------------------------------------------------
#' Predictions for ml_lm objects
#'
#' @param object A fitted `ml_lm` object.
#' @param newdata Optional data frame for out-of-sample predictions.
#' @param type Character string indicating what to predict. See **Details**.
#' @param se.fit Logical. If `TRUE`, also return standard errors (delta method).
#' @param vcov Optional user-supplied variance-covariance matrix.
#' @param vcov.type Type of variance-covariance matrix. See [vcov.mlmodel()].
#' @param cl_var Clustering variable (name or vector).
#' @param repetitions Number of bootstrap replications when `vcov.type = "boot"`.
#' @param seed Random seed for reproducibility.
#' @param progress Logical. Show bootstrap progress bar? Default is `FALSE` in
#' higher-level functions.
#' @param ... Not currently used.
#'
#' @return If `se.fit = FALSE` (default), a numeric vector of predictions.
#' If `se.fit = TRUE`, a list with components `fit` and `se.fit`.
#'
#' @details
#' The `type` argument controls what quantity is returned. Behavior differs
#' depending on whether the outcome was modeled in logs (log(y)).
#'
#' ### Prediction types
#'
#' | Type | Normal (linear) case | Lognormal case (log(y)) | Notes |
#' |---------------------|---------------------------------------|------------------------------------------------------|-------|
#' | `link` | Linear predictor for scale (zd) | Linear predictor on log scale (μ_log) | Scale equation |
#' | `fitted` | xb (mean predictor) | xb (original log-scale predictor) | Mean equation |
#' | `response`, `mean` | xb (E\code{[y]}) | E\code{[y]} = exp(μ_log + σ²/2) - shift | Proper expected value on original scale |
#' | `median` | xb (same as mean) | exp(μ_log) - shift | Median of y |
#' | `sigma`, `sd` | σ (sd of ε) | σ (sd of log(y)) | On log scale |
#' | `sigma_y`, `sd_y` | same as `sigma` | sd(y) | Only meaningful in lognormal case |
#' | `variance` | σ² | σ² (variance of log(y)) | On log scale |
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
                  "link" = zd,
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
  full_vcov <- get_vcov(object,
                        vcov = vcov,
                        vcov.type = vcov.type,
                        cl_var = cl_var,
                        repetitions = repetitions,
                        seed = seed,
                        progress = progress)
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

  se_fit <- sqrt(rowSums(g * (g %*% full_vcov)))

  # Align in-sample SEs
  if (is.null(newdata) && any(!sample_idx)) {
    full_se <- rep(NA_real_, n_orig)
    full_se[sample_idx] <- se_fit
    se_fit <- full_se
  }
  list(fit = out, se.fit = se_fit)
}

## FITTED VALUES ---------------------------------------------------------------
#' Extract Fitted Values from ml_lm objects
#'
#' Returns the fitted values (predicted conditional mean) from the value equation.
#'
#' @param object A fitted `ml_lm` object.
#' @param ... Not currently used.
#'
#' @return A numeric vector of fitted values with length equal to the number of
#'   observations used in estimation.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @export
fitted.ml_lm <- function(object, ...)
{
  if (!inherits(object, "ml_lm"))
    cli::cli_abort("`object` must be a model of class 'ml_lm' (from ml_lm).",
                   call = NULL)


  if (!is.null(object$model$fitted.values))
    return(object$model$fitted.values)

  # Fallback to predict if not stored
  predict(object, type = "response")
}

## RESIDUALS -------------------------------------------------------------------
#' Extract Residuals from ml_lm objects
#'
#' Returns the residuals (observed minus fitted values) from the value equation.
#'
#' @param object A fitted `ml_lm` object.
#' @param ... Not currently used.
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

## SUMMARY ---------------------------------------------------------------------
#' Summary for ml_lm objects
#'
#' @param object A fitted model object of class `"ml_lm"`.
#' @param correlation Logical. Should the correlation matrix of the estimated
#'   parameters be included in the output? Default is `FALSE`.
#' @param vcov Optional user-supplied variance-covariance matrix. If provided,
#'   it will be used instead of computing one internally.
#' @param vcov.type Character string specifying the type of variance-covariance
#'   matrix to use. One of `"oim"` (default), `"robust"`, `"opg"`, `"cluster"`,
#'   or `"boot"`. See [vcov.mlmodel()] for details.
#' @param cl_var Character string or vector. Name of the clustering variable
#'   or the vector itself. Only used when `vcov.type = "cluster"` or when
#'   `vcov.type = "boot"` with clustering.
#' @param repetitions Integer. Number of bootstrap replications when
#'   `vcov.type = "boot"`. Default is 999.
#' @param seed Integer. Random seed for reproducibility when `vcov.type = "boot"`.
#'   If `NULL`, a random seed is generated.
#' @param progress Logical. Should a progress bar be displayed during
#'   bootstrapping? Default is `FALSE` (silent) when called from `summary()`.
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class `"summary.ml_lm"` containing:
#'   - `call`: the original call
#'   - `coefficients`: coefficient table with Estimate, Std. Error, z value, Pr(>|z|)
#'   - `vcov.type`: type of variance-covariance matrix used
#'   - `vcov.cluster`: clustering information (if applicable)
#'   - `logLik`, `AIC`, `BIC`, `r.squared`, `adj.r.squared` (if converged)
#'   - `sigma` (if homoskedastic)
#'   - `significance`: joint Wald tests for overall, mean, and scale equations
#'
#' @author Alfonso Sanchez-Penalver
#'
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

  converged <- object$code %in% c(1L, 2L, 8L)
  # Get variance-covariance matrix once
  vcov_mat <- get_vcov(object,
                       vcov = vcov,
                       vcov.type   = vcov.type,
                       cl_var      = cl_var,
                       repetitions = repetitions,
                       seed        = seed,
                       progress    = progress)

  n <- object$model$n_used
  k_total <- length(coef(object))
  k_scale <- ncol(object$model$scale$predictors)
  is_heteroskedastic <- !is.null(object$model$scale_formula)
  k_mean <- k_total - k_scale

  # Start building the summary object
  s <- list()

  # Basic information
  s$call           <- object$call                    # ← Now using root-level call
  s$formula        <- object$model$formula
  s$scale_formula  <- object$model$scale_formula
  s$nobs           <- n
  s$df.residual    <- n - k_mean
  s$converged      <- converged
  s$vcov.type      <- vcov.type
  s$is_heteroskedastic <- is_heteroskedastic

  if (vcov.type == "cluster" && !is.null(cl_var))
    s$vcov.cluster <- vcov_cluster_info(object, cl_var)

  # Coefficient table
  se <- se.mlmodel(object,
                   vcov = vcov_mat)

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

    ll <- as.numeric(logLik(object))
    s$logLik         <- ll
    s$AIC            <- -2 * ll + 2 * k_total
    s$BIC            <- -2 * ll + log(n) * k_total

    if (!is_heteroskedastic) {
      s$sigma <- exp(coef(object)[k_mean + 1])
    }

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
  } else {
    s$r.squared <- s$adj.r.squared <- s$AIC <- s$BIC <- s$sigma <- s$significance <- NULL
  }

  s$model_type <- if (is_heteroskedastic) {
    "Heteroskedastic Gaussian Linear Model"
  } else {
    "Homoskedastic Gaussian Linear Model"
  }

  if (correlation) {
    s$correlation <- cov2cor(vcov_mat)
  }

  class(s) <- c("summary.ml_lm", "summary")
  s
}

## PRINT SUMMARY ---------------------------------------------------------------
#' Print the summary stastics from an `ml_lm` estimation.
#'
#' @param x An object of class `summary.ml_lm`, usually from [summary.ml_lm()][mlmodels::summary.ml_lm()].
#'
#' @param digits A numeric scalar with the number of decimal places to use in
#'    the statistics
#'
#' @param ... Currently not used.
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
  }

  # Log-Likelihood + Joint tests
  if (x$converged) {
    cat("Log-Likelihood:", format(x$logLik, nsmall = 2, digits = digits + 1), "\n\n")
    cat("Joint significance tests:\n")
    w <- x$significance$all
    p_str <- if (isTRUE(w$singular)) "<singular>" else
      if (w$pval < 1e-8) "< 1e-8" else sprintf("%.4f", w$pval)
    cat(sprintf(" Overall: Chisq(%d) = %.3f, Pr(>Chisq) = %s\n",
                w$df, w$waldstat, p_str))

    if (!is.null(x$significance$mean)) {
      w <- x$significance$mean
      p_str <- if (isTRUE(w$singular)) "<singular>" else
        if (w$pval < 1e-8) "< 1e-8" else sprintf("%.4f", w$pval)
      cat(sprintf(" Mean: Chisq(%d) = %.3f, Pr(>Chisq) = %s\n",
                  w$df, w$waldstat, p_str))

      w <- x$significance$scale
      p_str <- if (isTRUE(w$singular)) "<singular>" else
        if (w$pval < 1e-8) "< 1e-8" else sprintf("%.4f", w$pval)
      cat(sprintf(" Scale: Chisq(%d) = %.3f, Pr(>Chisq) = %s\n",
                  w$df, w$waldstat, p_str))
    }
  }

  cat("\nVariance type:", x$vcov.type)
  if (!is.null(x$vcov.cluster)) {
    cat(" | Clusters:", x$vcov.cluster$n_cluster)
    if (!is.null(x$vcov.cluster$var_name))
      cat(" (", x$vcov.cluster$var_name, ")", sep = "")
  }
  cat("\n---------------------------------------\n")

  old_pen <- getOption("scipen")
  options(scipen = 2)

  # ── Split coefficients into Value and Scale equations ─────────────────────
  value_coefficients <- x$coefficients[grepl("^value::", rownames(x$coefficients)), , drop = FALSE]
  scale_coefficients <- x$coefficients[grepl("^scale::", rownames(x$coefficients)), , drop = FALSE]

  # Strip prefixes for clean display
  # rownames(value_coefficients) <- sub("^value::", "", rownames(value_coefficients))
  # rownames(scale_coefficients) <- sub("^scale::", "", rownames(scale_coefficients))

  # Print Value equation
  cat("Value equation:\n")
  cat("----------------\n")
  cat("  ")
  cat(capture.output(printCoefmat(value_coefficients,
                                  digits = digits,
                                  signif.legend = TRUE)),
      sep = "\n  ")
  # printCoefmat(value_coefficients, digits = digits, signif.legend = TRUE)

  # Print Scale equation
  cat("\nScale equation:\n")
  cat("----------------\n")
  cat("  ")
  cat(capture.output(printCoefmat(scale_coefficients,
                                  digits = digits,
                                  signif.legend = TRUE)),
      sep = "\n  ")

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

# FUNCTIONS TO MAKE ML_LM WORK WITH SANDWICH -----------------------------------

#' @export
terms.ml_lm <- function(x, ...) {
  # Return NULL or the terms from the value equation if available
  if (!is.null(x$model$value$terms)) {
    x$model$value$terms
  } else {
    NULL
  }
}

#' Update an ml_lm model
#'
#' This method allows updating the formulas or other arguments of a fitted
#' `ml_lm` model. It is primarily used by bootstrap methods such as
#' `sandwich::vcovBS()`.
#'
#' @param object An object of class `"ml_lm"`.
#' @param formula. An optional updated formula for the **value** (mean) equation.
#' @param scale. An optional updated formula for the **scale** equation.
#'   Use `NULL` to remove the scale equation (i.e., make the model homoskedastic).
#' @param data Optional new data frame to use.
#' @param ... Other arguments to be passed to `ml_lm()` (e.g. `data`, `weights`,
#'   `noint_value`, `noint_scale`, etc.).
#' @param evaluate Logical. If `TRUE` (default), the updated call is evaluated.
#'   If `FALSE`, the updated call is returned as a language object.
#'
#' @return
#' If `evaluate = TRUE`, returns a new fitted `ml_lm` object.
#' If `evaluate = FALSE`, returns the updated call as a language object.
#'
#' @export
update.ml_lm <- function(object,
                         formula.,
                         scale. = NULL,
                         weights = NULL,
                         ...,
                         evaluate = TRUE)
{
  if (is.null(call <- object$call))
    cli::cli_abort("need an object with call component", call = NULL)

  if (!missing(formula.)) {
    call$value <- update.formula(formula(object), formula.)
  }

  if (!is.null(scale.)) {
    call$scale <- if (is.null(call$scale)) scale. else update.formula(call$scale, scale.)
  }

  # Forward weights if supplied
  if (!is.null(weights)) {
    call$weights <- weights
  }

  # Forward any other arguments
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0) {
    for (arg in names(extras)) {
      call[[arg]] <- extras[[arg]]
    }
  }

  if (evaluate) {
    eval(call, parent.frame())
  } else {
    call
  }
}
