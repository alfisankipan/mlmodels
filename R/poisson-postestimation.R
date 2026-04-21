## PREDICT =====================================================================
# Remember that the parameter documentation is done at the generic predict for
# mlmodel class.

#' @details
#' ### ml_poisson prediction types
#'
#' The `type` argument controls what quantity is returned. In addition to
#' standard types, Poisson models support flexible probability requests
#' using the `P(...)` syntax.
#'
#' | Type          | Description                              | Notes |
#' |---------------|------------------------------------------|-------|
#' | `"link"`      | Linear predictor ( xb )                  | log-mean |
#' | `"response"`  | Expected count ( \code{mu} = \code{exp(xb)} )          | Default |
#' | `"mean"`      | Alias for `"response"`                   | - |
#' | `"mu"`        | Alias for `"response"`                   | - |
#' | `"fitted"`    | Alias for `"response"`                   | - |
#' | `P(k)`        | P(Y = k)                                 | Exact probability, k integer ≥ 0 |
#' | `P(,k)`       | P(Y ≤ k)                                 | Cumulative (lower tail) |
#' | `P(k,)`       | P(Y ≥ k)                                 | Survival (upper tail) |
#' | `P(a,b)`      | P(a ≤ Y ≤ b)                             | Interval probability, a ≤ b, a ≥ 0 |
#'
#' When `se.fit = TRUE`, standard errors are computed using the delta method
#' for all supported types.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @rdname predict.mlmodel
#' @export
predict.ml_poisson <- function(object,
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
  if (!inherits(object, "ml_poisson"))
    cli::cli_abort("`object` must be of class 'ml_poisson'.")
  
  # Send type to the parser, and then check others.
  
  parsed_type <- .predict_types_parsing(type)
  
  if(parsed_type$base_type != "prob")
  {
    parsed_type$base_type <- rlang::arg_match(type,
                                              c("response", "fitted", "mean", "mu", "link"))
  }
  
  # ── Prepare predictors (using hardhat) ─────────────────────────────
  if (is.null(newdata)) {
    val_predictors <- object$model$value$predictors
    sample_idx <- object$model$sample
    n_orig <- length(sample_idx)
  } else {
    val_bp <- object$model$value$blueprint
    forged <- hardhat::forge(newdata, blueprint = val_bp, outcomes = FALSE)
    val_predictors <- forged$predictors
    sample_idx <- rep(TRUE, nrow(val_predictors))
    n_orig <- nrow(val_predictors)
  }
  
  # ── Extract coefficients and compute linear predictors ─────────────
  beta <- coef(object)
  X <- as.matrix(val_predictors)
  xb <- as.vector(X %*% beta)
  mu <- exp(xb)
  if (parsed_type$base_type == "prob") {
    
    # Safety checks
    if (parsed_type$prob_type == "exact" && parsed_type$lower < 0) {
      cli::cli_abort("P(k) requires k >= 0.", call = NULL)
    }
    
    if (parsed_type$prob_type == "geq" && parsed_type$lower < 0) {
      cli::cli_abort("P(k,) requires k >= 0.", call = NULL)
    }
    
    if (parsed_type$prob_type == "interval") {
      if (parsed_type$lower < 0) {
        cli::cli_abort("Lower bound in P(a,b) must be >= 0.", call = NULL)
      }
      if (parsed_type$upper <= parsed_type$lower) {
        cli::cli_abort("Invalid interval P({parsed_type$lower},{parsed_type$upper}). Upper must be greater than lower.", 
                       call = NULL)
      }
    }
    
    out <- switch(parsed_type$prob_type,
                  "exact"    = dpois(parsed_type$lower, lambda = mu),
                  "leq"      = ppois(parsed_type$upper, lambda = mu),
                  "geq"      = ppois(parsed_type$lower - 1, lambda = mu, lower.tail = FALSE),
                  "interval" = ppois(parsed_type$upper, lambda = mu) - 
                    ppois(parsed_type$lower - 1, lambda = mu),
                  cli::cli_abort("Unknown probability type.", call = NULL)
    )
    
  } else {
    # Standard types (link, response, mean)
    out <- switch(parsed_type$base_type,
                  "link"     = xb,
                  "mean"     = ,
                  "fitted"   = ,
                  "mu"       = ,
                  "response" = mu,
                  cli::cli_abort("Unknown prediction type '{parsed_type$base_type}'.", 
                                 call = NULL))
  }
  
  # ── Align in-sample predictions to original data length ────────────
  if (is.null(newdata) && any(!sample_idx)) {
    full_out <- rep(NA_real_, n_orig)
    full_out[sample_idx] <- out
    out <- full_out
  }
  
  if (!se.fit)
  {
    res <- list(
      fit = out,
      se.fit = NULL
    )
    class(res) <- c("predict.ml_poisson", "predict.mlmodel")
    return(res)
  }
  
  n_obs <- length(mu)
  n_beta <- length(beta)
  if (parsed_type$base_type == "prob") {
    
    # Derivative of the probability w.r.t. μ
    dP_dmu <- switch(parsed_type$prob_type,
                     "exact" = dpois(parsed_type$lower, mu) * (parsed_type$lower / mu - 1),
                     "leq"   = - dpois(parsed_type$upper, mu),
                     "geq"   = if (parsed_type$lower == 0) rep(0, n_obs) 
                               else dpois(parsed_type$lower - 1L, mu),
                     "interval" = dpois(parsed_type$lower - 1L, mu) - 
                       dpois(parsed_type$upper, mu),
                     cli::cli_abort("Unknown prediction type '{parsed_type$base_type}'.", 
                                    call = NULL)
    )
    
    # Chain rule: d(prob)/dβ = d(prob)/dμ * dμ/dβ = dP_dmu * mu * X
    g <- dP_dmu * mu * X
    
  } else {
    g <- switch(parsed_type$base_type,
                "link"     = X,
                "mean"     = ,
                "fitted"   = ,
                "mu"       = ,
                "response" = mu * X,
                cli::cli_abort("Unknown prediction type '{parsed_type$base_type}'.", 
                               call = NULL))
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
  res <- list(
    fit = out,
    se.fit = se_fit
  )
  class(res) <- c("predict.ml_poisson", "predict.mlmodel")
  return(res)
}

## PRINT SUMMARY ===============================================================
#' @export
print.summary.ml_poisson <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!inherits(x, "summary.ml_poisson"))
    cli::cli_abort("`x` needs to be a `summary.ml_poisson` object.")
  
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
  
  cat("\nVariance type:", x$var_description)
  cat("\n---------------------------------------\n")
  
  old_pen <- getOption("scipen")
  options(scipen = .mlmodels_get_default("scipen"))
  
  # Determining the number of leading zeroes in the estimates and standard errors.
  format_coef <- x$coefficients
  checkvals <- abs(format_coef[, 1:2])
  checkvals <- checkvals[checkvals > 0 & !is.na(checkvals)]
  # - log10() gives us the number of leading zeroes (and decimal). Floor then
  # tells us where the first nonzero value is. Then max() gets us the maximum, so
  # using that + 2 guarantees that the number with more leading zeroes has two
  # nonzero decimals. Finally, max(digits, -) ensures that at least we have digits
  # decimal places (for estimations with few decimal places)
  num_digits <- if(length(checkvals) > 0) max(digits, max(floor(-log10(checkvals))) + 2)
  else digits
  # Rounding the estimate and standard errors.
  format_coef[, 1:2] <- round(format_coef[, 1:2], num_digits)
  
  # Capture the whole output of printCoefmat into a vector of strings.
  captured <- capture.output(printCoefmat(format_coef,
                                          digits = digits,
                                          signif.legend = TRUE))
  
  # Get the number of coefficients in each equation.
  k1 <- sum(grepl("^value::", rownames(x$coefficients)))
  
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
  # It seems as if there is an empty line between the coefficients and the legend
  # so we need to add one more line and start at k1+k2+3, to avoid that empty one.
  cat("---------------------------------------",
      captured[(k1+3):length(captured)],
      sep = "\n")
  
  options(scipen = old_pen)
  
  if (x$converged) {
    cat("---\n")
    cat("Number of observations:", x$nobs, 
        " Deg. of freedom: ", x$df.residual, "\n", sep = "")
    cat("Pseudo R-squared - Cor.Sq.: ",
        format(x$r.squared$cor, digits = digits),
        " McFadden: ", format(x$r.squared$mcfadden, digits = digits),
        "\n",
        sep = "")
    cat("AIC:", format(x$AIC, nsmall = 2, digits = digits + 1),
        " BIC:", format(x$BIC, nsmall = 2, digits = digits + 1), "\n")
    cat("\nCount Diagnostics:\n")
    cat("  Dispersion Ratio (Pearson):", format(x$ov, nsmall = 2, digits = digits + 1), "\n")
    cat("  Zeros - Observed:", x$zero$count, "Predicted:", round(x$zero$pred, 2), "\n")
  } else {
    cat("\nGoodness-of-fit statistics not available (model did not converge).\n")
  }
  
  invisible(x)
}

## SUMMARY =====================================================================
#' @rdname summary.mlmodel
#' @export
summary.ml_poisson <- function(object,
                               correlation = FALSE,
                               vcov = NULL,           # User-supplied variance matrix
                               vcov.type = "oim",
                               cl_var = NULL,
                               repetitions = 999,
                               seed = NULL,
                               progress = FALSE,
                               ...)
{
  if (!inherits(object, "ml_poisson"))
    cli::cli_abort("`object` must be a model of class 'ml_poisson'.")
  
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
  k_mean <- k_total
  
  # Start building the summary object
  s <- list()
  
  # Store the response's variable name
  s$response_name <- object$model$response_name
  
  # Call helper for variance type general description.
  s$var_description <- .vcov_description(vcov_mat)
  
  # Basic information
  s$logLik <- as.numeric(object$maximum %||% NA_real_)
  s$call           <- object$call                    # ← Now using root-level call
  s$formula        <- object$model$formula
  s$nobs           <- n
  s$df.residual    <- n - k_total
  s$converged      <- converged
  
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
    
    ll <- s$logLik
    s$AIC            <- -2 * ll + 2 * k_total
    s$BIC            <- -2 * ll + log(n) * k_total
    
    # overdispersion
    s$ov <- sum((y - yhat)^2 / yhat) / (n - k_total)
    
    # Null logLik
    y_bar <- mean(y)
    ll0 <- sum(y * log(y_bar) - y_bar - lfactorial(y))
    
    s$r.squared <- list(
      cor = cor(y, yhat)^2,
      mcfadden = 1 - ll / ll0
    )
    
    s$zero <- list(
      count = sum(y == 0),
      pred = sum(exp(-yhat))
    )
    
    
    if(usable_vcov)
    {
      # Joint significance test: only mean because no scale.
      idx_mean <- if (object$model$value$blueprint$intercept) 2:k_mean else 1:k_mean
      
      s$significance <- list(
        all  = waldtest(object, indices = idx_mean, vcov = vcov_mat),
        mean = NULL,
        scale = NULL
      )
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
    s$r.squared <- s$AIC <- s$BIC <- s$ov <- s$zero <- s$significance <- NULL
  }
  
  
  if(correlation && converged && usable_vcov)
    s$correlation <- cov2cor(vcov_mat)
  else
    s$correlation <- NULL
  
  s$model_type <- object$model$description
  class(s) <- c("summary.ml_poisson", "summary.mlmodel", "summary")
  s
}

## UPDATE ======================================================================
#' @rdname update.mlmodel
#' @export
update.ml_poisson <- function(object,
                              formula. = NULL,
                              data = NULL,
                              weights = NULL,
                              ...,
                              evaluate = TRUE)
{
  if (is.null(call <- object$call))
    cli::cli_abort("`object` does not contain a `call` component.", call = NULL)
  
  # Update value formula if explicitly requested
  if (!is.null(formula.)) {
    call$value <- update.formula(formula(object), formula.)
  }
  
  # Update data if supplied (what vcovBS passes)
  if (!is.null(data)) {
    call$data <- data
  }
  
  # Update weights ONLY if explicitly supplied and non-NULL
  if (!is.null(weights)) {
    call$weights <- weights
  }
  
  # Forward any extra arguments
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0) {
    for (arg in names(extras)) {
      call[[arg]] <- extras[[arg]]
    }
  }
  
  # Evaluate or return the call
  if (evaluate) {
    eval(call, envir = parent.frame())
  } else {
    call
  }
}