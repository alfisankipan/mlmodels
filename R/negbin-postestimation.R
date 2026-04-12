## PREDICT =====================================================================
#' @author Alfonso Sanchez-Penalver
#'
#' @rdname predict.mlmodel
#' @export
predict.ml_negbin <- function(object,
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
  if (!inherits(object, "ml_negbin"))
    cli::cli_abort("`object` must be of class 'ml_negbin'.")
  
  # Send type to the parser, and then check others.
  
  parsed_type <- .predict_types_parsing(type)
  
  if(parsed_type$base_type != "prob")
  {
    parsed_type$base_type <- rlang::arg_match(type,
                                              c("response", "fitted", "mean", "link", "zd", "alpha", "variance", "sd", "sigma"))
  }
  
  # ── Prepare predictors (using hardhat) ─────────────────────────────
  is_heteroskedastic <- !is.null(object$model$scale_formula)
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
  xb <- as.vector(X %*% beta)
  Z <- as.matrix(scale_predictors)
  zd <- as.vector(Z %*% delta)
  mu <- exp(xb)
  alpha <- exp(zd)
  
  is_nb2 <- (object$model$dispersion == "NB2")
  
  var <- if(is_nb2) mu + alpha * mu^2 else mu + alpha * mu
  
  # Probability
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
                  "exact"    = if(is_nb2) dnbinom(parsed_type$lower,
                                                  size = 1 / alpha,
                                                  mu = mu)
                               else dnbinom(parsed_type$lower,
                                            size = mu / alpha,
                                            prob = 1 / (1 + alpha)),
                  "leq"      = if(is_nb2) pnbinom(parsed_type$upper,
                                                  size = 1 / alpha,
                                                  mu = mu)
                               else pnbinom(parsed_type$upper,
                                            size = mu / alpha,
                                            prob = 1 / (1 + alpha)),
                  "geq"      = if(is_nb2) pnbinom(parsed_type$lower - 1,
                                                  size = 1 / alpha,
                                                  mu = mu,
                                                  lower.tail = FALSE)
                               else pnbinom(parsed_type$lower - 1,
                                            size = mu / alpha,
                                            prob = 1 / (1 + alpha),
                                            lower.tail = FALSE),
                  "interval" = if(is_nb2) pnbinom(parsed_type$upper,
                                                  size = 1 / alpha,
                                                  mu = mu) -
                                          pnbinom(parsed_type$lower - 1,
                                                  size = 1 / alpha,
                                                  mu = mu)
                               else pnbinom(parsed_type$upper,
                                            size = mu / alpha,
                                            prob = 1 / (1 + alpha)) -
                                    pnbinom(parsed_type$lower - 1,
                                            size = mu / alpha,
                                            prob = 1 / (1 + alpha)),
                  cli::cli_abort("Unknown probability type.", call = NULL)
    )
    
  } else {
    # Standard types (link, response, mean)
    out <- switch(parsed_type$base_type,
                  "link"     = xb,
                  "mean"     = ,
                  "fitted"   = ,
                  "response" = mu,
                  "alpha"    = alpha,
                  "zd"       = zd,
                  "variance" = var,
                  "sd"       = ,
                  "sigma"    = sqrt(var),
                  cli::cli_abort("Unknown prediction type '{parsed_type$base_type}'.", 
                                 call = NULL))
  }
  
  # ── Align in-sample predictions to original data length ────────────
  if (is.null(newdata) && any(!sample_idx)) {
    full_out <- rep(NA_real_, n_orig)
    full_out[sample_idx] <- out
    out <- full_out
  }
  
  if (!se.fit) return(out)
  
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
  list(fit = out, se.fit = se_fit)
}



## PRINT SUMMARY ===============================================================
#' @export
print.summary.ml_negbin <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!inherits(x, "summary.ml_negbin"))
    cli::cli_abort("`x` needs to be a `summary.ml_negbin` object.")
  
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
  cat("Scale (log(alpha)):")
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
#' Summary for ml_negbin objects
#'
#' @param object A fitted model object of class `"ml_negbin"`.
#' @param correlation Logical. Should the correlation matrix of the estimated
#'   parameters be included in the output? Default is `FALSE`. If `TRUE` the
#'   correlation matrix will be computed, and stored in the `'summary.ml_lm'`
#'   object the function returns.
#' @param vcov Optional user-supplied variance-covariance matrix. If provided,
#'   it will be used instead of computing one internally.
#' @param vcov.type Character string specifying the type of variance-covariance
#'   matrix to use. See [vcov][mlmodels::vcov.mlmodel].
#' @param cl_var Character string or vector. Name of the clustering variable
#'   or the vector itself. See [vcov][mlmodels::vcov.mlmodel].
#' @param repetitions Integer. Number of bootstrap replications when
#'   `vcov.type = "boot"`. Default is 999.
#' @param seed Integer. Random seed for reproducibility when `vcov.type = "boot"`.
#'   If `NULL`, a random seed is generated.
#' @param progress Logical. Should a progress bar be displayed during
#'   bootstrapping or jackknifing? Default is `FALSE` (silent).
#' @param ... Further arguments passed to methods.
#'
#' @details
#' Coefficient names in the fitted object use the prefixes `value::` for
#' consistency with our other estimators that model the scale as well.
#' 
#' @return An object of class `c("summary.ml_negbin", "summary.mlmodel", "summary")`.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @export
summary.ml_negbin <- function(object,
                              correlation = FALSE,
                              vcov = NULL,           # User-supplied variance matrix
                              vcov.type = "oim",
                              cl_var = NULL,
                              repetitions = 999,
                              seed = NULL,
                              progress = FALSE,
                              ...)
{
  if (!inherits(object, "ml_negbin"))
    cli::cli_abort("`object` must be a model of class 'ml_negbin'.")
  
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
  s$df.residual    <- n - k_total
  s$converged      <- converged
  s$is_heteroskedastic <- is_heteroskedastic
  s$dispersion     <- object$model$dispersion
  
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
    alphahat <- object$model$fitted.alpha
    
    ll <- s$logLik
    s$AIC            <- -2 * ll + 2 * k_total
    s$BIC            <- -2 * ll + log(n) * k_total
    
    # overdispersion
    s$ov <- sum((y - yhat)^2 / yhat) / (n - k_total)
    
    # Null logLik
    y_bar <- mean(y)
    ll0 <- object$model$ll0
    
    s$r.squared <- list(
      cor = cor(y, yhat)^2,
      mcfadden = 1 - ll / ll0
    )
    
    pred_zero <- if(object$model$dispersion == "NB2")
      sum((1 + alphahat * yhat)^(-1 / alphahat))
    else
      sum((1 / (1 + alphahat))^(yhat / alphahat))
    
    s$zero <- list(
      count = sum(y == 0),
      pred = pred_zero
    )
    
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
    s$r.squared <- s$AIC <- s$BIC <- s$ov <- s$zero <- s$significance <- NULL
  }
  
  
  if(correlation && converged && usable_vcov)
    s$correlation <- cov2cor(vcov_mat)
  else
    s$correlation <- NULL
  
  
  s$model_type <- if (is_heteroskedastic) {
    paste0("Heteroskedastic Negative Binomial (", object$model$dispersion, ") Model")
  } else {
    paste0("Homoskedastic Negative Binomial (", object$model$dispersion, ") Model")
  }
  
  class(s) <- c("summary.ml_negbin", "summary.mlmodel", "summary")
  s
}


## UPDATE ======================================================================
#' Update method for ml_negbin objects
#' 
#' @param object An `ml_negbin` estimation object.
#' @param formula. The formula of the value equation (optional).
#' @param scale. The formula of the scale equation (optional).
#' @param data A data.frame with the data to do the estimation (optional).
#' @param weights A vector with the weights (optional).
#' @param ... Currently not implemented.
#' @param evaluate Should the updated call be evaluated? Defaults to `TRUE`.
#'
#' @details
#' Re-evaluates the original call with new data/weights while preserving
#' the original scale formula, constraints, control, etc.
#'
#' **Note on weights**: If the original model was weighted, sandwich
#' usually passes weights = NULL. In that case we keep the original weights.
#' This means the bootstrap may not be properly re-weighted. For accurate
#' weighted bootstrap use our own [vcov()][mlmodels::vcov.mlmodel] instead.
#' 
#' **Note on sandwich**: This package's functions does not work reliably with 
#' `ml_negbin` objects.  We, therefore, built our own bootstrap and jackknife
#' implementations. We strongly recommend using [vcov()][mlmodels::vcov.mlmodel] instead.
#'
#' @export
update.ml_negbin <- function(object,
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