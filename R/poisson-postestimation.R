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
#' | `P(k)`        | P(Y = k)                                 | Exact probability, k integer >= 0 |
#' | `P(,k)`       | P(Y <= k)                                 | Cumulative (lower tail) |
#' | `P(k,)`       | P(Y >= k)                                 | Survival (upper tail) |
#' | `P(a,b)`      | P(a <= Y <= b)                             | Interval probability, a <= b, a >= 0 |
#'
#' When `se.fit = TRUE`, standard errors are computed using the delta method
#' for all supported types.
#'
#' @rdname predict.mlmodel
#' @method predict ml_poisson
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
  
  # -- Prepare predictors (using hardhat) --------------------------------------
  predictors <- .prepare_prediction_data(object, newdata = newdata)
  X <- predictors$X
  
  # -- Extract coefficients and compute linear predictors ----------------------
  beta <- coef(object)
  
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
  
  # -- Align in-sample predictions to original data length ---------------------
  if (is.null(newdata)) {
    out <- .predict_align_estimates(object, out)
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
    
    # Derivative of the probability w.r.t. mu
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
  
  # -- Check for unusable variance matrix --------------------------------------
  if (any(!is.finite(full_vcov)) || any(is.na(full_vcov))) {
    cli::cli_warn(
      c("Variance matrix is unusable (contains NAs or non-finite values).",
        "i" = "This usually happens with bootstrap when constraints are present.",
        "i" = "Standard errors will be returned as NA.")
    )
    se_fit <- rep(NA_real_, length(out))
  } else {
    # -- Delta-method standard errors ------------------------------------------
    se_fit <- sqrt(rowSums(g * (g %*% full_vcov)))
  }
  
  # Align in-sample SEs
  if (is.null(newdata)) {
    se_fit <- .predict_align_estimates(object, se_fit)
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
    # Use helper function to print the weight/loglikelihood part of the header
    .print_weight_loglik(x, digits = digits)
    
    # Wald Tests
    .print_wald_tests(x, digits = digits)
  }
  
  cat("\nVariance type:", x$var_description)
  cat("\n---------------------------------------\n")
  
  old_pen <- getOption("scipen")
  options(scipen = .mlmodels_get_default("scipen"))
  
  # Determining the number of leading zeroes in the estimates and standard errors.
  format_coef <- .format_coef_matrix(x$coefficients, digits = digits)
  
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
    
    cat("Observations:\n")
    labels <- c("Res. Deg. of Freedom:", "Sample:")
    if(x$weight_info$is_weighted)
    {
      labels <- c(labels,
                  "Effective (Sum Wts):")
      width <- max(nchar(labels)) + 1
      cat(sprintf("  %-*s %d", width, "Sample:", x$nobs),
          sprintf("  %-*s %d", width, "Effective (Sum Wts):", x$weight_info$sum_weights),
          sprintf("  %-*s %d", width, "Res. Deg. of Freedom:", x$df.residual),
          sep = "\n")
    }
    else
    {
      width <- max(nchar(labels)) + 1
      cat(sprintf("  %-*s %d", width, "Sample:", x$nobs),
          sprintf("  %-*s %d", width, "Res. Deg. of Freedom:", x$df.residual),
          sep = "\n")
    }
    
    cat("\nGoodness of Fit:",
        "  Pseudo R-Squared:", sep = "\n")
    
    labels <- c("Cor.Sq.:", "McFadden:")
    width <- max(nchar(labels)) + 1
    cat(sprintf("    %-*s %.4f", width, "Cor.Sq.:", x$r.squared$cor),
        sprintf("    %-*s %.4f", width, "McFadden:", x$r.squared$mcfadden),
        sep = "\n")
    
    # Call helper to print the AIC and BIC with or without scaling
    .print_information_criteria(x, digits)
    
    cat("\nCount Diagnostics:\n")
    labels <- c("Observed:",
                "Predicted:")
    
    width <- max(nchar(labels)) + 1
    cat(sprintf("  Dispersion Ratio (Pearson):  %.2f", x$ov),
        "  Zeros:",
        sprintf("    %-*s %d", width, "Observed:", x$zero$count),
        sprintf("    %-*s %.2f", width, "Predicted:", x$zero$pred),
        sep = "\n")
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
  s$call           <- object$call                    # <- Now using root-level call
  s$formula        <- object$model$formula
  s$weight_info    <- .generate_weight_info(object)
  s$nobs           <- n
  if(s$weight_info$is_weighted)
    n_w <- s$weight_info$sum_weights
  else
    n_w <- n
  s$df.residual    <- n_w - k_total
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
    
    s$AIC <- AIC(object, scaled = FALSE)
    s$BIC <- BIC(object, scaled = FALSE)
    
    # Weights vector: real weights or 1s for unweighted case
    w <- object$model$weights %||% rep(1, n)   # n = actual number of obs
    
    # overdispersion
    s$ov <- sum(w * (y - yhat)^2 / yhat, na.rm = TRUE) / (n_w - k_total)
    
    s$zero <- list(
      count = sum(w[y == 0], na.rm = TRUE),
      pred = sum(w * exp(-yhat), na.rm = TRUE)
    )
    
    # Null logLik
    ll <- s$logLik
    y_bar <- mean(y)
    ll0 <- sum(w * (y * log(y_bar) - y_bar - lfactorial(y)), na.rm = TRUE)
    
    s$r.squared <- list(
      cor = cor(y, yhat)^2,
      mcfadden = 1 - ll / ll0
    )
    
    # Weight Information (from helper)
    s$weight_info <- .generate_weight_info(object)
    
    if(usable_vcov)
    {
      # Joint significance test: only mean because no scale.
      idx_mean <- if (object$model$value$blueprint$intercept) 2:k_mean else 1:k_mean
      
      s$significance <- list(
        overall  = waldtest(object, indices = idx_mean, vcov = vcov_mat),
        mean = NULL,
        scale = NULL
      )
    }
    else
    {
      s$significance <- list(
        overall  = NULL,
        mean = NULL,
        scale = NULL
      )
    }
    
  } else {
    s$r.squared <- s$AIC <- s$BIC <- s$ov <- s$zero <- s$significance <- s$weight_info <- NULL
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
                              subset = NULL,
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
  
  # Process Subset (sandwich uses a vector fo indices in it to resample data)
  if (!is.null(subset) && is.numeric(subset) && !all(subset %in% c(0L, 1L))) {
    
    # This is a vector of indices of the rows from the estimated data that
    # have to be in the new estimation data (resampling).
    
    # We pull the original dataset and index it to reduce it to the observations
    # used in the original estimation.
    est_data <- object$model$data[object$model$sample, , drop = FALSE]
    
    if (min(subset) < 1 || max(subset) > nrow(est_data)) {
      cli::cli_abort("Invalid indices in `subset` vector.")
    }
    
    # If bootstrapping the vector will have a length equal to the number of
    # observations in est_data, but if jackknifing it will be the number of
    # observations minus one (the one dropped out).
    est_n <- nrow(est_data)
    
    if (length(subset) == est_n) {
      # Bootstrapping. We just index est_data with the vector to form the new
      # dataframe.
      new_data <- est_data[subset, , drop = FALSE]
      
    } else if (length(subset) == est_n - 1) {
      # Jackknifing, we form a logical vector with length equal to est_n that will
      # end up having FALSE in the observation that was dropped.
      full_idx <- rep(FALSE, est_n)
      full_idx[subset] <- TRUE        # subset = indices to KEEP
      # And now we index est_data to return all observations except the FALSE one.
      new_data <- est_data[full_idx, , drop = FALSE]
    } else {
      cli::cli_abort("`subset` vector length ({length(subset)}) is invalid for bootstrap or jackknife.")
    }
    
    call$data   <- new_data      # Pass the resampled data
    # Since we are passing a dataframe with only data used in the original estimation
    # we must set the subset argument in the call to NULL, in case the original
    # estimation was done with a subset condition, because the dataset has already
    # been subset.
    call$subset <- NULL
  }
  else if (!is.null(subset)) {
    # If subset is not a vector of integers different from 0 and 1, then it must
    # be a true subset condition: we modify the argument in the call
    call$subset <- subset
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