# SUMMARY ----------------------------------------------------------------------
#' Summary for ml_logit objects
#'
#' @param object A fitted model object of class `"ml_logit"`.
#' @param correlation Logical. Should the correlation matrix of the estimated
#'   parameters be computed and stored in the returned summary object as
#'   `$correlation`? Default is `FALSE`. This is provided for compatibility
#'   with `maxLik::summary`. The matrix is **not** printed automatically.
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
#' @details
#' Coefficient names in the fitted object use the prefixes `value::` and
#' `scale::` to identify to which equation they belong to, and to avoid
#' confusion when the same variable(s) appear(s) in both the value and scale
#' equations.
#'
#' For binary logit models, two pseudo-R-squared measures are computed and
#' stored in the `$r.squared` component:
#' - `cor`: Squared correlation between observed and fitted values.
#' - `mcfadden`: McFadden's pseudo-R².
#'
#' @return An object of class `c("summary.ml_logit", "summary.mlmodel", "summary")`.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @export
summary.ml_logit <- function(object,
                             correlation = FALSE,
                             vcov = NULL,           # User-supplied variance matrix
                             vcov.type = "oim",
                             cl_var = NULL,
                             repetitions = 999,
                             seed = NULL,
                             progress = FALSE,
                             ...)
{
  if (!inherits(object, "ml_logit"))
    cli::cli_abort("`object` must be a model of class 'ml_logit'.")
  
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
  
  # Extract observations, and number of parameters.
  n <- object$model$n_used
  k_total <- length(coef(object))
  k_scale <- if (!is.null(object$model$scale)) ncol(object$model$scale$predictors) else 0L
  is_heteroskedastic <- !is.null(object$model$scale)
  k_mean <- k_total - k_scale
  
  # Start building the summary object
  s <- list()
  
  s$response_name <- object$model$response_name
  s$call          <- object$call
  s$formula       <- object$model$formula
  s$scale_formula <- object$model$scale_formula
  s$nobs          <- n
  s$converged     <- converged
  s$is_heteroskedastic <- is_heteroskedastic
  
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
    
    ll <- as.numeric(logLik(object))
    s$logLik         <- ll
    s$AIC            <- -2 * ll + 2 * k_total
    s$BIC            <- -2 * ll + log(n) * k_total
    
    p_bar <- mean(y)
    ll0   <- sum(y) * log(p_bar) + sum(1 - y) * log(1 - p_bar)
    
    s$r.squared <- list(
      cor = cor(y, yhat)^2,
      mcfadden = 1 - ll / ll0
    )
    
    if(usable_vcov)
    {
      # Joint significance tests (reuse vcov_mat)
      idx_mean <- if (object$model$value$blueprint$intercept) 2:k_mean else 1:k_mean
      
      if (is_heteroskedastic) {
        idx_scale <- (k_mean + 1):k_total
        
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
  
  
  s$model_type <- if (is_heteroskedastic) {
    "Heteroskedastic Binary Logit"
  } else {
    "Homoskedastic Binary Logit"
  }
  
  class(s) <- c("summary.ml_logit", "summary.mlmodel", "summary")
  s
}

# PRINT SUMMARY ----------------------------------------------------------------
#' @export
print.summary.ml_logit <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!inherits(x, "summary.ml_logit"))
    cli::cli_abort("`x` needs to be a `summary.ml_logit` object.")
  
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
    cat("Joint significance tests:\n")
    
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
                         "boot" = if(is.null(x$vcov.cluster)) "Bootstrapped" else "Cluster Bootstrapped",
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
  options(scipen = 2)
  
  # Capture the whole output of printCoefmat into a vector of strings.
  captured <- capture.output(printCoefmat(x$coefficients,
                                          digits = digits,
                                          signif.legend = TRUE))
  
  # Get the number of coefficients in each equation.
  k1 <- sum(grepl("^value::", rownames(x$coefficients)))
  if(x$is_heteroskedastic)
    k2 <- sum(grepl("^scale::", rownames(x$coefficients)))
  else
    k2 <- 0
  
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
  if(x$is_heteroskedastic)
  {
    cat("Scale (log(sigma)):")
    cat("  ", captured[(k1+2):(k1+k2+1)],
        sep = "\n  ")
  }
  # It seems as if there is an empty line between the coefficients and the legend
  # so we need to add one more line and start at k1+k2+3, to avoid that empty one.
  cat("---------------------------------------",
      captured[(k1+k2+3):length(captured)],
      sep = "\n")
  
  options(scipen = old_pen)
  
  if (x$converged) {
    cat("---\n")
    cat("Number of observations:", x$nobs, "\n")
    cat("Pseudo R-squared - Cor.Sq.: ",
        format(x$r.squared$cor, digits = digits),
        " McFadden: ", format(x$r.squared$mcfadden, digits = digits),
        "\n",
        sep = "")
    cat("AIC:", format(x$AIC, nsmall = 2, digits = digits + 1),
        " BIC:", format(x$BIC, nsmall = 2, digits = digits + 1), "\n")
  } else {
    cat("\nGoodness-of-fit statistics not available (model did not converge).\n")
  }
  
  invisible(x)
}