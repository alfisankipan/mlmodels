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
    cli::cli_warn(
      c("Variance matrix is not usable (contains NAs or non-finite values).",
        "i" = "This usually happens with bootstrap when constraints are present.",
        "i" = "Joint significance tests and correlation matrix will be skipped.")
    )
  } else {
    eig <- eigen(vcov_mat, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig < sqrt(.Machine$double.eps))) {
      usable_vcov <- FALSE
      cli::cli_warn(
        c("Variance matrix is singular or nearly singular.",
          "i" = "Joint significance tests and correlation matrix could not be computed.",
          "i" = "Consider using `vcov.type = 'robust'` instead.")
      )
    }
  }
  
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