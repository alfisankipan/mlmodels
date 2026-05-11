## ML_TRUNC_LM =================================================================
#' Fit Truncated (log)normal Model by Maximum Likelihood
#'
#' @param value Formula for the conditional mean (value) equation.
#' @param scale Formula for log(sigma) (optional). If `NULL`, a homoskedastic
#'   model is fitted.
#' @param left Left censoring value, variable, or vector. Defaults to `-Inf`.
#' @param right Right censoring value, variable, or vector. Defaults to `Inf`.
#' @param weights Optional weights variable. It can be either the name of the
#'   variable in `data`, or a vector with the weights.
#' @param data Data frame.
#' @param subset Optional subset expression. Only observations for which this
#'   expression evaluates to `TRUE` are used in the estimation. This can be
#'   a logical vector or an expression (e.g. `subset = age > 30`).
#' @param noint_value Logical. Should the value equation omit the intercept?
#'   Default is `FALSE`.
#' @param noint_scale Logical. Should the scale equation omit the intercept?
#'   Default is `FALSE`.
#' @param constraints Optional constraints on the parameters. Can be a character
#'   vector of string constraints, a named list of string constraints, or a raw
#'   maxLik constraints list. See **Details**.
#' @param start Numeric vector of starting values for the coefficients. Required
#'   if constraints are being supplied. If supplied without constraints they
#'   will be ignored. See **Details**.
#' @param method A string with the method used for optimization. See
#'   [maxLik][maxLik::maxLik()] for options, and see **Details**.
#' @param start Numeric vector of starting values for the coefficients. Required
#'   if constraints are being supplied. If supplied without constraints they
#'   will be ignored. See **Details**.
#' @param control A list of control parameters passed to [maxLik][maxLik::maxLik].
#'   If `NULL` (default), a sensible set of options is chosen automatically
#'   depending on whether constraints are used. See [maxControl][maxLik::maxControl].
#' @param ... Additional arguments passed to [maxLik][maxLik::maxLik].
#'
#' @details
#' **Important:** Do not use the usual R syntax to remove the intercept in the
#' formula (`- 1` or `+ 0`) for the value or scale equations. Use the dedicated
#' arguments `noint_value` and `noint_scale` instead.
#'
#' Coefficient names in the fitted object use the prefixes `value::` and
#' `scale::` to clearly identify to which equation each coefficient belongs to,
#' and to avoid confusion when the same variable(s) appear(s) in both the value
#' and scale equations.
#'
#' Either inequality or equality linear constraints are accepted, but not both.
#' A constraint cannot have a linear combination of more than two coefficients.
#'
#' **Important**: When `constraints` are supplied, `start` cannot be `NULL`.
#' You **must** provide initial values that yield a feasible log-likelihood.
#' If no constraints are used, any supplied `start` is ignored.
#'
#' When constraints are used, `ml_lm` automatically chooses the optimizer:
#' - Equality constraints => Nelder-Mead (`"NM"`)
#' - Inequality constraints => BFGS (`"BFGS"`)
#'
#' In these cases your supplied `method` argument (if any) is ignored.
#'
#' @return An object of class `ml_trunc_lm` that extends `mlmodel`.
#' 
#' @examples
#' 
#' # Homoskedastic linear model
#' data(mroz)
#' fit_lin <- ml_lm(faminc ~ age + I(age^2) + huswage + educ + unem, 
#'                  data = mroz)
#' summary(fit_lin, vcov.type = "robust")
#' 
#' # Heteroskedastic linear model
#' fit_het <- ml_lm(faminc ~ age + I(age^2) + huswage + educ + unem,
#'                  scale = ~ educ + exper,
#'                  data = mroz)
#' summary(fit_het, vcov.type = "robust")
#' 
#' # Lognormal (log-linear) model
#' fit_log <- ml_lm(log(faminc) ~ age + I(age^2) + huswage + educ + unem,
#'                  data = mroz)
#' summary(fit_log, vcov.type = "robust")
#' 
#' # Different predict types
#' head(predict(fit_log, type = "response")$fit)      # Expected value E[y]
#' head(predict(fit_log, type = "median")$fit)        # Median of y
#' head(predict(fit_log, type = "variance_y")$fit)    # Variance of y
#' head(predict(fit_log, type = "var")$fit)           # Variance of log(y)
#' 
#' # Fitted values and residuals
#' head(fitted(fit_lin))
#' head(residuals(fit_lin))
#' head(residuals(fit_lin, type = "pearson"))
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @export
ml_trunc_lm <- function(value,
                        scale = NULL,
                        weights = NULL,
                        left = -Inf,
                        right = Inf,
                        data,
                        subset = NULL,
                        noint_value = FALSE,
                        noint_scale = FALSE,
                        constraints = NULL,
                        start = NULL,
                        method = "NR",
                        control = NULL,
                        ...)
{
  cl <- match.call()
  # We call the helper to form the basic frame
  mf <- .model_frame_mlmodel(value = value,
                             scale = scale,
                             weights = rlang::enquo(weights),
                             data = data,
                             subset = rlang::enquo(subset),
                             cl = cl)
  
  # Pull out different values that we use.
  data <- mf$data
  wts <- mf$weights
  w_name <- mf$w_name
  sample <- mf$sample
  n_orig <- mf$n_orig
  
  # -- 5. Detect log transformation ------------------------------
  # We evaluate on the full data so invalid_idx has length = n_orig
  log_info <- .detect_log_transformations(list(value = value),
                                          data)
  
  # Count how many *additional* observations are invalid due to the log
  # (only among those that survived subset + NAs)
  n_invalid_log <- sum(sample & log_info$value$invalid_idx)
  
  if (log_info$value$is_log && n_invalid_log > 0) {
    cli::cli_alert_info(
      "Outcome is log-transformed. Dropped {n_invalid_log} observation(s) \\
       because they would produce invalid values (<= 0)."
    )
    
    # Final reduction of sample
    sample <- sample & !log_info$value$invalid_idx
  }
  
  # Add the actual dropped count to log_info for later use
  log_info$value$n_invalid_log <- n_invalid_log
  
  # Expansion
  mf$sample <- sample
  left  <- .expand_trunc_point(left,  mf = mf, name = "left")
  right <- .expand_trunc_point(right, mf = mf, name = "right")
  
  # -- Truncation point validation and log-transform -----------------------
  is_lognormal <- log_info$value$is_log
  
  # 1. Universal check: left must be strictly less than right everywhere
  if (any(left >= right - 1e-8, na.rm = TRUE)) {
    cli::cli_abort(c(
      "x" = "All left truncation points must be strictly less than right truncation points.",
      "i" = "At least one observation has left >= right."
    ))
  }
  
  # 2. Lognormal-specific checks (on original scale)
  if (is_lognormal) {
    # Finite left values must be strictly positive
    if (any(is.finite(left) & left <= 0, na.rm = TRUE)) {
      cli::cli_abort(c(
        "x" = "Left truncation point must be strictly positive when the outcome is log-transformed.",
        "i" = "Some `left` values are ≤ 0."
      ))
    }
    
    # Finite right values must be strictly positive
    if (any(is.finite(right) & right <= 0, na.rm = TRUE)) {
      cli::cli_abort(c(
        "x" = "Right truncation point must be strictly positive when the outcome is log-transformed.",
        "i" = "Some `right` values are ≤ 0."
      ))
    }
    
    left  <- ifelse(is.finite(left),  log(left),  left)
    right <- ifelse(is.finite(right), log(right), right)
  }
  
  if (all(is.infinite(left)) && all(is.infinite(right))) {
    cli::cli_abort(c(
      "x" = "`left` and `right` indicate no truncation.",
      "i" = "Use {.fn ml_lm} for a standard (non-truncated) Gaussian model."
    ))
  }
  
  # Check that observations fall between left and right
  lhs_var <- rlang::f_lhs(value)
  y_vec <- rlang::eval_tidy(lhs_var, data = data)
  y_vec <- y_vec[sample]
  valid_left <- y_vec >= left
  valid_right <- y_vec <= right
  valid <- valid_left & valid_right
  trunc_info <- list(
    left = left,
    right = right,
    n_before = sum(sample),
    tr_left = sum(!valid_left),
    tr_right = sum(!valid_right),
    tr_total = sum(!valid),
    left_is_constant  = length(unique(left))  == 1L,
    right_is_constant = length(unique(right)) == 1L
  )
  sample <- sample & valid
  
  # -- 6. Create clean dataset for modeling ----------------------
  data_clean <- data[sample, , drop = FALSE]
  wts_clean <- if (!is.null(wts)) wts[sample] else rep(1, sum(sample))
  left <- left[sample]
  right <- right[sample]
  
  # Add this safety check:
  if (length(wts_clean) != sum(sample) || any(is.na(wts_clean))) {
    cli::cli_abort("Final weights vector has wrong length or contains NAs.", call = NULL)
  }
  
  model_value <- hardhat::mold(value,
                               data_clean,
                               blueprint = hardhat::default_formula_blueprint(intercept = !noint_value))
  
  molds <- list(
    value = model_value
  )
  
  y <- model_value$outcomes[[1]]
  x <- as.matrix(model_value$predictors)
  
  if(!is.null(scale))
  {
    model_scale <- hardhat::mold(scale,
                                 data_clean,
                                 blueprint = hardhat::default_formula_blueprint(intercept = !noint_scale))
    
    molds$scale <- model_scale
    
    z <- as.matrix(model_scale$predictors)
  }
  else
  {
    z <- matrix(1, nrow = nrow(x), ncol = 1,
                dimnames = list(NULL, "lnsigma"))
    model_scale <- list(
      predictors = tibble::tibble(lnsigma = as.vector(z)),
      blueprint = NULL
    )
  }
  
  # -- 7. Map factor variables in relevant equations ---------
  factor_mapping <- .build_factor_mapping(molds)
  
  # -- 8. Managing control and constraints -------------------
  # Default control lists
  default_NR <- list(tol = -1,
                     reltol = 1e-12,
                     gradtol = 1e-12,
                     lambdatol = 1e-20)
  
  default_BFGS <- list(reltol = 1e-8)
  default_NM <- list(reltol = 1e-8,
                     iterlim = 1000)
  
  # Parse constraints (if any)
  if (!is.null(constraints)) {
    if (is.null(start))
      cli::cli_abort("Constrained optimization requires a vector of initial values (`start`).",
                     call = NULL)
    if (any(!is.numeric(start)))
      cli::cli_abort("Initial values must be numeric.", call = NULL)
    if (length(start) != (ncol(x) + ncol(z)))
      cli::cli_abort("The vector of initial values has the wrong dimension. It requires {.val {ncol(x) + ncol(z)}} values.")
    
    coef_names <- c(paste0("value::", colnames(x)), paste0("scale::", colnames(z)))
    parsed_constraints <- .parse_constraints(constraints, coef_names)
    
    if (!is.null(parsed_constraints$maxLik$eqA))
    {
      cli::cli_alert_info("Equality constraints detected => using Nelder-Mead optimizer.")
      method <- "NM"
      if (is.null(control))
        control <- default_NM
    }
    else
    {
      method = "BFGS"
      cli::cli_alert_info("Inequality constraints detected => using BFGS optimizer.")
      if (is.null(control))
        control <- default_BFGS
    }
  } else {
    parsed_constraints <- list(names = NULL, strings = NULL, maxLik = NULL)
    start <- NULL
    
    # Unconstrained optimization
    if (is.null(control)) {
      if (method %in% c("NR", "BHHH")) {
        control <- default_NR
      } else if(method == "NM") {
        control <- default_NM
      } else {
        control <- default_BFGS
      }
    }
  }
  
  # -- 9. Fitting the model with maxLik ----------------------------------------
  
  # -- 9a. Scaling the weights to ease optimization ----------------------------
  sc_factor <- sum(wts_clean)
  w_scaled <- wts_clean / sc_factor
  
  #-- 9b. Calling the fit function with scaled weights -------------------------
  ml <- .ml_trunc_lm.fit(y = y,
                         x = x,
                         z = z,
                         left = left,
                         right = right,
                         w = w_scaled,
                         lognormal = is_lognormal,
                         constraints = parsed_constraints$maxLik,
                         start = start,
                         method = method,
                         control = control,
                         ...)
  
  # -- 9c. Scaling the log-likelihood, scores and hessian back -----------------
  ml$hessian <- ml$hessian * sc_factor
  ml$gradient <- ml$gradient * sc_factor
  ml$gradientObs <- ml$gradientObs * sc_factor
  ml$maximum <- ml$maximum * sc_factor
  
  # -- 10. Forming the dataset name ------------------------------
  # Safely get a readable name for the dataset (for printing/storage)
  d_name <- tryCatch(
    deparse(substitute(data)),
    error = function(e) "<unknown data>"
  )
  
  if (length(d_name) > 1 || d_name == "NULL" || grepl("^\\s*\\(", d_name)) {
    d_name <- "<unknown data>"
  }
  
  # -- 11. Internal safety check(s) (for development/testing) --------
  # They should be the same value, since we never indexed sample (that i can remember),
  # so if we get an alert, we must check the code.
  if (length(sample) != n_orig) {
    cli::cli_alert_danger(
      "Internal error: length of 'sample' ({length(sample)}) does not match n_orig ({n_orig})."
    )
  }
  
  # -- 12. Forming the the model list --------------------------------
  
  # -- 12.a. The functions list --------------------------------------
  
  functions <- list(
    # predict        = predict.ml_lm,
    # gradientObs    = .ml_lm_gradientObs,
    # hessianObs     = .ml_lm_hessianObs,
    # loglikeObs     = .ml_lm_loglikeObs,
    # loglik         = .ml_lm_ll,
    # fit            = .ml_lm.fit
  )
  
  # -- 12.b. The common structure --------------------------------------
  if(!is.null(scale))
  {
    description <- if(log_info$value$is_log)
      "Heteroskedastic Lognormal Model"
    else
      "Heteroskedastic Linear Model"
  }
  else
  {
    description <- if(log_info$value$is_log)
      "Homoskedastic Lognormal Model"
    else
      "Homoskedastic Linear Model"
  }
  
  model_list <- list(
    description   = description,
    value         = model_value,
    scale         = model_scale,
    factor_mapping = factor_mapping,
    formula       = model_value$blueprint$formula,
    scale_formula = if(!is.null(scale)) model_scale$blueprint$formula else NULL,
    weights       = wts_clean,
    w_name        = w_name,
    sample        = sample,
    subset_sample = mf$keep,
    usable_sample = mf$usable_obs,
    data          = data,
    data_name     = d_name,
    functions     = functions,
    response_name = names(model_value$outcomes)[1],
    n_used        = sum(sample),
    n_orig        = n_orig,
    log_info      = log_info,
    control       = control,
    constraints   = parsed_constraints,
    start         = start,
    method        = method,
    trunc_info    = trunc_info
  )
  
  if (!(ml$code %in% c(0, 1, 2, 8))) {
    cli::cli_alert_warning(
      "Estimation did not converge (code {.strong {ml$code}}).\nMessage: {ml$message}"
    )
    cli::cli_alert_info(
      "Returning model without fitted/residual values. Use coef() to inspect parameters."
    )
    model_list$fitted.values <- NULL
    model_list$residuals     <- NULL
    ml$model <- model_list
    return(ml)
  }
  
  # -- 12.c. The fitted values and residuals ------------------------
  # Converged: compute fitted/residuals
  coefs <- coef(ml)
  
  beta <- coefs[1:ncol(x)]
  yhat <- as.vector(x %*% beta)
  delta <- coefs[(ncol(x) + 1):length(coefs)]
  sigma <- as.vector(exp(z %*% delta))
  
  model_list$fitted.values <- yhat
  model_list$residuals     <- y - yhat
  model_list$sigma         <- sigma
  
  # -- 13. Add the model to the maxLik object ----------------------
  ml$model <- model_list
  ml$call <- mf$cl
  
  # -- 14. Call the function to create the class and return  ----------
  new_ml_trunc_lm(ml)
}

# Hidden function to create the class and return the object.
new_ml_trunc_lm <- function(object, ...) {
  # object is the result from maxLik::maxLik()
  structure(
    object,
    class = unique(c("ml_trunc_lm", "mlmodel", class(object)))
  )
}

# ML_TRUNC_LM FIT --------------------------------------------------------------
#' @keywords internal
.ml_trunc_lm.fit <- function(y, x, z, w = NULL,
                             left = -Inf,
                             right = Inf,
                             method = "NR",
                             start = NULL,
                             constraints = NULL,
                             lognormal = FALSE,
                             control = list(tol = -1,
                                            reltol = 1e-12,
                                            gradtol = 1e-12,
                                            lambdatol = 1e-20),
                             ...)
{
  # At this point start has been checked and if supplied is numeric and
  # has the right dimension.
  if(!is.null(start))
  {
    ll <- .ml_trunc_lm_ll(start,y = y, x = x, z = z,
                          left = left, right = right, w = w, lognormal = lognormal)
    
    if(any(!is.finite(ll)))
      cli::cli_abort("Infeasible log-likelihood value at supplied `start` vector.",
                     call = NULL)
    
    names(start) <- c(paste0("value::", colnames(x)),
                      paste0("scale::", colnames(z)))
  }
  else
  {
    # Beta starting values using low-level lm.fit()
    fit_beta <- .lm.fit(x, y)
    b0 <- fit_beta$coefficients
    
    names(b0) <- paste0("value::", colnames(x))
    
    # Residuals
    resid <- y - x %*% b0
    # Scale starting values
    if (ncol(z) == 1 && colnames(z)[1] == "lnsigma") {
      # Homoskedastic case => exact MLE
      n <- length(resid)
      rss <- sum(resid^2)
      g0 <- 0.5 * (log(rss) - log(n))
    } else {
      # Heteroskedastic case => auxiliary regression on log(resid^2)
      ln_res2 <- log(resid^2)
      fit_aux <- .lm.fit(z, ln_res2)
      g0 <- fit_aux$coefficients / 2
    }
    # Apply scale:: prefix to all scale coefficients
    names(g0) <- paste0("scale::", colnames(z))
    
    start <- c(b0, g0)
    
    start <- .initial_values.mlmodel(.ml_trunc_lm_ll, start,
                                     y = y, x = x, z = z, w = w, left = left,
                                     right = right,
                                     lognormal = lognormal)
    if(isFALSE(attr(start, "feasible")))
      cli::cli_abort("Couldn't find feasible initial values.",
                     call = NULL)
  }
  maxLik::maxLik(.ml_trunc_lm_ll,
                 start = start,
                 y = y,
                 x = x,
                 z = z,
                 left = left,
                 right = right,
                 w = w,
                 lognormal = lognormal,
                 constraints = constraints,
                 method = method,
                 control = control,
                 ...)
}