#' Fit linear model by Maximum Likelihood
#'
#' @param value Formula for the conditional mean (value) equation.
#' @param scale Formula for log(sigma) (optional). If NULL, a homoskedastic
#'   model is fitted.
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
#' @param control A list of control parameters passed to [maxLik::maxLik()].
#'   The default values are chosen to work well with difficult likelihoods.
#'   See [maxLik::maxLik()] for details.
#' @param ... Additional arguments passed to [maxLik::maxLik()].
#'
#' @details
#' **Important:** Do not use the usual R syntax to remove the intercept in the
#' formulas (`- 1` or `+ 0`). Use the dedicated arguments `noint_value` and
#' `noint_scale` instead.
#'
#' @return An object of class `ml_lm` that extends `mlmodel`.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @examples
#' if (requireNamespace("wooldridge", quietly = TRUE)) {
#'
#'   # Homoskedastic model
#'   linm <- ml_lm(cigs ~ log(income) + age + I(age^2) + educ + restaurn,
#'                 data = wooldridge::smoke)
#'   summary(linm)
#'
#'   # Heteroskedastic model
#'   linm_het <- ml_lm(value = cigs ~ log(income) + age + I(age^2) + educ + restaurn,
#'                     scale = ~ age + I(age^2) + restaurn,
#'                     data = wooldridge::smoke)
#'   summary(linm_het, vcov.type = "robust")
#' }
#'
#' @export
ml_lm <- function(value,
                  scale = NULL,
                  weights = NULL,
                  data,
                  subset = NULL,
                  noint_value = FALSE,
                  noint_scale = FALSE,
                  control = list(tol = -1,
                                 reltol = 1e-12,
                                 gradtol = 1e-12,
                                 lambdatol = 1e-20,
                                 qac = "marquardt"),
                  ...)
{
  # -- Basic input validation ------------------------------------------
  if (!rlang::is_formula(value, lhs = TRUE)) {
    cli::cli_abort("`value` must be a two-sided formula with an outcome variable on the left-hand side.",
                   call = NULL)
  }

  if (!is.null(scale) && !rlang::is_formula(scale, lhs = FALSE)) {
    cli::cli_abort("`scale` must be a one-sided formula (no outcome on the left-hand side).",
                   call = NULL)
  }

  cl <- match.call()

  # -- Save original data dimensions and create keep vector ------------
  n_orig <- nrow(data)
  keep <- rep(TRUE, n_orig)          # Start with all observations kept

  # -- 1. Handle subset argument --------------------------------------
  if (!is.null(subset)) {
    if (is.logical(subset) && length(subset) == n_orig) {
      subset_idx <- subset
    } else {
      subset_expr <- rlang::enquo(subset)
      subset_idx <- rlang::eval_tidy(subset_expr, data = data)
      if (!is.logical(subset_idx)) {
        subset_idx <- as.logical(subset_idx)
      }
    }

    if (length(subset_idx) != n_orig) {
      cli::cli_abort("`subset` must have the same length as `data`.", call = NULL)
    }

    keep <- keep & subset_idx
  }

  # Always store subset in the call (even if NULL)
  cl$subset <- subset

  # -- 2. Weights handling ------------------------------------------
  w_expr <- rlang::enquo(weights)
  if (!rlang::quo_is_null(w_expr)) {
    # Try to see if it looks like a column name
    w_name <- tryCatch(rlang::as_name(w_expr), error = function(e) NULL)
    wts    <- rlang::eval_tidy(w_expr, data)
  } else {
    w_name <- NULL
    wts    <- NULL
  }

  # Safety: if user passed a vector, w_name should not be treated as a column
  if (!is.null(w_name) && !w_name %in% names(data)) {
    w_name <- NULL   # it was a direct vector, not a column
  }

  # -- 3. Identify usable observations on full data ----------------
  if (is.null(scale)) {
    v_vars <- all.vars(value)
  } else {
    v_vars <- unique(c(all.vars(value), all.vars(scale)))
  }

  cols_to_check <- v_vars
  if (!is.null(w_name)) {
    cols_to_check <- unique(c(cols_to_check, w_name))
  }

  usable_obs <- complete.cases(data[, cols_to_check, drop = FALSE])

  # Extra check if user passed a weights vector directly
  if (!is.null(wts) && is.null(w_name)) {
    usable_obs <- usable_obs & complete.cases(wts)
  }

  # Count observations dropped due to missing values *within the subset*
  nas_dropped <- sum(keep & !usable_obs)
  if (nas_dropped > 0) {
    cli::cli_alert_info("Dropped {nas_dropped} observations due to missing values.")
  }

  # -- 4. Modify sample = subset ∩ complete cases ----------------
  sample <- keep & usable_obs

  # -- 5. Detect log transformation ------------------------------
  # We evaluate on the full data so invalid_idx has length = n_orig
  log_info <- detect_log_transformations(list(value = value),
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

  # -- 6. Create clean dataset for modeling ----------------------
  data_clean <- data[sample, , drop = FALSE]
  wts_clean <- if (!is.null(wts)) wts[sample] else rep(1, sum(sample))

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

  # -- 7. Map factor variables in relevant equations ------------------
  factor_mapping <- build_factor_mapping(molds)

  # -- 8. Fitting the model with maxLik ----------------------
  ml <- .ml_lm.fit(y = y,
                   x = x,
                   z = z,
                   w = wts_clean,
                   control = control,
                   ...)

  # -- 9. Forming the dataset name ------------------------------
  # Safely get a readable name for the dataset (for printing/storage)
  d_name <- tryCatch(
    deparse(substitute(data)),
    error = function(e) "<unknown data>"
  )

  if (length(d_name) > 1 || d_name == "NULL" || grepl("^\\s*\\(", d_name)) {
    d_name <- "<unknown data>"
  }

  # -- 10. Internal safety check(s) (for development/testing) --------
  # They should be the same value, since we never indexed sample (that i can remember),
  # so if we get an alert, we must check the code.
  if (length(sample) != n_orig) {
    cli::cli_alert_danger(
      "Internal error: length of 'sample' ({length(sample)}) does not match n_orig ({n_orig})."
    )
  }

  # -- 11. Forming the the model list --------------------------------

  # -- 11.a. The functions list --------------------------------------

  functions <- list(
    predict        = predict.ml_lm,
    hessianObs     = ml_lm_hessianObs,
    update         = update.ml_lm,
    loglik         = ml_lm_ll,
    fit            = .ml_lm.fit
  )

  # Common model list structure
  model_list <- list(
    value         = model_value,
    scale         = model_scale,
    factor_mapping = factor_mapping,
    formula       = model_value$blueprint$formula,
    scale_formula = if(!is.null(scale)) model_scale$blueprint$formula else NULL,
    weights       = wts_clean,
    w_name        = w_name,
    sample        = sample,
    subset_sample = keep,
    usable_sample = usable_obs,
    data          = data,
    data_name     = d_name,
    functions     = functions,
    response_name = names(model_value$outcomes)[1],
    n_used        = sum(sample),
    n_orig        = n_orig,
    log_info      = log_info
  )

  if (!(ml$code %in% c(1, 2, 8))) {
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

  # Converged: compute fitted/residuals
  beta <- coef(ml)[1:ncol(x)]
  yhat <- as.vector(x %*% beta)

  model_list$fitted.values <- yhat
  model_list$residuals     <- y - yhat

  # -- 11. Add the model to the maxLik object ----------------------
  ml$model <- model_list
  ml$call <- cl

  # -- 12. Call the function to create tge class and return  ----------
  new_ml_lm(ml)
}

# Hidden function to create the class and return the object.
new_ml_lm <- function(object, ...) {
  # object is the result from maxLik::maxLik()
  structure(
    object,
    class = unique(c("ml_lm", "mlmodel", class(object)))
  )
}

# ML_LM FIT --------------------------------------------------------------
.ml_lm.fit <- function(y, x, z, w,
                       control = list(tol = -1,
                                      reltol = 1e-12,
                                      gradtol = 1e-12,
                                      lambdatol = 1e-20,
                                      qac = "marquardt"),
                       ...)
{
  # Beta starting values using low-level lm.fit()
  fit_beta <- .lm.fit(x, y)
  b0 <- fit_beta$coefficients
  names(b0) <- colnames(x)

  # Residuals
  resid <- y - x %*% b0

  # Scale starting values
  if (ncol(z) == 1 && colnames(z)[1] == "lnsigma") {
    # Homoskedastic case → exact MLE
    n  <- length(resid)
    rss <- sum(resid^2)
    g0  <- 0.5 * (log(rss) - log(n))
    names(g0) <- "lnsigma"
  } else {
    # Heteroskedastic case → auxiliary regression on log(resid²)
    ln_res2 <- log(resid^2)
    fit_aux <- .lm.fit(z, ln_res2)
    g0 <- fit_aux$coefficients / 2
    names(g0) <- paste0("scale::", colnames(z))
  }

  # Final estimation
  maxLik::maxLik(ml_lm_ll,
                 start = c(b0, g0),
                 y = y,
                 x = x,
                 z = z,
                 w = w,
                 control = control,
                 ...)
}
