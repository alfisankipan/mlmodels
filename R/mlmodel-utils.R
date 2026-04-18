# mlmodels: Exported utility functions for all mlmodel objects

## AIC =========================================================================
# --- General ------------------------------------------------------------------
#' Extract AIC from mlmodel objects
#'
#' @param object An object of class `"mlmodel"` or `"summary.mlmodel"`.
#' @param k Numeric. The penalty per parameter to be used. Default is `k = 2`,
#'   which gives the standard AIC. See [stats::AIC()] for details.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' For fitted `mlmodel` objects, AIC is computed as `-2 * logLik + k * npar`,
#' where `npar` is the total number of coefficients.
#'
#' For `summary.mlmodel` objects, the pre-computed AIC (calculated with `k = 2`)
#' is returned. The `k` argument is accepted for compatibility but is ignored.
#'
#' @return A numeric value with the AIC.
#' 
#' @importFrom stats AIC
#' @export
AIC <- function(object, ..., k = 2) UseMethod("AIC")

# --- mlmodel ------------------------------------------------------------------
#' @rdname AIC
#' @export
AIC.mlmodel <- function(object, ..., k = 2)
{
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must inherit from class 'mlmodel'.", call = NULL)
  
  if (!(object$code %in% c(0L, 1L, 2L, 8L))) {
    cli::cli_abort("AIC is not available (model did not converge).", call = NULL)
  }
  
  ll <- logLik(object)
  npar <- attr(ll, "df") %||% length(coef(object))
  
  -2 * as.numeric(ll) + k * npar
}

# --- summary.mlmodel ----------------------------------------------------------
#' @rdname AIC
#' @export
AIC.summary.mlmodel <- function(object, ..., k = 2)
{
  if (!inherits(object, "summary.mlmodel"))
    cli::cli_abort("`object` must inherit from class 'summary.mlmodel'.", call = NULL)
  
  if (is.null(object$AIC) || !isTRUE(object$converged))
    cli::cli_abort("AIC is not available (model did not converge).", call = NULL)
  
  # Note: we ignore the `k` argument for summary objects because we already computed AIC with k=2
  object$AIC
}

## BIC =========================================================================
# --- General ------------------------------------------------------------------
#' Extract BIC
#'
#' @param object An object of class `"mlmodel"` or `"summary.mlmodel"`.
#' @param ... Further arguments passed to methods.
#'
#' @importFrom stats BIC
#' @export
BIC <- function(object, ...) UseMethod("BIC")

# --- mlmodel ------------------------------------------------------------------
#' @rdname BIC
#' @export
BIC.mlmodel <- function(object, ...)
{
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must inherit from class 'mlmodel'.", call = NULL)
  
  if (!(object$code %in% c(0L, 1L, 2L, 8L))) {
    cli::cli_abort("AIC is not available (model did not converge).", call = NULL)
  }
  
  ll <- logLik(object)
  npar <- attr(ll, "df") %||% length(coef(object))
  nobs <- attr(ll, "nobs") %||% object$model$n_used
  
  -2 * as.numeric(ll) + log(nobs) * npar
}

# --- summary.mlmodel ----------------------------------------------------------
#' @rdname BIC
#' @export
BIC.summary.mlmodel <- function(object, ...)
{
  if (!inherits(object, "summary.mlmodel"))
    cli::cli_abort("`object` must inherit from class 'summary.mlmodel'.", call = NULL)
  
  if (is.null(object$BIC))
    cli::cli_abort("BIC is not available (model did not converge).", call = NULL)
  
  object$BIC
}

## COEFFICIENTS ================================================================
#' Gets the coefficients from an mlmodel estimation.
#'
#' @param object An mlmodel from an estimator, that you want the coefficients for.
#'
#' @param ... Currently not used.
#'
#' @returns A vector with the coefficients from the estimation.
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @method coef mlmodel
#' @export
coef.mlmodel <- function(object, ...) {
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must be a model of class 'mlmodel'.",
                   call = NULL)
  object$estimate
}

## CONFINT =====================================================================
#' Gets the confidence interval for an estimated model.
#' 
#' @param object Any model estimated with one of our estimators.
#' @param parm A vector with the names or indices of the parameters you want the
#'    interval for. If missing all parameters are considered.
#' @param level The confidence level for the interval.
#' @param vcov Optional user-supplied variance-covariance matrix.
#' @param vcov.type Type of variance-covariance matrix. See [vcov][vcov.mlmodel].
#' @param cl_var Clustering variable (name or vector).
#' @param repetitions Number of bootstrap replications when `vcov.type = "boot"`.
#' @param seed Random seed for the boostrapping, for reproducibility.
#' @param progress Logical. Show bootstrap/jackknife progress bar? Default is
#'   `FALSE` in higher-level functions.
#' @param ... Additional arguments for methods.
#' 
#' @method confint mlmodel
#' @export
confint.mlmodel <- function(object,
                            parm,
                            level = 0.95,
                            vcov = NULL,
                            vcov.type = "oim",
                            cl_var = NULL,
                            repetitions = 999,
                            seed = NULL,
                            progress = FALSE,
                            ...)
{
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must be a model of class 'mlmodel'.", call = NULL)
  
  coefs <- coef(object)
  
  if(missing(parm))
    cfs_idx <- rep(TRUE, length(coefs))
  else
  {
    if (is.numeric(parm))
    {
      # indices
      if (any(is.na(parm)) || any(parm > length(coefs)) || any(parm < 1) ||
          any(parm != floor(parm)))
        cli::cli_abort("Invalid indices in `parm`.", call = NULL)
      cfs_idx <- seq_along(coefs) %in% parm
    }
    else
    {
      if(!is.character(parm))
        cli::cli_abort("Invalid format for `parm`.", call = NULL)
      cfs_idx <- names(coefs) %in% parm
    }
  }
  # This could happen if the names in parm are not in the coefficient names, but
  # it is safer to do it out of the if statement in case I missed something in
  # the is.numeric part.
  if (!any(cfs_idx))
    cli::cli_abort("Invalid parameters in `parm`.", call = NULL)
  
  # Preliminary values.
  a <- (1 - level) / 2
  crit <- qnorm(c(a, 1 - a))
  parm_names <- names(coefs)[cfs_idx]
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
  
  # Process the vcov
  vcov <- .process_vcov(object,
                        vcov = vcov,
                        vcov.type = vcov.type,
                        cl_var = cl_var,
                        repetitions = repetitions,
                        seed = seed,
                        progress = progress)
  
  
  # ── Check for unusable variance matrix ─────────────────────────────
  if (any(!is.finite(vcov)) || any(is.na(vcov))) {
    cli::cli_warn(
      c("Variance matrix is unusable (contains NAs or non-finite values).",
        "i" = "This usually happens with bootstrap when constraints are present.",
        "i" = "Standard errors will be returned as NA.")
    )
    ci_matrix <- matrix(NA, nrow = length(parm_names), ncol = 2)
  }
  else
  {
    ses <- sqrt(diag(vcov))
    ci_matrix <- coefs + ses %o% crit
    ci_matrix <- ci_matrix[cfs_idx, , drop = FALSE]
  }
  rownames(ci_matrix) <- parm_names
  colnames(ci_matrix) <- pct
  
  ci_matrix
}

## FITTED ======================================================================
# --- mlmodel ------------------------------------------------------------------
#' Returns the fitted values from an object of class `mlmodel`.
#' 
#' @param object Object of class `mlmodel` estmated with one of the models in
#'    package.
#' 
#' @param ... Additional arguments passed to methods.
#' 
#' @returns A vector with the fitted values of the model.
#' 
#' @importFrom stats fitted
#' @export
fitted.mlmodel <- function(object, ...)
{
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must be a model of class 'mlmodel'.",
                   call = NULL)
  
  
  if (!is.null(object$model$fitted.values))
    return(object$model$fitted.values)
  
  # Fallback to predict if not stored
  predict(object, type = "response")
}

#' @rdname fitted.mlmodel
#' @importFrom stats fitted.values
#' @export
fitted.values.mlmodel <- fitted.mlmodel

## LOGLIK ======================================================================
# --- General ------------------------------------------------------------------
#' Return the Log-likelihood value.
#'
#' @param object Object of class `mlmodel` or `summary.mlmodel`, estimated with
#' one of the estimators in this package.
#'
#' @param ... Additional arguments to methods.
#'
#' @returns A scalar numeric: the log-likelihood of the model. It includes the
#' number of observations as the attribute 'nobs'.
#'
#' @importFrom stats logLik
#' @export
logLik <- function(object, ...) UseMethod("logLik")

# --- mlmodel ------------------------------------------------------------------
#' @rdname logLik
#' @export
logLik.mlmodel <- function(object, ...)
{
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must inherit from class 'mlmodel'.",
                   call = NULL)
  
  ll <- object$maximum
  
  # Attach useful attributes
  attr(ll, "nobs") <- object$model$n_used
  attr(ll, "df")   <- length(coef(object))   # number of free parameters
  
  ll
}

# --- summary.mlmodel ----------------------------------------------------------
#' @rdname logLik
#' @export
logLik.summary.mlmodel <- function(object, ...)
{
  if(!inherits(object, "summary.mlmodel"))
    cli::cli_abort("`object` must inerit from class `summary.mlmodel`.",
                   call = NULL)
  if (is.null(object$logLik) || is.na(object$logLik)) {
    cli::cli_warn("Log-likelihood value is not available.")
    return(NA_real_)
  }
  
  ll <- object$logLik
  
  attr(ll, "nobs") <- object$nobs
  attr(ll, "df") <- nrow(object$coefficients)
  
  ll
}

## MLE FUNCTIONS ===============================================================
# -- gradientObs ----------------------------------------------------------------
#' Gradient by observation
#' 
#' Get the gradients (scores), evaluated at the optimum from an object of class
#' `mlmodel`.
#' 
#' @param object An object of class `mlmodel`, i.e. an model estimated with one
#'               of the estimators in the package.
#' 
#' @returns Matrix with each observation's gradient (score) in a row.
#' 
#' @author Alfonso Sanchez-Penalver
#' 
#' @export
gradientObs <- function(object)
  UseMethod("gradientObs")

#' @rdname gradientObs
#' @export
gradientObs.mlmodel <- function(object)
{
  if(!inherits(object, "mlmodel"))
    cli::cli_abort("`object` needs to be of class `'mlmodel'`")
  object$model$functions$gradientObs(object)
}

# -- hessianObs ----------------------------------------------------------------
#' Hessian matrix by observation
#' 
#' Get the hessian matrix, evaluated at the optimum, by observation from an
#' object of class `mlmodel`.
#' 
#' @param object An object of class `mlmodel`, i.e. an model estimated with one
#'               of the estimators in the package.
#' 
#' @returns Matrix with the evaluated Hessian matrices per observation.
#' 
#' @details
#' If there are N observations in the estimation, and the model has K parameters,
#' this method returns a (N * K) by K matrix. The Hessian for each observation
#' stacked on top of each other.
#' 
#' @author Alfonso Sanchez-Penalver
#' 
#' @export
hessianObs <- function(object)
  UseMethod("hessianObs")

#' @rdname hessianObs
#' @export
hessianObs.mlmodel <- function(object)
{
  if(!inherits(object, "mlmodel"))
    cli::cli_abort("`object` needs to be of class `'mlmodel'`")
  object$model$functions$hessianObs(object)
}

# -- loglikeObs ----------------------------------------------------------------
#' Log-Likelihood by observation
#' 
#' Get the log-likelihood, evaluated at the optimum from an object of class
#' `mlmodel`.
#' 
#' @param object An object of class `mlmodel`, i.e. an model estimated with one
#'               of the estimators in the package.
#' 
#' @returns Vector with the evaluated log-likelihoods per observation.
#' 
#' @author Alfonso Sanchez-Penalver
#' 
#' @export
loglikeObs <- function(object)
  UseMethod("loglikeObs")

#' @rdname loglikeObs
#' @export
loglikeObs.mlmodel <- function(object)
{
  if(!inherits(object, "mlmodel"))
    cli::cli_abort("`object` needs to be of class `'mlmodel'`")
  object$model$functions$loglikeObs(object)
}

# NOBS =========================================================================
#' Returns the number of observations used in an estimation of an `mlmodel` model.
#' 
#' @param object An `mlmodel` estimation model.
#' @param ... Not currently implemented.
#' 
#' @returns A value with the number of observations used in the estimation.
#' 
#' @importFrom stats nobs
#' @export
nobs.mlmodel <- function(object, ...) {
  if(!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must be of class `'mlmodel'`",
                   call = NULL)
  n <- object$model$n_used %||% 
    object$nobs %||% 
    length(object$model$sample) %||% 
    length(object$model$usable_sample)
  
  if (is.null(n) || n == 0) {
    cli::cli_warn("Could not determine number of observations.")
    return(NA_integer_)
  }
  
  as.integer(n)
}

## NULL DEFAULT ================================================================
#' Null default operator
#'
#' @description
#' Provides a convenient infix operator for replacing `NULL` values.
#' `x \%||\% y` is equivalent to `if (is.null(x)) y else x`.
#'
#' @param x,y Any R objects.
#'
#' @name null-default
#' @usage x \%||\% y
#' @export
#' @importFrom rlang %||%
#' @examples
#' NULL %||% "fallback"
#' list(a = 1) %||% list(b = 2)
`%||%` <- rlang::`%||%`

## PREDICT GENERIC =============================================================
#' Predictions from mlmodel models
#'
#' Generic method for computing predictions from models fitted with the
#' `mlmodels` package.
#' 
#' @param object An object from an estimation with one of our models.
#' @param newdata Optional data frame for out-of-sample predictions.
#' @param type Character string indicating what to predict. See **Details**.
#' @param se.fit Logical. If `TRUE`, also return standard errors (delta method).
#' @param vcov Optional user-supplied variance-covariance matrix.
#' @param vcov.type Type of variance-covariance matrix. See [vcov][vcov.mlmodel].
#' @param cl_var Clustering variable (name or vector).
#' @param repetitions Number of bootstrap replications when `vcov.type = "boot"`.
#' @param seed Random seed for bootstrapping, for reproducibility.
#' @param progress Logical. Show bootstrap/jackknife progress bar? Default is
#'    `FALSE` in higher-level functions.
#' @param ... Additional arguments passed to methods.
#'
#' @returns An object that inherits from `predict.mlmodel` and has two elements:
#' \describe{
#'    \item{fit}{Vector with the predictions.}
#'    \item{se.fit}{If `se.fit` is `TRUE` a vector with the delta-method standard
#'    errors, using analytical gradients. If `se.fit` is `FALSE`, it is set to
#'    `NULL`.}
#' }
#' 
#' @method predict mlmodel
#' @export
predict.mlmodel <- function(object, ...) {
  UseMethod("predict")
}

## SE ==========================================================================
# --- Generic ------------------------------------------------------------------
#' Extracts standard errors for an `mlmodel` object.
#'
#' @param object An object of class `"mlmodel"` or any model that inherits
#'   from it (e.g. `"ml_lm"`).
#' @param vcov An optional user-supplied variance-covariance matrix.
#' @param vcov.type Character string specifying the type of variance-covariance
#'   matrix to use. One of `"oim"` (default), `"robust"`, `"opg"`, `"cluster"`,
#'   or `"boot"`.
#' @param cl_var Character string or vector. Name of the clustering variable
#'   or the vector itself. Only used when `vcov.type = "cluster"`.
#' @param repetitions Integer. Number of bootstrap replications to use when
#'   `vcov.type = "boot"`. Default is 999.
#' @param seed Integer. Random seed for reproducibility when `vcov.type = "boot"`.
#'   If `NULL`, a random seed is generated.
#' @param progress Logical. Should a progress bar be displayed during
#'   bootstrapping? Default is `FALSE` (silent) when called from `summary()`.
#' @param ... Not currently used.
#'
#' @return A named numeric vector of standard errors.
#'
#' @seealso [vcov.mlmodel]
#'
#' @author Alfonso Sanchez-Penalver
#' 
#' @export
se <- function(object, ...) {
  UseMethod("se")
}

# --- mlmodel ------------------------------------------------------------------
#' Extracts standard errors for an `mlmodel` object.
#' 
#' @rdname se
#' @export
se.mlmodel <- function(object,
                       vcov = NULL,
                       vcov.type = "oim",
                       cl_var = NULL,
                       repetitions = 999,
                       seed = NULL,
                       progress = FALSE,
                       ...)
{
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must be of 'mlmodel' class.")

  var <- .process_vcov(object,
                       vcov = vcov,
                       vcov.type   = vcov.type,
                       cl_var      = cl_var,
                       repetitions = repetitions,
                       seed        = seed,
                       progress    = progress)
  se <- sqrt(diag(var))
  names(se) <- colnames(var)
  return(se)
}

## SUMMARY GENERIC =============================================================
#' Summary for ml_lm objects
#'
#' @param object A fitted model object of class `"ml_lm"`.
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
#' Coefficient names in the fitted object use the prefixes `value::` and
#' `scale::` to identify to which equation they belong to, and to avoid
#' confusion when the same variable(s) appear(s) in both the value and scale
#' equations.
#' 
#' @author Alfonso Sanchez-Penalver
#' 
#' @method summary mlmodel
#' @export
summary.mlmodel <- function(object,
                            correlation = FALSE,
                            vcov = NULL,           # User-supplied variance matrix
                            vcov.type = "oim",
                            cl_var = NULL,
                            repetitions = 999,
                            seed = NULL,
                            progress = FALSE,
                            ...)
{
  UseMethod("summary")
}

# TERMS ========================================================================
#' Extract terms from mlmodel objects
#'
#' Returns the terms object for the value (mean/probability) equation.
#' For heteroskedastic models, the scale equation is **not** included in the
#' returned terms object.
#' 
#' @param x An `mlmodel` object from an estimation of one of our models.
#' @param ... Not currently implemented.
#'
#' @note The `terms()` method is provided for minimal compatibility with other
#' packages, but it is incomplete for heteroskedastic models. For reliable
#' bootstrap standard errors, use `vcov(object, type = "boot")` instead of
#' functions from the sandwich package.
#'
#' @export
terms.mlmodel <- function(x, ...) {
  if (!is.null(x$model$value$terms)) {
    x$model$value$terms
  } else {
    NULL
  }
}

## UPDATE GENERIC ==============================================================
#' Update method for ml_lm objects
#' 
#' @param object An `mlmodel` estimation object.
#' @param formula. The formula of the value equation (optional).
#' @param scale. The formula of the scale equation (optional).
#' @param data A data.frame with the data to do the estimation (optional).
#' @param weights A vector with the weights (optional).
#' @param ... Currently not implemented.
#' @param evaluate Should the updated call be evaluated? Defaults to `TRUE`.
#'
#' @details
#' Used in bootstrapping of tests. Re-evaluates the original call with arguments
#' passed into the function.
#' 
#' **Note on sandwich::vcovBS()**: This function does not work reliably with 
#' `mlmodel` objects, even in simple homoskedastic cases.  We, therefore, built
#' our own bootstrap implementation. We strongly recommend using
#' `vcov(object, type = "boot")` instead. See [vcov][mlmodels::vcov.mlmodel].
#' 
#' @author Alfonso Sanchez-Penalver
#' @export
update.mlmodel <- function(object,
                           formula. = NULL,
                           scale. = NULL,
                           data = NULL,
                           weights = NULL,
                           ...,
                           evaluate = TRUE)
{
  UseMethod("update")
}

## VARIANCE ====================================================================
#' Variance-Covariance Matrix for mlmodel Objects
#'
#' Returns the variance-covariance matrix of the estimated parameters
#' using different methods.
#'
#' @param object An object of class `"mlmodel"` or a model inheriting from it
#'   (e.g. `"ml_lm"`).
#' @param type Character string specifying the type of variance-covariance
#'   matrix. One of `"oim"` (default), `"robust"`, `"opg"`, `"cluster"`,
#'   `"boot"`, and `"jack"` or `"jackknife"`.
#' @param cl_var Character string or vector. Name of the clustering variable
#'   in the data, or the vector itself.
#' @param repetitions Integer. Number of bootstrap replications to use when
#'   `type = "boot"`. Default is 999.
#' @param seed Integer. Random seed for reproducibility when `type = "boot"`.
#'   If `NULL`, a random seed is generated.
#' @param progress Logical. Should a progress bar be displayed? Default is
#'   `TRUE` when `type` is `"boot"` or `"jack"`/`"jackknife"`. Ignored for other
#'   types.
#' @param ... Further arguments passed to methods.
#'
#' @return A symmetric variance-covariance matrix with coefficient names
#'   on the rows and columns.
#' 
#' @details
#' Type `"cluster"` is an alias for `"robust"`, but it requires `cl_var` to be
#' passed. 
#' 
#' When `type` is set to `"robust"`, and `cl_var` is `NULL`, a standard
#' robust (sandwich) variance is returned.
#' 
#' When `type` is set to `"robust"` and `cl_var` is not `NULL`, a clustered-robust
#' variance is returned.
#' 
#' Type `"jackknife"` is an alias for type `"jack"`.
#' 
#' @author Alfonso Sanchez-Penalver
#'
#' @export
vcov.mlmodel <- function(object,
                         type = "oim",
                         cl_var = NULL,
                         repetitions = 999,
                         seed = NULL,
                         progress = TRUE,
                         ...)
{
  # 1. Inheritance check - must be first
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must be a model of class 'mlmodel'.", call = NULL)

  # 2. Validate and normalize type
  type <- rlang::arg_match(type, c("oim", "robust", "opg", "cluster", "boot", "jack", "jackknife"))

  if(type == "cluster" && is.null(cl_var))
    cli::cli_abort("`cl_var` cannot be null with 'cluster' `type`.",
                   call = NULL)

  # cluster is an alias for robust, since they must provide cl_var
  if (type == "cluster") type <- "robust"
  
  # jackknife is an alias for jack
  if (type %in% c("jack", "jackknife")) type <- "jack"

  # 3. Early validation for cl_var
  if (!is.null(cl_var) && !(type %in% c("robust", "boot", "jack")))
    cli::cli_abort(
      "`cl_var` can only be used when `type` is 'cluster', 'robust', 'boot', or 'jack'.",
      call = NULL
    )

  # 4. Process clustering variable if provided
  if (!is.null(cl_var)) {
    if (is.character(cl_var)) {
      # Retrieve data (new primary path first, then old d_name fallback)
      if (!is.null(object$model$data) && is.data.frame(object$model$data)) {
        d <- object$model$data
      } else if (!is.null(object$model$d_name) && object$model$d_name != "<unknown data>") {
        d <- tryCatch(get(object$model$d_name), error = function(e) {
          cli::cli_abort("Cannot retrieve the dataset to get the clustering variable.",
                         call = NULL)
        })
      } else {
        cli::cli_abort("Dataset and its name not stored; cannot retrieve clustering variable.",
                       call = NULL)
      }
      cl_var <- d[[cl_var]][object$model$sample]
    } else {
      # User passed a vector directly
      n_var  <- length(cl_var)
      n_orig <- object$model$n_orig
      n_used <- object$model$n_used
      if (n_var != n_orig && n_var != n_used)
        cli::cli_abort(
          c("The clustering vector has the wrong number of observations.",
            "It should have either {.val {n_orig}} (original) or {.val {n_used}} (used in estimation) observations."),
          call = NULL
        )
      if (n_var == n_orig)
        cl_var <- cl_var[object$model$sample]
    }

    # Check for unusable observations in clustering variable
    if (any(!complete.cases(cl_var))) {
      n_bad <- sum(!complete.cases(cl_var))
      cli::cli_abort(
        c("The clustering variable has unusable observations (NA or NaN).",
          "Found {.val {n_bad}} bad observation{?s} in the estimation sample."),
        call = NULL
      )
    }
  }

  # -- 5. Now, we're ready to select the method.
  if (type == "boot")
  {
    vcov_mat <- .vcov_boot(object,
                           repetitions = repetitions,
                           seed = seed,
                           cl_var = cl_var,
                           progress = progress,
                           ...)
    attr(vcov_mat, "vcov.type") <- type
    if (!is.null(cl_var)) {
      attr(vcov_mat, "clustered") <- TRUE
      attr(vcov_mat, "cluster.var") <- cl_var
    }
    attr(vcov_mat, "rep") <- repetitions
    return(vcov_mat)
  }
  
  if (type == "jack")
  {
    vcov_mat <- .vcov_jack(object,
                           cl_var = cl_var,
                           progress = progress,
                           ...)
    attr(vcov_mat, "vcov.type") <- type
    if (!is.null(cl_var)) {
      attr(vcov_mat, "clustered") <- TRUE
      attr(vcov_mat, "cluster.var") <- cl_var
    }
    return(vcov_mat)
  }

  # Regular (non-bootstrap or jackknife) variance types
  H <- object$hessian
  if (is.null(H))
    cli::cli_abort("Hessian is missing from the model object.", call = NULL)

  # Eigen decomposition
  eig <- eigen(-H, symmetric = TRUE)
  eigvals <- eig$values

  # Singularity warning for OIM
  if (type == "oim") {
    n_neg <- sum(eigvals < 0)
    if (n_neg > 0) {
      cli::cli_abort(c("Hessian is not negative semidefinite.",
                       "i" = "It has {.val {n_neg}} positive eigenvalue{?s}."),
                     call = NULL)
    }

    # Relative and absolute singularity checks
    rel_sing <- min(eigvals) / max(eigvals) < 1e-12
    abs_sing <- which(abs(eigvals) < 1e-6)

    if (rel_sing || length(abs_sing) > 0) {
      if (is.null(object$.oim_singularity_warned)) {
        pnames <- names(coef(object))
        eigvecs <- eig$vectors

        cli::cli_alert_warning("OIM variance may be unreliable due to singularity in the Hessian.")

        if (rel_sing) {
          max_idx <- which.max(eigvals)
          min_idx <- which.min(eigvals)
          cli::cli_alert_info(
            "Relative singularity: {.val {pnames[max_idx]}} has the largest eigenvalue ({round(eigvals[max_idx], 3)}), \\
             while {.val {pnames[min_idx]}} has the smallest ({round(eigvals[min_idx], 6)})."
          )
        }

        if (length(abs_sing) > 0) {
          cli::cli_alert_info("Near-zero eigenvalues detected. The following linear combinations are nearly flat:")

          for (i in abs_sing) {
            par_dim <- pnames[i]
            vec <- eigvecs[, i]
            large_idx <- which(abs(vec) > 0.2)

            if (length(large_idx) > 0) {
              terms <- paste0(round(vec[large_idx], 3), " * ", pnames[large_idx])
              combo <- paste(terms, collapse = " + ")
              cli::cli_alert_info(
                "  {.val {par_dim}} (eigenvalue = {.val {round(eigvals[i], 6)}}): {.code {combo}}"
              )
            }
          }
        }

        cli::cli_alert_info("Consider using `type = 'robust'` or revising the model specification.")
      }
    }
  }

  V <- chol2inv(chol(-H))
  G <- object$gradientObs

  # Compute variance
  if (!is.null(cl_var)) {
    # Cluster-robust variance
    clusters <- unique(cl_var)
    n_cl <- length(clusters)
    S <- matrix(0, nrow = n_cl, ncol = ncol(G))
    for (i in seq_along(clusters)) {
      idx <- cl_var == clusters[i]
      S[i, ] <- colSums(G[idx, , drop = FALSE])
    }
    n <- n_cl
    vcov_mat <- V %*% crossprod(S) %*% V * n / (n - 1)
  } else if (type == "opg") {
    vcov_mat <- chol2inv(chol(crossprod(G)))
  } else if (type == "robust") {
    n <- nrow(G)
    vcov_mat <- V %*% crossprod(G) %*% V * n / (n - 1)
  } else {  # oim
    vcov_mat <- V
  }

  attr(vcov_mat, "vcov.type") <- type
  if (!is.null(cl_var)) {
    attr(vcov_mat, "clustered") <- TRUE
    attr(vcov_mat, "cluster.var") <- cl_var
  }
  dimnames(vcov_mat) <- list(names(coef(object)), names(coef(object)))
  return(vcov_mat)
}