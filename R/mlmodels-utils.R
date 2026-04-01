# =============================================================================
# mlmodels: Exported utility functions for all mlmodel objects
# =============================================================================

## NULL DEFAULT ----------------------------------------------------------------
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

## PREDICT GENERIC -----------------------------------------------------------
#' Predictions for mlmodel objects
#'
#' Generic method for computing predictions from models fitted with the
#' `mlmodels` package.
#'
#' @param object A model object inheriting from `"mlmodel"`.
#' @param ... Further arguments passed to specific methods.
#'
#' @method predict mlmodel
#' @export
predict.mlmodel <- function(object, ...) {
  UseMethod("predict")
}

# COEFFICIENTS -----------------------------------------------------------------
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

# SE ---------------------------------------------------------------------------
#' Extract Standard Errors
#'
#' Computes standard errors for an `mlmodel` object.
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
#' @seealso [vcov.mlmodel()], [summary.ml_lm()]
#'
#' @author Alfonso Sanchez-Penalver
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

#' Extract Standard Errors
#'
#' Generic function to extract standard errors from a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Further arguments passed to methods.
#'
#' @export
se <- function(object, ...) {
  UseMethod("se")
}

# VARIANCE ---------------------------------------------------------------------
#' Variance-Covariance Matrix for mlmodel Objects
#'
#' Returns the variance-covariance matrix of the estimated parameters
#' using different methods.
#'
#' @param object An object of class `"mlmodel"` or a model inheriting from it
#'   (e.g. `"ml_lm"`).
#' @param type Character string specifying the type of variance-covariance
#'   matrix. One of `"oim"` (default), `"robust"`, `"opg"`, `"cluster"`,
#'   or `"boot"`.
#' @param cl_var Character string or vector. Name of the clustering variable
#'   in the data, or the vector itself.
#' @param repetitions Integer. Number of bootstrap replications to use when
#'   `type = "boot"`. Default is 999.
#' @param seed Integer. Random seed for reproducibility when `type = "boot"`.
#'   If `NULL`, a random seed is generated.
#' @param progress Logical. Should a progress bar be displayed during
#'   bootstrapping? Default is `TRUE`. Only relevant when `type = "boot"`.
#' @param ... Further arguments passed to methods.
#'
#' @return A symmetric variance-covariance matrix with coefficient names
#'   on the rows and columns.
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
  type <- rlang::arg_match(type, c("oim", "robust", "opg", "cluster", "boot"))

  if(type == "cluster" && is.null(cl_var))
    cli::cli_abort("`cl_var` cannot be null with 'cluster' `type`.",
                   call = NULL)

  # cluster is an alias for robust, since they must provide cl_var
  if (type == "cluster") type <- "robust"

  # 3. Early validation for cl_var
  if (!is.null(cl_var) && !(type %in% c("robust", "boot")))
    cli::cli_abort(
      "`cl_var` can only be used when `type` is 'cluster', 'robust' or 'boot'.",
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
    return(vcov_boot(object,
                     repetitions = repetitions,
                     seed = seed,
                     cl_var = cl_var,
                     progress = progress,
                     ...))

  # Regular (non-bootstrap) variance types
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

## LOGLIK ----------------------------------------------------------------------
#' Return the Log-likelihood value.
#'
#' @param object Object of class `mlmodel` or `summary.mlmodel`, estimated with
#' one of the estimators in this package.
#'
#' @param ... Additional arguments to methods.
#'
#' @returns A scalar numeric: the log-likelihood of the model. It has attribute(s)
#' 'df', the number of free parameters, and 'nobs' the number of observations if
#' `object` is of class `mlmodel`.
#'
#' @export
logLik.mlmodel <- function(object, ...)
{
  if (!inherits(object, c("mlmodel", "summary.mlmodel")))
    cli::cli_abort("`object` must be of either 'mlmodel' or 'summary.mlmodel' class.",
                   call = NULL)
  ll <- NextMethod("logLik", object, ...)
  if (inherits(object, "mlmodel"))
    if (!is.null(object$model$n_used))
      attr(ll, "nobs") <- object$model$n_used
  else
    attr(ll, "nobs") <- object$nobs
  return(ll)
}
