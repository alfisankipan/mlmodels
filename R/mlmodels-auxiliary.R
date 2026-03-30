## Private files that get called from several mediamodels top level functions to
## do a specific task.

## Function to pull the cluster info.
vcov_cluster_info <- function(object, cl_var)
{
  if(!inherits(object, "mlmodel"))
    cli::cli_abort("`object` needs to be of class `'mlmodel'`",
                   call = NULL)
  if (is.character(cl_var))
  {
    # User passed variable name
    d <- get(object$model$data)
    cl_vec <- d[[cl_var]][object$model$sample]
    vcov.cluster <- list(var_name = cl_var,
                         n_cluster = length(unique(cl_vec)),
                         ids = unique(cl_vec))
  }
  else
  {
    # User passed a vector directly
    cl_vec <- if (length(cl_var) == object$model$n_orig)
      cl_var[object$model$sample]
    else
      cl_var
    vcov.cluster <- list(var_name = NULL,
                         n_cluster = length(unique(cl_vec)),
                         ids = unique(cl_vec))
  }

  return(vcov.cluster)
}

#' Internal helper to obtain the variance-covariance matrix
#'
#' This helper is used by `summary.ml_lm()`, `waldtest.mlmodel()`,
#' `predict.ml_lm()` and other functions.
#'
#' If a matrix is supplied via the `vcov` argument, it is validated and
#' returned directly. Otherwise, it calls `vcov.mlmodel()`.
#'
#' @param object An `mlmodel` object.
#' @param vcov An optional user-supplied variance-covariance matrix.
#'   If provided, it will be used instead of computing a new one.
#' @param vcov.type Character string specifying the type of variance-covariance
#'   matrix to compute. One of `"oim"` (default), `"robust"`, `"opg"`,
#'   `"cluster"`, or `"boot"`.
#' @param cl_var Character string or vector. Clustering variable.
#'   Only used when `vcov.type = "cluster"`.
#' @param repetitions Integer. Number of bootstrap replications when
#'   `vcov.type = "boot"`. Default is 999.
#' @param seed Integer. Random seed for reproducibility when `vcov.type = "boot"`.
#'   If `NULL`, a random seed is generated.
#' @param ... Further arguments passed to [vcov.mlmodel()].
#' @param progress Logical. Should a progress bar be displayed during
#'   bootstrapping? Default is `FALSE`.
#'
#' @return A variance-covariance matrix.
#'
#' @keywords internal
get_vcov <- function(object,
                     vcov = NULL,           # User-supplied variance matrix
                     vcov.type = "oim",
                     cl_var = NULL,
                     repetitions = 999,
                     seed = NULL,
                     progress = FALSE,
                     ...)
{
  # Case 1: User passed a pre-computed variance matrix
  if (!is.null(vcov)) {
    if (!is.matrix(vcov))
      cli::cli_abort("`vcov` must be a matrix when supplied directly.",
                     call = NULL)

    n_coef <- length(coef(object))

    # Check dimensions
    if (!identical(dim(vcov), c(n_coef, n_coef)))
      cli::cli_abort(c("Supplied `vcov` matrix has incorrect dimensions ({nrow(vcov)} x {ncol(vcov)}).,
                       Expected {n_coef} x {n_coef} (matching number of coefficients)."),
                     call = NULL)

    # Optional: check row/column names
    coef_names <- names(coef(object))
    if (!is.null(rownames(vcov)) && !identical(rownames(vcov), coef_names))
      cli::cli_warn("Row names of supplied `vcov` do not match coefficient names.",
                    call = NULL)

    if (!is.null(colnames(vcov)) && !identical(colnames(vcov), coef_names))
      cli::cli_warn("Column names of supplied `vcov` do not match coefficient names.",
                    call = NULL)

    return(vcov)
  }

  # Case 2: Generate variance using vcov.mlmodel()
  vcov(object,
       type        = vcov.type,
       cl_var      = cl_var,
       repetitions = repetitions,
       seed        = seed,
       progress    = progress)
}


# VCOV_BOOT --------------------------------------------------------------------
# Generic
#' @keywords internal
vcov_boot <- function(object, ...) {
  UseMethod("vcov_boot")
}

#' Bootstrap Variance-Covariance Matrix
#'
#' Internal function to compute bootstrapped variance-covariance matrix.
#' Called by `vcov.mlmodel()` when `type = "boot"`.
#'
#' @param object An `mlmodel` object.
#' @param repetitions Number of bootstrap replications.
#' @param seed Random seed for reproducibility.
#' @param cl_var Clustering variable (if clustered bootstrap).
#' @param progress Logical. Whether to show progress bar.
#' @param ... Not currently used.
#'
#' @keywords internal
vcov_boot.mlmodel <- function(object,
                              repetitions = 999,
                              seed = NULL,
                              cl_var = NULL,
                              progress = TRUE,
                              ...)
{
  if(!inherits(object, "mlmodel"))
    cli::cli_abort("`object` needs to be of class 'mlmodel'.")
  if (is.null(seed)) seed <- sample.int(1e6, 1)
  set.seed(seed)

  # Recover original data
  original_data <- tryCatch({
    if (!is.null(object$call$data)) {
      eval(object$call$data, envir = parent.frame(2))
    } else if (!is.null(object$model$data) && object$model$data != "<unknown data>") {
      get(object$model$data, envir = .GlobalEnv)
    } else {
      stop("Could not recover original data", call. = FALSE)
    }
  }, error = function(e) {
    cli::cli_abort("Could not recover the original data for bootstrapping.", call = NULL)
  })

  # Use only the observations actually used in estimation
  used_data <- original_data[object$model$sample, , drop = FALSE]

  # Recover the final weights vector (always exists, length = n_used)
  w <- object$model$weights

  # ── Prepare for clustered bootstrap if requested ─────────────────────
  is_clustered <- !is.null(cl_var)
  if (is_clustered) {
    cluster_ids <- unique(cl_var[object$model$sample])
    n_cluster   <- length(cluster_ids)
  }

  # ── Bootstrap loop ───────────────────────────────────────────────────
  if (progress) {
    if (is_clustered) {
      cli::cli_alert_info("Clustered bootstrap with {.val {repetitions}} repetitions and {.val {n_cluster}} clusters.")
    } else {
      cli::cli_alert_info("Bootstrap with {.val {repetitions}} repetitions.")
    }
    cat(cli::col_blue(" 0"))
    for (i in seq(10, 50, by = 10)) cat(cli::col_blue(sprintf("%10d", i)))
    cat("\n")
    cat(cli::col_blue(strrep("=", 52), "\n"))
  }

  success     <- logical(repetitions)
  coef_matrix <- matrix(NA_real_, nrow = repetitions, ncol = length(coef(object)))

  for (i in seq_len(repetitions)) {
    if (progress && i %% 50 == 1 && i > 1) cat("\n ")
    else if(progress && i == 1) cat(" ")

    res <- tryCatch({
      if (is_clustered) {
        # Clustered bootstrap - sample whole clusters
        sampled_clusters <- sample(cluster_ids, size = n_cluster, replace = TRUE)
        boot_idx <- unlist(lapply(sampled_clusters, function(cid) {
          which(cl_var[object$model$sample] == cid)
        }))
      } else {
        # Regular (individual) bootstrap
        boot_idx <- sample(nrow(used_data), nrow(used_data), replace = TRUE)
      }

      boot_data <- used_data[boot_idx, , drop = FALSE]
      w_boot    <- w[boot_idx]                     # subset weights to bootstrap sample

      # Pass weights to update()
      updated <- update(object, data = boot_data, weights = w_boot, evaluate = TRUE)

      if (updated$code %in% c(1L, 2L, 8L)) {
        if (progress) cat(cli::col_green("."))
        success[i] <- TRUE
        coef(updated)
      } else {
        if (progress) cat(cli::col_red("x"))
        success[i] <- FALSE
        rep(NA_real_, length(coef(object)))
      }
    }, error = function(e) {
      if (progress) cat(cli::col_red("x"))
      success[i] <- FALSE
      rep(NA_real_, length(coef(object)))
    })

    coef_matrix[i, ] <- res
  }

  if (progress) {
    cat("\n")
    cat(cli::col_blue(strrep("=", 52), "\n"))
  }

  # Final summary
  success_rate <- mean(success) * 100
  cli::cli_text("Bootstrapping finished — {round(success_rate, 1)}% of replications converged.")

  if (success_rate < 70) {
    cli::cli_warn("Low convergence rate — bootstrap results may be unreliable.")
  }

  # Variance from successful replications only
  valid_rows <- complete.cases(coef_matrix)
  vcov_boot  <- var(coef_matrix[valid_rows, , drop = FALSE])

  dimnames(vcov_boot) <- list(names(coef(object)), names(coef(object)))
  vcov_boot
}
