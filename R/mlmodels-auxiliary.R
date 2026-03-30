# =============================================================================
# mlmodels: Internal helper functions for mlmodel objects
# =============================================================================

## Private files that get called from several mlmodels top level functions to
## do a specific task.

# ------------------------------------------------------------------------------
# Cluster information helper
# ------------------------------------------------------------------------------

#' Internal helper to extract clustering information
#'
#' Used by `vcov.mlmodel()` when `vcov.type = "cluster"` or `"cluster-boot"`.
#'
#' @param object An `mlmodel` object.
#' @param cl_var Character string or vector. The clustering variable.
#'
#' @return A list with `var_name`, `n_cluster`, and `ids`.
#'
#' @keywords internal
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

# ------------------------------------------------------------------------------
# get_vcov helper
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Bootstrap variance helpers
# ------------------------------------------------------------------------------

#' Generic for bootstrap variance-covariance matrix
#'
#' @keywords internal
vcov_boot <- function(object, ...) {
  UseMethod("vcov_boot")
}

#' Bootstrap Variance-Covariance Matrix (mlmodel method)
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


# POST-MOLD UTILITIES ==========================================================

## Factor Mapping --------------------------------------------------------------
#' Build factor mapping for a single equation
#'
#' Internal function. Creates a mapping of factor variables to their
#' main effect columns and interaction columns after `hardhat::mold()`.
#'
#' @param mold A mold object returned by `hardhat::mold()`.
#' @param equation_name Character string indicating the equation name
#'   ("value" or "scale").
#'
#' @return A named list containing the factor mapping for this equation.
#'
#' @keywords internal
build_factor_map <- function(mold, equation_name = "value")
{
  if (is.null(mold) || is.null(mold$predictors) || ncol(mold$predictors) == 0)
    return(list())

  actual_colnames <- colnames(mold$predictors)
  ptypes <- mold$blueprint$ptypes$predictors

  # Initialize as empty list - this is the key fix
  mapping <- list()

  # === Dangerous transformation check ===
  actual_colnames <- colnames(mold$predictors)
  if (length(actual_colnames) == 0) return(list())

  ptypes <- mold$blueprint$ptypes$predictors

  # Find columns with dangerous transformations
  dangerous_names <- actual_colnames[grepl("I\\(|poly\\(", actual_colnames)]
  if (length(dangerous_names) > 0) {

    # Check variables to identify factor variable.
    for (pt_name in names(ptypes)) {
      if (!is.factor(ptypes[[pt_name]])) next

      # Which dangerous columns mention this factor?
      involved <- grepl(pt_name, dangerous_names, fixed = TRUE)
      if (!any(involved)) next

      # Check if any other variable names also contain this factor name
      other_names <- names(ptypes)[names(ptypes) != pt_name]
      other_involved <- grepl(pt_name, other_names, fixed = TRUE)

      if (any(other_involved)) {
        names_to_check <- other_names[other_involved]

        count_factor <- sapply(dangerous_names, function(col) {
          length(gregexpr(pt_name, col, fixed = TRUE)[[1]])
        })

        # Get the count factor to be the same dimensions as the l_matrix
        count_factor <- as.matrix(count_factor)

        # Matrix: rows = dangerous_names, columns = other variables
        l_matrix <- sapply(names_to_check, function(on) {
          grepl(on, dangerous_names, fixed = TRUE)
        })

        f_matrix <- !l_matrix & count_factor < 2

        if(length(f_matrix) == 1)
          caused_by_factor <- f_matrix
        else
          caused_by_factor <- rowSums(f_matrix) > 0

        if (any(caused_by_factor))
        {
          bad_cols <- dangerous_names[caused_by_factor]

          cli::cli_abort(c("Complex transformation involving factor '{pt_name}' in the {.strong {equation_name}} equation is not supported.",
                           "i" = "Found: {.val {paste(bad_cols, collapse = ', ')}}",
                           "i" = "Please use only main effects or simple interactions (* or :) with factors.",
                           "x" = "Mixed interactions like I((factor * continuous)^2) aren't allowed.",
                           "i" = "The proper syntax would be factor * I(continuous)^2.")
          )
        }
      }
      else
      {
        # All dangerous columns are caused by this factor
        cli::cli_abort(c("Complex transformation involving factor '{pt_name}' in the {.strong {equation_name}} is not supported.",
                         "i" = "Found: {.val {paste(bad_cols, collapse = ', ')}}",
                         "i" = "Please use only main effects or simple interactions (* or :) with factors.")
        )
      }
    }
  }

  # === 2. Build mapping for main effects and interactions ===
  for (var_name in names(ptypes)) {
    if (!is.factor(ptypes[[var_name]]))
    {
      mapping[[var_name]] <- NULL
      next
    }

    levels <- levels(ptypes[[var_name]])

    main_effect_cols <- character(0)
    interaction_cols <- character(0)

    for (lvl in levels) {
      # Main effect: exact match "varnamelevel"
      main_pattern <- paste0("^", var_name, lvl, "$")
      main_matches <- actual_colnames[grepl(main_pattern, actual_colnames)]
      main_effect_cols <- c(main_effect_cols, main_matches)

      # Interaction: contains "varnamelevel:" or ":varnamelevel"
      inter_pattern <- paste0(var_name, lvl, ":|:", var_name, lvl)
      inter_matches <- actual_colnames[grepl(inter_pattern, actual_colnames)]
      interaction_cols <- c(interaction_cols, inter_matches)
    }

    dummy_cols <- unique(c(main_effect_cols, interaction_cols))

    # Determine base level (NULL if all levels are present = no intercept case)
    possible_main_dummies <- paste0(var_name, levels)
    missing_levels <- setdiff(possible_main_dummies, main_effect_cols)

    base_level <- if (length(missing_levels) == 1) {
      sub(var_name, "", missing_levels[1])
    } else if (length(missing_levels) == 0) {
      NULL   # All levels present → no intercept
    } else {
      cli::cli_alert_warning(
        "Unexpected number of missing levels for factor '{var_name}' in equation '{equation_name}'."
      )
      NULL
    }
    mapping[[var_name]] <- list(
      equation         = equation_name,
      var_name         = var_name,
      levels           = levels,
      main_effect_cols = main_effect_cols,     # ← separate
      interaction_cols = interaction_cols,     # ← separate
      dummy_cols       = dummy_cols,           # combined for convenience
      base_level       = base_level,
      n_dummies        = length(dummy_cols),
      has_interaction  = length(interaction_cols) > 0
    )
  }

  mapping
}

## Build factor mappings for all equations -------------------------------------
#' Build factor mappings for all equations
#'
#' Internal function. Takes a named list of mold objects and returns
#' factor mappings for each equation. Also performs safety checks
#' for all-NA and invalid predictor columns.
#'
#' @param molds A named list of mold objects (e.g. `list(value = ..., scale = ...)`).
#'
#' @return A named list with the same names as `molds`, each containing
#'   the factor mapping for that equation.
#'
#' @keywords internal
build_factor_mapping <- function(molds)
{
  if (!is.list(molds) || is.null(names(molds))) {
    cli::cli_abort("`molds` must be a named list of mold objects.",
                   call = NULL)
  }

  mapping <- list()

  for (eq_name in names(molds)) {
    mold <- molds[[eq_name]]

    # Check for all-NA columns first (fail fast)
    check_for_invalid_predictors(mold, equation_name = eq_name)

    # Only map factors if we passed the check
    mapping[[eq_name]] <- build_factor_map(mold, equation_name = eq_name)
  }

  mapping
}

## Invalid Predictors ----------------------------------------------------------
#' Check for invalid predictor columns after molding
#'
#' Internal helper that aborts with a clear error if any predictor
#' column contains invalid values (NA, NaN, Inf, or -Inf).
#'
#' @param mold A mold object returned by `hardhat::mold()`.
#' @param equation_name Character string indicating the equation name
#'   ("value" or "scale").
#'
#' @return Invisibly returns `TRUE` if the check passes.
#'
#' @keywords internal
check_for_invalid_predictors <- function(mold, equation_name = "value")
{
  if (is.null(mold) || is.null(mold$predictors) || ncol(mold$predictors) == 0) {
    return(invisible(TRUE))
  }

  predictors <- mold$predictors

  # Check each column for any invalid value
  has_invalid <- vapply(predictors, function(col) {
    any(is.na(col) | is.nan(col) | is.infinite(col))
  }, logical(1))

  if (any(has_invalid))
  {
    bad_cols <- colnames(predictors)[has_invalid]
    cli::cli_abort(c("Column(s) in the {.strong {equation_name}} equation contain invalid values (NA, NaN, or Inf).",
                     "i" = "Found: {.val {paste(bad_cols, collapse = ', ')}}",
                     "x" = "This usually happens with invalid transformations involving factors.",
                     "i" = "Please revise your formula.")
    )
  }

  invisible(TRUE)
}


## Log Transformations ---------------------------------------------------------
#' Detect log transformation on the outcome variable
#'
#' Returns a consistent list structure whether a log transformation is detected or not.
#'
#' @param formula a formula with a lhs expression to check for.
#' @param data A data.frame that holds the variables in the estimation.
#'
#' @returns A list with all the information of a log transformation.
#'
#' @keywords internal
detect_log_transformation <- function(formula, data) {

  lhs_expr <- rlang::f_lhs(formula)

  # No left-hand side → not a log transformation
  if (is.null(lhs_expr)) {
    return(list(
      is_log        = FALSE,
      log_fun       = NULL,
      shift         = 0,
      multiplier    = 1,
      complex_inner = FALSE,
      var_names     = character(0),
      n_invalid     = 0L,
      invalid_idx   = rep(FALSE, nrow(data)),
      raw_expr      = NULL
    ))
  }

  # Get the expression on the LHS
  expr <- rlang::quo_get_expr(rlang::quo(!!lhs_expr))

  # Not a log function → not a log transformation
  if (!is.call(expr) || !as.character(expr[[1]]) %in% c("log", "log10", "log1p")) {
    return(list(
      is_log        = FALSE,
      log_fun       = NULL,
      shift         = 0,
      multiplier    = 1,
      complex_inner = FALSE,
      var_names     = character(0),
      n_invalid     = 0L,
      invalid_idx   = rep(FALSE, nrow(data)),
      raw_expr      = NULL
    ))
  }

  # ── Log transformation detected ─────────────────────────────────────
  fun_name <- as.character(expr[[1]])
  inner_expr <- expr[[2]]

  shift <- if (fun_name == "log1p") 1 else 0
  multiplier <- 1
  complex_inner <- FALSE
  var_names <- character(0)

  # Check for operations inside the log: y + c, c * y, y + z, y * z, etc.
  if (is.call(inner_expr) && length(inner_expr) == 3) {
    op <- as.character(inner_expr[[1]])
    left  <- inner_expr[[2]]
    right <- inner_expr[[3]]

    left_is_var  <- is.symbol(left)
    right_is_var <- is.symbol(right)

    if (left_is_var && right_is_var) {
      # Complex case: two variables (e.g. log(y + z))
      complex_inner <- TRUE
      var_names <- c(as.character(left), as.character(right))
      shift <- 0
      multiplier <- 1
    } else if (op == "+") {
      # Simple shift: log(y + c)
      if (is.numeric(right)) shift <- right
      else if (is.numeric(left)) shift <- left
    } else if (op == "*") {
      # Simple multiplier: log(c * y)
      if (is.numeric(right)) multiplier <- right
      else if (is.numeric(left)) multiplier <- left
    }
  }

  # Evaluate inner expression to count invalid observations
  inner_val <- tryCatch(
    rlang::eval_tidy(rlang::quo(!!inner_expr), data = data),
    error = function(e) rep(NaN, nrow(data))
  )

  if (fun_name %in% c("log", "log10")) {
    invalid <- inner_val <= 0
  } else if (fun_name == "log1p") {
    invalid <- inner_val <= -1
  } else {
    invalid <- rep(FALSE, nrow(data))
  }

  n_invalid <- sum(invalid, na.rm = TRUE)

  if (complex_inner) {
    cli::cli_alert_warning(
      "Complex transformation inside {fun_name}() detected: {deparse(inner_expr)}. \\
       Full back-transformation support may not be available yet."
    )
  }

  list(
    is_log        = TRUE,
    log_fun       = fun_name,
    shift         = shift,
    multiplier    = multiplier,
    complex_inner = complex_inner,
    var_names     = var_names,
    n_invalid     = n_invalid,
    invalid_idx   = invalid,
    raw_expr      = deparse(inner_expr)
  )
}

#' Detect log transformations across multiple formulas
#'
#' Returns a named list with detection results for each formula.
#' Useful for future multi-equation models.
#'
#' @param formulas A list with the different formulas for the different
#' equations that need checking.
#' @param data A data.frame that holds all the data in the estimation.
#'
#' @returns A list with the same number of elements as in `formulas_list`, each
#' with an inner list with all the information about the log-transformation for
#' the respective equation.
#'
#' @keywords internal
detect_log_transformations <- function(formulas_list, data) {
  if (!is.list(formulas_list) || is.null(names(formulas_list))) {
    cli::cli_abort("`formulas_list` must be a named list of formulas.",
                   call = NULL)
  }

  lapply(formulas_list, function(f) {
    if (is.null(f)) return(list(is_log = FALSE))
    detect_log_transformation(f, data)
  })
}

