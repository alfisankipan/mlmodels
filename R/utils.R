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

## Invalid Predictors ----------------------------------------------------------
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

#' Detect log transformation on the outcome variable
#'
#' Returns a consistent list structure whether a log transformation is detected or not.
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
