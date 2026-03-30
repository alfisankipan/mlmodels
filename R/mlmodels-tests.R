# =============================================================================
# mlmodels: Hypothesis testing functions (exported)
# =============================================================================

## LR TEST ---------------------------------------------------------------------

#' Likelihood Ratio Test
#'
#' Performs a likelihood ratio test for comparing two nested models fitted
#' with the same estimator (e.g., `ml_lm`).
#'
#' @param object_1 A fitted model object inheriting from `"mlmodel"`.
#'   This is typically the restricted (smaller) model.
#' @param object_2 A fitted model object inheriting from `"mlmodel"`.
#'   This is typically the unrestricted (larger) model.
#'   Note: The order of `object_1` and `object_2` does not matter —
#'   the function automatically determines which is the restricted model.
#'
#' @param ... Further arguments passed to methods (currently not used).
#'
#' @details
#' The likelihood ratio test statistic is calculated as:
#' \deqn{LR = 2 \times | \log L_{full} - \log L_{restricted} |}
#' which follows a \eqn{\chi^2} distribution with degrees of freedom equal to
#' the difference in the number of parameters between the two models.
#'
#' **Important**: The restricted model must have a lower (or equal) log-likelihood
#' than the unrestricted model. If this condition is violated, the test is invalid
#' and an error will be thrown.
#'
#' @return An object of class `"lrtest.mlmodel"`.
#'
#' @seealso [waldtest()], [IMtest()]
#'
#' @examples
#' \dontrun{
#' fit_small <- ml_lm(mpg ~ wt, data = mtcars)
#' fit_large <- ml_lm(mpg ~ wt + hp + qsec, data = mtcars)
#' lrtest(fit_small, fit_large)
#' }
#'
#' @export
lrtest <- function(object_1, object_2, ...) {
  UseMethod("lrtest")
}

#' @rdname lrtest
#' @export
lrtest.mlmodel <- function(object_1, object_2, ...)
{
  if (!inherits(object_1, "mlmodel") || !inherits(object_2, "mlmodel"))
    cli::cli_abort("Both `object_1` and `object_2` must inherit from class 'mlmodel'.",
                   call = NULL)

  if (class(object_1)[1] != class(object_2)[1])
    cli::cli_abort("`object_1` and `object_2` must be fitted with the same estimator type.",
                   call = NULL)

  k1 <- length(coef(object_1))
  k2 <- length(coef(object_2))

  if (k1 == k2) {
    cli::cli_abort(c("The two models have the same number of parameters ({.val {k1}}).",
                     "A likelihood ratio test is not meaningful in this case."),
                   call = NULL)
  }

  ll1 <- as.numeric(logLik(object_1))
  ll2 <- as.numeric(logLik(object_2))

  # Determine which is the restricted (smaller) model
  if (k1 < k2) {
    restricted <- object_1
    full       <- object_2
    ll_r       <- ll1
    ll_f       <- ll2
    cli::cli_alert_info("`object_1` is the restricted model (nested in `object_2`).")
  } else {
    restricted <- object_2
    full       <- object_1
    ll_r       <- ll2
    ll_f       <- ll1
    cli::cli_alert_info("`object_2` is the restricted model (nested in `object_1`).")
  }

  # Validity check
  if (ll_r > ll_f + 1e-8) {
    cli::cli_abort(c("Invalid LR test: The restricted model has a higher log-likelihood than the full model.",
                      "i" = "This usually indicates the models are not properly nested or there was an optimization issue."),
                    call = NULL)
  }

  lrstat <- 2 * (ll_f - ll_r)
  df     <- abs(k2 - k1)

  res <- list(
    chisq    = lrstat,
    df       = df,
    pval     = pchisq(lrstat, df, lower.tail = FALSE),
    logLik_r = ll_r,
    logLik_f = ll_f,
    npar_r   = k1,
    npar_f   = k2
  )

  class(res) <- "lrtest.mlmodel"
  res
}

#' @export
print.lrtest.mlmodel <- function(x, digits = 3, ...)
{
  if (!inherits(x, "lrtest.mlmodel"))
    cli::cli_abort("`x` must be an object of class 'lrtest.mlmodel'.",
                   call = NULL)

  cat("Likelihood Ratio Test\n")
  cat("--------------------------------------------\n")

  p_str <- if (x$pval < 1e-8) "< 1e-8" else sprintf("%.4f", x$pval)

  cat(sprintf("Chisq(%d) = %.3f    Pr(>Chisq) = %s\n",
              x$df, x$chisq, p_str))

  cat("--------------------------------------------\n")
  cat(sprintf("LogLik (restricted)   : %.3f   (df = %d)\n", x$logLik_r, x$npar_r))
  cat(sprintf("LogLik (unrestricted) : %.3f   (df = %d)\n", x$logLik_f, x$npar_f))

  invisible(x)
}



## IM TEST ---------------------------------------------------------------------
#' Information Matrix Test
#'
#' Performs the Information Matrix (IM) test for misspecification on models
#' fitted with `ml_lm()` and other `ml_*` models.
#'
#' @param object A fitted model object inheriting from `"mlmodel"`.
#' @param method Character. One of `"opg"` (default), `"quad"`, `"boot_opg"`,
#'   or `"boot_quad"`.
#' @param repetitions Integer. Number of bootstrap replications
#'   (only used for bootstrap methods). Default is 999.
#' @param seed Integer. Random seed for bootstrap.
#' @param ... Further arguments passed to methods.
#'
#' @export
IMtest <- function(object, ...) {
  UseMethod("IMtest")
}


#' @rdname IMtest
#' @export
IMtest.mlmodel <- function(object,
                           method = "opg",
                           repetitions = 999,
                           seed = 1234L,
                           ...)
{
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must inherit from class 'mlmodel'.", call = NULL)

  method <- rlang::arg_match(method, c("opg", "quad", "boot_opg", "boot_quad"))

  if (is.null(object$gradientObs))
    cli::cli_abort("Model object does not contain `gradientObs`.", call = NULL)

  if (is.null(object$model$functions$hessianObs))
    cli::cli_abort("No `hessianObs` function found in model$functions.", call = NULL)

  S <- object$gradientObs
  n <- nrow(S)
  k <- ncol(S)
  m <- k * (k + 1) / 2

  # Compute per-observation Hessians
  H_per_obs <- object$model$functions$hessianObs(object)

  # Build IM indicators
  ID <- matrix(0, nrow = n, ncol = m)
  for (i in seq_len(n)) {
    start <- (i-1)*k + 1
    end   <- i*k
    si <- S[i, , drop = FALSE]
    Hi <- H_per_obs[start:end, , drop = FALSE]
    ID[i, ] <- matrixcalc::vech(Hi + crossprod(si))
  }

  # OPG-based methods (Chesher/Lancaster)
  if (method %in% c("opg", "boot_opg")) {
    X <- cbind(S, ID)
    y <- rep(1, n)
    reg <- lm(y ~ X - 1, singular.ok = TRUE)
    tstat <- sum(reg$fitted.values^2)

    res <- list(
      tstat = tstat,
      df = m,
      pval = list(analytical = pchisq(tstat, m, lower.tail = FALSE)),   # analytical is primary
      version = list(description = "Chesher/Lancaster OPG", method = method)
    )

    # Bootstrap version
    if (method == "boot_opg") {
      set.seed(seed)
      boot_stats <- numeric(repetitions)
      n_success  <- 0

      w <- object$model$weights

      if (!is.null(object$model$data) && is.data.frame(object$model$data)) {
        orig_data <- object$model$data
      } else if (!is.null(object$model$d_name) && object$model$d_name != "<unknown data>") {
        orig_data <- tryCatch(get(object$model$d_name), error = function(e) {
          cli::cli_abort("Cannot retrieve the dataset to get the clustering variable.",
                         call = NULL)
        })
      } else {
        cli::cli_abort("Dataset and its name not stored; cannot retrieve clustering variable.",
                       call = NULL)
      }

      for (r in seq_len(repetitions)) {
        idx <- sample.int(n, n, replace = TRUE)

        boot_obj <- tryCatch({
          object$model$functions$update(
            object,
            data   = orig_data[idx, , drop = FALSE],
            weights = w[idx]
          )
        }, error = function(e) NULL)

        if (is.null(boot_obj) || !(boot_obj$code %in% c(1L, 2L, 8L))) {
          boot_stats[r] <- NA
          next
        }

        S_r <- boot_obj$gradientObs
        H_r <- boot_obj$model$functions$hessianObs(boot_obj)

        ID_r <- matrix(0, nrow = n, ncol = m)
        for (i in seq_len(n)) {
          start <- (i-1)*k + 1
          end   <- i*k
          si <- S_r[i, , drop = FALSE]
          Hi <- H_r[start:end, , drop = FALSE]
          ID_r[i, ] <- matrixcalc::vech(Hi + crossprod(si))
        }

        Xr <- cbind(S_r, ID_r)
        reg_r <- lm(y ~ Xr - 1, singular.ok = TRUE)
        boot_stats[r] <- sum(reg_r$fitted.values^2)
        n_success <- n_success + 1
      }

      res$pval$bootstrapped <- mean(boot_stats >= tstat, na.rm = TRUE)
      res$repetitions <- list(total = repetitions, valid = n_success)
      res$version$description <- "Chesher/Lancaster OPG + model-based bootstrap"
    }
  }
  else
  {
    # Quadratic methods.
    M <- ID
    G <- colSums(M)

    XS <- crossprod(S)
    if (rcond(XS) < 1e-12)
      cli::cli_alert_warning("Score matrix is nearly singular; IM test may be unreliable.")

    proj_coeff <- S %*% solve(XS, crossprod(S, M))
    R_mat      <- M - proj_coeff
    W          <- crossprod(R_mat)

    if (method == "quad") {
      tstat <- as.numeric(t(G) %*% MASS::ginv(W) %*% G)
      res <- list(
        tstat   = tstat,
        df      = m,
        pval    = pchisq(tstat, df = m, lower.tail = FALSE),
        version = list(description = "Orthogonalized Quadratic Form", method = method)
      )
    } else { # boot_quad
      set.seed(seed)
      boot_D <- matrix(0, nrow = repetitions, ncol = m)

      for (r in seq_len(repetitions)) {
        idx   <- sample.int(n, n, replace = TRUE)
        M_r   <- M[idx, , drop = FALSE]
        S_r   <- S[idx, , drop = FALSE]
        proj_r <- S_r %*% solve(crossprod(S_r), crossprod(S_r, M_r))
        boot_D[r, ] <- colMeans(M_r - proj_r)
      }

      V_boot <- cov(boot_D)
      tstat  <- as.numeric(t(G/n) %*% MASS::ginv(V_boot) %*% (G/n))

      res <- list(
        tstat   = tstat,
        df      = m,
        pval    = pchisq(tstat, df = m, lower.tail = FALSE),
        version = list(description = "Orthogonalized Quadratic Form with bootstrapped variance",
                       method = method),
        repetitions = repetitions
      )
    }
  }

  class(res) <- "IMtest"
  res
}


#' @export
print.IMtest <- function(x, digits = 3, ...)
{
  if (!inherits(x, "IMtest"))
    cli::cli_abort("`x` must be an object of class 'IMtest'.")

  cat("Information Matrix Test\n")
  cat(" Method:", x$version$description, "\n")
  cat("--------------------------------------------\n")

  if (x$version$method %in% c("boot_opg", "boot_quad")) {
    cat(" Repetitions:", x$repetitions$total, "\n")
  }

  cat(sprintf(" Chisq(%i) = %.3f", x$df, x$tstat))

  if (x$version$method == "boot_opg") {
    cat("\n P(>Chisq): Analytical   =", sprintf("%.4f", x$pval$analytical),
        "\n            Bootstrapped =", sprintf("%.4f", x$pval$bootstrapped))
  } else {
    cat("    Pr(>Chisq) =", sprintf("%.4f", x$pval))
  }

  cat("\n--------------------------------------------\n")
  invisible(x)
}

## WALDTEST --------------------------------------------------------------------
#' Wald Test for mlmodel objects
#'
#' Performs a Wald test of linear restrictions on the parameters.
#'
#' @param object An object of class `"mlmodel"` or any model that inherits
#'   from it (e.g. `"ml_lm"`).
#' @param indices Integer vector. Positions of the coefficients to test.
#' @param coef_names Character vector. Names of the coefficients to test.
#' @param rest_matrix Numeric matrix. A q × k restriction matrix (advanced).
#' @param rhs Numeric. Value(s) the linear combination(s) should equal.
#'   Default is 0.
#' @param vcov Optional user-supplied variance-covariance matrix.
#'   If provided, it will be used instead of computing one internally.
#'   Must be a square matrix with dimensions matching the number of
#'   coefficients in the model. Useful when you have pre-computed a
#'   bootstrap or clustered variance and want to reuse it.
#' @param vcov.type Character string specifying the type of variance-covariance
#'   matrix. One of `"oim"` (default), `"robust"`, `"opg"`, `"cluster"`,
#'   or `"boot"`. See [vcov.mlmodel()] for details.
#' @param cl_var Character string or vector. Clustering variable.
#'   Only used when `vcov.type = "cluster"` or `vcov.type = "boot"`.
#' @param repetitions Integer. Number of bootstrap replications when
#'   `vcov.type = "boot"`. Default is 999.
#' @param seed Integer. Random seed for reproducibility when `vcov.type = "boot"`.
#' @param progress Logical. Should a progress bar be displayed during
#'   bootstrapping? Default is `FALSE`.
#' @param ... Further arguments passed to methods.
#'
#' @export
waldtest <- function(object,
                     ...)
{
  UseMethod("waldtest")
}


#' @rdname waldtest
#' @export
waldtest.mlmodel <- function(object,
                             indices = NULL,
                             coef_names = NULL,
                             rest_matrix = NULL,
                             rhs = 0,
                             vcov = NULL,        # User-supplied variance matrix
                             vcov.type = "oim",
                             cl_var = NULL,
                             repetitions = 999,
                             seed = NULL,
                             progress = FALSE,
                             ...)
{
  if (!inherits(object, "mlmodel"))
    cli::cli_abort("`object` must be a model of class 'mlmodel'.",
                   call = NULL)

  # Check that exactly one way of specifying restrictions is used
  inputs <- sum(!is.null(indices), !is.null(coef_names), !is.null(rest_matrix))
  if (inputs == 0)
    cli::cli_abort("You must provide one of: `indices`, `coef_names`, or `rest_matrix`.",
                   call = NULL)
  if (inputs > 1)
    cli::cli_abort("You can only specify one of `indices`, `coef_names`, or `rest_matrix`.",
                   call = NULL)

  b <- coef(object)
  k <- length(b)

  # Get variance-covariance matrix using our existing vcov method
  V <- get_vcov(object,
                vcov = vcov,
                vcov.type   = vcov.type,
                repetitions = repetitions,
                seed        = seed,
                cl_var      = cl_var,
                progress    = progress,
              )

  # Build restriction matrix R
  if (!is.null(indices))
  {
    if (any(indices < 1 | indices > k))
      cli::cli_abort("`indices` must be between 1 and {k}.",
                     call = NULL)
    R <- diag(k)[indices, , drop = FALSE]

  }
  else if (!is.null(coef_names))
  {
    if (!all(coef_names %in% names(b)))
      cli::cli_abort("Some coefficient names were not found.",
                     call = NULL)
    idx <- match(coef_names, names(b))
    R <- diag(k)[idx, , drop = FALSE]

  }
  else
  {  # rest_matrix
    if (ncol(rest_matrix) != k)
      cli::cli_abort("`rest_matrix` must have {k} columns (one per parameter).",
                     call = NULL)
    if (nrow(rest_matrix) >= k)
      cli::cli_abort("Number of restrictions cannot be greater than or equal to number of parameters.")
    R <- rest_matrix
  }

  # Handle rhs
  q <- nrow(R)
  if (length(rhs) == 1) rhs <- rep(rhs, q)
  if (length(rhs) != q)
    cli::cli_abort("`rhs` must be a single number or a vector of length equal to the number of restrictions.",
                   call = NULL)

  # Build readable restriction strings
  restrictions <- character(q)
  b_names <- names(b)

  for (i in seq_len(q))
  {
    row <- R[i, ]
    non_zero <- which(abs(row) > 1e-10)

    if (length(non_zero) == 0)
    {
      restrictions[i] <- paste("0 =", rhs[i])
      next
    }

    terms <- character(length(non_zero))
    for (m in seq_along(non_zero))
    {
      j <- non_zero[m]
      coef_val <- row[j]
      name <- b_names[j]
      if (abs(coef_val) == 1)
        terms[m] <- ifelse(coef_val > 0, name, paste0("-", name))
      else
        terms[m] <- paste0(round(coef_val, 4), " * ", name)
    }
    lhs <- paste(terms, collapse = " + ")
    lhs <- gsub(" \\+ -", " - ", lhs)
    restrictions[i] <- paste(lhs, "=", rhs[i])
  }

  # Compute Wald statistic
  # === Compute Wald statistic with proper error handling ===
  wald_result <- tryCatch(
    {
      diff <- R %*% b - rhs
      restricted_var <- R %*% V %*% t(R)
      W <- as.vector(t(diff) %*% chol2inv(chol(restricted_var)) %*% diff)
      list(success = TRUE, W = W)
    },
    error = function(e) {
      list(success = FALSE,
           message = conditionMessage(e))
    }
  )

  if (!wald_result$success) {
    cli::cli_alert_warning(
      "Wald test could not be computed: the restricted variance matrix is singular."
    )
    cli::cli_alert_info(
      "This often happens with cluster-robust standard errors when there are few clusters \\
       or low within-cluster variation. Consider using {.code vcov.type = 'robust'} instead."
    )

    res <- list(
      waldstat     = NA_real_,
      df           = q,
      pval         = NA_real_,
      restrictions = restrictions,
      vcov.type    = vcov.type,
      singular     = TRUE,
      error_msg    = wald_result$message
    )
    class(res) <- "waldtest.mlmodel"
    return(res)
  }

  # Normal successful case
  res <- list(
    waldstat     = wald_result$W,
    df           = q,
    pval         = pchisq(wald_result$W, q, lower.tail = FALSE),
    restrictions = restrictions,
    vcov.type    = vcov.type,
    singular     = FALSE
  )

  # Store cluster information (if applicable)
  if (vcov.type == "cluster" && !is.null(cl_var))
    res$vcov.cluster <- vcov_cluster_info(object, cl_var)
  else
    res$vcov.cluster <- NULL

  class(res) <- "waldtest.mlmodel"

  return(res)
}

## PRINT WALDTEST --------------------------------------------------------------
#' @export
print.waldtest.mlmodel <- function(x, digits = 3, ...)
{
  if (!inherits(x, "waldtest.mlmodel"))
    cli::cli_abort("`x` must be an object of class 'waldtest.mlmodel'.")

  cat("\nWald Test of Linear Restrictions\n")

  # --- Smart variance type printing ---
  cat(" ",
    switch (x$vcov.type,
            "oim" = "Variance type: OIM",
            "robust" = "Variance type: Robust (sandwich)",
            "opg" = "Variance type: Outer Product of Gradients (OPG)",
            "cluster" = paste("Variance type: Cluster-robust   |  Clusters:",
                              x$vcov.cluster.n_cluster)    )
  )
  cat("\n")

  cat("--------------------------------------------\n")

  cat("Restrictions:\n")
  for (i in seq_along(x$restrictions)) {
    cat(sprintf("  %d: %s\n", i, x$restrictions[i]))
  }

  if (isTRUE(x$singular)) {
    cat("Wald statistic could not be computed (singular matrix).\n")
    cat("See warning message for details.\n")
  } else {
    p_str <- if (x$pval < 1e-8) "< 1e-8" else format.pval(x$pval, digits = 4)

    cat(sprintf("Chisq(%d) = %.3f    Pr(>Chisq) = %s\n",
                x$df, x$waldstat, p_str))
  }

  cat("--------------------------------------------\n\n")

  invisible(x)
}
