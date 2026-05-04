# =============================================================================
# mlmodels: Hypothesis testing functions (exported)
# =============================================================================

## LR TEST ---------------------------------------------------------------------
#' Likelihood Ratio Test for Nested mlmodel Objects
#'
#' Performs a likelihood ratio test comparing two nested models fitted with
#' the same estimator (e.g. `ml_lm`, `ml_logit`, `ml_negbin`, etc.).
#'
#' @param object_1 A fitted model object inheriting from `"mlmodel"`.
#'   Typically the restricted (smaller) model.
#' @param object_2 A fitted model object inheriting from `"mlmodel"`.
#'   Typically the unrestricted (larger) model.
#'   The order of `object_1` and `object_2` does not matter — the function
#'   automatically determines which is the restricted model.
#' @param ... Further arguments passed to methods (currently not used).
#'
#' @details
#' The likelihood ratio test statistic is calculated as:
#' \deqn{LR = 2 \times (\log L_{\text{unrestricted}} - \log L_{\text{restricted}})}
#' 
#' Under the null hypothesis that the restricted model is correct, `LR` follows
#' a \eqn{\chi^2} distribution with degrees of freedom equal to the difference
#' in the number of parameters between the two models.
#'
#' **Important:** The two models must be nested (the restricted model must be
#' a special case of the unrestricted one) and fitted on exactly the same sample.
#' The restricted model must have a lower (or equal) log-likelihood.
#'
#' @return An object of class `"lrtest.mlmodel"` with the test statistic,
#'   degrees of freedom, and p-value.
#'
#' @seealso [waldtest()], [IMtest()], [vuongtest()]
#'
#' @examples
#' 
#' # Linear model example
#' data(mroz)
#' mroz$incthou <- mroz$faminc / 1000
#' 
#' fit_small <- ml_lm(incthou ~ age + huswage, data = mroz)
#' fit_large <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
#'                    data = mroz)
#' 
#' lrtest(fit_small, fit_large)
#' 
#' # You can also reverse the order — the function detects the restricted model
#' lrtest(fit_large, fit_small)
#' 
#' @author Alfonso Sanchez-Penalver
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
  
  # -- 1. Check that models were estimated on compatible samples/weights -------
  .compare_estimation_samples(object_1, object_2)
  
  ll1 <- as.numeric(logLik(object_1, scaled = TRUE))
  ll2 <- as.numeric(logLik(object_2, scaled = TRUE))
  
  ll1_us <- as.numeric(logLik(object_1))
  
  # If they are NOT the same (difference > epsilon), it's a weighted model
  if(abs(ll1 - ll1_us) > 1e-8) {
    cli::cli_warn(c(
      "!" = "Estimations were done using weights. Likelihood-ratio test may be inappropriate.",
      "i" = "Using log-likelihoods scaled to the sample size (N = {.val {nobs(object_1)}}).",
      "i" = "We strongly recommend using `waldtest()` instead with robust standard errors."
    ))
  }

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
#' Information Matrix Test for Model Misspecification
#'
#' Performs the Information Matrix (IM) test for misspecification on models
#' fitted with the `mlmodels` package.
#'
#' @param object A fitted model object inheriting from `"mlmodel"`.
#' @param method Character string. Specifies the version of the test:
#'   * `"quad"` (default): Quadratic form of the Information Matrix test (most common).
#'   * `"opg"`: Outer Product of Gradients version (Chesher-Lancaster).
#'   * `"boot_quad"`: Analytical chi-square and p-value, plus bootstrap p-value for the quadratic form.
#'   * `"boot_opg"`: Analytical chi-square and p-value, plus bootstrap p-value for the OPG version.
#' @param repetitions Integer. Number of bootstrap replications when using a 
#'   bootstrap method. Default is 999.
#' @param seed Integer. Random seed for reproducibility in bootstrap methods.
#'   If `NULL`, a random seed is generated.
#' @param ... Further arguments passed to methods (currently not used).
#'
#' @details
#' The Information Matrix test checks whether the model is correctly specified
#' by testing the equality between the Hessian and the outer product of the 
#' gradient (information matrix equality). Rejection of the null hypothesis
#' indicates model misspecification (e.g., incorrect functional form, 
#' heteroskedasticity not properly modeled, omitted variables, etc.).
#'
#' Two main versions are implemented:
#' - **Quadratic form** (`"quad"`): Generally preferred for its better finite-sample properties.
#' - **OPG version** (`"opg"`): Chesher and Lancaster (1983) version.
#'
#' Bootstrap versions (`"boot_quad"` and `"boot_opg"`) provide p-values based on 
#' the empirical distribution of the test statistic and are useful when asymptotic 
#' approximations may be unreliable.
#'
#' @return An object of class `"IMtest.mlmodel"` containing the analytical test
#'   statistic, degrees of freedom and p-value, plus the bootstrapped p-value 
#'   (if a bootstrap method was selected).
#'
#' @seealso [waldtest()], [lrtest()], [vuongtest()]
#'
#' @examples
#' 
#' # Linear model example
#' data(mroz)
#' mroz$incthou <- mroz$faminc / 1000
#' 
#' fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
#'              data = mroz)
#' 
#' # Default quadratic form test
#' IMtest(fit)
#' 
#' # OPG version
#' IMtest(fit, method = "opg")
#' 
#' # Bootstrap p-value (quadratic form)
#' IMtest(fit, method = "boot_quad", repetitions = 100, seed = 123)
#' 
#' # Heteroskedastic model
#' fit_het <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem,
#'                  scale = ~ educ, data = mroz)
#' IMtest(fit_het)
#' 
#' @references
#' Chesher, A. (1983). The information matrix test: Simplified calculation 
#' via a score test interpretation. Economics Letters, 13(1), 45-48.
#' 
#' Lancaster, T. (1984). The covariance matrix of the information matrix test. 
#' Econometrica, 52(4), 1051-1053.
#' 
#' White, H. (1982). Maximum likelihood estimation of misspecified models. 
#' Econometrica, 50(1), 1-25.
#' 
#' @author Alfonso Sanchez-Penalver
#'
#' @export
IMtest <- function(object, ...) {
  UseMethod("IMtest")
}


#' @rdname IMtest
#' @export
IMtest.mlmodel <- function(object,
                           method = "quad",
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

  # Get weights for scaling inside the loop
  weights <- object$model$weights %||% rep(1, nobs(object))
  sc_factor <- length(weights) / sum(weights, na.rm = TRUE)
  
  S <- object$gradientObs * sc_factor # Overall scaling for later use
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
    Hi <- H_per_obs[start:end, , drop = FALSE] * sc_factor
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

      if (!is.null(object$model$data) && is.data.frame(object$model$data)) {
        orig_data <- object$model$data[object$model$sample, , drop = FALSE]
      } else if (!is.null(object$model$d_name) && object$model$d_name != "<unknown data>") {
        orig_data <- tryCatch(get(object$model$d_name), error = function(e) {
          cli::cli_abort("Cannot retrieve the dataset to get the clustering variable.",
                         call = NULL)
        })
        orig_data <- orig_data[object$model$sample, , drop = FALSE]
      } else {
        cli::cli_abort("Dataset and its name not stored; cannot retrieve clustering variable.",
                       call = NULL)
      }

      for (r in seq_len(repetitions)) {
        idx <- sample.int(n, n, replace = TRUE)

        suppressMessages(
          boot_obj <- tryCatch({
            update(
              object,
              data   = orig_data[idx, , drop = FALSE],
              weights = weights[idx]
            )
          }, error = function(e) NULL)
        )

        if (is.null(boot_obj) || !(boot_obj$code %in% c(0L, 1L, 2L, 8L))) {
          boot_stats[r] <- NA_real_
          next
        }
        
        w_r <- boot_obj$model$weights
        sc_r <- length(w_r) / sum(w_r, na.rm = TRUE)
        
        S_r <- boot_obj$gradientObs * sc_r
        H_r <- boot_obj$model$functions$hessianObs(boot_obj)

        ID_r <- matrix(0, nrow = n, ncol = m)
        for (i in seq_len(n)) {
          start <- (i-1)*k + 1
          end   <- i*k
          si <- S_r[i, , drop = FALSE]
          Hi <- H_r[start:end, , drop = FALSE] * sc_r
          ID_r[i, ] <- matrixcalc::vech(Hi + crossprod(si))
        }

        Xr <- cbind(S_r, ID_r)
        
        reg_r <- lm(y ~ Xr - 1, singular.ok = TRUE)
        boot_stats[r] <- sum(reg_r$fitted.values^2)
        n_success <- n_success + 1
      }

      res$pval$bootstrapped <- mean(boot_stats >= tstat, na.rm = TRUE)
      res$repetitions <- list(total = repetitions, valid = n_success)
      res$version$description <- "Chesher/Lancaster OPG + Model-based bootstrap"
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
    W          <- crossprod(R_mat) * (n / (n - k))  # <- small sample adjustment.
    
    tstat <- as.numeric(t(G) %*% MASS::ginv(W) %*% G)
    res <- list(
      tstat   = tstat,
      df      = m,
      pval = list(analytical = pchisq(tstat, m, lower.tail = FALSE)),   # analytical is primary
      version = list(description = "Orthogonalized Quadratic Form", method = method)
    )

    if (method == "boot_quad") {
      set.seed(seed)
      boot_stats <- numeric(repetitions)
      n_success  <- 0
      
      if (!is.null(object$model$data) && is.data.frame(object$model$data)) {
        orig_data <- object$model$data[object$model$sample, , drop = FALSE]
      } else if (!is.null(object$model$d_name) && object$model$d_name != "<unknown data>") {
        orig_data <- tryCatch(get(object$model$d_name), error = function(e) {
          cli::cli_abort("Cannot retrieve the dataset to get the clustering variable.",
                         call = NULL)
        })
        orig_data <- orig_data[object$model$sample, , drop = FALSE]
      } else {
        cli::cli_abort("Dataset and its name not stored; cannot retrieve clustering variable.",
                       call = NULL)
      }
      
      # Loop, with re-estimation and storing the results.
      for (r in seq_len(repetitions)) {
        idx <- sample.int(n, n, replace = TRUE)
        
        suppressMessages(
          boot_obj <- tryCatch({
            update(
              object,
              data   = orig_data[idx, , drop = FALSE],
              weights = weights[idx]
            )
          }, error = function(e) NULL)
        )
        
        if (is.null(boot_obj) || !(boot_obj$code %in% c(0L, 1L, 2L, 8L))) {
          boot_stats[r] <- NA_real_
          next
        }
        
        w_r <- boot_obj$model$weights
        sc_r <- length(w_r) / sum(w_r, na.rm = TRUE)
        
        S_r <- boot_obj$gradientObs * sc_r
        H_r <- boot_obj$model$functions$hessianObs(boot_obj)
        
        ID_r <- matrix(0, nrow = n, ncol = m)
        for (i in seq_len(n)) {
          start <- (i-1)*k + 1
          end   <- i*k
          si <- S_r[i, , drop = FALSE]
          Hi <- H_r[start:end, , drop = FALSE] * sc_r
          ID_r[i, ] <- matrixcalc::vech(Hi + crossprod(si))
        }
        
        G_r <- colSums(ID_r)
        
        XS_r <- crossprod(S_r)
        
        proj_coeff_r <- tryCatch({
          S_r %*% solve(XS_r, crossprod(S_r, ID_r))
        }, error = function(e) NULL)
        
        if(is.null(proj_coeff_r))
        {
          boot_stats[r] = NA_real_
          next
        }
        
        R_mat_r      <- ID_r - proj_coeff_r
        W_r          <- crossprod(R_mat_r) * (n / (n - k))  # <- small sample adjustment.
        
        boot_stats[r] <- as.numeric(t(G_r) %*% MASS::ginv(W_r) %*% G_r)
        n_success <- n_success + 1
      }
      
      res$pval$bootstrapped <- mean(boot_stats >= tstat, na.rm = TRUE)
      res$repetitions <- list(total = repetitions, valid = n_success)
      res$version$description <- "Orthogonalized Quadratic Form + Model-based bootstrap"
    }
  }
  
  res$model <- object$model$description
  
  class(res) <- "IMtest.mlmodel"
  res
}

#' @export
print.IMtest.mlmodel <- function(x, digits = 4, ...)
{
  if (!inherits(x, "IMtest.mlmodel"))
    cli::cli_abort("`x` must be an object of class 'IMtest.mlmodel'.")

  cat("Information Matrix Test\n")
  cat(" Method:", x$version$description, "\n")
  cat(" Model: ", x$model, "\n")
  cat("--------------------------------------------\n")

  if (x$version$method %in% c("boot_opg", "boot_quad")) {
    cat(" Repetitions: Total", x$repetitions$total, "Successful", x$repetitions$valid,  "\n")
  }

  cat(sprintf(" Chisq(%i) = %.3f", x$df, x$tstat))

  if (x$version$method %in% c("boot_opg", "boot_quad")) {
    cat("\n P(>Chisq): Analytical   =", sprintf("%.4f", x$pval$analytical),
        "\n            Bootstrapped =", sprintf("%.4f", x$pval$bootstrapped))
  } else {
    cat("    Pr(>Chisq) =", sprintf("%.4f", x$pval$analytical))
  }

  cat("\n--------------------------------------------\n")
  invisible(x)
}

## VUONG's TEST ================================================================
#' Vuong's Test for Non-Nested Models
#'
#' Performs Vuong's (1989) test for comparing two non-nested models fitted
#' via maximum likelihood with the `mlmodels` package.
#'
#' @param object_1 A fitted model object inheriting from `"mlmodel"`.
#' @param object_2 A fitted model object inheriting from `"mlmodel"`.
#' @param boot Should bootstrapped p-values be calculated?
#' @param repetitions Number of iterations for the bootstrap method. Defaults to
#'                    999. Only relevant if `boot = TRUE`.
#' @param seed Integer with a seed to use for the random sampling for the
#'             bootstrapping. Only relevant if `boot = TRUE`. If none supplied
#'             a random one will be generated.
#' @param ... Further arguments passed to methods (currently not used).
#'
#' @details
#' Vuong's test compares two non-nested models by testing the null hypothesis
#' that the two models are equally close to the true data generating process.
#' 
#' The test statistic is based on the difference in the per-observation 
#' log-likelihood contributions between the two models. A positive significant 
#' value favors `object_1`, a negative significant value favors `object_2`, 
#' and a non-significant value leads to an "inconclusive" result.
#'
#' Both models must be estimated on exactly the same sample.
#'
#' @return An object of class `"vuongtest.mlmodel"` containing the test 
#'   statistic, p-value, and a conclusion (which model is preferred or 
#'   "inconclusive").
#'
#' @references
#' Vuong, Q. H. (1989). 'Likelihood Ratio Tests for Model Selection and Non-Nested 
#' Hypotheses.' *Econometrica*, 57(2), 307-333. 
#' \doi{10.2307/1912557}
#'
#' @seealso [lrtest()], [waldtest()], [IMtest()]
#'
#' @examples
#' 
#' # Linear models example (lognormal vs gamma)
#' data(mroz)
#' mroz$incthou <- mroz$faminc / 1000
#' 
#' fit_lognormal <- ml_lm(log(incthou) ~ age + I(age^2) + huswage + educ + unem,
#'                        data = mroz)
#' 
#' fit_gamma <- ml_gamma(incthou ~ age + I(age^2) + huswage + educ + unem,
#'                       data = mroz)
#' 
#' 
#' vuongtest(fit_lognormal, fit_gamma)
#' 
#' # Count models example
#' 
#' fit_poi <- ml_poisson(docvis ~ private + medicaid + age + I(age^2) + educyr +
#'                           actlim + totchr,
#'                       data = docvis)
#' 
#' fit_nb1 <- ml_negbin(docvis ~ private + medicaid + age + I(age^2) + educyr +
#'                           actlim + totchr,
#'                       data = docvis,
#'                       dispersion = "NB1")
#' 
#' fit_nb2 <- ml_negbin(docvis ~ private + medicaid + age + I(age^2) + educyr +
#'                           actlim + totchr,
#'                       data = docvis)
#'                       
#' # Poisson vs. NB1
#' vuongtest(fit_poi, fit_nb1)
#' 
#' # NB1 vs. NB2 (bootstrapped, low repetitions for speed)
#' vuongtest(fit_nb1, fit_nb2, boot = TRUE, repetitions = 50)
#' 
#' # Binary models example
#' data(smoke)
#' smoke$smokes <- smoke$cigs > 0
#' 
#' fit_logit <- ml_logit(smokes ~ cigpric + income + age, data = smoke)
#' fit_probit <- ml_probit(smokes ~ cigpric + income + age, data = smoke)
#' 
#' vuongtest(fit_logit, fit_probit)
#'
#' @export
vuongtest <- function(object_1, object_2, ...) UseMethod("vuongtest")

#' @rdname vuongtest
#' @export
vuongtest.mlmodel <- function(object_1, object_2,
                              boot = FALSE,
                              repetitions = 999,
                              seed = NULL,
                              ...)
{
  if (!inherits(object_1, "mlmodel") || !inherits(object_2, "mlmodel"))
    cli::cli_abort("Both `object_1` and `object_2` must inherit from 'mlmodel'.")
  
  # -- 1. Check that models were estimated on compatible samples/weights -------
  .compare_estimation_samples(object_1, object_2)
  
  # -- 2. Analytical Vuong Test ------------------------------------------------
  # Scale weights so test statistic is valid
  w1 <- object_1$model$weights %||% rep(1, nobs(object_1))
  w2 <- object_2$model$weights %||% rep(1, nobs(object_2))
  
  w1_scaled <- w1 / sum(w1) * length(w1)
  w2_scaled <- w2 / sum(w2) * length(w2)
  
  # Temporary objects with scaled weights
  temp1 <- object_1
  temp2 <- object_2
  temp1$model$weights <- w1_scaled
  temp2$model$weights <- w2_scaled
  
  ll1 <- object_1$model$functions$loglikeObs(temp1)
  ll2 <- object_2$model$functions$loglikeObs(temp2)
  
  diff <- ll1 - ll2
  n <- length(diff)
  
  v_stat <- (sqrt(n) * mean(diff)) / sd(diff)
  p_val <- 2 * pnorm(abs(v_stat), lower.tail = FALSE)
  
  res <- list(
    teststat = v_stat,
    pval = p_val,
    boot = NULL,
    models = c(object_1$model$description, object_2$model$description)
  )
  
  # --- 3. Bootstrap Version (if requested) -------------------------------
  if (boot) {
    if (is.null(seed)) seed <- sample.int(1e6, 1)
    set.seed(seed)
    
    data_orig <- object_1$model$data
    sample_idx <- object_1$model$sample
    n_used <- sum(sample_idx)
    used_data <- data_orig[sample_idx, , drop = FALSE]
    w_used <- (object_1$model$weights %||% rep(1, nobs(object_1)))
    
    boot_stats <- numeric(repetitions)
    n_success <- 0
    
    for (r in seq_len(repetitions)) {
      idx <- sample.int(n_used, n_used, replace = TRUE)
      d_boot <- used_data[idx, , drop = FALSE]
      w_boot <- w_used[idx]
      
      # Scale weights for this bootstrap sample
      w_boot <- w_boot / sum(w_boot) * length(w_boot)
      
      suppressMessages({
        boot1 <- tryCatch(update(object_1, data = d_boot, weights = w_boot), error = function(e) NULL)
        boot2 <- tryCatch(update(object_2, data = d_boot, weights = w_boot), error = function(e) NULL)
      })
      
      if (is.null(boot1) || is.null(boot2) || 
          !boot1$code %in% c(0L, 1L, 2L, 8L) || 
          !boot2$code %in% c(0L, 1L, 2L, 8L)) {
        boot_stats[r] <- NA_real_
        next
      }
      
      # Scaled log-likelihoods for this bootstrap sample
      ll1b <- object_1$model$functions$loglikeObs(boot1)
      ll2b <- object_2$model$functions$loglikeObs(boot2)
      
      diffb <- ll1b - ll2b
      boot_stats[r] <- (sqrt(n_used) * mean(diffb)) / sd(diffb)
      n_success <- n_success + 1
    }
    
    boot_pval <- mean(boot_stats >= v_stat, na.rm = TRUE)
    
    res$boot <- list(
      pval = boot_pval,
      repetitions = list(total = repetitions, successful = n_success)
    )
  }
  
  class(res) <- "vuongtest.mlmodel"
  res
}

#' @export
print.vuongtest.mlmodel <- function(x, digits = 4, ...)
{
  if(!inherits(x, "vuongtest.mlmodel"))
    cli::cli_abort("`x` needs to be an object of class `vuongtest.mlmodel`")
  
  cat("\nVuong's (1989) Test\n")
  cat("--------------------------------------------------\n")
  cat("  Model 1:", x$models[1], "\n")
  cat("  Model 2:", x$models[2], "\n")
  cat("--------------------------------------------------\n")
  
  cat("Analyitical Results:",
      sprintf("  z-stat:  %.3f", x$teststat),
      sprintf("  p-value: %.4f", x$pval),
      sep = "\n")
  if(! is.null(x$boot))
  {
    cat(sprintf("Bootstrap (%d/%d repetitions):",
                x$boot$repetitions$successful,x$boot$repetitions$total),
        sprintf("  p-value: %.4f", x$boot$pval),
        sep = "\n")
  }
  
  cat("--------------------------------------------------\n")
  # Decision Logic
  if (!is.null(x$boot)) {
    anal_sig <- x$pval < 0.10
    boot_sig  <- x$boot$pval < 0.10
    
    if (anal_sig && boot_sig) {
      # Strong agreement
      winner <- if (x$teststat > 0) x$models[1] else x$models[2]
      cat(" ", winner, "is preferred.\n")
    } 
    else if (!anal_sig && !boot_sig) {
      # Strong agreement on no difference
      cat("  Inconclusive test: neither model is clearly preferred.\n")
    } 
    else if (anal_sig && !boot_sig) {
      # Most common interesting case
      winner <- if (x$teststat > 0) x$models[1] else x$models[2]
      cat(" ", winner, "appears better overall (analytical test),\n")
      cat("  but the bootstrap p-value indicates this preference is highly sensitive\n")
      cat("  to sample variation, and that choice is not robust to it.\n")
    } 
    else {
      # Rare case: bootstrap sees difference, analytical does not
      cat("  The bootstrap p-value contradicts the analytical test, suggesting that\n")
      cat("  the analytical version may be losing power (possibly due to high variance\n")
      cat("  or outliers in the likelihood ratios). The bootstrap result is likely\n")
      cat("  more robust in this sample.\n")
    }
  } else {
    # Only analytical test available
    if (x$pval < 0.05) {
      winner <- if (x$teststat > 0) x$models[1] else x$models[2]
      cat(" ", winner, "seems to be preferred.\n")
    } else if (x$pval < 0.10) {
      cat("  Weak evidence:", ifelse(x$teststat > 0, x$models[1], x$models[2]), 
          "seems to be preferred.\n")
    } else {
      cat("  Inconclusive test: neither model is clearly preferred.\n")
    }
  }
  
  invisible(x)
}

## WALDTEST --------------------------------------------------------------------
#' Wald Test for Linear Restrictions
#'
#' Performs a Wald test of linear restrictions on the parameters of an 
#' `mlmodel` object.
#'
#' @param object An object of class `"mlmodel"`.
#' @param constraints Specification of the linear constraints to test. 
#'   Can be one of the following:
#'   * An **integer vector** with the positions (indices) of the coefficients 
#'     to test (e.g. `2:5` or `c(1, 3, 7)`).
#'   * A **character vector** with coefficient names or linear combinations 
#'     on the left-hand side (e.g. `c("value::age", "value::educ + value::huswage")`).
#'   * A **numeric matrix** containing the full restriction matrix `R` 
#'     (advanced use).
#' @param rhs Numeric vector. Value(s) the linear combination(s) should equal.
#'   Default is 0.
#' @param vcov Optional user-supplied variance-covariance matrix.
#' @param vcov.type Character string. Type of variance-covariance matrix to use.
#'   One of `"oim"` (default), `"opg"`, `"robust"`, `"boot"`, or `"jack"`.
#'   See [vcov.mlmodel()] for details.
#' @param cl_var Character string or vector. Clustering variable when 
#'   `vcov.type = "robust"` or `"boot"`.
#' @param repetitions Integer. Number of bootstrap replications when 
#'   `vcov.type = "boot"`. Default is 999.
#' @param seed Integer. Random seed for bootstrap.
#' @param progress Logical. Show progress bar during bootstrapping? Default `FALSE`.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' The Wald test evaluates linear restrictions of the form \eqn{R\beta = r}.
#'
#' The `constraints` argument specifies the left-hand side of the linear
#' restrictions (without the equality sign). The right-hand side values are
#' controlled separately with the `rhs` argument.
#' 
#' `rhs` can be:
#' * A single number — all constraints are tested against this value (default is `0`).
#' * A vector of the same length as the number of constraints — each constraint
#'   is tested against its corresponding value.
#'
#' The test statistic follows a \eqn{\chi^2} distribution with degrees of 
#' freedom equal to the number of restrictions under the null hypothesis.
#'
#' Internally, the test always computes a chi-squared statistic. However, when
#' there is only **one restriction** (`df = 1`), the printed output shows the
#' equivalent **z-statistic** (`z = \sqrt{\text{Chisq(1)}}`) instead. 
#' 
#' This is because a \eqn{\chi^2(1)} random variable is the square of a standard
#' normal (\eqn{z}) random variable. Reporting the z-statistic in this case is
#' conventional, but both distributions are equivalent in that case.
#' 
#' @return An object of class `"waldtest.mlmodel"`.
#'
#' @seealso [lrtest()], [IMtest()], [confint.mlmodel()]
#'
#' @examples
#' 
#' data(mroz)
#' mroz$incthou <- mroz$faminc / 1000
#' 
#' fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
#'              data = mroz)
#' 
#' # 1. Joint sginificance using positions (OIM variance)
#' waldtest(fit, constraints = c(2, 5))
#' 
#' # 2. Different magnitudes for different coefficients using names (robust variance)
#' waldtest(fit, constraints = c("value::educ", "value::unem"),
#'          rhs = c(1,0), vcov.type = "robust")
#'          
#' # 3. Linear combination: educ + huswage = 3
#' waldtest(fit, constraints = "value::educ + value::huswage", rhs = 3,
#'          vcov.type = "robust")
#' 
#' # 4. Same test using the restriction matrix.
#' R <- matrix(c(0, 0, 0, 1, 1, 0, 0), nrow = 1)
#' waldtest(fit, constraints = R, rhs = 3, vcov.type = "robust")
#' 
#' @author Alfonso Sanchez-Penalver
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
                             # indices = NULL,
                             # coef_names = NULL,
                             # rest_matrix = NULL,
                             constraints = NULL,
                             rhs = 0,
                             vcov = NULL,
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
  if (is.null(constraints))
    cli::cli_abort("`constraints` cannot be null.",
                   call = NULL)
  
  if (!(is.matrix(constraints) || is.numeric(constraints) || is.character(constraints))) {
    cli::cli_abort(c(
      "Invalid {.arg constraints} argument.",
      "i" = "It must be one of:",
      "*" = "Character vector (coefficient names or linear combinations);",
      "*" = "Integer vector (indices of coefficients);",
      "*" = "Matrix (full restriction matrix).",
      " " = " ",
      "x" = "You supplied an object of class {.cls {class(constraints)[1]}}."
    ))
  }

  b <- coef(object)
  k <- length(b)

  # Get variance-covariance matrix using our existing vcov method
  V <- .process_vcov(object,
                     vcov = vcov,
                     vcov.type   = vcov.type,
                     cl_var      = cl_var,
                     repetitions = repetitions,
                     seed        = seed,
                     progress    = progress)

  # -- Check for unusable variance ---------------------------------------------
  if (any(!is.finite(V)) || any(is.na(V))) {
    cli::cli_abort(
      c("Cannot perform Wald test: variance matrix is unusable.",
        "i" = "This usually happens when using bootstrap variance with constraints.",
        "i" = "Consider using `vcov.type = 'robust'` or `vcov.type = 'oim'` instead."),
      call = NULL
    )
  }
  
  # Pass the matrix to the description helper to create the string
  var_description <- .vcov_description(V)
  
  # Check fractional response inference for oim or opg (the helper checks if
  # they're logit or probit, so no need to do it here)
  .fractional_response_inference_alert(object, V)
  
  # Build restriction matrix R
  if (is.matrix(constraints))
  {
    # Restriction matrix.
    if (ncol(constraints) != k) {
      cli::cli_abort(c(
        "Invalid restriction matrix.",
        "x" = "It has {.val {ncol(constraints)}} columns.",
        "i" = "It must have exactly {.val {k}} columns (one per parameter)."
      ))
    }
    
    if (nrow(constraints) > k - 1) {
      cli::cli_abort(c(
        "Too many restrictions.",
        "x" = "You provided {.val {nrow(constraints)}} restrictions.",
        "i" = "The maximum allowed is {.val {k-1}} (otherwise the restriction matrix is singular)."
      ))
    }
    
    R <- constraints
  }
  else if (is.character(constraints))
  {
    if (length(constraints) > k - 1) {
      cli::cli_abort(c(
        "Too many restrictions.",
        "x" = "You provided {.val {nrow(constraints)}} restrictions.",
        "i" = "The maximum allowed is {.val {k-1}} (otherwise the restriction matrix is singular)."
      ))
    }
    
    R <- .make_restriction_matrix(object, constraints)
  }
  else
  {
    # Indices
    if (!all(floor(constraints) == constraints)) {
      cli::cli_abort("Coefficient indices must be integers.", call = NULL)
    }
    
    if (any(constraints < 1 | constraints > k)) {
      cli::cli_abort(c(
        "Coefficient indices out of range.",
        "i" = "They must be between {.val 1} and {.val {k}}."
      ))
    }
    
    if (length(constraints) > k - 1) {
      cli::cli_abort(c(
        "Too many restrictions.",
        "x" = "You provided {.val {length(constraints)}} constraints.",
        "i" = "The maximum allowed is {.val {k-1}}."
      ))
    }
    
    R <- diag(k)[constraints, , drop = FALSE]
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
      var_description = var_description,
      # var_cluster = vcov_cluster,
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
    var_description = var_description,
    # var_cluster = vcov_cluster,
    singular     = FALSE
  )

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
  cat("\nVariance type:", x$var_description)
  cat("\n---------------------------------------\n")

  cat("Restrictions:\n")
  for (i in seq_along(x$restrictions)) {
    cat(sprintf("  %d: %s\n", i, x$restrictions[i]))
  }

  if (isTRUE(x$singular)) {
    cat("Wald statistic could not be computed (singular matrix).\n")
    cat("See warning message for details.\n")
  } else {
    
    chisq <- x$waldstat
    df    <- x$df
    pval  <- x$pval
    
    if(df == 1)
    {
      # Single restriction z-statistic
      zstat <- sqrt(chisq)
      cat(sprintf("z = %.3f    Pr(>|z|) = %s\n",
                  zstat, format.pval(pval, digits = 4, eps = 1e-8)))
    }
    else
      cat(sprintf("Chisq(%d) = %.3f    Pr(>Chisq) = %s\n",
                  df, chisq, format.pval(pval, digits = 4, eps = 1e-8)))
  }

  cat("--------------------------------------------\n\n")

  invisible(x)
}
