## GOFtest =====================================================================
#' Goodness of Fit Test for Count Models
#' 
#' Performs the Manjón y Martínez (2014) Score Test of Overdispersion.
#' 
#' @param object An object of class `mlmodel.count`, that has been estimated with
#'   one of our estimators for count data (ml_poisson, ml_negbin, ...)
#'
#' @param bins A vector of integers that sets the boundaries of the bins to
#'   analyze. Defaults to `0:5`.
#' 
#' @returns An object of class `GOFtest` with the following elements:
#' \describe{
#'    \item{model}{A string with the estimation model type for which the test
#'    is being performed.}
#'    \item{matrix}{A matrix with the different statistics to display for the
#'    different bins.}
#'    \item{test}{A list with the values of the actual test. The elements are:
#'        \describe{
#'            \item{teststat}{The value of the test statistic.}
#'            \item{df}{The number of degrees of freedom in the test.}
#'            \item{pval}{The p-value of the test.}
#'        }
#'    }
#' }
#' 
#' @details
#' While performing the test, it also prepares a table to display collecting
#' the frequency, proportion of cases, average predicted probability, absolute
#' difference between these two last values, and the pearson statistic for each
#' of the bins implied by the vector of numbers in `bins`.
#' 
#' @examples
#' set.seed(123)
#' d <- data.frame(x = rnorm(100))
#' d$y <- rpois(100, lambda = exp(d$x * 0.5))
#' fit <- ml_poisson(y ~ x, data = d)
#' GOFtest(fit, bins = 1:3)
#' 
#' 
#' @author Alfonso Sanchez-Penalver
#' 
#' @references 
#' Manjón, M., & Martínez, O. (2014). The chi-squared goodness-of-fit 
#' test for count-data models. \emph{The Stata Journal}, 14(4), 798-816. 
#' \doi{10.1177/1536867X1401400406}
#' 
#' 
#' @export
GOFtest <- function(object, bins = 0:5)
{
  if(!inherits(object, "mlmodel.count"))
    cli::cli_abort("`object` needs to be of class 'mlmodel.count'", call = NULL)
  
  # 1. Type validation
  if (!is.numeric(bins)) {
    cli::cli_abort("{.arg bins} must be a numeric vector.", call = NULL)
  }
  
  # 2. Structure cleaning
  bins <- sort(unique(floor(bins)))
  
  # 3. Domain validation
  if (any(bins < 0)) {
    cli::cli_abort("All values in {.arg bins} must be non-negative.", call = NULL)
  }
  
  # 4. Check the estimator for two things:
  #    - a. Set the name of the model.
  #    - b. Add zero to the bins if not there, for the appropriate estimators
  if(inherits(object, "ml_poisson"))
  {
    model <- "Poisson"
    if(!(0 %in% bins))
      bins <- c(0, bins)
  }
  else if (inherits(object, "ml_negbin"))
  {
    if(is.null(object$model$scale_formula))
      model <- paste0("Homoskedastic Negative Binomial (", object$model$dispersion, ")")
    else
      model <- paste0("Heteroskedastic Negative Binomial (", object$model$dispersion, ")")
    if(!(0 %in% bins))
      bins <- c(0, bins)
  }
  else
    model <- NULL  # For later models we will check here.
  
  y <- as.vector(object$model$value$outcomes[[1]])
  n_obs <- length(y)
  n_bins <- length(bins)
  
  # Regression matrix: rows: # observations, cols: # bins.
  R <- matrix(0, nrow = n_obs, ncol = n_bins)
  
  # Display matrix: rows # bins + 1, 
  D <- matrix(0, nrow = n_bins + 1, ncol = 5)
  
  row_names <- character(n_bins + 1)
  # need to loop through each value in bins.
  for(i in seq_len(length(bins)))
  {
    # now the bin we want will be formed by the current value and the next one
    # except for the last bin
    if (i == length(bins))
    {
      row_names[i] <- paste(bins[i], "-", bins[i])
      # We're on the last bin so it's probability equal to that value. 
      d_j <- (y == bins[i])
      p_j <- predict(object, type = paste0("P(", bins[i], ")"))
    }
    else
    {
      low <- bins[i]
      high <- bins[i+1] - 1
      row_names[i] <- paste(low, "-", high)
      if (low == high)
      {
        # exact probability
        d_j <- (y == low)
        p_j <- predict(object, type = paste0("P(", low, ")"))
      }
      else
      {
        # interval
        d_j <- (y >= low & y <= high)
        p_j <- predict(object, type = paste0("P(", low, ",", high, ")"))
      }
    }
    # Extract the actual prediction from the new class object
    p_j <- p_j$fit
    R[, i] <- d_j - p_j
    D[i, ] <- rbind(
      sum(d_j),                                # Frequency
      mean(d_j),                               # Proportion
      mean(p_j),                               # Probability
      abs(mean(d_j) - mean(p_j)),              # |Difference|
      n_obs * ((mean(d_j)) - mean(p_j))^2 / mean(p_j)  # Pearson Chisq.
    )
  }
  
  # We need to add the last one to the display, which is the probability of
  # greater than the last value.
  k <- bins[n_bins] + 1
  d_j <- y >= k
  p_j <- predict(object, type = paste0("P(", k, ",)"))$fit
  D[nrow(D), ] <- rbind(
    sum(d_j),                                # Frequency
    mean(d_j),                               # Proportion
    mean(p_j),                               # Probability
    abs(mean(d_j) - mean(p_j)),              # |Difference|
    n_obs * ((mean(d_j)) - mean(p_j))^2 / mean(p_j)  # Pearson Chisq.
  )
  row_names[nrow(D)] <- paste(k,"+")
  
  colnames(D) <- c("Frequency",
                   "Proportion",
                   "Probability",
                   "|Difference|",
                   "Pearson")
  rownames(D) <- row_names
  
  # Regression.
  ones <- rep(1, length(y))
  x <- cbind(R, object$gradientObs)
  ols <- .lm.fit(x, ones)
  teststat <- n_obs - sum(ols$residuals^2)
  df <- n_bins
  pval <- pchisq(teststat, df, lower.tail = FALSE)
  
  out <- list(
    model = model,
    matrix = D,
    test = list(
      teststat = teststat,
      df = df,
      pval = pval
    )
  )
  class(out) <- "GOFtest"
  return(out)
}


#' @export
print.GOFtest <- function(x, ...) {
  if(!inherits(x, "GOFtest"))
    cli::cli_abort("`x` must be of class `GOFtest`.",
                   call = NULL)
  
  cat("\nGoodness-of-fit test for count models\n")
  if(!is.null(x$model))
    cat("   Model:", x$model)
  cat("\n--------------------------------------------------\n")
  cat("Manjón & Martínez (2014) Score Test\n\n")
  
  # Print the Comparison Table
  print(round(x$matrix, 4))
  
  cat("\n--------------------------------------------------\n")
  cat(sprintf("  Chisq(%d):             %0.4f\n", x$test$df, x$test$teststat))
  cat(sprintf("  p-value:               %0.4f\n", x$test$pval))
  cat("--------------------------------------------------\n")
  
  invisible(x)
}

## OVDtest =====================================================================
#' Overdispersion tests for count outcomes.
#' 
#' Performs Cameron and Trivedi's test for overdispersion.
#' 
#' @param object An object of class `'mlmodel.count'`, i.e estimated by one of
#'   our count variable models (ml_poisson, ml_negbin,...)
#' 
#'
#' @returns A list with the results of overdispersion tests of a poisson model
#'   against both a NB1 and NB2 negative binomial models.
#'   
#' @references 
#' Cameron, A. C. and Trivedi, P. K. (1990). Regression-based tests for overdispersion in 
#' the Poisson model. *Journal of Econometrics*, 46(3), 347-364.
#' 
#' Cameron, A. C. and Trivedi, P. K. (2013). *Regression Analysis of Count Data*. 
#' 2nd Edition. Cambridge University Press.
#' 
#' @seealso \code{\link{IMtest}} for general misspecification testing.
#'
#' @export
OVDtest <- function(object)
{
  if(!inherits(object, "mlmodel.count"))
    cli::cli_abort("`object` needs to be of class 'mlmodel.count'", call = NULL)
  
  # No matter what model, we need to pull y from the value object
  y <- as.vector(object$model$value$outcomes[[1]])
  if(inherits(object, "ml_poisson"))
  {
    # Poisson already estimated fitted values and residuals.
    
    u <- object$model$residuals
    mu <- object$model$fitted.values
  }
  else
  {
    # NB1 or NB2. Need to estimate a poisson. Need x and w.
    x <- as.matrix(object$model$value$predictors)
    w <- object$model$weights
    
    pois_fit <- suppressWarnings({
      .ml_poisson.fit(y,x,w)
    })
    
    # this returns the maxLik object.
    beta <- pois_fit$estimate
    mu <- as.vector(exp(x %*% beta))
    u <- y - mu
  }
  
  index <- (u^2 - y) / mu
  df <- length(index) - 1
  
  nb2_reg <- .lm.fit(as.matrix(mu), index)
  nb1_reg <- .lm.fit(matrix(1, nrow = length(index)), index)
  
  # okay have to do the t-tests myself. 
  nb2_alpha <- nb2_reg$coefficients
  nb2_resid <- nb2_reg$residuals
  nb2_mse <- sum(nb2_resid^2) / df
  nb2_se  <- sqrt(nb2_mse / sum(mu^2))
  nb2_t <- nb2_alpha / nb2_se
  nb2_pval <- pt(nb2_t, df = df, lower.tail = FALSE)
  
  nb1_alpha <- nb1_reg$coefficients
  nb1_resid <- nb1_reg$residuals
  nb1_mse <- sum(nb1_resid^2) / df
  nb1_se  <- sqrt(nb1_mse / length(index))
  nb1_t <- nb1_alpha / nb1_se
  nb1_pval <- pt(nb1_t, df = df, lower.tail = FALSE)
  
  test <- list(nb2 = list(
    alpha = nb2_alpha,
    se = nb2_se,
    tstat = nb2_t,
    pval = nb2_pval
  ), nb1 = list(
    alpha = nb1_alpha,
    se = nb1_se,
    tstat = nb1_t,
    pval = nb1_pval
  ), nobs = length(y))
  
  class(test) <- "OVDtest"
  return(test)
}

#' @export
print.OVDtest <- function(x, digits = 4, ...)
{
  if(!inherits(x, "OVDtest"))
    cli::cli_abort("`x` needs to be of class `'OVDtest'`")
  
  cat("\nCameron and Trivedi Overdispersion Test:",
      "--------------------------------------",
      "  H0: Poisson (alpha = 0)",
      "  H1: Overdispersion (alpha > 0)",
      "--------------------------------------",
      sep = "\n")
  
  # NB2
  res <- x$nb2
  
  stars <- cut(res$pval,
               breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), 
               labels = c("***", "**", "*", ".", " "))
  
  nb2_df <- data.frame(
    Estimate = format(round(res$alpha, digits), nsmall = digits),
    `t-stat`  = format(round(res$tstat, digits), nsmall = digits),
    `p-value` = format(round(res$pval, digits), nsmall = digits),
    ` `      = as.character(stars),
    check.names = FALSE
  )
  
  res <- x$nb1
  
  stars <- cut(res$pval,
               breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), 
               labels = c("***", "**", "*", ".", " "))
  
  nb1_df <- data.frame(
    Estimate = format(round(res$alpha, digits), nsmall = digits),
    `t-stat`  = format(round(res$tstat, digits), nsmall = digits),
    `p-value` = format(round(res$pval, digits), nsmall = digits),
    ` `      = as.character(stars),
    check.names = FALSE
  )
  
  out <- rbind(nb2_df, nb1_df)
  rownames(out) <- c("NB2", "NB1")
  
  print(out, quote = FALSE, right = TRUE)
  
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat("Observations:", x$nobs, "\n")
  cat("Note: P-values are based on a one-tailed t-test (Right Tail).\n\n")
  
  invisible(x)
}