#
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
  print(nb2_alpha)
  print(nb2_se)
  print(nb2_t)
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