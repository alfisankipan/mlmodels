# Summary for mlmodel objects

Summary for mlmodel objects

## Usage

``` r
# S3 method for class 'ml_beta'
summary(
  object,
  correlation = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_gamma'
summary(
  object,
  correlation = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_lm'
summary(
  object,
  correlation = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_logit'
summary(
  object,
  correlation = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'mlmodel'
summary(
  object,
  correlation = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_negbin'
summary(
  object,
  correlation = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_poisson'
summary(
  object,
  correlation = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_probit'
summary(
  object,
  correlation = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted model object of class `"mlmodel"`.

- correlation:

  Logical. Should the correlation matrix of the estimated parameters be
  included in the output? Default is `FALSE`. If `TRUE` the correlation
  matrix will be computed, and stored in the `'summary.mlmodel'` object
  the function returns.

- vcov:

  Optional user-supplied variance-covariance matrix. If provided, it
  will be used instead of computing one internally.

- vcov.type:

  Character string specifying the type of variance-covariance matrix to
  use. See
  [vcov](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md).

- cl_var:

  Character string or vector. Name of the clustering variable or the
  vector itself. See
  [vcov](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md).

- repetitions:

  Integer. Number of bootstrap replications when `vcov.type = "boot"`.
  Default is 999.

- seed:

  Integer. Random seed for reproducibility when `vcov.type = "boot"`. If
  `NULL`, a random seed is generated.

- progress:

  Logical. Should a progress bar be displayed during bootstrapping or
  jackknifing? Default is `FALSE` (silent).

- ...:

  Further arguments passed to methods.

## Details

Coefficient names in the fitted object use the prefixes `value::` and
`scale::` to identify to which equation they belong to, and to avoid
confusion when the same variable(s) appear(s) in both the value and
scale equations.

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

data(mroz)
mroz$incthou <- mroz$faminc / 1000

fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

# Default: observed information matrix
summary(fit)
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = incthou ~ age + I(age^2) + huswage + educ + unem, 
#>     data = mroz)
#> 
#> Log-Likelihood: -2638.68 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 972.786, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -29.4778     8.6369  -3.413 0.000643 ***
#>   value::age           1.2623     0.3996   3.159 0.001585 ** 
#>   value::I(age^2)     -0.0136     0.0046  -2.943 0.003250 ** 
#>   value::huswage       1.9566     0.0733  26.691  < 2e-16 ***
#>   value::educ          0.9636     0.1356   7.106 1.19e-12 ***
#>   value::unem         -0.2538     0.0957  -2.651 0.008028 ** 
#> Scale (log(sigma)):  
#>   scale::lnsigma       2.0853     0.0258  80.924  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 5291.36  BIC: 5323.72 
#> Residual standard error (sigma): 8.047 

# Different variance types
summary(fit, vcov.type = "opg")                    # Outer product of gradients
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = incthou ~ age + I(age^2) + huswage + educ + unem, 
#>     data = mroz)
#> 
#> Log-Likelihood: -2638.68 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 2553.948, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Outer Product of Gradients (BHHH)
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -29.4778     9.2817  -3.176  0.00149 ** 
#>   value::age           1.2623     0.4318   2.924  0.00346 ** 
#>   value::I(age^2)     -0.0136     0.0050  -2.712  0.00670 ** 
#>   value::huswage       1.9566     0.0443  44.161  < 2e-16 ***
#>   value::educ          0.9636     0.1167   8.254  < 2e-16 ***
#>   value::unem         -0.2538     0.0968  -2.622  0.00874 ** 
#> Scale (log(sigma)):  
#>   scale::lnsigma       2.0853     0.0144 144.847  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 5291.36  BIC: 5323.72 
#> Residual standard error (sigma): 8.047 
summary(fit, vcov.type = "robust")                 # Robust/sandwich estimator
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = incthou ~ age + I(age^2) + huswage + educ + unem, 
#>     data = mroz)
#> 
#> Log-Likelihood: -2638.68 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 325.926, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -29.4778     8.4113  -3.505 0.000457 ***
#>   value::age           1.2623     0.3823   3.302 0.000960 ***
#>   value::I(age^2)     -0.0136     0.0044  -3.089 0.002011 ** 
#>   value::huswage       1.9566     0.1360  14.391  < 2e-16 ***
#>   value::educ          0.9636     0.1652   5.834 5.42e-09 ***
#>   value::unem         -0.2538     0.0987  -2.572 0.010109 *  
#> Scale (log(sigma)):  
#>   scale::lnsigma       2.0853     0.0556  37.538  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 5291.36  BIC: 5323.72 
#> Residual standard error (sigma): 8.047 

# Clustered robust standard errors
summary(fit, vcov.type = "robust", cl_var = "age")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = incthou ~ age + I(age^2) + huswage + educ + unem, 
#>     data = mroz)
#> 
#> Log-Likelihood: -2638.68 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 393.650, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust | Clusters: 31 (age)
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -29.4778     7.2386  -4.072 4.66e-05 ***
#>   value::age           1.2623     0.3257   3.876 0.000106 ***
#>   value::I(age^2)     -0.0136     0.0038  -3.538 0.000403 ***
#>   value::huswage       1.9566     0.1445  13.541  < 2e-16 ***
#>   value::educ          0.9636     0.1402   6.873 6.31e-12 ***
#>   value::unem         -0.2538     0.0804  -3.156 0.001598 ** 
#> Scale (log(sigma)):  
#>   scale::lnsigma       2.0853     0.0592  35.253  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 5291.36  BIC: 5323.72 
#> Residual standard error (sigma): 8.047 

# Using a pre-computed variance matrix (e.g. bootstrap)
v_boot <- vcov(fit, type = "boot", repetitions = 100, seed = 123)
#> ℹ Bootstrap with 100 repetitions.
#>  0        10        20        30        40        50
#> ====================================================
#>  ..................................................
#>  ..................................................
#> ====================================================
#> 
#> Bootstrapping finished - 100% of replications converged.
summary(fit, vcov = v_boot)
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = incthou ~ age + I(age^2) + huswage + educ + unem, 
#>     data = mroz)
#> 
#> Log-Likelihood: -2638.68 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 359.486, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Bootstrap (100/100 reps. - 100.00% rate)
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -29.4778     7.4135  -3.976 7.00e-05 ***
#>   value::age           1.2623     0.3413   3.699 0.000217 ***
#>   value::I(age^2)     -0.0136     0.0040  -3.420 0.000626 ***
#>   value::huswage       1.9566     0.1509  12.967  < 2e-16 ***
#>   value::educ          0.9636     0.1512   6.375 1.83e-10 ***
#>   value::unem         -0.2538     0.1052  -2.413 0.015827 *  
#> Scale (log(sigma)):  
#>   scale::lnsigma       2.0853     0.0526  39.622  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 5291.36  BIC: 5323.72 
#> Residual standard error (sigma): 8.047 
```
