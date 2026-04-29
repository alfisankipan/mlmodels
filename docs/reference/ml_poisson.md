# Fit Poisson model by Maximum Likelihood

Fit Poisson model by Maximum Likelihood

## Usage

``` r
ml_poisson(
  value,
  weights = NULL,
  data,
  subset = NULL,
  noint_value = FALSE,
  constraints = NULL,
  start = NULL,
  method = "NR",
  control = NULL,
  ...
)
```

## Arguments

- value:

  Formula for the conditional mean (value) equation.

- weights:

  Optional weights variable. It can be either the name of the variable
  in `data`, or a vector with the weights.

- data:

  Data frame.

- subset:

  Optional subset expression. Only observations for which this
  expression evaluates to `TRUE` are used in the estimation. This can be
  a logical vector or an expression (e.g. `subset = age > 30`).

- noint_value:

  Logical. Should the value equation omit the intercept? Default is
  `FALSE`.

- constraints:

  Optional constraints on the parameters. Can be a character vector of
  string constraints, a named list of string constraints, or a raw
  maxLik constraints list. See **Details**.

- start:

  Numeric vector of starting values for the coefficients. Required if
  constraints are being supplied. If supplied without constraints they
  will be ignored. See **Details**.

- method:

  A string with the method used for optimization. See
  [maxLik](https://rdrr.io/pkg/maxLik/man/maxLik.html) for options, and
  see **Details**.

- control:

  A list of control parameters passed to
  [maxLik](https://rdrr.io/pkg/maxLik/man/maxLik.html). If `NULL`
  (default), a sensible set of options is chosen automatically depending
  on whether constraints are used. See
  [maxControl](https://rdrr.io/pkg/maxLik/man/maxControl.html).

- ...:

  Additional arguments passed to
  [maxLik](https://rdrr.io/pkg/maxLik/man/maxLik.html).

## Value

An object of class `ml_poisson` that extends `mlmodel.count` and
`mlmodel`.

## Details

**Important:** Do not use the usual R syntax to remove the intercept in
the formula (`- 1` or `+ 0`) for the value equation. Use the dedicated
argument `noint_value` instead.

Coefficient names in the fitted object use the prefixes `value::`. This
is for consistency with other `mlmodel` estimators that model the scale
(dispersion) as well.

Either inequality or equality linear constraints are accepted, but not
both. A constraint cannot have a linear combination of more than two
coefficients.

**Important**: When `constraints` are supplied, `start` cannot be
`NULL`. You **must** provide initial values that yield a feasible
log-likelihood. If no constraints are used, any supplied `start` is
ignored.

When constraints are used, `ml_lm` automatically chooses the optimizer:

- Equality constraints =\> Nelder-Mead (`"NM"`)

- Inequality constraints =\> BFGS (`"BFGS"`)

In these cases your supplied `method` argument (if any) is ignored.

The Poisson model assumes equidispersion (mean = variance). When the
data show overdispersion (as is common), consider using
[ml_negbin](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md)
instead.

## See also

[ml_negbin](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Poisson model
data(docvis)
fit_pois <- ml_poisson(docvis ~ age + educyr + totchr, 
                       data = docvis)

summary(fit_pois, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Poisson 
#> ---------------------------------------
#> Call:
#> ml_poisson(value = docvis ~ age + educyr + totchr, data = docvis)
#> 
#> Log-Likelihood: -15200.67 
#> 
#> Wald significance tests:
#>  all: Chisq(3) = 610.354, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept)   0.6733     0.2035   3.309 0.000936 ***
#>   value::age           0.0046     0.0025   1.815 0.069518 .  
#>   value::educyr        0.0286     0.0041   6.896 5.33e-12 ***
#>   value::totchr        0.2749     0.0114  24.206  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3673
#> Pseudo R-squared - Cor.Sq.: 0.1394 McFadden: 0.1192
#> AIC: 30409.34  BIC: 30434.18 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 6.5225 
#>   Zeros - Observed: 401 Predicted: 23.95 

# Different predict types
head(predict(fit_pois, type = "response")$fit)   # Expected count
#> [1] 10.169550  6.266430  6.795981  9.255944  5.471075  5.002632
head(predict(fit_pois, type = "P(2,)")$fit)      # Probability of at least 2
#> [1] 0.9995720 0.9862011 0.9912821 0.9990201 0.9727781 0.9596609
head(predict(fit_pois, type = "P(3)")$fit)       # Probability of exactly 3
#> [1] 0.006716988 0.077881334 0.058498950 0.012627153 0.114817735 0.140226107

# Fitted values and residuals
head(fitted(fit_pois))
#> [1] 10.169550  6.266430  6.795981  9.255944  5.471075  5.002632
head(residuals(fit_pois))
#> [1] -6.1695499 -0.2664299 -4.7959814  1.7440558 -2.4710755 -3.0026324
head(residuals(fit_pois, type = "pearson"))
#> [1] -1.9346509 -0.1064322 -1.8397186  0.5732578 -1.0564517 -1.3424647
```
