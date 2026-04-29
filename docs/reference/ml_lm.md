# Fit linear model by Maximum Likelihood

Fit linear model by Maximum Likelihood

## Usage

``` r
ml_lm(
  value,
  scale = NULL,
  weights = NULL,
  data,
  subset = NULL,
  noint_value = FALSE,
  noint_scale = FALSE,
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

- scale:

  Formula for log(sigma) (optional). If `NULL`, a homoskedastic model is
  fitted.

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

- noint_scale:

  Logical. Should the scale equation omit the intercept? Default is
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

An object of class `ml_lm` that extends `mlmodel`.

## Details

**Important:** Do not use the usual R syntax to remove the intercept in
the formula (`- 1` or `+ 0`) for the value or scale equations. Use the
dedicated arguments `noint_value` and `noint_scale` instead.

Coefficient names in the fitted object use the prefixes `value::` and
`scale::` to clearly identify to which equation each coefficient belongs
to, and to avoid confusion when the same variable(s) appear(s) in both
the value and scale equations.

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

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Homoskedastic linear model
data(mroz)
fit_lin <- ml_lm(faminc ~ age + I(age^2) + huswage + educ + unem, 
                 data = mroz)
summary(fit_lin, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = faminc ~ age + I(age^2) + huswage + educ + unem, 
#>     data = mroz)
#> 
#> Log-Likelihood: -7840.22 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 325.926, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                          Estimate  Std. Error z value Pr(>|z|)     
#> Value (faminc):  
#>   value::(Intercept) -29477.7659   8411.2714  -3.505 0.000457 ***
#>   value::age           1262.3145    382.2848   3.302 0.000960 ***
#>   value::I(age^2)       -13.5649      4.3919  -3.089 0.002011 ** 
#>   value::huswage       1956.5535    135.9529  14.391  < 2e-16 ***
#>   value::educ           963.6067    165.1769   5.834 5.42e-09 ***
#>   value::unem          -253.8058     98.6769  -2.572 0.010109 *  
#> Scale (log(sigma)):  
#>   scale::lnsigma          8.9930      0.0556 161.889  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 15694.44  BIC: 15726.80 
#> Residual standard error (sigma): 8047 

# Heteroskedastic linear model
fit_het <- ml_lm(faminc ~ age + I(age^2) + huswage + educ + unem,
                 scale = ~ educ + exper,
                 data = mroz)
summary(fit_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = faminc ~ age + I(age^2) + huswage + educ + unem, 
#>     scale = ~educ + exper, data = mroz)
#> 
#> Log-Likelihood: -7823.32 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 455.971, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(5) = 437.442, Pr(>Chisq) = < 1e-8
#>  Scale: Chisq(2) = 10.428, Pr(>Chisq) = 0.0054
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                          Estimate  Std. Error z value Pr(>|z|)     
#> Value (faminc):  
#>   value::(Intercept) -28864.2463   7794.2680  -3.703 0.000213 ***
#>   value::age           1282.8031    364.8667   3.516 0.000438 ***
#>   value::I(age^2)       -13.9022      4.1614  -3.341 0.000836 ***
#>   value::huswage       1933.1470    129.2254  14.960  < 2e-16 ***
#>   value::educ           879.5471    129.8162   6.775 1.24e-11 ***
#>   value::unem          -211.4286     92.0444  -2.297 0.021617 *  
#> Scale (log(sigma)):  
#>   scale::(Intercept)      8.4627      0.3170  26.694  < 2e-16 ***
#>   scale::educ             0.0503      0.0229   2.200 0.027814 *  
#>   scale::exper           -0.0103      0.0063  -1.652 0.098598 .  
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5634 Adjusted R-squared: 0.5605
#> AIC: 15664.64  BIC: 15706.25 
#> 
#> Distribution of Std. Deviation (sigma):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    4368    7251    7898    7942    8479   11016 
#> 

# Lognormal (log-linear) model
fit_log <- ml_lm(log(faminc) ~ age + I(age^2) + huswage + educ + unem,
                 data = mroz)
summary(fit_log, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Lognormal Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = log(faminc) ~ age + I(age^2) + huswage + educ + 
#>     unem, data = mroz)
#> 
#> Log-Likelihood: -7774.72 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 332.611, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (log(faminc)):  
#>   value::(Intercept)  7.55088    0.38075  19.832  < 2e-16 ***
#>   value::age          0.05826    0.01732   3.365 0.000766 ***
#>   value::I(age^2)    -0.00062    0.00020  -3.137 0.001707 ** 
#>   value::huswage      0.07577    0.00637  11.892  < 2e-16 ***
#>   value::educ         0.04799    0.00669   7.176 7.18e-13 ***
#>   value::unem        -0.01141    0.00463  -2.465 0.013702 *  
#> Scale (log(sigma)):  
#>   scale::lnsigma     -1.01459    0.04204 -24.132  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5069 Adjusted R-squared: 0.5036
#> AIC: 15563.44  BIC: 15595.81 
#> Residual standard error (sigma): 0.3626 

# Different predict types
head(predict(fit_log, type = "response")$fit)      # Expected value E[y]
#> [1] 15799.37 19816.98 16049.73 15759.18 25483.20 20879.54
head(predict(fit_log, type = "median")$fit)        # Median of y
#> [1] 14794.40 18556.45 15028.82 14756.75 23862.25 19551.42
head(predict(fit_log, type = "variance_y")$fit)    # Variance of y
#> [1] 35065084 55165838 36185150 34886880 91222728 61240248
head(predict(fit_log, type = "var")$fit)           # Variance of log(y)
#> [1] 0.1314437 0.1314437 0.1314437 0.1314437 0.1314437 0.1314437

# Fitted values and residuals
head(fitted(fit_lin))
#> [1] 15202.65 21471.12 15386.32 14983.68 27262.99 21921.35
head(residuals(fit_lin))
#> [1]  1107.34544   328.88277  5653.67868 -7683.67928    37.01175 -2426.35278
head(residuals(fit_lin, type = "pearson"))
#> [1]  0.137612157  0.040870957  0.702594595 -0.954867058  0.004599528
#> [6] -0.301527987
```
