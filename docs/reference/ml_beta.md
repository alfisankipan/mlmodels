# Fit Beta Model by Maximum Likelihood

Fit Beta Model by Maximum Likelihood

## Usage

``` r
ml_beta(
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

  Formula for the conditional log(mean) equation.

- scale:

  Formula for log(phi) equation (precision parameter - optional). If
  `NULL`, a homoskedastic (constant precision) model is fitted.

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

An object of class `ml_beta` that extends `mlmodel`.

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

The Beta model is only defined for a strictly fractional response
variable in the open interval (0, 1). Observations where the response is
`y <= 0` or `y >= 1` are automatically dropped with a warning. If your
data contains boundary values (0 or 1), consider using
[`ml_logit()`](https://alfisankipan.github.io/mlmodels/reference/ml_logit.md)
instead.

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Homoskedastic beta regression (fractional response)
data(pw401k)

# Beta regression requires 0 < y < 1. 
# Observations at the boundaries (0 or 1) are automatically dropped.
fit_beta <- ml_beta(prate ~ mrate + I(mrate^2) + log(totemp) + 
                    I(log(totemp)^2) + age + I(age^2) + sole, 
                    data = pw401k, 
                    subset = prate < 1)   # drop y = 1
#> ℹ Improving initial values by scaling (factor = 0.5).
#> ℹ Initial log-likelihood: -311.974
#> ℹ Final scaled log-likelihood: 79.945

summary(fit_beta, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Beta Model 
#> ---------------------------------------
#> Call:
#> ml_beta(value = prate ~ mrate + I(mrate^2) + log(totemp) + I(log(totemp)^2) + 
#>     age + I(age^2) + sole, data = pw401k, subset = prate < 1)
#> 
#> Log-Likelihood: 1677.81 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 351.410, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                            Estimate Std. Error z value Pr(>|z|)     
#> Value (prate):  
#>   value::(Intercept)       2.99043    0.32638   9.162  < 2e-16 ***
#>   value::mrate             0.90356    0.09600   9.412  < 2e-16 ***
#>   value::I(mrate^2)       -0.21881    0.02716  -8.057 7.80e-16 ***
#>   value::log(totemp)      -0.57702    0.08733  -6.608 3.91e-11 ***
#>   value::I(log(totemp)^2)  0.03076    0.00564   5.455 4.89e-08 ***
#>   value::age               0.04639    0.00557   8.332  < 2e-16 ***
#>   value::I(age^2)         -0.00065    0.00012  -5.229 1.70e-07 ***
#>   value::sole             -0.10352    0.04001  -2.587  0.00968 ** 
#> Scale (log(phi)):  
#>   scale::lnphi             1.87118    0.03502  53.429  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:2711 Deg. of freedom: 2703
#> Pseudo R-squared - Cor.Sq.: 0.1051
#> AIC: -3337.62  BIC: -3284.48 
#> Precision Param.: 6.50

# Heteroskedastic beta regression
fit_beta_het <- ml_beta(prate ~ mrate + I(mrate^2) + log(totemp) + 
                        I(log(totemp)^2) + age + I(age^2) + sole,
                        scale = ~ totemp + sole,
                        data = pw401k, 
                        subset = prate < 1)
#> ℹ Improving initial values by scaling (factor = 0.5).
#> ℹ Initial log-likelihood: -311.974
#> ℹ Final scaled log-likelihood: 79.945

summary(fit_beta_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Beta Model 
#> ---------------------------------------
#> Call:
#> ml_beta(value = prate ~ mrate + I(mrate^2) + log(totemp) + I(log(totemp)^2) + 
#>     age + I(age^2) + sole, scale = ~totemp + sole, data = pw401k, 
#>     subset = prate < 1)
#> 
#> Log-Likelihood: 1693.08 
#> 
#> Wald significance tests:
#>  all: Chisq(9) = 371.864, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(7) = 340.291, Pr(>Chisq) = < 1e-8
#>  Scale: Chisq(2) = 18.676, Pr(>Chisq) = 0.0001
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                              Estimate Std. Error z value Pr(>|z|)     
#> Value (prate):  
#>   value::(Intercept)       2.7863712  0.3490165   7.983 1.42e-15 ***
#>   value::mrate             0.8779909  0.0919348   9.550  < 2e-16 ***
#>   value::I(mrate^2)       -0.2112053  0.0258366  -8.175 2.97e-16 ***
#>   value::log(totemp)      -0.5103614  0.0952286  -5.359 8.35e-08 ***
#>   value::I(log(totemp)^2)  0.0261633  0.0063336   4.131 3.61e-05 ***
#>   value::age               0.0468373  0.0064764   7.232 4.76e-13 ***
#>   value::I(age^2)         -0.0006665  0.0001534  -4.344 1.40e-05 ***
#>   value::sole             -0.1624332  0.0411977  -3.943 8.05e-05 ***
#> Scale (log(phi)):  
#>   scale::(Intercept)       2.0001168  0.0435566  45.920  < 2e-16 ***
#>   scale::totemp           -0.0000054  0.0000028  -1.903   0.0571 .  
#>   scale::sole             -0.2804443  0.0684349  -4.098 4.17e-05 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:2711 Deg. of freedom: 2703
#> Pseudo R-squared - Cor.Sq.: 0.106
#> AIC: -3364.16  BIC: -3299.21 
#> 
#> Distribution of Predicted Precision (phi):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    0.67    5.58    7.23    6.64    7.36    7.39 
#> 


# Note: All predictions (including those from predict(), fitted(), and 
# residuals()) return values aligned to the original data, with NA 
# for observations dropped due to subset or boundary values.

# Different predict types
head(predict(fit_beta, type = "response")$fit)   # Expected value E[y]
#> [1] 0.8034086 0.7705266        NA 0.7500891 0.7293741        NA
head(predict(fit_beta, type = "variance")$fit)   # Variance of y
#> [1] 0.02107037 0.02358800         NA 0.02500743 0.02633242         NA
head(predict(fit_beta, type = "phi")$fit)        # Precision parameter
#> [1] 6.495988 6.495988       NA 6.495988 6.495988       NA

# Fitted values and residuals
head(fitted(fit_beta))
#> [1] 0.8034086 0.7705266        NA 0.7500891 0.7293741        NA
head(residuals(fit_beta))
#> [1] -0.14462481  0.07282318          NA  0.19091386  0.10077445          NA
head(residuals(fit_beta, type = "pearson"))
#> [1] -0.9963381  0.4741591         NA  1.2072658  0.6210193         NA
```
