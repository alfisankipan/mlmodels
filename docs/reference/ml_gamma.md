# Fit Gamma Model by Maximum Likelihood

Fit Gamma Model by Maximum Likelihood

## Usage

``` r
ml_gamma(
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

  Formula for log(nu) equation (shape parameter - optional). If `NULL`,
  a homoskedastic (constant shape) model is fitted.

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

An object of class `ml_gamma` that extends `mlmodel`.

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

The Gamma model requires a strictly positive response variable
(`y > 0`). Observations where `y <= 0` are automatically dropped with a
warning.

If your data contains zeros or non-positive values, consider using
[`ml_poisson()`](https://alfisankipan.github.io/mlmodels/reference/ml_poisson.md)
or
[`ml_negbin()`](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md)
instead, as they are frequently applied to continuous non-negative
outcomes.

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Homoskedastic gamma regression
data(mroz)
fit_gamma <- ml_gamma(faminc ~ hours + hushrs + age + educ, 
                      data = mroz)

summary(fit_gamma, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Gamma Model 
#> ---------------------------------------
#> Call:
#> ml_gamma(value = faminc ~ hours + hushrs + age + educ, data = mroz)
#> 
#> Log-Likelihood: -7957.39 
#> 
#> Wald significance tests:
#>  all: Chisq(4) = 130.046, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value  Pr(>|z|)     
#> Value (faminc):  
#>   value::(Intercept) 8.383279   0.157845  53.111   < 2e-16 ***
#>   value::hours       0.000073   0.000019   3.931 0.0000846 ***
#>   value::hushrs      0.000123   0.000035   3.478  0.000506 ***
#>   value::age         0.007990   0.002116   3.776  0.000159 ***
#>   value::educ        0.078817   0.008697   9.062   < 2e-16 ***
#> Scale (log(nu)):  
#>   scale::lnnu        1.610431   0.069391  23.208   < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:753 Deg. of freedom: 748
#> Pseudo R-squared - Cor.Sq.: 0.1653 McFadden: 0.009728
#> AIC: 15926.77  BIC: 15954.52 
#> Shape Param.: 5.00  - Coef.Var.: 0.45 

# Heteroskedastic gamma regression
fit_gamma_het <- ml_gamma(faminc ~ hours + hushrs + age + educ,
                          scale = ~ kidslt6,
                          data = mroz)

summary(fit_gamma_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Gamma Model 
#> ---------------------------------------
#> Call:
#> ml_gamma(value = faminc ~ hours + hushrs + age + educ, scale = ~kidslt6, 
#>     data = mroz)
#> 
#> Log-Likelihood: -7954.63 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 131.538, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(4) = 130.023, Pr(>Chisq) = < 1e-8
#>  Scale: Chisq(1) = 3.847, Pr(>Chisq) = 0.0498
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                        Estimate Std. Error z value  Pr(>|z|)     
#> Value (faminc):  
#>   value::(Intercept)  8.403495   0.157455  53.371   < 2e-16 ***
#>   value::hours        0.000074   0.000019   3.975 0.0000705 ***
#>   value::hushrs       0.000113   0.000034   3.309  0.000936 ***
#>   value::age          0.007908   0.002086   3.790  0.000150 ***
#>   value::educ         0.079100   0.008598   9.200   < 2e-16 ***
#> Scale (log(nu)):  
#>   scale::(Intercept)  1.667046   0.078418  21.258   < 2e-16 ***
#>   scale::kidslt6     -0.207617   0.105847  -1.961  0.049824 *  
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:753 Deg. of freedom: 748
#> Pseudo R-squared - Cor.Sq.: 0.1661 McFadden: 0.01007
#> AIC: 15923.26  BIC: 15955.63 
#> 
#> Distribution of Shape Related Params.:
#> ---------------------------------------
#>            Min. 1st Qu. Median Mean 3rd Qu. Max.
#> Shape (nu) 2.84    5.30   5.30 5.07    5.30 5.30
#> Coef. Var. 0.43    0.43   0.43 0.45    0.43 0.59
#> 

# Different predict types
head(predict(fit_gamma, type = "response")$fit)   # Expected value E[y]
#> [1] 22810.01 21452.12 25102.54 19332.81 24215.43 22864.27
head(predict(fit_gamma, type = "variance")$fit)   # Variance of y
#> [1] 103955935  91947317 125902327  74677310 117160917 104451145

# Fitted values and residuals
head(fitted(fit_gamma))
#> [1] 22810.01 21452.12 25102.54 19332.81 24215.43 22864.27
head(residuals(fit_gamma))
#> [1]  -6500.006    347.877  -4062.538 -12032.814   3084.573  -3369.271
head(residuals(fit_gamma, type = "pearson"))
#> [1] -0.63751306  0.03627908 -0.36205997 -1.39242843  0.28497299 -0.32966989

# Comparison with lognormal model (often very similar mean predictions)
fit_lognorm <- ml_lm(log(faminc) ~ hours + hushrs + age + educ, 
                     data = mroz)

head(predict(fit_gamma, type = "response")$fit)
#> [1] 22810.01 21452.12 25102.54 19332.81 24215.43 22864.27
head(predict(fit_lognorm, type = "response")$fit)
#> [1] 23359.88 22222.96 25712.64 19503.64 25342.53 24298.59
```
