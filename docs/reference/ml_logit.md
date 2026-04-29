# Fit Binary Logit Model by Maximum Likelihood

Fit Binary Logit Model by Maximum Likelihood

## Usage

``` r
ml_logit(
  value,
  scale = NULL,
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

  Two-sided formula for the probability equation.

- scale:

  Optional one-sided formula for heteroskedasticity.

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

An object of class `ml_logit` that extends `mlmodel`.

## Details

**Important:** Do not use the usual R syntax to remove the intercept in
the formula (`- 1` or `+ 0`) for the value equation. Use the dedicated
argument `noint_value` instead. For the scale equation (if modeling
heteroskedasticity), the formula must contain only the predictors
(right-hand side).

`ml_logit()` handles both strictly binary (`0/1`), and fractional
response (`0 <= y <= 1`) outcomes. When using fractional responses, it
is recommended to use robust standard errors (`vcov.type = "robust"`).

Coefficient names in the fitted object use the prefixes `value::` and
`scale::` (when heteroskedasticity is modeled) to clearly identify which
equation each coefficient belongs to.

Either inequality or equality linear constraints are accepted, but not
both. A constraint cannot have a linear combination of more than two
coefficients.

**Important**: When `constraints` are supplied, `start` cannot be
`NULL`. You **must** provide initial values that yield a feasible
log-likelihood. If no constraints are used, any supplied `start` is
ignored.

When constraints are used, `ml_logit` automatically chooses `method`:

- Equality constraints =\> Nelder-Mead (`"NM"`)

- Inequality constraints =\> BFGS (`"BFGS"`)

In these cases your supplied `method` argument (if any) is ignored.

## See also

[ml_probit](https://alfisankipan.github.io/mlmodels/reference/ml_probit.md)
[ml_beta](https://alfisankipan.github.io/mlmodels/reference/ml_beta.md)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Homoskedastic binary logit model
data(smoke)
smoke$smokes <- smoke$cigs > 0

fit_logit <- ml_logit(smokes ~ cigpric + income + age, 
                      data = smoke)

summary(fit_logit, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Binary Logit 
#> ---------------------------------------
#> Call:
#> ml_logit(value = smokes ~ cigpric + income + age, data = smoke)
#> 
#> Log-Likelihood: -533.27 
#> 
#> Wald significance tests:
#>  all: Chisq(3) = 9.010, Pr(>Chisq) = 0.0292
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                         Estimate Std. Error z value Pr(>|z|)    
#> Value (smokes):  
#>   value::(Intercept)  0.6885299  0.9527423   0.723  0.46988   
#>   value::cigpric     -0.0104209  0.0154879  -0.673  0.50105   
#>   value::income      -0.0000019  0.0000080  -0.243  0.80786   
#>   value::age         -0.0121257  0.0041786  -2.902  0.00371 **
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:807 (Successes: 310, Failures: 497)
#> Pseudo R-squared - Cor.Sq.: 0.009815 McFadden: 0.007888
#> AIC: 1074.53  BIC: 1093.31 

# Heteroskedastic binary logit model
fit_logit_het <- ml_logit(smokes ~ cigpric + income + age,
                          scale = ~ educ,
                          data = smoke)
#> ℹ Improving initial values by scaling (factor = 0.5).
#> ℹ Initial log-likelihood: -555.841
#> ℹ Final scaled log-likelihood: -537.617

summary(fit_logit_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Binary Logit 
#> ---------------------------------------
#> Call:
#> ml_logit(value = smokes ~ cigpric + income + age, scale = ~educ, 
#>     data = smoke)
#> 
#> Log-Likelihood: -522.93 
#> 
#> Wald significance tests:
#>  all: Chisq(4) = 215.327, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(3) = 0.761, Pr(>Chisq) = 0.8588
#>  Scale: Chisq(1) = 13.936, Pr(>Chisq) = 0.0002
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                         Estimate Std. Error z value Pr(>|z|)     
#> Value (smokes):  
#>   value::(Intercept) -6.334e-03  3.278e-02  -0.193 0.846792    
#>   value::cigpric      8.513e-05  5.654e-04   0.151 0.880325    
#>   value::income      -9.200e-09  2.932e-07  -0.031 0.974960    
#>   value::age         -4.163e-04  4.950e-04  -0.841 0.400312    
#> Scale (log(sigma)):  
#>   scale::educ        -2.439e-01  6.534e-02  -3.733 0.000189 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:807 (Successes: 310, Failures: 497)
#> Pseudo R-squared - Cor.Sq.: 0.03261 McFadden: 0.02712
#> AIC: 1055.86  BIC: 1079.32 
#> 
#> Distribution of Std. Deviation (sigma):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.012   0.037   0.054   0.063   0.087   0.231 
#> 

# Different predict types
head(predict(fit_logit, type = "response")$fit)   # Predicted probability
#> [1] 0.3684765 0.3874225 0.3375703 0.4213065 0.4590092 0.2718363
head(predict(fit_logit, type = "link")$fit)       # Linear predictor (log-odds)
#> [1] -0.5387582 -0.4581597 -0.6741409 -0.3174125 -0.1643319 -0.9853257

# Fitted values and residuals
head(fitted(fit_logit))
#> [1] 0.3684765 0.3874225 0.3375703 0.4213065 0.4590092 0.2718363
head(residuals(fit_logit))
#> [1] -0.3684765 -0.3874225  0.6624297 -0.4213065 -0.4590092 -0.2718363
head(residuals(fit_logit, type = "pearson"))
#> [1] -0.7638536 -0.7952650  1.4008377 -0.8532470 -0.9211191 -0.6109972
```
