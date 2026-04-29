# Fit Binary Probit Model by Maximum Likelihood

Fit Binary Probit Model by Maximum Likelihood

## Usage

``` r
ml_probit(
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

An object of class `ml_probit` that extends `mlmodel`.

## Details

**Important:** Do not use the usual R syntax to remove the intercept in
the formula (`- 1` or `+ 0`) for the value equation. Use the dedicated
argument `noint_value` instead. For the scale equation (if modeling
heteroskedasticity), the formula must contain only the predictors
(right-hand side).

The dependent variable must be binary (0/1 or TRUE/FALSE).

Coefficient names in the fitted object use the prefixes `value::` and
`scale::` (when heteroskedasticity is modeled) to clearly identify which
equation each coefficient belongs to.

`ml_probit()` handles both strictly binary (`0/1`), and fractional
response (`0 <= y <= 1`) outcomes. When using fractional responses, it
is recommended to use robust standard errors (`vcov.type = "robust"`).

Either inequality or equality linear constraints are accepted, but not
both. A constraint cannot have a linear combination of more than two
coefficients.

**Important**: When `constraints` are supplied, `start` cannot be
`NULL`. You **must** provide initial values that yield a feasible
log-likelihood. If no constraints are used, any supplied `start` is
ignored.

When constraints are used, `ml_probit` automatically chooses `method`:

- Equality constraints =\> Nelder-Mead (`"NM"`)

- Inequality constraints =\> BFGS (`"BFGS"`)

In these cases your supplied `method` argument (if any) is ignored.

## See also

[ml_logit](https://alfisankipan.github.io/mlmodels/reference/ml_logit.md)
[ml_beta](https://alfisankipan.github.io/mlmodels/reference/ml_beta.md)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Probit model (strictly binary outcome)
data(smoke)
smoke$smokes <- smoke$cigs > 0

fit_probit <- ml_probit(smokes ~ cigpric + income + age, 
                        data = smoke)

summary(fit_probit, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Binary Probit 
#> ---------------------------------------
#> Call:
#> ml_probit(value = smokes ~ cigpric + income + age, data = smoke)
#> 
#> Log-Likelihood: -533.21 
#> 
#> Wald significance tests:
#>  all: Chisq(3) = 9.071, Pr(>Chisq) = 0.0284
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                         Estimate Std. Error z value Pr(>|z|)    
#> Value (smokes):  
#>   value::(Intercept)  0.4259399  0.5879224   0.724  0.46877   
#>   value::cigpric     -0.0064471  0.0095486  -0.675  0.49956   
#>   value::income      -0.0000011  0.0000049  -0.228  0.81999   
#>   value::age         -0.0075935  0.0026078  -2.912  0.00359 **
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:807 (Successes: 310, Failures: 497)
#> Pseudo R-squared - Cor.Sq.: 0.009922 McKelvey & Zavoina: 0.01753
#> AIC: 1074.42  BIC: 1093.20 

# Heteroskedastic probit model
fit_probit_het <- ml_probit(smokes ~ cigpric + income + age,
                            scale = ~ educ,
                            data = smoke)

summary(fit_probit_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Binary Probit 
#> ---------------------------------------
#> Call:
#> ml_probit(value = smokes ~ cigpric + income + age, scale = ~educ, 
#>     data = smoke)
#> 
#> Log-Likelihood: -522.88 
#> 
#> Wald significance tests:
#>  all: Chisq(4) = 224.005, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(3) = 0.791, Pr(>Chisq) = 0.8516
#>  Scale: Chisq(1) = 14.732, Pr(>Chisq) = 0.0001
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                         Estimate Std. Error z value Pr(>|z|)     
#> Value (smokes):  
#>   value::(Intercept) -2.060e-03  2.307e-02  -0.089 0.928856    
#>   value::cigpric      1.701e-05  3.907e-04   0.044 0.965268    
#>   value::income      -2.300e-09  2.020e-07  -0.011 0.990848    
#>   value::age         -2.953e-04  3.393e-04  -0.870 0.384057    
#> Scale (log(sigma)):  
#>   scale::educ        -2.340e-01  6.097e-02  -3.838 0.000124 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:807 (Successes: 310, Failures: 497)
#> Pseudo R-squared - Cor.Sq.: 0.03257 McKelvey & Zavoina: 0.06789
#> AIC: 1055.77  BIC: 1079.24 
#> 
#> Distribution of Std. Deviation (sigma):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.015   0.042   0.060   0.070   0.096   0.246 
#> 

# Different predict types
head(predict(fit_probit, type = "response")$fit)   # Probability of success
#> [1] 0.3684478 0.3879168 0.3372339 0.4216946 0.4595126 0.2686198
head(predict(fit_probit, type = "prob0")$fit)      # Probability of failure
#> [1] 0.6315522 0.6120832 0.6627661 0.5783054 0.5404874 0.7313802
head(predict(fit_probit, type = "link")$fit)       # Linear predictor (z)
#> [1] -0.3359671 -0.2847528 -0.4200241 -0.1975601 -0.1016618 -0.6169926

# Fitted values and residuals
head(fitted(fit_probit))
#> [1] 0.3684478 0.3879168 0.3372339 0.4216946 0.4595126 0.2686198
head(residuals(fit_probit))
#> [1] -0.3684478 -0.3879168  0.6627661 -0.4216946 -0.4595126 -0.2686198
head(residuals(fit_probit, type = "pearson"))
#> [1] -0.7638065 -0.7960934  1.4018918 -0.8539264 -0.9220530 -0.6060346
```
