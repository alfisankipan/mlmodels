# Fit negative binomial models by Maximum Likelihood

Fit negative binomial models by Maximum Likelihood

## Usage

``` r
ml_negbin(
  value,
  scale = NULL,
  dispersion = "NB2",
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

  Formula for the dispersion parameter log(alpha) (optional). If `NULL`,
  a constant alpha is used for all observations.

- dispersion:

  Either NB1 (proportional to mean variance), or NB2 (quadratic to mean
  variance). Defaults to NB2.

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

An object of class `ml_negbin` that extends `mlmodel.count` and
`mlmodel`.

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

## See also

[ml_poisson](https://alfisankipan.github.io/mlmodels/reference/ml_poisson.md)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Homoskedastic NB2 model (default dispersion)
data(docvis)
fit_nb2 <- ml_negbin(docvis ~ age + educyr + totchr, 
                     data = docvis)

summary(fit_nb2, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Negative Binomial (NB2) Model 
#> ---------------------------------------
#> Call:
#> ml_negbin(value = docvis ~ age + educyr + totchr, data = docvis)
#> 
#> Log-Likelihood: -10625.13 
#> 
#> Wald significance tests:
#>  all: Chisq(3) = 611.258, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept)   0.5124     0.2158   2.374   0.0176 *  
#>   value::age           0.0060     0.0027   2.182   0.0291 *  
#>   value::educyr        0.0289     0.0042   6.891 5.55e-12 ***
#>   value::totchr        0.3022     0.0125  24.200  < 2e-16 ***
#> Scale (log(alpha)):  
#>   scale::lnalpha      -0.4194     0.0382 -10.978  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3672
#> Pseudo R-squared - Cor.Sq.: 0.1366 McFadden: 0.03197
#> AIC: 21260.26  BIC: 21291.30 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 1.2816 
#>   Zeros - Observed: 401 Predicted: 338.95 

# Homoskedastic NB1 model
fit_nb1 <- ml_negbin(docvis ~ age + educyr + totchr, 
                     dispersion = "NB1",
                     data = docvis)

summary(fit_nb1, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Negative Binomial (NB1) Model 
#> ---------------------------------------
#> Call:
#> ml_negbin(value = docvis ~ age + educyr + totchr, dispersion = "NB1", 
#>     data = docvis)
#> 
#> Log-Likelihood: -10563.04 
#> 
#> Wald significance tests:
#>  all: Chisq(3) = 813.461, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept)   0.6716     0.1618   4.151 3.32e-05 ***
#>   value::age           0.0051     0.0020   2.533   0.0113 *  
#>   value::educyr        0.0268     0.0036   7.512 5.81e-14 ***
#>   value::totchr        0.2693     0.0097  27.856  < 2e-16 ***
#> Scale (log(alpha)):  
#>   scale::lnalpha       1.4792     0.0451  32.766  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3672
#> Pseudo R-squared - Cor.Sq.: 0.1397 McFadden: 0.03763
#> AIC: 21136.08  BIC: 21167.13 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 1.2087 
#>   Zeros - Observed: 401 Predicted: 395.03 

# Heteroskedastic NB2 model
fit_nb2_het <- ml_negbin(docvis ~ age + educyr + totchr,
                         scale = ~ female + bh,
                         data = docvis)

summary(fit_nb2_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Negative Binomial (NB2) Model 
#> ---------------------------------------
#> Call:
#> ml_negbin(value = docvis ~ age + educyr + totchr, scale = ~female + 
#>     bh, data = docvis)
#> 
#> Log-Likelihood: -10611.44 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 659.260, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(3) = 587.792, Pr(>Chisq) = < 1e-8
#>  Scale: Chisq(2) = 19.635, Pr(>Chisq) = 0.0001
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept)   0.5789     0.2104   2.752  0.00592 ** 
#>   value::age           0.0052     0.0026   1.967  0.04921 *  
#>   value::educyr        0.0288     0.0042   6.791 1.11e-11 ***
#>   value::totchr        0.2996     0.0126  23.866  < 2e-16 ***
#> Scale (log(alpha)):  
#>   scale::(Intercept)  -0.3848     0.0583  -6.606 3.96e-11 ***
#>   scale::female       -0.1917     0.0732  -2.621  0.00878 ** 
#>   scale::bh            0.3145     0.0985   3.193  0.00141 ** 
#> ---------------------------------------
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3670
#> Pseudo R-squared - Cor.Sq.: 0.137 McFadden: 0.03322
#> AIC: 21236.88  BIC: 21280.35 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 1.2531 
#>   Zeros - Observed: 401 Predicted: 342.61 
#> 
#> Distribution of Dispersion (alpha):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    0.56    0.56    0.68    0.67    0.77    0.93 
#> 

# Different predict types
head(predict(fit_nb2, type = "response")$fit)   # Expected count
#> [1] 10.580216  6.312773  6.843628  9.521127  5.323215  4.818987
head(predict(fit_nb2, type = "var")$fit)       # Variance
#> [1] 84.17289 32.51184 37.63424 69.11783 23.95239 20.08610
head(predict(fit_nb2, type = "alpha")$fit)     # Dispersion parameter
#> [1] 0.657424 0.657424 0.657424 0.657424 0.657424 0.657424

# Fitted values and residuals
head(fitted(fit_nb2))
#> [1] 10.580216  6.312773  6.843628  9.521127  5.323215  4.818987
head(residuals(fit_nb2))
#> [1] -6.5802161 -0.3127731 -4.8436275  1.4788730 -2.3232151 -2.8189871
head(residuals(fit_nb2, type = "pearson"))
#> [1] -0.71722270 -0.05485403 -0.78954919  0.17788356 -0.47469539 -0.62899215
```
