# Information Matrix Test for Model Misspecification

Performs the Information Matrix (IM) test for misspecification on models
fitted with the `mlmodels` package.

## Usage

``` r
IMtest(object, ...)

# S3 method for class 'mlmodel'
IMtest(object, method = "quad", repetitions = 999, seed = 1234L, ...)
```

## Arguments

- object:

  A fitted model object inheriting from `"mlmodel"`.

- ...:

  Further arguments passed to methods (currently not used).

- method:

  Character string. Specifies the version of the test:

  - `"quad"` (default): Quadratic form of the Information Matrix test
    (most common).

  - `"opg"`: Outer Product of Gradients version (Chesher-Lancaster).

  - `"boot_quad"`: Analytical chi-square and p-value, plus bootstrap
    p-value for the quadratic form.

  - `"boot_opg"`: Analytical chi-square and p-value, plus bootstrap
    p-value for the OPG version.

- repetitions:

  Integer. Number of bootstrap replications when using a bootstrap
  method. Default is 999.

- seed:

  Integer. Random seed for reproducibility in bootstrap methods. If
  `NULL`, a random seed is generated.

## Value

An object of class `"IMtest.mlmodel"` containing the analytical test
statistic, degrees of freedom and p-value, plus the bootstrapped p-value
(if a bootstrap method was selected).

## Details

The Information Matrix test checks whether the model is correctly
specified by testing the equality between the Hessian and the outer
product of the gradient (information matrix equality). Rejection of the
null hypothesis indicates model misspecification (e.g., incorrect
functional form, heteroskedasticity not properly modeled, omitted
variables, etc.).

Two main versions are implemented:

- **Quadratic form** (`"quad"`): Generally preferred for its better
  finite-sample properties.

- **OPG version** (`"opg"`): Chesher and Lancaster (1983) version.

Bootstrap versions (`"boot_quad"` and `"boot_opg"`) provide p-values
based on the empirical distribution of the test statistic and are useful
when asymptotic approximations may be unreliable.

## References

Chesher, A. (1983). The information matrix test: Simplified calculation
via a score test interpretation. Economics Letters, 13(1), 45-48.

Lancaster, T. (1984). The covariance matrix of the information matrix
test. Econometrica, 52(4), 1051-1053.

White, H. (1982). Maximum likelihood estimation of misspecified models.
Econometrica, 50(1), 1-25.

## See also

[`waldtest()`](https://alfisankipan.github.io/mlmodels/reference/waldtest.md),
[`lrtest()`](https://alfisankipan.github.io/mlmodels/reference/lrtest.md),
[`vuongtest()`](https://alfisankipan.github.io/mlmodels/reference/vuongtest.md)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Linear model example
data(mroz)
mroz$incthou <- mroz$faminc / 1000

fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

# Default quadratic form test
IMtest(fit)
#> Information Matrix Test
#>  Method: Orthogonalized Quadratic Form 
#>  Model:  Homoskedastic Linear Model 
#> --------------------------------------------
#>  Chisq(28) = 122.732    Pr(>Chisq) = 0.0000
#> --------------------------------------------

# OPG version
IMtest(fit, method = "opg")
#> Information Matrix Test
#>  Method: Chesher/Lancaster OPG 
#>  Model:  Homoskedastic Linear Model 
#> --------------------------------------------
#>  Chisq(28) = 298.563    Pr(>Chisq) = 0.0000
#> --------------------------------------------

# Bootstrap p-value (quadratic form)
IMtest(fit, method = "boot_quad", repetitions = 100, seed = 123)
#> Information Matrix Test
#>  Method: Orthogonalized Quadratic Form + Model-based bootstrap 
#>  Model:  Homoskedastic Linear Model 
#> --------------------------------------------
#>  Repetitions: Total 100 Successful 100 
#>  Chisq(28) = 122.732
#>  P(>Chisq): Analytical   = 0.0000 
#>             Bootstrapped = 0.7700
#> --------------------------------------------

# Heteroskedastic model
fit_het <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem,
                 scale = ~ educ, data = mroz)
IMtest(fit_het)
#> Information Matrix Test
#>  Method: Orthogonalized Quadratic Form 
#>  Model:  Heteroskedastic Linear Model 
#> --------------------------------------------
#>  Chisq(36) = 163.110    Pr(>Chisq) = 0.0000
#> --------------------------------------------
```
