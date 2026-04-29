# Wald Test for Linear Restrictions

Performs a Wald test of linear restrictions on the parameters of an
`mlmodel` object.

## Usage

``` r
waldtest(object, ...)

# S3 method for class 'mlmodel'
waldtest(
  object,
  indices = NULL,
  coef_names = NULL,
  rest_matrix = NULL,
  rhs = 0,
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

  An object of class `"mlmodel"`.

- ...:

  Further arguments passed to methods.

- indices:

  Integer vector. Positions of the coefficients to be tested.

- coef_names:

  Character vector. Names of the coefficients to test.

- rest_matrix:

  Numeric matrix. A q × k restriction matrix (advanced use).

- rhs:

  Numeric vector. Value(s) the linear combination(s) should equal.
  Default is 0.

- vcov:

  Optional user-supplied variance-covariance matrix.

- vcov.type:

  Character string. Type of variance-covariance matrix to use. One of
  `"oim"` (default), `"opg"`, `"robust"`, `"boot"`, or `"jack"`. See
  [`vcov.mlmodel()`](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md)
  for details.

- cl_var:

  Character string or vector. Clustering variable when
  `vcov.type = "robust"` or `"boot"`.

- repetitions:

  Integer. Number of bootstrap replications when `vcov.type = "boot"`.
  Default is 999.

- seed:

  Integer. Random seed for bootstrap.

- progress:

  Logical. Show progress bar during bootstrapping? Default `FALSE`.

## Value

An object of class `"waldtest.mlmodel"`.

## Details

The Wald test evaluates linear restrictions of the form \\R\beta = r\\.

Three convenient interfaces are provided:

- `indices` or `coef_names`: Test individual coefficients (or groups of
  coefficients) against the value(s) in `rhs` (defaults to 0, which is
  useful for joint significance tests).

- `rest_matrix` + `rhs`: Test general linear combinations of
  coefficients (advanced use).

The test statistic follows a \\\chi^2\\ distribution with degrees of
freedom equal to the number of restrictions under the null hypothesis.

## See also

[`lrtest()`](https://alfisankipan.github.io/mlmodels/reference/lrtest.md),
[`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md),
[`confint.mlmodel()`](https://alfisankipan.github.io/mlmodels/reference/confint.mlmodel.md)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

data(mroz)
mroz$incthou <- mroz$faminc / 1000

fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

# 1. Test single coefficients using indices (default OIM)
waldtest(fit, indices = c(2, 5))
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#> Restrictions:
#>   1: value::age = 0
#>   2: value::educ = 0
#> Chisq(2) = 59.026    Pr(>Chisq) = < 1e-8
#> --------------------------------------------
#> 

# 2. Test using coefficient names and robust standard errors
waldtest(fit, coef_names = c("value::educ", "value::unem"), 
         vcov.type = "robust")
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Robust
#> ---------------------------------------
#> Restrictions:
#>   1: value::educ = 0
#>   2: value::unem = 0
#> Chisq(2) = 43.248    Pr(>Chisq) = < 1e-8
#> --------------------------------------------
#> 
         
# 3. Test explicit constraints
waldtest(fit, coef_names = "value::educ", rhs = 1, vcov.type = "robust")
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Robust
#> ---------------------------------------
#> Restrictions:
#>   1: value::educ = 1
#> Chisq(1) = 0.049    Pr(>Chisq) = 0.8256
#> --------------------------------------------
#> 

# 4. Test a linear combination of two coefficients using a restriction matrix
# H0: educ + huswage = 3
R <- matrix(c(0, 0, 0, 1, 1, 0, 0), nrow = 1)
waldtest(fit, rest_matrix = R, rhs = 3, vcov.type = "boot", 
         repetitions = 100, seed = 123)
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Bootstrap (100/100 reps. - 100.00% rate)
#> ---------------------------------------
#> Restrictions:
#>   1: value::huswage + value::educ = 3
#> Chisq(1) = 0.209    Pr(>Chisq) = 0.6473
#> --------------------------------------------
#> 
```
