# Variance-Covariance Matrix for mlmodel Objects

Returns the variance-covariance matrix of the estimated parameters using
different methods.

## Usage

``` r
# S3 method for class 'mlmodel'
vcov(
  object,
  type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = TRUE,
  ...
)
```

## Arguments

- object:

  An object of class `"mlmodel"` or a model inheriting from it (e.g.
  `"ml_lm"`).

- type:

  Character string specifying the type of variance-covariance matrix.
  One of `"oim"` (default), `"robust"`, `"opg"`, `"cluster"`, `"boot"`,
  and `"jack"` or `"jackknife"`.

- cl_var:

  Character string or vector. Name of the clustering variable in the
  data, or the vector itself.

- repetitions:

  Integer. Number of bootstrap replications to use when `type = "boot"`.
  Default is 999.

- seed:

  Integer. Random seed for reproducibility when `type = "boot"`. If
  `NULL`, a random seed is generated.

- progress:

  Logical. Should a progress bar be displayed? Default is `TRUE` when
  `type` is `"boot"` or `"jack"`/`"jackknife"`. Ignored for other types.

- ...:

  Further arguments passed to methods.

## Value

A symmetric variance-covariance matrix with coefficient names on the
rows and columns.

## Details

The package provides several variance-covariance estimators through the
`type` argument:

- `"oim"` - Observed Information Matrix (default)

- `"opg"` - Outer Product of Gradients (BHHH)

- `"robust"` - Robust (sandwich) estimator

- `"cluster"` - Alias for `"robust"` when clustering the variance but
  requires `cl_var` to be set.

- `"boot"` - Bootstrap (with optional clustering)

- `"jack"` - Jackknife (with optional clustering)

- `"jackknife"` - alias for `"jack"`

Clustered standard errors are obtained by setting `cl_var` when using
`"robust"`/`"cluster"`, `"boot"`, or `"jack"`.

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

data(mroz)
mroz$incthou <- mroz$faminc / 1000

fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

# Different variance-covariance estimators
v_oim   <- vcov(fit, type = "oim")      # Observed Information Matrix (default)
v_opg   <- vcov(fit, type = "opg")      # Outer Product of Gradients (BHHH)
v_robust <- vcov(fit, type = "robust")   # Robust / Sandwich estimator

# Clustered robust standard errors
v_clust <- vcov(fit, type = "robust", cl_var = "age")

# Bootstrap variance-covariance matrix
v_boot  <- vcov(fit, type = "boot", repetitions = 100, seed = 123)
#> ℹ Bootstrap with 100 repetitions.
#>  0        10        20        30        40        50
#> ====================================================
#>  ..................................................
#>  ..................................................
#> ====================================================
#> 
#> Bootstrapping finished - 100% of replications converged.

# Jackknife variance-covariance matrix
v_jack  <- vcov(fit, type = "jack")
#> ℹ Jackknife variance.
#>  0        10        20        30        40        50
#> ====================================================
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ...
#> ====================================================
#> 
#> Jackknife finished - 100% of replications converged.

# Compare standard errors across methods
sterrors <- data.frame(
  oim = sqrt(diag(v_oim)),
  opg = sqrt(diag(v_opg)),
  robust = sqrt(diag(v_robust)),
  cluster = sqrt(diag(v_clust)),
  bootstrap = sqrt(diag(v_boot)),
  jackknife = sqrt(diag(v_jack))
)

sterrors
#>                            oim         opg      robust    cluster   bootstrap
#> value::(Intercept) 8.636925227 9.281669188 8.411271366 7.23859394 7.413482653
#> value::age         0.399627579 0.431781719 0.382284787 0.32571570 0.341291878
#> value::I(age^2)    0.004609191 0.005002672 0.004391868 0.00383376 0.003966148
#> value::huswage     0.073304149 0.044304744 0.135952916 0.14449057 0.150892743
#> value::educ        0.135599149 0.116741374 0.165176906 0.14021038 0.151153737
#> value::unem        0.095744322 0.096801039 0.098676875 0.08041205 0.105188018
#> scale::lnsigma     0.025768404 0.014396434 0.055550631 0.05915267 0.052628926
#>                      jackknife
#> value::(Intercept) 8.507232115
#> value::age         0.386935880
#> value::I(age^2)    0.004446869
#> value::huswage     0.144060112
#> value::educ        0.167837527
#> value::unem        0.100966601
#> scale::lnsigma     0.057408605
```
