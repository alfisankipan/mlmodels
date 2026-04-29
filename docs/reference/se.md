# Extract Standard Errors from mlmodel Objects

Extract Standard Errors from mlmodel Objects

Extracts standard errors for an `mlmodel` object.

## Usage

``` r
se(object, ...)

# S3 method for class 'mlmodel'
se(
  object,
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

  Further arguments passed to methods (currently not used).

- vcov:

  An optional user-supplied variance-covariance matrix.

- vcov.type:

  Character string specifying the type of variance-covariance matrix to
  use. One of `"oim"` (default), `"opg"`, `"robust"`, `"boot"`, or
  `"jack"`. See
  [`vcov.mlmodel()`](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md)
  for details.

- cl_var:

  Character string or vector. Name of the clustering variable (or the
  vector itself) when `vcov.type = "robust"` (or its alias `"cluster"`).

- repetitions:

  Integer. Number of bootstrap replications when `vcov.type = "boot"`.
  Default is 999.

- seed:

  Integer. Random seed for reproducibility when bootstrapping. If
  `NULL`, a random seed is generated internally.

- progress:

  Logical. Should a progress bar be shown during bootstrapping or
  jackknifing? Default is `FALSE`.

## Value

A named numeric vector of standard errors, with the same names as
`coef(object)`.

## See also

[vcov.mlmodel](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md),
[summary.mlmodel](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md),
[confint.mlmodel](https://alfisankipan.github.io/mlmodels/reference/confint.mlmodel.md)
