# Update an mlmodel Call

Update an mlmodel Call

## Usage

``` r
# S3 method for class 'mlmodel'
update(
  object,
  formula. = NULL,
  scale. = NULL,
  data = NULL,
  weights = NULL,
  subset = NULL,
  ...,
  evaluate = TRUE
)

# S3 method for class 'ml_poisson'
update(
  object,
  formula. = NULL,
  data = NULL,
  weights = NULL,
  subset = NULL,
  ...,
  evaluate = TRUE
)
```

## Arguments

- object:

  An `mlmodel` object.

- formula.:

  An updated formula for the value (location/mean) equation.

- scale.:

  An updated formula for the scale equation (if supported by the model).

- data:

  A data frame to be used when re-fitting the model.

- weights:

  Optional case weights.

- subset:

  An expression or logical vector to subset, or a vector of indices to
  resample (`sandwich` uses it for that).

- ...:

  Further arguments passed to methods (currently ignored).

- evaluate:

  Logical. If `TRUE` (default), the updated call is evaluated and the
  new fitted model is returned. If `FALSE`, the updated call (as a
  language object) is returned without evaluation.

## Details

This method re-evaluates the original model call after modifying
selected arguments. It is used internally by
[`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md)
and serves as a fallback mechanism in bootstrap and jackknife variance
estimation when a model-specific implementation is not available.

**`sandwich` package compatibility**

The functions
[`sandwich::vcovBS()`](https://sandwich.R-Forge.R-project.org/reference/vcovBS.html)
and
[`sandwich::vcovJK()`](https://sandwich.R-Forge.R-project.org/reference/vcovJK.html)
are supported through this
[`update()`](https://rdrr.io/r/stats/update.html) method. They produce
numerically equivalent results to our own `vcov(object, type = "boot")`
and `vcov(object, type = "jack")` when all bootstrap/jackknife
replications converge, taking longer to compute them.

**Important difference**: When some replications fail to converge,
`sandwich` includes those failed iterations in the variance calculation,
while our [`vcov()`](https://rdrr.io/r/stats/vcov.html) implementation
uses **only successful replications**. The latter is statistically more
appropriate.

We therefore strongly recommend using the native
[`vcov()`](https://rdrr.io/r/stats/vcov.html) methods provided by
**mlmodels** for bootstrap and jackknife variance-covariance matrices.
