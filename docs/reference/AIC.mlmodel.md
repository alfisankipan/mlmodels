# Extract AIC from mlmodel objects

Extract AIC from mlmodel objects

## Usage

``` r
# S3 method for class 'mlmodel'
AIC(object, ..., k = 2)

# S3 method for class 'summary.mlmodel'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  An object of class `"mlmodel"` or `"summary.mlmodel"`.

- ...:

  Further arguments passed to methods.

- k:

  Numeric. The penalty per parameter. Default is `k = 2` (standard AIC).
  See [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) for details.

## Value

A numeric value with the AIC.

## Details

For `mlmodel` objects, AIC is computed as
`-2 * logLik(object) + k * npar`.

For `summary.mlmodel` objects, the pre-computed AIC (with `k = 2`) is
returned; the `k` argument is accepted for compatibility but ignored.
