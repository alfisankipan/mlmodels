# Extract BIC from mlmodel objects

Extract BIC from mlmodel objects

## Usage

``` r
# S3 method for class 'mlmodel'
BIC(object, ...)

# S3 method for class 'summary.mlmodel'
BIC(object, ...)
```

## Arguments

- object:

  An object of class `"mlmodel"` or `"summary.mlmodel"`.

- ...:

  Further arguments passed to methods.

## Value

A numeric value with the BIC.

## Details

BIC is computed as `-2 * logLik(object) + log(nobs) * npar`.
