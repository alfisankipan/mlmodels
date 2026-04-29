# Extract Log-Likelihood from mlmodel objects

Extract Log-Likelihood from mlmodel objects

## Usage

``` r
# S3 method for class 'mlmodel'
logLik(object, ...)

# S3 method for class 'summary.mlmodel'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `mlmodel` or `summary.mlmodel`.

- ...:

  Additional arguments passed to methods.

## Value

An object of class `"logLik"` with the log-likelihood value and the
attributes `nobs` and `df`.

## Details

The returned object is of class `"logLik"` and has two important
attributes:

- `nobs`: number of observations used in estimation.

- `df`: number of estimated parameters (usually called *K*), computed as
  `length(coef(object))`. This includes coefficients from both the
  location (mean/value) and scale equations when present.
