# Log-Likelihood by Observation

Extract the per-observation log-likelihood contributions from an
`mlmodel` object.

## Usage

``` r
loglikeObs(object)

# S3 method for class 'mlmodel'
loglikeObs(object)
```

## Arguments

- object:

  An `mlmodel` object.

## Value

A numeric vector of length `nobs(object)` containing the log-likelihood
contribution of each observation.

## Details

These individual contributions are useful for Vuong tests, robust
variance estimation, or custom model diagnostics.
