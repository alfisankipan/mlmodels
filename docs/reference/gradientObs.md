# Gradient (Score) by Observation

Extract the per-observation gradients (scores) evaluated at the
estimated parameters from an `mlmodel` object.

## Usage

``` r
gradientObs(object)

# S3 method for class 'mlmodel'
gradientObs(object)
```

## Arguments

- object:

  An `mlmodel` object.

## Value

A numeric matrix with one row per observation and one column per
parameter. Each row contains the gradient of the log-likelihood for that
observation.

## Details

These are the individual contributions to the score vector. They are
mainly useful for advanced users who want to implement custom tests or
diagnostics.
