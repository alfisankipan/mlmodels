# Hessian by Observation

Extract the per-observation Hessian matrices evaluated at the estimated
parameters from an `mlmodel` object.

## Usage

``` r
hessianObs(object)

# S3 method for class 'mlmodel'
hessianObs(object)
```

## Arguments

- object:

  An `mlmodel` object.

## Value

A numeric matrix of dimension `(N*K) x K`, where `N` is the number of
observations and `K` is the number of parameters. The Hessian for each
observation is stacked vertically.

## Details

This is mainly intended for advanced use (e.g., custom diagnostics or
information matrix tests). For most users, the functions `IMtest` or
`vcov` are more convenient.
