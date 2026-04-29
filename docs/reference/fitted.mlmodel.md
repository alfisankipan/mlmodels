# Extract Fitted Values from mlmodel

Extract Fitted Values from mlmodel

## Usage

``` r
# S3 method for class 'mlmodel'
fitted(object, ...)

# S3 method for class 'values.mlmodel'
fitted(object, ...)
```

## Arguments

- object:

  An `mlmodel` object.

- ...:

  Further arguments passed to methods (currently ignored).

## Value

A numeric vector of fitted values, aligned to the original data. Dropped
observations (due to `NA`s or `subset`) return `NA`.
