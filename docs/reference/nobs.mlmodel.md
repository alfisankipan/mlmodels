# Extract the Number of Observations from an mlmodel

Extract the Number of Observations from an mlmodel

## Usage

``` r
# S3 method for class 'mlmodel'
nobs(object, ...)
```

## Arguments

- object:

  An object of class `"mlmodel"`.

- ...:

  Further arguments passed to methods (currently not used).

## Value

An integer giving the number of observations used in the estimation
(after removing missing values and applying any `subset`).
