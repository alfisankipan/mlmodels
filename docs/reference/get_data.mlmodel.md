# Extract data used to fit the model (for insight/marginaleffects compatibility)

Extract data used to fit the model (for insight/marginaleffects
compatibility)

## Usage

``` r
# S3 method for class 'mlmodel'
get_data(x, ...)

get_modeldata.mlmodel(x, ...)
```

## Arguments

- x:

  An object of class `"mlmodel"`

- ...:

  Further arguments (currently ignored)

## Value

The original data frame used when fitting the model
