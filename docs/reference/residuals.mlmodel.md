# Extract Model Residuals

Extract Model Residuals

## Usage

``` r
# S3 method for class 'mlmodel'
residuals(object, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  An `mlmodel` object.

- type:

  Character string. Type of residuals to return. Currently supported:
  `"response"` (default) or `"pearson"`.

- ...:

  Further arguments passed to methods (currently not used).

## Value

A numeric vector of residuals aligned to the original data frame.
Observations dropped during estimation (due to `NA`s or `subset`) return
`NA`.

## Details

`"response"` residuals are the raw residuals: observed minus fitted
values.

`"pearson"` residuals are standardized by the model-implied standard
deviation: \\(y - \hat{y}) / \sqrt{\text{Var}(y)}\\. For Poisson models
they use \\\sqrt{\hat{\mu}}\\, for binary models
\\\sqrt{\hat{p}(1-\hat{p})}\\, and for other models the appropriate
variance from `predict(object, type = "var")` or `type = "var_y"`.
