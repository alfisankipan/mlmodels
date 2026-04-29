# Predictions for mlmodel models

Methods for computing predictions from models fitted with the `mlmodels`
package.

## Usage

``` r
# S3 method for class 'ml_beta'
predict(
  object,
  newdata = NULL,
  type = "response",
  se.fit = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_gamma'
predict(
  object,
  newdata = NULL,
  type = "response",
  se.fit = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_lm'
predict(
  object,
  newdata = NULL,
  type = "response",
  se.fit = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_logit'
predict(
  object,
  newdata = NULL,
  type = "response",
  se.fit = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'mlmodel'
predict(object, ...)

# S3 method for class 'ml_negbin'
predict(
  object,
  newdata = NULL,
  type = "response",
  se.fit = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_poisson'
predict(
  object,
  newdata = NULL,
  type = "response",
  se.fit = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)

# S3 method for class 'ml_probit'
predict(
  object,
  newdata = NULL,
  type = "response",
  se.fit = FALSE,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- object:

  An object from an estimation with one of our models.

- newdata:

  Optional data frame for out-of-sample predictions.

- type:

  Character string indicating what to predict. See **Details**.

- se.fit:

  Logical. If `TRUE`, also return standard errors (delta method).

- vcov:

  Optional user-supplied variance-covariance matrix.

- vcov.type:

  Type of variance-covariance matrix. See
  [vcov](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md).

- cl_var:

  Clustering variable (name or vector).

- repetitions:

  Number of bootstrap replications when `vcov.type = "boot"`.

- seed:

  Random seed for bootstrapping, for reproducibility.

- progress:

  Logical. Show bootstrap/jackknife progress bar? Default is `FALSE` in
  higher-level functions.

- ...:

  Additional arguments passed to methods.

## Value

An object that inherits from `predict.mlmodel` and has two elements:

- fit:

  Vector with the predictions.

- se.fit:

  If `se.fit` is `TRUE` a vector with the delta-method standard errors,
  using analytical gradients. If `se.fit` is `FALSE`, it is set to
  `NULL`.

## Details

### ml_beta prediction types

The `type` argument controls what quantity is returned.

|  |  |  |
|----|----|----|
| Type | Description | Notes |
| `"link"` | Linear mean predictor ( xb ) | logit-mean |
| `"response"` | Expected proportion (outcome) | Default |
| `"mean"` | Alias for `"response"` | \- |
| `"fitted"` | Alias for `"response"` | \- |
| `"odds"` | Odds ratio | exp(xb) |
| `"zd"` | Linear precision predictor | log-phi |
| `"phi"` | Dispersion parameter | \- |
| `"shape1"` | Shape parameter of the beta distribution | mu \* phi |
| `"shape2"` | Shape parameter of the beta distribution | (1 - mu) \* phi |
| `"mode"` | Mode prediction (See below) | (shape1 - 1) / (shape1 + shape2 - 2) |
| `"variance"` | Variance of the outcome variable | mu \* (1 - mu) / (1 + phi) |
| `"var"` | Alias for `"variance"` | \- |
| `"sigma"` | Standard deviation of outcome variable | sqrt(`"variance"`) |
| `"sd"` | Alias for `"sigma"` | \- |

When `se.fit = TRUE`, standard errors are computed using the delta
method for all supported types.

**Mode Indeterminations**

The mode is only defined if `shape1 > 1` **and** `shape2 > 1` **and**
`shape1 + shape2 != 2`. If these conditions are not met the prediction
and standard error will be `NA`.

### ml_gamma prediction types

The `type` argument controls what quantity is returned.

|              |                                        |                    |
|--------------|----------------------------------------|--------------------|
| Type         | Description                            | Notes              |
| `"link"`     | Linear mean predictor ( xb )           | log-mean           |
| `"response"` | Expected outcome                       | Default            |
| `"mean"`     | Alias for `"response"`                 | \-                 |
| `"fitted"`   | Alias for `"response"`                 | \-                 |
| `"zd"`       | Linear shape predictor                 | log-nu             |
| `"nu"`       | Shape parameter                        | \-                 |
| `"variance"` | Variance of the outcome variable       | \-                 |
| `"var"`      | Alias for `"variance"`                 | \-                 |
| `"sigma"`    | Standard deviation of outcome variable | sqrt(`"variance"`) |
| `"sd"`       | Alias for `"sigma"`                    | \-                 |

When `se.fit = TRUE`, standard errors are computed using the delta
method for all supported types.

### ml_lm prediction types

The `type` argument controls what quantity is returned. Behavior differs
depending on whether the outcome was modeled in logs (`log(y)`).

|  |  |  |  |
|----|----|----|----|
| Type | Normal (linear) case | Lognormal case (`log(y)`) | Notes |
| `link` | Linear predictor for scale (zd) | Linear predictor on log scale (mu-log) | Scale equation |
| `fitted` | xb (mean predictor) | xb (original log-scale predictor) | Mean equation |
| `response`, `mean`, `mu` | xb (E`[y]`) | E`[y]` = exp(mu-log + sigma^2/2) - shift | Proper expected value on original scale |
| `median` | xb (same as mean) | exp(mu-log) - shift | Median of y |
| `sigma`, `sd` | sd of y | sd of `log(y)` | On log scale |
| `sigma_y`, `sd_y` | same as `sigma` | sd of y | Only meaningful in lognormal case |
| `variance`, `var` | sigma^2 | sigma^2 (variance of `log(y)`) | On log scale |
| `variance_y`, `var_y` | same as `variance` | Var(y) = exp(2 mu-log + sigma^2)(exp(sigma^2) - 1) | Only meaningful in lognormal case |
| `zd` | Linear predictor for scale (zd) | Linear predictor for scale (zd) | Alias for `link` |

When the outcome is log-transformed, `response` (or `mean`) returns the
correct lognormal expected value on the original scale of y. The
`median` is the simple exponential back-transform.

### ml_logit prediction types

The `type` argument controls what quantity is returned. Behavior differs
depending on whether the model is homoskedastic or heteroskedastic.

|  |  |  |  |
|----|----|----|----|
| Type | Homoskedastic case | Heteroskedastic case | Notes |
| `"xb"` | Linear predictor xb | Linear predictor xb | Linear predictor for value |
| `"response"` | P(y=1 \| x) | P(y=1 \| x) | Prob. of success (default) |
| `"prob"` | Alias for `"response"` | Alias for `"response"` | \- |
| `"fitted"` | Alias for `"response"` | Alias for `"response"` | \- |
| `"prob0"` | P(y=0 \| x) | P(y=0 \| x) | Prob. of failure |
| `"link"` | Linear predictor xb | xb / exp(zd) | Log-odds |
| `"odds"` | Odds = exp(xb) | Odds = exp(xb / exp(zd)) | \- |
| `"sigma"` | 1 (constant) | Std. Deviation: exp(zd) | Only available if heteroskedastic |
| `"variance"` | 1 (constant) | Variance: exp(2\*zd) | Only available if heteroskedastic |
| `"zd"` | 0 (constant) | Linear predictor zd | Linear predictor for scale |

In binary logit models, the **overall scale** of the latent error term
is not identified and is normalized to 1. In the homoskedastic case
there is no scale equation, so sigma is fixed at 1. In the
heteroskedastic case, the scale equation has no intercept. Therefore,
the predicted `"sigma"` and `"variance"` represent **individual-level
deviations** from the normalized overall scale, not the absolute
standard deviation or variance.

When `se.fit = TRUE`, standard errors are computed using the delta
method. Standard errors are not available (and will return `NA`) for
`"sigma"`, `"variance"`, and `"zd"` in homoskedastic models.

### ml_negbin prediction types

The `type` argument controls what quantity is returned. In addition to
standard types, Negative Binomial models support flexible probability
requests using the `P(...)` syntax.

|  |  |  |
|----|----|----|
| Type | Description | Notes |
| `"link"` | Linear mean predictor ( xb ) | log-mean |
| `"response"` | Expected count ( `mu` = `exp(xb)` ) | Default |
| `"mean"` | Alias for `"response"` | \- |
| `"fitted"` | Alias for `"response"` | \- |
| `"zd"` | Linear dispersion predictor | log-alpha |
| `"alpha"` | Dispersion parameter | \- |
| `"variance"` | Variance of the outcome variable | \- |
| `"var"` | Alias for `"variance"` | \- |
| `"sigma"` | Standard deviation of outcome variable | sqrt(`"variance"`) |
| `"sd"` | Alias for `"sigma"` | \- |
| `P(k)` | P(Y = k) | Exact probability, k integer \>= 0 |
| `P(,k)` | P(Y \<= k) | Cumulative (lower tail) |
| `P(k,)` | P(Y \>= k) | Survival (upper tail) |
| `P(a,b)` | P(a \<= Y \<= b) | Interval probability, a \<= b, a \>= 0 |

When `se.fit = TRUE`, standard errors are computed using the delta
method for all supported types.

### ml_poisson prediction types

The `type` argument controls what quantity is returned. In addition to
standard types, Poisson models support flexible probability requests
using the `P(...)` syntax.

|  |  |  |
|----|----|----|
| Type | Description | Notes |
| `"link"` | Linear predictor ( xb ) | log-mean |
| `"response"` | Expected count ( `mu` = `exp(xb)` ) | Default |
| `"mean"` | Alias for `"response"` | \- |
| `"mu"` | Alias for `"response"` | \- |
| `"fitted"` | Alias for `"response"` | \- |
| `P(k)` | P(Y = k) | Exact probability, k integer \>= 0 |
| `P(,k)` | P(Y \<= k) | Cumulative (lower tail) |
| `P(k,)` | P(Y \>= k) | Survival (upper tail) |
| `P(a,b)` | P(a \<= Y \<= b) | Interval probability, a \<= b, a \>= 0 |

When `se.fit = TRUE`, standard errors are computed using the delta
method for all supported types.

### ml_probit prediction types

The `type` argument controls what quantity is returned. Behavior differs
depending on whether the model is homoskedastic or heteroskedastic.

|  |  |  |  |
|----|----|----|----|
| Type | Homoskedastic case | Heteroskedastic case | Notes |
| `"xb"` | Linear predictor xb | Linear predictor xb | Linear predictor for value |
| `"response"` | P(y=1 \| x) | P(y=1 \| x) | Prob. of success (default) |
| `"prob"` | Alias for `"response"` | Alias for `"response"` | \- |
| `"fitted"` | Alias for `"response"` | Alias for `"response"` | \- |
| `"prob0"` | P(y=0 \| x) | P(y=0 \| x) | Prob. of failure |
| `"link"` | Linear predictor xb | xb / exp(zd) | Probit index |
| `"odds"` | Odds = prob / prob0 | Odds = prob / prob0. | \- |
| `"sigma"` | 1 (constant) | Std. Deviation: exp(zd) | Only available if heteroskedastic |
| `"variance"` | 1 (constant) | Variance: exp(2\*zd) | Only available if heteroskedastic |
| `"zd"` | 0 (constant) | Linear predictor zd | Linear predictor for scale |

In binary probit models, the **overall scale** of the latent error term
is not identified and is normalized to 1. In the homoskedastic case
there is no scale equation, so sigma is fixed at 1. In the
heteroskedastic case, the scale equation has no intercept. Therefore,
the predicted `"sigma"` and `"variance"` represent **individual-level
deviations** from the normalized overall scale, not the absolute
standard deviation or variance.

The `"link"` type returns the value on the **probit scale**, which is
the inverse of the standard normal cumulative distribution function (p =
Phi^(-1)(p)). This is the linear prediction (p = xb) ih homoskedastic
models, and the standardized linear predictor (p = xb / sigma) in
heteroskedastic models.

When `se.fit = TRUE`, standard errors are computed using the delta
method. Standard errors are not available (and will return `NA`) for
`"sigma"`, `"variance"`, and `"zd"` in homoskedastic models.

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Basic usage and different predict types
data(docvis)
fit_pois <- ml_poisson(docvis ~ age + educyr + totchr, data = docvis)

head(predict(fit_pois, type = "response")$fit)     # Expected count
#> [1] 10.169550  6.266430  6.795981  9.255944  5.471075  5.002632
head(predict(fit_pois, type = "P(3)")$fit)         # Prob of exactly 3
#> [1] 0.006716988 0.077881334 0.058498950 0.012627153 0.114817735 0.140226107

# Prediction at the mean (typical case)
typical <- data.frame(age = mean(docvis$age), 
                      educyr = mean(docvis$educyr), 
                      totchr = mean(docvis$totchr))
predict(fit_pois, newdata = typical, type = "response")
#> $fit
#> [1] 6.31286
#> 
#> $se.fit
#> NULL
#> 
#> attr(,"class")
#> [1] "predict.ml_poisson" "predict.mlmodel"   

# In-sample vs full-data prediction with subset / boundary dropping
data(pw401k)
fit_beta <- ml_beta(prate ~ mrate + I(mrate^2) + log(totemp) + 
                    I(log(totemp)^2) + age + I(age^2) + sole,
                    data = pw401k, 
                    subset = prate < 1)
#> ℹ Improving initial values by scaling (factor = 0.5).
#> ℹ Initial log-likelihood: -311.974
#> ℹ Final scaled log-likelihood: 79.945

# In-sample prediction (NAs for dropped observations)
head(predict(fit_beta, type = "response")$fit)
#> [1] 0.8034086 0.7705266        NA 0.7500891 0.7293741        NA

# Full-data prediction (predicts for all rows, including dropped ones)
head(predict(fit_beta, newdata = pw401k, type = "response")$fit)
#> [1] 0.8034086 0.7705266 0.8783125 0.7500891 0.7293741 0.8903587
```
