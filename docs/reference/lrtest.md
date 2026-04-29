# Likelihood Ratio Test for Nested mlmodel Objects

Performs a likelihood ratio test comparing two nested models fitted with
the same estimator (e.g. `ml_lm`, `ml_logit`, `ml_negbin`, etc.).

## Usage

``` r
lrtest(object_1, object_2, ...)

# S3 method for class 'mlmodel'
lrtest(object_1, object_2, ...)
```

## Arguments

- object_1:

  A fitted model object inheriting from `"mlmodel"`. Typically the
  restricted (smaller) model.

- object_2:

  A fitted model object inheriting from `"mlmodel"`. Typically the
  unrestricted (larger) model. The order of `object_1` and `object_2`
  does not matter — the function automatically determines which is the
  restricted model.

- ...:

  Further arguments passed to methods (currently not used).

## Value

An object of class `"lrtest.mlmodel"` with the test statistic, degrees
of freedom, and p-value.

## Details

The likelihood ratio test statistic is calculated as: \$\$LR = 2 \times
(\log L\_{\text{unrestricted}} - \log L\_{\text{restricted}})\$\$

Under the null hypothesis that the restricted model is correct, `LR`
follows a \\\chi^2\\ distribution with degrees of freedom equal to the
difference in the number of parameters between the two models.

**Important:** The two models must be nested (the restricted model must
be a special case of the unrestricted one) and fitted on exactly the
same sample. The restricted model must have a lower (or equal)
log-likelihood.

## See also

[`waldtest()`](https://alfisankipan.github.io/mlmodels/reference/waldtest.md),
[`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md),
[`vuongtest()`](https://alfisankipan.github.io/mlmodels/reference/vuongtest.md)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Linear model example
data(mroz)
mroz$incthou <- mroz$faminc / 1000

fit_small <- ml_lm(incthou ~ age + huswage, data = mroz)
fit_large <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
                   data = mroz)

lrtest(fit_small, fit_large)
#> ℹ `object_1` is the restricted model (nested in `object_2`).
#> Likelihood Ratio Test
#> --------------------------------------------
#> Chisq(3) = 61.194    Pr(>Chisq) = < 1e-8
#> --------------------------------------------
#> LogLik (restricted)   : -2669.275   (df = 4)
#> LogLik (unrestricted) : -2638.678   (df = 7)

# You can also reverse the order — the function detects the restricted model
lrtest(fit_large, fit_small)
#> ℹ `object_2` is the restricted model (nested in `object_1`).
#> Likelihood Ratio Test
#> --------------------------------------------
#> Chisq(3) = 61.194    Pr(>Chisq) = < 1e-8
#> --------------------------------------------
#> LogLik (restricted)   : -2669.275   (df = 7)
#> LogLik (unrestricted) : -2638.678   (df = 4)
```
