# Vuong's Test for Non-Nested Models

Performs Vuong's (1989) test for comparing two non-nested models fitted
via maximum likelihood with the `mlmodels` package.

## Usage

``` r
vuongtest(object_1, object_2, ...)

# S3 method for class 'mlmodel'
vuongtest(object_1, object_2, ...)
```

## Arguments

- object_1:

  A fitted model object inheriting from `"mlmodel"`.

- object_2:

  A fitted model object inheriting from `"mlmodel"`.

- ...:

  Further arguments passed to methods (currently not used).

## Value

An object of class `"vuongtest.mlmodel"` containing the test statistic,
p-value, and a conclusion (which model is preferred or "inconclusive").

## Details

Vuong's test compares two non-nested models by testing the null
hypothesis that the two models are equally close to the true data
generating process.

The test statistic is based on the difference in the per-observation
log-likelihood contributions between the two models. A positive
significant value favors `object_1`, a negative significant value favors
`object_2`, and a non-significant value leads to an "inconclusive"
result.

Both models must be estimated on exactly the same sample.

## References

Vuong, Q. H. (1989). 'Likelihood Ratio Tests for Model Selection and
Non-Nested Hypotheses.' *Econometrica*, 57(2), 307-333.
[doi:10.2307/1912557](https://doi.org/10.2307/1912557)

## See also

[`lrtest()`](https://alfisankipan.github.io/mlmodels/reference/lrtest.md),
[`waldtest()`](https://alfisankipan.github.io/mlmodels/reference/waldtest.md),
[`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md)

## Examples

``` r

# Linear models example (lognormal vs gamma)
data(mroz)
mroz$incthou <- mroz$faminc / 1000

fit_lognormal <- ml_lm(log(incthou) ~ age + I(age^2) + huswage + educ + unem,
                       data = mroz)

fit_gamma <- ml_gamma(incthou ~ age + I(age^2) + huswage + educ + unem,
                      data = mroz)


vuongtest(fit_lognormal, fit_gamma)
#> 
#> Vuong's (1989) Test
#> --------------------------------------------------
#>   Model 1: Homoskedastic Lognormal Model 
#>   Model 2: Homoskedastic Gamma Model 
#> --------------------------------------------------
#>   z-stat:  -1.592
#>   p-value: 0.1114
#> --------------------------------------------------
#>  Inconclusive test: neither model is preferred.

# Count models example

fit_poi <- ml_poisson(docvis ~ private + medicaid + age + I(age^2) + educyr +
                          actlim + totchr,
                      data = docvis)

fit_nb1 <- ml_negbin(docvis ~ private + medicaid + age + I(age^2) + educyr +
                          actlim + totchr,
                      data = docvis,
                      dispersion = "NB1")

fit_nb2 <- ml_negbin(docvis ~ private + medicaid + age + I(age^2) + educyr +
                          actlim + totchr,
                      data = docvis)
                      
# Poisson vs. NB1
vuongtest(fit_poi, fit_nb1)
#> 
#> Vuong's (1989) Test
#> --------------------------------------------------
#>   Model 1: Poisson 
#>   Model 2: Homoskedastic Negative Binomial (NB1) Model 
#> --------------------------------------------------
#>   z-stat:  -12.512
#>   p-value: 0.0000
#> --------------------------------------------------
#>   Homoskedastic Negative Binomial (NB1) Model seems to be preferred.

# NB1 vs. NB2
vuongtest(fit_nb1, fit_nb2)
#> 
#> Vuong's (1989) Test
#> --------------------------------------------------
#>   Model 1: Homoskedastic Negative Binomial (NB1) Model 
#>   Model 2: Homoskedastic Negative Binomial (NB2) Model 
#> --------------------------------------------------
#>   z-stat:  3.638
#>   p-value: 0.0003
#> --------------------------------------------------
#>   Homoskedastic Negative Binomial (NB1) Model seems to be preferred.

# Binary models example
data(smoke)
smoke$smokes <- smoke$cigs > 0

fit_logit <- ml_logit(smokes ~ cigpric + income + age, data = smoke)
fit_probit <- ml_probit(smokes ~ cigpric + income + age, data = smoke)

vuongtest(fit_logit, fit_probit)
#> 
#> Vuong's (1989) Test
#> --------------------------------------------------
#>   Model 1: Homoskedastic Binary Logit 
#>   Model 2: Homoskedastic Binary Probit 
#> --------------------------------------------------
#>   z-stat:  -1.158
#>   p-value: 0.2467
#> --------------------------------------------------
#>  Inconclusive test: neither model is preferred.
```
