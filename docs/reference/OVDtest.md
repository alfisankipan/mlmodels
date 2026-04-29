# Overdispersion Tests for Count Models

Performs Cameron and Trivedi's (1990) regression-based tests for
overdispersion in count models.

## Usage

``` r
OVDtest(object)
```

## Arguments

- object:

  An object of class `"mlmodel.count"` (fitted with
  [`ml_poisson()`](https://alfisankipan.github.io/mlmodels/reference/ml_poisson.md)
  or
  [`ml_negbin()`](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md)).

## Value

A list containing the results of the tests against NB1 and NB2
alternatives, with coefficient estimates, t-statistics, and p-values.

## Details

These tests evaluate the null hypothesis that the conditional variance
equals the conditional mean (the Poisson assumption). Rejection
indicates overdispersion and suggests that a negative binomial model may
be more appropriate.

When the input object is not a Poisson model, a Poisson regression is
fitted internally using the value (mean) equation specification from
`object` in order to perform the test.

## References

Cameron, A. C., & Trivedi, P. K. (1990). 'Regression-based tests for
overdispersion in the Poisson model.' *Journal of Econometrics*, 46(3),
347-364.
[doi:10.1016/0304-4076(90)90014-K](https://doi.org/10.1016/0304-4076%2890%2990014-K)

Cameron, A. C., & Trivedi, P. K. (2013). *Regression Analysis of Count
Data* (2nd ed.). Cambridge University Press.
[doi:10.1017/CBO9781139013567](https://doi.org/10.1017/CBO9781139013567)

## See also

[`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md),
[`GOFtest()`](https://alfisankipan.github.io/mlmodels/reference/GOFtest.md)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Poisson model
fit_pois <- ml_poisson(docvis ~ private + medicaid + age + I(age^2) + 
                       educyr + actlim + totchr, data = docvis)
OVDtest(fit_pois)
#> 
#> Cameron and Trivedi (1990) Overdispersion Test:
#> --------------------------------------
#>   H0: Poisson (alpha = 0)
#>   H1: Overdispersion (alpha > 0)
#> --------------------------------------
#>     Estimate t-stat p-value    
#> NB2   0.7047 6.8029  0.0000 ***
#> NB1   5.3043 6.8383  0.0000 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Observations: 3677 
#> Note: P-values are based on a one-tailed t-test (Right Tail).
#> 

# Negative binomial model (the test still fits a Poisson internally using 
# only the value equation, so results are identical)
fit_nb2 <- ml_negbin(docvis ~ private + medicaid + age + I(age^2) + 
                     educyr + actlim + totchr, data = docvis)
OVDtest(fit_nb2)
#> 
#> Cameron and Trivedi (1990) Overdispersion Test:
#> --------------------------------------
#>   H0: Poisson (alpha = 0)
#>   H1: Overdispersion (alpha > 0)
#> --------------------------------------
#>     Estimate t-stat p-value    
#> NB2   0.7047 6.8029  0.0000 ***
#> NB1   5.3043 6.8383  0.0000 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Observations: 3677 
#> Note: P-values are based on a one-tailed t-test (Right Tail).
#> 
```
