# Goodness-of-Fit Test for Count Models

Performs the Manjon and Martinez (2014) chi-squared goodness-of-fit test
for count data models.

## Usage

``` r
GOFtest(object, bins = 0:5)

# S3 method for class 'mlmodel'
GOFtest(object, bins = 0:5)
```

## Arguments

- object:

  An object of class `"mlmodel.count"` (typically from
  [`ml_poisson()`](https://alfisankipan.github.io/mlmodels/reference/ml_poisson.md)
  or
  [`ml_negbin()`](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md)).

- bins:

  Integer vector. Defines the boundaries of the bins used to group
  counts. Default is `0:5`.

## Value

An object of class `"GOFtest.mlmodel"` with components:

- model:

  Description of the fitted model.

- matrix:

  A table with observed and predicted frequencies, proportions, absolute
  differences, and Pearson contributions per bin.

- test:

  A list containing `teststat`, `df`, and `pval` for the overall
  goodness-of-fit test.

## Details

The test compares the observed frequencies with the expected frequencies
predicted by the model across different count bins. It produces both a
binned comparison table and an overall regression-based chi-squared test
statistic.

A low p-value indicates that the model's predicted probabilities do not
adequately match the observed count distribution (model
misspecification).

## References

Manjon, M., & Martinez, O. (2014). 'The chi-squared goodness-of-fit test
for count-data models.' *The Stata Journal*, 14(4), 798-816.
[doi:10.1177/1536867X1401400406](https://doi.org/10.1177/1536867X1401400406)

## Author

Alfonso Sanchez-Penalver

## Examples

``` r

# Poisson model
fit_pois <- ml_poisson(docvis ~ private + medicaid + age + I(age^2) + 
                       educyr + actlim + totchr, data = docvis)

GOFtest(fit_pois, bins = 0:5)
#> 
#> Goodness-of-fit test for count models
#>    Model: Poisson
#> --------------------------------------------------
#> Manjon & Martinez (2014) Score Test
#> 
#>       Frequency Proportion Probability |Difference|   Pearson
#> 0 - 0       401     0.1091      0.0074       0.1017 5168.2331
#> 1 - 1       314     0.0854      0.0296       0.0558  387.8678
#> 2 - 2       358     0.0974      0.0630       0.0344   69.0000
#> 3 - 3       334     0.0908      0.0954       0.0045    0.7889
#> 4 - 4       339     0.0922      0.1159       0.0237   17.8613
#> 5 - 5       266     0.0723      0.1212       0.0489   72.4408
#> 6 +        1665     0.4528      0.5676       0.1148   85.3671
#> 
#> --------------------------------------------------
#>   Chisq(6):             1097.8023
#>   p-value:               0.0000
#> --------------------------------------------------

# Negative binomial model
fit_nb2 <- ml_negbin(docvis ~ private + medicaid + age + I(age^2) + 
                     educyr + actlim + totchr, data = docvis)

GOFtest(fit_nb2)
#> 
#> Goodness-of-fit test for count models
#>    Model: Homoskedastic Negative Binomial (NB2) Model
#> --------------------------------------------------
#> Manjon & Martinez (2014) Score Test
#> 
#>       Frequency Proportion Probability |Difference| Pearson
#> 0 - 0       401     0.1091      0.0913       0.0178 12.7078
#> 1 - 1       314     0.0854      0.1079       0.0225 17.2884
#> 2 - 2       358     0.0974      0.1054       0.0081  2.2695
#> 3 - 3       334     0.0908      0.0962       0.0053  1.0863
#> 4 - 4       339     0.0922      0.0849       0.0073  2.3328
#> 5 - 5       266     0.0723      0.0735       0.0012  0.0720
#> 6 +        1665     0.4528      0.4408       0.0120  1.2058
#> 
#> --------------------------------------------------
#>   Chisq(6):             41.8194
#>   p-value:               0.0000
#> --------------------------------------------------
```
