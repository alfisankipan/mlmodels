# Introduction to Count Data

``` r

library(mlmodels)
```

## Introduction

Count data — non-negative integers such as the number of doctor visits,
patents filed, or accidents — are common in economics, health, and
social sciences.

The `mlmodels` package provides a consistent, flexible framework for
Poisson and negative binomial models. A key advantage is the ability to
model heteroskedasticity in the dispersion parameter for both NB1 and
NB2 specifications.

To illustrate the functionality of our modeling functions for these type
of data, we take on reproducing the examples in chapter 20 of
*Microeconometrics Using Stata 2nd ed.* (Cameron and Trivedi, 2022),
while demonstrating additional capabilities, most notably
heteroskedastic NB1 estimation and unified post-estimation tools.

## The Data

The data come from a cross-sectional sample of the U.S. Medical
Expenditure Panel Survey (MEPS) for the year 2003. This dataset is
supplied with the book and is available on the [Stata Press
website](https://www.stata-press.com/data/mus2.html). We have included
it in the `mlmodels` package under the name `docvis`.

``` r

  data("docvis", package = "mlmodels")
```

## Poisson Model

The first model we estimate is the Poisson regression. All count models
in this vignette share the same mean specification, so we store the
formula once for convenience.

``` r

  mean_for <- docvis ~ private + medicaid + age + I(age^2) + educyr + actlim + totchr
  pois <- ml_poisson(mean_for, data = docvis)
```

The first estimation results in the chapter are the estimates of this
Poisson model using the default Original Information Matrix standard
errors.

``` r

  # oim standard errors are the default in summary
  summary(pois)
#> 
#> Maximum Likelihood Model
#>  Type: Poisson 
#> ---------------------------------------
#> Call:
#> ml_poisson(value = docvis ~ private + medicaid + age + I(age^2) + 
#>     educyr + actlim + totchr, data = docvis)
#> 
#> Log-Likelihood: -15019.64 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 4686.272, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#>                        Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept) -10.18221    0.97201 -10.475  < 2e-16 ***
#>   value::private       0.14223    0.01433   9.925  < 2e-16 ***
#>   value::medicaid      0.09700    0.01893   5.124 2.99e-07 ***
#>   value::age           0.29367    0.02596  11.314  < 2e-16 ***
#>   value::I(age^2)     -0.00193    0.00017 -11.199  < 2e-16 ***
#>   value::educyr        0.02956    0.00188  15.705  < 2e-16 ***
#>   value::actlim        0.18642    0.01457  12.798  < 2e-16 ***
#>   value::totchr        0.24839    0.00464  53.478  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3669
#> Pseudo R-squared - Cor.Sq.: 0.1531 McFadden: 0.1297
#> AIC: 30055.28  BIC: 30104.96 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 6.3089 
#>   Zeros - Observed: 401 Predicted: 27.06
```

A useful set of statistics appears under the header **Count
Diagnostics**. These include the Pearson index of dispersion and a
comparison between the actual and predicted number of zeros.

For a model that correctly captures dispersion, the Pearson index should
be close to 1. Here we see a value of approximately 6.3, which clearly
indicates that the Poisson model is not handling the overdispersion
present in the data.

The second statistic evaluates the model’s ability to predict the number
of zeros. The Poisson model severely underpredicts the number of zeros
(401 observed vs. only about 27 predicted), highlighting its poor fit in
this dimension as well.

Notice also the two pseudo R-squared measures reported by
[`summary()`](https://rdrr.io/r/base/summary.html):

- **Cor.Sq.** (squared correlation between predicted and observed
  values) is provided consistently across all models in `mlmodels`. This
  is the only R-squared-like measure that is directly comparable between
  linear and non-linear models.

- **McFadden** is the more traditional pseudo R-squared for discrete
  choice and count models.

Cameron and Trivedi had to calculate the correlated squared R-squared by
hand. We provide it automatically for every model.

Because the Poisson model shows clear signs of overdispersion, Cameron
and Trivedi recommend using robust standard errors with Poisson
regressions.

In `mlmodels` you do not need to re-estimate the model. Simply request
the desired variance type when calling
[`summary()`](https://rdrr.io/r/base/summary.html):

``` r

summary(pois, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Poisson 
#> ---------------------------------------
#> Call:
#> ml_poisson(value = docvis ~ private + medicaid + age + I(age^2) + 
#>     educyr + actlim + totchr, data = docvis)
#> 
#> Log-Likelihood: -15019.64 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 720.431, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                        Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept) -10.18221    2.36921  -4.298 1.73e-05 ***
#>   value::private       0.14223    0.03636   3.912 9.15e-05 ***
#>   value::medicaid      0.09700    0.05683   1.707   0.0878 .  
#>   value::age           0.29367    0.06298   4.663 3.11e-06 ***
#>   value::I(age^2)     -0.00193    0.00042  -4.636 3.55e-06 ***
#>   value::educyr        0.02956    0.00485   6.100 1.06e-09 ***
#>   value::actlim        0.18642    0.03966   4.701 2.59e-06 ***
#>   value::totchr        0.24839    0.01258  19.747  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3669
#> Pseudo R-squared - Cor.Sq.: 0.1531 McFadden: 0.1297
#> AIC: 30055.28  BIC: 30104.96 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 6.3089 
#>   Zeros - Observed: 401 Predicted: 27.06
```

Notice that switching to robust standard errors not only changes the
standard errors, z-statistics, and p-values of the individual
coefficients, but also updates the Wald joint significance tests at the
top of the table to use the robust variance-covariance matrix, or
whichever variance you select to use.

This makes it extremely easy to compare different variance-covariance
estimators (Original Information Matrix, robust, bootstrap, jackknife,
clustered, etc.) without refitting the model. See the vignette on
variance types for more details.

## Test of Overdispersion

One of the key concerns with the Poisson model is overdispersion.
Cameron and Trivedi present a regression-based test of overdispersion
(against the NB2 parameterization) and provide code to implement it.

In `mlmodels` we have implemented this test in a convenient function
that compares the Poisson model against **both** NB1 and NB2 forms of
overdispersion:

``` r

  OVDtest(pois)
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
```

The output shows strong evidence of overdispersion in both
specifications (highly significant t-statistics for both NB1 and NB2).

This confirms two important points:

- If you are only interested in inference about the **mean** parameters,
  using robust standard errors with the Poisson model is usually
  sufficient.
- However, if you also care about correctly modeling the **variance** or
  making accurate **probability predictions** (e.g., P(Y=0), P(Y=k),
  etc.), the Poisson model is inadequate. In such cases, moving to a
  negative binomial model is more appropriate.

## Average Marginal Effects

To produce average marginal effects, we have made `mlmodels` compatible
with the leading package for these estimates:
[`marginaleffects`](https://marginaleffects.com/).

Thanks to this integration, you can now compute marginal effects and
average marginal effects directly on any model estimated with
`mlmodels`. For example:

``` r

  library(marginaleffects)
  avg_slopes(pois, vcov = "robust")
#> 
#>      Term Contrast Estimate Std. Error     z Pr(>|z|)     S    2.5 % 97.5 %
#>  actlim      1 - 0   1.2959     0.2851  4.55   <0.001  17.5  0.73724 1.8546
#>  age         dY/dX   0.0386     0.0172  2.24   0.0249   5.3  0.00486 0.0723
#>  educyr      dY/dX   0.2017     0.0338  5.97   <0.001  28.6  0.13544 0.2679
#>  medicaid    1 - 0   0.6831     0.4153  1.64   0.1000   3.3 -0.13096 1.4971
#>  private     1 - 0   0.9702     0.2473  3.92   <0.001  13.5  0.48546 1.4549
#>  totchr      dY/dX   1.6947     0.0909 18.65   <0.001 255.3  1.51655 1.8728
#> 
#> Type: response
```

Even though the order of the variables differs from the book (they are
sorted alphabetically by `marginaleffects`), the estimates match
exactly.

You can also see how `marginaleffects` automatically identifies discrete
variables (showing contrast `1 - 0`) versus continuous variables
(showing contrast `dy/dx`).

For more information about predictions and our integration with
`marginaleffects`, see the dedicated vignette: **[Predictions and
marginaleffects
integration](https://alfisankipan.github.io/mlmodels/articles/mlmodels-predictions.md)**.

## Negative Binomial Estimation

Cameron and Trivedi focus primarily on the NB2 parameterization, which
is the most commonly used negative binomial model in the literature.

In `mlmodels` it is very easy to estimate both NB1 and NB2 models using
the same mean specification. This allows us to directly compare how well
each form of overdispersion fits the data.

``` r

# NB1 (linear variance)
nb1 <- ml_negbin(mean_for, data = docvis, dispersion = "NB1")
summary(nb1, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Negative Binomial (NB1) Model 
#> ---------------------------------------
#> Call:
#> ml_negbin(value = docvis ~ private + medicaid + age + I(age^2) + 
#>     educyr + actlim + totchr, dispersion = "NB1", data = docvis)
#> 
#> Log-Likelihood: -10531.05 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 936.865, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept) -9.30261    1.92387  -4.835 1.33e-06 ***
#>   value::private      0.15492    0.02928   5.291 1.22e-07 ***
#>   value::medicaid     0.06203    0.04120   1.505    0.132    
#>   value::age          0.27062    0.05137   5.268 1.38e-07 ***
#>   value::I(age^2)    -0.00177    0.00034  -5.179 2.23e-07 ***
#>   value::educyr       0.02499    0.00389   6.427 1.30e-10 ***
#>   value::actlim       0.13043    0.03109   4.194 2.73e-05 ***
#>   value::totchr       0.25150    0.01022  24.604  < 2e-16 ***
#> Scale (log(alpha)):  
#>   scale::lnalpha      1.45654    0.04542  32.069  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3668
#> Pseudo R-squared - Cor.Sq.: 0.1521 McFadden: 0.04054
#> AIC: 21080.11  BIC: 21136.00 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 1.2033 
#>   Zeros - Observed: 401 Predicted: 395.77

# NB2 (quadratic variance — default)
nb2 <- ml_negbin(mean_for, data = docvis)
summary(nb2, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Negative Binomial (NB2) Model 
#> ---------------------------------------
#> Call:
#> ml_negbin(value = docvis ~ private + medicaid + age + I(age^2) + 
#>     educyr + actlim + totchr, data = docvis)
#> 
#> Log-Likelihood: -10589.34 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 725.694, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                        Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept) -10.29748    2.42413  -4.248 2.16e-05 ***
#>   value::private       0.16409    0.03689   4.449 8.64e-06 ***
#>   value::medicaid      0.10034    0.05674   1.768    0.077 .  
#>   value::age           0.29413    0.06463   4.551 5.34e-06 ***
#>   value::I(age^2)     -0.00193    0.00043  -4.501 6.75e-06 ***
#>   value::educyr        0.02869    0.00492   5.831 5.51e-09 ***
#>   value::actlim        0.18954    0.03939   4.812 1.50e-06 ***
#>   value::totchr        0.27764    0.01324  20.964  < 2e-16 ***
#> Scale (log(alpha)):  
#>   scale::lnalpha      -0.44528    0.03780 -11.780  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3668
#> Pseudo R-squared - Cor.Sq.: 0.1498 McFadden: 0.03523
#> AIC: 21196.68  BIC: 21252.57 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 1.2544 
#>   Zeros - Observed: 401 Predicted: 335.69
```

## Predictions of the Mean: Poisson vs Negative Binomial

Before comparing the two negative binomial models, it is useful to
examine whether the Poisson model already provides reasonable
predictions for the **mean** (the quantity of primary interest in many
applications).

``` r

mpoi <- predict(pois, vcov.type = "robust", se.fit = TRUE)
mnb1 <- predict(nb1, vcov.type = "robust", se.fit = TRUE)
mnb2 <- predict(nb2, vcov.type = "robust", se.fit = TRUE)

comp <- data.frame(
  poi_mean = mpoi$fit,
  poi_se   = mpoi$se.fit,
  nb1_mean = mnb1$fit,
  nb1_se   = mnb1$se.fit,
  nb2_mean = mnb2$fit,
  nb2_se   = mnb2$se.fit
)

comp[1:10, ]
#>    poi_mean    poi_se nb1_mean    nb1_se nb2_mean    nb2_se
#> 1  8.657887 0.4280228 9.055702 0.3801113 9.025250 0.4516771
#> 2  6.103619 0.2330346 6.083600 0.2132264 6.089498 0.2452112
#> 3  6.253173 0.3403346 6.255303 0.2839657 6.218107 0.3274187
#> 4  9.513269 0.3336518 9.786868 0.2728545 9.860059 0.3446714
#> 5  5.973006 0.1883343 6.068180 0.1727429 5.828602 0.1917382
#> 6  4.438282 0.1576532 4.478689 0.1364806 4.218450 0.1501005
#> 7  3.136027 0.1574681 3.165791 0.1279083 2.887349 0.1446395
#> 8  8.280863 0.3279198 8.021330 0.2711635 8.323788 0.3305824
#> 9  6.300009 0.3214674 6.253927 0.2597699 6.235221 0.3138449
#> 10 5.363114 0.1913217 5.478609 0.1656597 5.257852 0.1847931
```

The individual predicted means and their standard errors are very close
across the three models.

We can also look at the average predicted value across the sample:

``` r

avg_predictions(pois, vcov = "robust")
#> 
#>  Estimate Std. Error    z Pr(>|z|)   S 2.5 % 97.5 %
#>      6.82      0.112 60.8   <0.001 Inf   6.6   7.04
#> 
#> Type: response
avg_predictions(nb1, vcov = "robust")
#> 
#>  Estimate Std. Error    z Pr(>|z|)   S 2.5 % 97.5 %
#>      6.82      0.112 60.7   <0.001 Inf   6.6   7.04
#> 
#> Type: response
avg_predictions(nb2, vcov = "robust")
#> 
#>  Estimate Std. Error    z Pr(>|z|)   S 2.5 % 97.5 %
#>      6.89      0.115 59.7   <0.001 Inf  6.66   7.12
#> 
#> Type: response
```

The average predictions (and their standard errors) are even closer to
each other.

Even the **average marginal effects** on the mean are very similar
across the three models:

``` r

avg_slopes(pois, vcov = "robust")
#> 
#>      Term Contrast Estimate Std. Error     z Pr(>|z|)     S    2.5 % 97.5 %
#>  actlim      1 - 0   1.2959     0.2851  4.55   <0.001  17.5  0.73724 1.8546
#>  age         dY/dX   0.0386     0.0172  2.24   0.0249   5.3  0.00486 0.0723
#>  educyr      dY/dX   0.2017     0.0338  5.97   <0.001  28.6  0.13544 0.2679
#>  medicaid    1 - 0   0.6831     0.4153  1.64   0.1000   3.3 -0.13096 1.4971
#>  private     1 - 0   0.9702     0.2473  3.92   <0.001  13.5  0.48546 1.4549
#>  totchr      dY/dX   1.6947     0.0909 18.65   <0.001 255.3  1.51655 1.8728
#> 
#> Type: response
avg_slopes(nb1, vcov = "robust")
#> 
#>      Term Contrast Estimate Std. Error     z Pr(>|z|)     S   2.5 % 97.5 %
#>  actlim      1 - 0   0.9020     0.2194  4.11   <0.001  14.6  0.4719 1.3320
#>  age         dY/dX   0.0457     0.0131  3.50   <0.001  11.1  0.0201 0.0713
#>  educyr      dY/dX   0.1705     0.0268  6.35   <0.001  32.1  0.1179 0.2231
#>  medicaid    1 - 0   0.4319     0.2933  1.47    0.141   2.8 -0.1430 1.0068
#>  private     1 - 0   1.0562     0.1996  5.29   <0.001  23.0  0.6650 1.4474
#>  totchr      dY/dX   1.7159     0.0756 22.68   <0.001 376.0  1.5677 1.8642
#> 
#> Type: response
avg_slopes(nb2, vcov = "robust")
#> 
#>      Term Contrast Estimate Std. Error     z Pr(>|z|)     S   2.5 % 97.5 %
#>  actlim      1 - 0   1.3288     0.2859  4.65  < 0.001  18.2  0.7685  1.889
#>  age         dY/dX   0.0438     0.0169  2.59  0.00963   6.7  0.0106  0.077
#>  educyr      dY/dX   0.1977     0.0346  5.72  < 0.001  26.4  0.1299  0.266
#>  medicaid    1 - 0   0.7142     0.4196  1.70  0.08872   3.5 -0.1081  1.537
#>  private     1 - 0   1.1304     0.2540  4.45  < 0.001  16.8  0.6325  1.628
#>  totchr      dY/dX   1.9130     0.0997 19.19  < 0.001 270.2  1.7176  2.108
#> 
#> Type: response
```

This illustrates the important practical point we made earlier: **if you
are primarily interested in inference about the conditional mean**, the
Poisson model with robust standard errors often performs quite well,
even when overdispersion is present.

The negative binomial models become more valuable when you also care
about correctly modeling the **variance** or making accurate
**probability predictions**, as we shall see.

## Model Comparison: NB1 versus NB2

Looking at the estimation results, the NB1 model appears to fit these
data better than the NB2 model. The log-likelihood is higher for NB1
(−10,531.05 versus −10,589.34), which also produces lower values of both
AIC and BIC.

In addition, the Pearson dispersion ratio is closer to 1 under NB1, the
predicted number of zeros is much closer to the observed count, and both
pseudo R-squared measures (Cor.Sq. and McFadden) are higher for the NB1
specification.

To formally test whether this difference in fit is statistically
significant, we use Vuong’s (1989) test for non-nested models:

``` r

vuongtest(nb1, nb2)
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
```

The test confirms that the NB1 model is preferred over the NB2 model for
this particular specification.

## Comparing Predicted Probabilities

Cameron and Trivedi compare predicted probabilities across models using
a couple of community-supplied packages: `countfit` (Long and Freese
2014) and `chi2gof` (Manjón and Martínez 2014).

We have implemented Manjón and Martínez’s approach, which allows
flexibility in the choice of bins and performs an asymptotically correct
goodness-of-fit test. In addition, we also report the **Pearson
contribution of each bin** — a feature present in `countfit` but not in
the original `chi2gof`. This extra information makes it much easier to
identify which specific bins are contributing most to any lack of fit.

``` r

# Default bins are 0 to 5
GOFtest(pois)
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
GOFtest(nb1)
#> 
#> Goodness-of-fit test for count models
#>    Model: Homoskedastic Negative Binomial (NB1) Model
#> --------------------------------------------------
#> Manjon & Martinez (2014) Score Test
#> 
#>       Frequency Proportion Probability |Difference| Pearson
#> 0 - 0       401     0.1091      0.1076       0.0014  0.0690
#> 1 - 1       314     0.0854      0.1024       0.0170 10.4096
#> 2 - 2       358     0.0974      0.0951       0.0023  0.2063
#> 3 - 3       334     0.0908      0.0865       0.0043  0.7862
#> 4 - 4       339     0.0922      0.0778       0.0144  9.8359
#> 5 - 5       266     0.0723      0.0693       0.0031  0.5047
#> 6 +        1665     0.4528      0.4613       0.0085  0.5773
#> 
#> --------------------------------------------------
#>   Chisq(6):             29.4554
#>   p-value:               0.0000
#> --------------------------------------------------
GOFtest(nb2)
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

A quick look at the results shows that all three models fail the
goodness-of-fit test at conventional significance levels. However, the
test statistic is smallest for the NB1 model, suggesting it fits the
data relatively better than the Poisson or NB2 models.

Beyond the overall test statistic, the table we provide makes it easy to
identify **which specific counts** each model struggles to predict.

- The **Poisson model** has difficulty predicting almost all counts,
  with the notable exception of 3. This pattern is typical when
  overdispersion is present and was one of the main reasons we
  recommended moving to a negative binomial specification.
- The **NB1 model** performs poorly mainly at counts of 1 and 4.
- The **NB2 model** struggles particularly with counts of 0 and 1.
  Although its Pearson contributions are generally higher than those of
  NB1 (indicating worse fit), it does slightly better at counts of 4 and
  5.

These detailed diagnostics help us understand not just *whether* a model
fits poorly, but *where* the misfit occurs.

## Heteroskedastic Negative Binomial Models

The last Maximum Likelihood specification in Cameron and Trivedi’s
chapter is a heteroskedastic NB2 model, where they model the natural log
of the dispersion parameter as a function of two binary variables:
`female` and `bh` (black and Hispanic).

In `mlmodels` we can easily estimate heteroskedastic versions of
**both** NB1 and NB2 models using the same syntax:

``` r

nb1_het <- ml_negbin(mean_for, 
                     scale = ~ female + bh, 
                     dispersion = "NB1", 
                     data = docvis)

nb2_het <- ml_negbin(mean_for, 
                     scale = ~ female + bh, 
                     data = docvis)

summary(nb1_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Negative Binomial (NB1) Model 
#> ---------------------------------------
#> Call:
#> ml_negbin(value = docvis ~ private + medicaid + age + I(age^2) + 
#>     educyr + actlim + totchr, scale = ~female + bh, dispersion = "NB1", 
#>     data = docvis)
#> 
#> Log-Likelihood: -10522.97 
#> 
#> Wald significance tests:
#>  all: Chisq(9) = 938.480, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(7) = 900.628, Pr(>Chisq) = < 1e-8
#>  Scale: Chisq(2) = 12.093, Pr(>Chisq) = 0.0024
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept) -9.40807    1.91860  -4.904 9.41e-07 ***
#>   value::private      0.14790    0.02951   5.012 5.40e-07 ***
#>   value::medicaid     0.06261    0.04132   1.515   0.1297    
#>   value::age          0.27514    0.05126   5.367 7.99e-08 ***
#>   value::I(age^2)    -0.00180    0.00034  -5.294 1.20e-07 ***
#>   value::educyr       0.02335    0.00388   6.018 1.77e-09 ***
#>   value::actlim       0.13000    0.03115   4.173 3.00e-05 ***
#>   value::totchr       0.24957    0.01033  24.151  < 2e-16 ***
#> Scale (log(alpha)):  
#>   scale::(Intercept)  1.48924    0.06011  24.773  < 2e-16 ***
#>   scale::female      -0.14900    0.07195  -2.071   0.0384 *  
#>   scale::bh           0.23478    0.09779   2.401   0.0164 *  
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3666
#> Pseudo R-squared - Cor.Sq.: 0.1522 McFadden: 0.04128
#> AIC: 21067.94  BIC: 21136.24 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 1.1863 
#>   Zeros - Observed: 401 Predicted: 396.78 
#> 
#> Distribution of Dispersion (alpha):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     3.8     3.8     4.4     4.3     4.8     5.6
summary(nb2_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Negative Binomial (NB2) Model 
#> ---------------------------------------
#> Call:
#> ml_negbin(value = docvis ~ private + medicaid + age + I(age^2) + 
#>     educyr + actlim + totchr, scale = ~female + bh, data = docvis)
#> 
#> Log-Likelihood: -10576.26 
#> 
#> Wald significance tests:
#>  all: Chisq(9) = 774.748, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(7) = 703.038, Pr(>Chisq) = < 1e-8
#>  Scale: Chisq(2) = 19.792, Pr(>Chisq) = 0.0001
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                        Estimate Std. Error z value Pr(>|z|)     
#> Value (docvis):  
#>   value::(Intercept) -10.54756    2.39369  -4.406 1.05e-05 ***
#>   value::private       0.15718    0.03674   4.278 1.88e-05 ***
#>   value::medicaid      0.08602    0.05394   1.595  0.11078    
#>   value::age           0.30188    0.06383   4.729 2.25e-06 ***
#>   value::I(age^2)     -0.00198    0.00042  -4.687 2.77e-06 ***
#>   value::educyr        0.02848    0.00493   5.771 7.87e-09 ***
#>   value::actlim        0.18754    0.03874   4.842 1.29e-06 ***
#>   value::totchr        0.27615    0.01328  20.800  < 2e-16 ***
#> Scale (log(alpha)):  
#>   scale::(Intercept)  -0.41191    0.05820  -7.077 1.47e-12 ***
#>   scale::female       -0.18719    0.07330  -2.554  0.01065 *  
#>   scale::bh            0.31031    0.09498   3.267  0.00109 ** 
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:3677 Deg. of freedom: 3666
#> Pseudo R-squared - Cor.Sq.: 0.1501 McFadden: 0.03642
#> AIC: 21174.52  BIC: 21242.83 
#> 
#> Count Diagnostics:
#>   Dispersion Ratio (Pearson): 1.2316 
#>   Zeros - Observed: 401 Predicted: 339.24 
#> 
#> Distribution of Dispersion (alpha):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    0.55    0.55    0.66    0.65    0.75    0.90
```

**Important note:** The dispersion parameters are not directly
comparable across NB1 and NB2 specifications because they represent
different functional forms. However, the statistical significance of the
coefficients in the scale equation is directly comparable.

We see that both `female` and `bh` are statistically significant (at
least at the 5% level) in the scale equations of both models. The Wald
tests for the scale parameters (shown at the top of the coefficient
tables) confirm that the two variables are jointly significant in
explaining dispersion in both specifications.

``` r

vuongtest(nb1_het, nb2_het)
#> 
#> Vuong's (1989) Test
#> --------------------------------------------------
#>   Model 1: Heteroskedastic Negative Binomial (NB1) Model 
#>   Model 2: Heteroskedastic Negative Binomial (NB2) Model 
#> --------------------------------------------------
#>   z-stat:  3.391
#>   p-value: 0.0007
#> --------------------------------------------------
#>   Heteroskedastic Negative Binomial (NB1) Model seems to be preferred.
GOFtest(nb1_het)
#> 
#> Goodness-of-fit test for count models
#>    Model: Heteroskedastic Negative Binomial (NB1) Model
#> --------------------------------------------------
#> Manjon & Martinez (2014) Score Test
#> 
#>       Frequency Proportion Probability |Difference| Pearson
#> 0 - 0       401     0.1091      0.1079       0.0011  0.0448
#> 1 - 1       314     0.0854      0.1015       0.0161  9.3910
#> 2 - 2       358     0.0974      0.0944       0.0030  0.3487
#> 3 - 3       334     0.0908      0.0862       0.0047  0.9323
#> 4 - 4       339     0.0922      0.0776       0.0146 10.0353
#> 5 - 5       266     0.0723      0.0693       0.0031  0.4953
#> 6 +        1665     0.4528      0.4631       0.0103  0.8461
#> 
#> --------------------------------------------------
#>   Chisq(6):             26.2436
#>   p-value:               0.0002
#> --------------------------------------------------
GOFtest(nb2_het)
#> 
#> Goodness-of-fit test for count models
#>    Model: Heteroskedastic Negative Binomial (NB2) Model
#> --------------------------------------------------
#> Manjon & Martinez (2014) Score Test
#> 
#>       Frequency Proportion Probability |Difference| Pearson
#> 0 - 0       401     0.1091      0.0923       0.0168 11.2454
#> 1 - 1       314     0.0854      0.1068       0.0214 15.7195
#> 2 - 2       358     0.0974      0.1043       0.0069  1.6807
#> 3 - 3       334     0.0908      0.0954       0.0046  0.8136
#> 4 - 4       339     0.0922      0.0845       0.0077  2.5548
#> 5 - 5       266     0.0723      0.0735       0.0012  0.0680
#> 6 +        1665     0.4528      0.4432       0.0096  0.7592
#> 
#> --------------------------------------------------
#>   Chisq(6):             36.2833
#>   p-value:               0.0000
#> --------------------------------------------------
```

Even after allowing for heteroskedasticity in the dispersion parameter,
the different tests continue to favor the NB1 model. Interestingly,
modeling dispersion has a larger impact on improving the goodness-of-fit
for the NB2 model, which is consistent with `bh` being more significant
in its scale equation.

## Concluding Remarks

This vignette has illustrated the core count data models available in
`mlmodels`: Poisson, NB1, and NB2, both in their homoskedastic and
heteroskedastic forms.

A consistent finding across diagnostics (log-likelihood, AIC, BIC,
Pearson dispersion, goodness-of-fit tests, and Vuong tests) is that the
**NB1 parameterization** provides a better fit for the `docvis` data
than the more commonly used NB2 model. Allowing the dispersion to vary
with covariates further improves model fit, particularly for NB2, but
the overall preference for NB1 remains.

The `mlmodels` package offers several practical advantages:

- A unified, consistent interface across all models
- Straightforward modeling of heteroskedasticity in the dispersion
  parameter for both NB1 and NB2
- Easy switching between variance estimators without refitting the model
- Built-in overdispersion and detailed goodness-of-fit tests with
  per-bin diagnostics
- Seamless integration with `marginaleffects` for post-estimation
  analysis

We hope this introduction helps you explore count data modeling more
flexibly and transparently in R.

## References

Cameron, A. C., & Trivedi, P. K. (2022). *Microeconometrics Using Stata:
Volumes I and II* (2nd ed.). Stata Press.

Long, J. S., & Freese, J. (2014). *Regression Models for Categorical
Dependent Variables Using Stata* (3rd ed.). Stata Press.

Manjón, M., & Martínez, O. (2014). ‘The chi-squared goodness-of-fit test
for count-data models.’ *The Stata Journal*, 14(4), 798–816.
<https://doi.org/10.1177/1536867X1401400406>

Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
non-nested hypotheses. *Econometrica*, 57(2), 307–333.
<https://doi.org/10.2307/1912557>
