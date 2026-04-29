# Maximum Likelihood Models in R

``` r

library(mlmodels)
library(marginaleffects)
```

The `mlmodels` package provides a consistent and flexible framework for
estimating a wide range of maximum likelihood models in R.

Unlike [`glm()`](https://rdrr.io/r/stats/glm.html), which is limited to
a small set of distributions and does not allow direct modeling of the
scale parameter (variance, dispersion, precision, or shape), `mlmodels`
lets you explicitly model both the **location/mean** and the **scale**
of the outcome in a unified way.

Even for homoskedastic models, maximum likelihood estimation is often
statistically more efficient than the quasi-likelihood approach used by
[`glm()`](https://rdrr.io/r/stats/glm.html). With `mlmodels`, you get
this efficiency **plus** the ability to model heteroskedasticity,
dispersion, precision, shape or other parameters particular to the
model.

## Core Features

- Consistent S3 interface across all models.
- Support for modeling scale parameters alongside the mean.
- Rich [`predict()`](https://rdrr.io/r/stats/predict.html) method with
  many useful output types and integration with `marginaleffects`.
- Multiple variance-covariance estimators (`oim`, `opg`, `robust`,
  `boot`, `jack`, clustered)
- Standard R methods that “just work”:
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html), etc.
- Comprehensive suite of hypothesis tests (`waldtest`, `lrtest`,
  `IMtest`, `vuongtest`, `OVDtest`, `GOFtest`)

All models share the same syntax and workflow, making it easy to compare
specifications and switch between model families. The package is
designed with practitioners in mind, focusing on the information and
flexibility researchers need in their work.

The package builds on the excellent
[`maxLik`](https://cran.r-project.org/package=maxLik) package for
maximum likelihood optimization (Henningsen and Toomet, 2011). This
foundation allows `mlmodels` to deliver precise, numerically stable
estimations while providing a consistent and intuitive user interface
across all model families.

## Quick Start

Let’s fit a simple linear model to get a feel for the package:

``` r

# Load example data
data(mroz)
mroz$incthou <- mroz$faminc / 1000

# Fit a homoskedastic linear model
fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

# View the results with robust standard errors
summary(fit, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = incthou ~ age + I(age^2) + huswage + educ + unem, 
#>     data = mroz)
#> 
#> Log-Likelihood: -2638.68 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 325.926, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -29.4778     8.4113  -3.505 0.000457 ***
#>   value::age           1.2623     0.3823   3.302 0.000960 ***
#>   value::I(age^2)     -0.0136     0.0044  -3.089 0.002011 ** 
#>   value::huswage       1.9566     0.1360  14.391  < 2e-16 ***
#>   value::educ          0.9636     0.1652   5.834 5.42e-09 ***
#>   value::unem         -0.2538     0.0987  -2.572 0.010109 *  
#> Scale (log(sigma)):  
#>   scale::lnsigma       2.0853     0.0556  37.538  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 5291.36  BIC: 5323.72 
#> Residual standard error (sigma): 8.047
```

To get the summary statistics from the model, you use the familiar
function [`summary()`](https://rdrr.io/r/base/summary.html). Through the
use of the `vcov.type` argument, you select what type of standard errors
you want to use for inference. In this case, `robust` (sandwich
estimator).

The presentation is clean, and consistent for all models. It has a
header with the type of model you fitted. It, then, presents the call,
the log-likelihood and joint significance test(s). In this case, since
the estimation was homoskedastic, the overall significance test is the
only one presented.

After that, you are presented with the type of standard errors used for
inference and, immediately after, the coefficients table separated into
two clear groups: value and scale. The value section holds the
parameters of the conditional mean, and the scale section the parameters
of the standard deviation, in this case.

All our models add the prefixes `value::` and `scale::` to the names of
the coefficients in the corresponding equation. The purpose is twofold:

- It avoids naming the coefficients equally if you entered the same
  variable as a predictor in both equations.

- It points immediately to which equation the coefficient belongs to,
  giving you an easy way to select just the coefficients of a given
  equation, in case you’d like to separate them.

As we mentioned, we use the name scale consistently for that equation
across models. The following is a table with what the scale equation
models for each of the estimators that model scale in the package:

| Model | Scale Parameter | Interpretation |
|----|----|----|
| Normal / Lognormal | `log(sigma)` | Log of the standard deviation |
| Logit / Probit | `log(sigma)` | Log of the deviations from the unidentified overall standard deviation |
| Negative Binomial | `log(alpha)` | Log of the dispersion parameter |
| Gamma | `log(nu)` | Log of the shape parameter |
| Beta | `log(phi)` | Log of the precision parameter |

The purpose of naming the equation scale consistently across models,
even though it may reflect different type of parameters across models,
is that you can use the same format to access the estimates of the
equation, independent of what model you used. That way you don’t have to
remember exactly which parameter you are modeling in each case to access
its coefficients and standard errors. Nevertheless, in the header for
the scale equation, you are reminded of what that equation actually
models (`log(sigma)` in this case).

Finally, after the coefficients’ table, we present information about the
number of observations, goodness of fit measures, and information about
the predicted values of the scale equation. Since this estimation was of
a homoskedastic normal model, the prediction is if `sigma`.

## Making Predictions

The [`predict()`](https://rdrr.io/r/stats/predict.html) method is very
flexible:

``` r

# Basic prediction (expected value)
head(predict(fit, type = "response")$fit)
#> [1] 15.20265 21.47112 15.38632 14.98368 27.26299 21.92135

# Variance of the outcome
head(predict(fit, type = "variance")$fit)
#> [1] 64.75192 64.75192 64.75192 64.75192 64.75192 64.75192

# Prediction with standard errors using robust variance
pred <- predict(fit, type = "response", vcov.type = "robust", se.fit = TRUE)
head(pred$fit)      # Fitted values / expected values
#> [1] 15.20265 21.47112 15.38632 14.98368 27.26299 21.92135
head(pred$se.fit)   # Standard errors
#> [1] 0.5987687 0.7375031 0.5413427 0.5557155 0.7622242 0.5022675
```

## Modeling Heteroskedasticity

One of the most powerful features of `mlmodels` is the ability to
explicitly model the scale parameter (variance/dispersion) alongside the
mean. To do that we add the `scale` argument to the call, with a formula
with just right hand side variables: the predictors for that equation.

``` r

# Heteroskedastic linear model
fit_het <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem,
                 scale = ~ educ + exper, 
                 data = mroz)

summary(fit_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = incthou ~ age + I(age^2) + huswage + educ + unem, 
#>     scale = ~educ + exper, data = mroz)
#> 
#> Log-Likelihood: -2621.78 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 455.971, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(5) = 437.442, Pr(>Chisq) = < 1e-8
#>  Scale: Chisq(2) = 10.428, Pr(>Chisq) = 0.0054
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -28.8642     7.7943  -3.703 0.000213 ***
#>   value::age           1.2828     0.3649   3.516 0.000438 ***
#>   value::I(age^2)     -0.0139     0.0042  -3.341 0.000836 ***
#>   value::huswage       1.9331     0.1292  14.960  < 2e-16 ***
#>   value::educ          0.8795     0.1298   6.775 1.24e-11 ***
#>   value::unem         -0.2114     0.0920  -2.297 0.021617 *  
#> Scale (log(sigma)):  
#>   scale::(Intercept)   1.5549     0.3170   4.905 9.36e-07 ***
#>   scale::educ          0.0503     0.0229   2.200 0.027814 *  
#>   scale::exper        -0.0103     0.0063  -1.652 0.098598 .  
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5634 Adjusted R-squared: 0.5605
#> AIC: 5261.56  BIC: 5303.17 
#> 
#> Distribution of Std. Deviation (sigma):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     4.4     7.3     7.9     7.9     8.5    11.0
```

Notice how the presentation follows the pattern we described for the
homoskedastic case. The main differences are that now:

- The joint significance tests include one for the parameters of the
  value equation, and another for the parameters of the scale equation,
  in addition to the overall one.
- At the bottom we present summary statistics of the predicted standard
  deviation, instead of the single value.

This allows you to immediately identify if there is heteroskedasticity
present, and how much the standard deviation varies within the sample
used in the estimation.

The real benefit becomes clear when you make predictions:

``` r

  head(predict(fit_het, type = "variance")$fit)
#> [1] 56.10345 67.57323 54.95579 66.19094 79.28079 37.88290
```

Compare this to the homoskedastic model earlier, where we predicted that
the variance was the same for every observation. Modeling the scale
parameter allows the uncertainty to vary realistically across
individuals.

To get the partial effects on the standard deviation, you can use the
functions from the `marginaleffects` package:

``` r

# AMEs for the standard deviation using robust standard errors
avg_slopes(fit_het, variables = c("educ", "exper"), type = "sigma", vcov = "robust")
#> 
#>   Term Estimate Std. Error     z Pr(>|z|)   S   2.5 % 97.5 %
#>  educ    0.3993     0.1826  2.19   0.0287 5.1  0.0415 0.7571
#>  exper  -0.0821     0.0506 -1.62   0.1049 3.3 -0.1813 0.0171
#> 
#> Type: sigma
#> Comparison: dY/dX
```

The only model that does not accept a scale argument is
[`ml_poisson()`](https://alfisankipan.github.io/mlmodels/reference/ml_poisson.md),
because its probability mass function is fully determined by its mean,
as is its variance.

For a deeper look at predictions, see the [Predictions with
`mlmodels`](https://alfisankipan.github.io/mlmodels/articles/mlmodels-predictions.md)
vignette.

## Controlling Intercepts

All `mlmodels` estimating functions allow you to explicitly control
whether an intercept is included in the `value` (mean) equation and/or
the `scale` equation. This is done through the dedicated arguments
`noint_value` and `noint_scale`. Trying to do it through the formulas in
the `value` and/or `scale` arguments will throw an error and abort the
estimation.

Let’s see how:

``` r

# No intercept in the value function (wrong way)
fit_no_int <- tryCatch({
  ml_lm(incthou ~ 0 + age + I(age^2) + huswage + educ + unem,
                 scale = ~ educ + exper, 
                 data = mroz)
}, error = function(e) {
  print(e)
  invisible(NULL)
})
#> <error/rlang_error>
#> Error in `hardhat::mold()`:
#> ! `formula` must not contain the intercept removal term: `+ 0` or `0 +`.
#> ---
#> Backtrace:
#>     ▆
#>  1. ├─base::tryCatch(...)
#>  2. │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#>  3. │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#>  4. │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#>  5. └─mlmodels::ml_lm(...)
#>  6.   ├─hardhat::mold(value, data_clean, blueprint = hardhat::default_formula_blueprint(intercept = !noint_value))
#>  7.   └─hardhat:::mold.formula(value, data_clean, blueprint = hardhat::default_formula_blueprint(intercept = !noint_value))

# No intercept in the value function (right way)
fit_no_int <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem,
                 scale = ~ educ + exper, 
                 noint_value = TRUE,
                 data = mroz)
summary(fit_no_int)
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Linear Model 
#> ---------------------------------------
#> Call:
#> ml_lm(value = incthou ~ age + I(age^2) + huswage + educ + unem, 
#>     scale = ~educ + exper, data = mroz, noint_value = TRUE)
#> 
#> Log-Likelihood: -2627.74 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 6994.167, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(5) = 6919.276, Pr(>Chisq) = < 1e-8
#>  Scale: Chisq(2) = 33.037, Pr(>Chisq) = 0.0000
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::age         -0.00878    0.08750  -0.100  0.92009    
#>   value::I(age^2)     0.00066    0.00125   0.528  0.59763    
#>   value::huswage      1.95042    0.07461  26.142  < 2e-16 ***
#>   value::educ         0.78254    0.12799   6.114 9.70e-10 ***
#>   value::unem        -0.22951    0.09294  -2.469  0.01354 *  
#> Scale (log(sigma)):  
#>   scale::(Intercept)  1.54220    0.11733  13.144  < 2e-16 ***
#>   scale::educ         0.05110    0.00943   5.421 5.94e-08 ***
#>   scale::exper       -0.00935    0.00323  -2.896  0.00378 ** 
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 748 
#> Multiple R-squared: 0.5568 Adjusted R-squared: 0.5545
#> AIC: 5271.49  BIC: 5308.48 
#> 
#> Distribution of Std. Deviation (sigma):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     4.5     7.3     7.9     8.0     8.6    11.0
```

In the above code we catch the error and reproduce the stop message you
would have received from
[`hardhat::mold()`](https://hardhat.tidymodels.org/reference/mold.html),
so that the rest of the code could be executed, and you can see the
wording you will receive if you try that.

You can see how in the second estimation’s results there is no intercept
estimate for the value equation, which is what we wanted.

**Important Notes**

- Using formula syntax like `~0 + ...` or `~-1 + ...` will not work as
  expected. It will trigger errors from
  [`hardhat::mold()`](https://hardhat.tidymodels.org/reference/mold.html).
- Always use the dedicated arguments:
  - `noint_value = TRUE` → removes intercept from the mean/value
    equation
  - `noint_scale = TRUE` → removes intercept from the scale equation
    (when present)
- [`ml_poisson()`](https://alfisankipan.github.io/mlmodels/reference/ml_poisson.md)
  only has noint_value (its proability doesn’t have a scale parameter).

**Why this design?**

Our models use `hardhat` blueprints for consistent and safe formula
handling. The `noint_*` arguments give you explicit, reliable control
over intercepts while maintaining compatibility with the internal
structure.

## Binary Outcomes: Logit

The same consistent syntax applies to all our models. Here we use the
`inlf` variable from the `mroz` dataset (whether the wife is in the
labor force) to estimate a logit model, because it is already binary
(0/1).

``` r

# Homoskedastic logit model
fit_logit <- ml_logit(inlf ~ age + I(age^2) + huswage + educ + unem, 
                      data = mroz)
summary(fit_logit, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Binary Logit 
#> ---------------------------------------
#> Call:
#> ml_logit(value = inlf ~ age + I(age^2) + huswage + educ + unem, 
#>     data = mroz)
#> 
#> Log-Likelihood: -490.46 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 43.750, Pr(>Chisq) = 0.0000
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (inlf):  
#>   value::(Intercept)  -6.5743     2.2340  -2.943  0.00325 ** 
#>   value::age           0.2484     0.1033   2.405  0.01616 *  
#>   value::I(age^2)     -0.0030     0.0012  -2.538  0.01114 *  
#>   value::huswage      -0.0683     0.0192  -3.554  0.00038 ***
#>   value::educ          0.2126     0.0379   5.615 1.96e-08 ***
#>   value::unem         -0.0156     0.0248  -0.628  0.52992    
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:753 (Successes: 428, Failures: 325)
#> Pseudo R-squared - Cor.Sq.: 0.06239 McFadden: 0.04742
#> AIC: 992.91  BIC: 1020.66

# Heteroskedastic logit model
fit_logit_het <- ml_logit(inlf ~ age + I(age^2) + huswage + educ + unem,
                          scale = ~ educ + exper,
                          data = mroz)

summary(fit_logit_het, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Binary Logit 
#> ---------------------------------------
#> Call:
#> ml_logit(value = inlf ~ age + I(age^2) + huswage + educ + unem, 
#>     scale = ~educ + exper, data = mroz)
#> 
#> Log-Likelihood: -467.67 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 1834.641, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(5) = 1.769, Pr(>Chisq) = 0.8800
#>  Scale: Chisq(2) = 72.448, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                        Estimate Std. Error z value Pr(>|z|)     
#> Value (inlf):  
#>   value::(Intercept)  0.135979   0.108547   1.253   0.2103    
#>   value::age         -0.004567   0.003716  -1.229   0.2190    
#>   value::I(age^2)     0.000037   0.000031   1.184   0.2364    
#>   value::huswage     -0.000011   0.000313  -0.036   0.9714    
#>   value::educ         0.000150   0.000373   0.403   0.6869    
#>   value::unem         0.000430   0.000487   0.883   0.3773    
#> Scale (log(sigma)):  
#>   scale::educ        -0.117931   0.061091  -1.930   0.0536 .  
#>   scale::exper       -0.162190   0.020692  -7.838 4.57e-15 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:753 (Successes: 428, Failures: 325)
#> Pseudo R-squared - Cor.Sq.: 0.1413 McFadden: 0.09168
#> AIC: 951.34  BIC: 988.33 
#> 
#> Distribution of Std. Deviation (sigma):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00018 0.01923 0.04869 0.08073 0.11794 0.55452
```

Notice how the modeling of the homoskedastic and heteroskedastic models
are identical in structure to the linear case. The same pattern works
for
[`ml_probit()`](https://alfisankipan.github.io/mlmodels/reference/ml_probit.md),
[`ml_negbin()`](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md),
[`ml_gamma()`](https://alfisankipan.github.io/mlmodels/reference/ml_gamma.md),
and
[`ml_beta()`](https://alfisankipan.github.io/mlmodels/reference/ml_beta.md).

For a deeper discussion of the different models see the dedicated
vignettes to the topics that involve them.

## Where to Go Next

You have now seen the core philosophy of the `mlmodels` package: a
consistent interface that works the same way whether you are estimating
a linear model, a count model, a binary choice model, or a fractional
response model.

To go deeper, explore the following vignettes:

- **[Predictions with
  `mlmodels`](https://alfisankipan.github.io/mlmodels/articles/mlmodels-predictions.md)**
  – Detailed guide to the unified
  [`predict()`](https://rdrr.io/r/stats/predict.html) method and how to
  compute predictions, marginal effects, and bootstrap the predictions
  with the `marginaleffects` package.
- **[Variance-Covariance Estimation in
  `mlmodels`](https://alfisankipan.github.io/mlmodels/articles/mlmodels-variance.md)**
  – Learn about the different variance-covariance estimators available
  (`oim`, `opg`, `robust`, `cluster`, `boot`, and `jack`), when to use
  each one, how they relate to one another, and why we implemented our
  own [`vcov()`](https://rdrr.io/r/stats/vcov.html) method instead of
  relying solely on external packages.
- **[Introduction to Count
  Data](https://alfisankipan.github.io/mlmodels/articles/mlmodels-countintro.md)**
  – We illustrate how to work with our count data models while
  reproducing classic examples from *Microeconometrics Using Stata*
  (Cameron and Trivedi 2022), including overdispersion and
  goodness-of-fit tests. We also demonstrate additional functionality,
  most notably **heteroskedastic NB1 estimation**, allowing direct
  comparison between homoskedastic and heteroskedastic versions of both
  NB1 and NB2 models.
- **[Gamma versus
  Lognormal](https://alfisankipan.github.io/mlmodels/articles/mlmodels-gamma-lognormal.md)**
  – Comparing two strong estimators for positive, right-skewed outcomes.
  We look at estimation, predicted moments of the outcome (mean,
  variance, etc.), and practical differences between the two approaches.
- **[Fractional Response
  Outcomes](https://alfisankipan.github.io/mlmodels/articles/mlmodels-fractional.md)**
  – We consider QMLE estimation with logit or probit, and MLE estimation
  with the Beta model. We go over when either type is suitable, and how,
  as well as look at an application from the seminal article, Papke and
  Wooldridge (1996), to compare the estimators.
- **[Diagnostic Tools in
  `mlmodels`](https://alfisankipan.github.io/mlmodels/articles/mlmodels-diagnostics.md)**
  – Learn about the diagnostic tests implemented in the package
  (`IMtest`, `waldtest`, `OVDtest`, etc.), how to use them, and how to
  interpret the results.

We hope this vignette has given you a solid foundation for getting
started with maximum likelihood modeling in R.

Happy modeling!

## References

Cameron, A. C., & Trivedi, P. K. (2022). *Microeconometrics Using Stata:
Volumes I and II* (2nd ed.). Stata Press.

Henningsen, A., & Toomet, O. (2011). “maxLik: A package for maximum
likelihood estimation in R.” *Computational Statistics*, 26(3), 443–458.
<https://doi.org/10.1007/s00180-010-0217-1> (Also available as the R
package `maxLik` on CRAN: <https://cran.r-project.org/package=maxLik>)

Papke, L. E., & Wooldridge, J. M. (1996). “Econometric methods for
fractional response variables with an application to 401(k) plan
participation rates.” *Journal of Applied Econometrics*, 11(6), 619–632.
<https://doi.org/10.1002/(SICI)1099-1255(199611)11:6%3C619::AID-JAE418%3E3.0.CO;2-1>
