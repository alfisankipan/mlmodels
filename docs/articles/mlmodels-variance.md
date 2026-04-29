# Variance-Covariance Estimation in \`mlmodels\`

``` r

library(mlmodels)
```

## Introduction

Accurate variance-covariance estimation is essential for statistical
inference in maximum likelihood models. The variance-covariance matrix
allows us to:

- Compute **standard errors** for parameter estimates
- Construct **confidence intervals**
- Perform **hypothesis tests** (Wald, z-tests, etc.)
- Calculate **p-values** and assess statistical significance

Without a reliable variance-covariance matrix, we cannot make valid
inferences about the parameters of our model.

The `mlmodels` package provides a comprehensive set of
variance-covariance estimators, through its
[`vcov()`](https://rdrr.io/r/stats/vcov.html) implementation, ranging
from classical information matrix methods to robust and resampling-based
approaches. This vignette explains each method, when to use them, and
how they relate to each other.

### Why Different Variance Estimators Matter

Different estimators make different assumptions about the model:

- **Classical methods** (`oim`, `opg`) assume the model is correctly
  specified (the Information Matrix Identity holds).
- **Robust methods** relax this assumption and are valid under more
  general conditions.
- **Resampling methods** (bootstrap, jackknife) make even fewer
  parametric assumptions and are often the most robust, at the cost of
  higher computation time.

Choosing the right variance estimator is therefore an important part of
responsible statistical practice.

## Classical Information Matrix Estimators

The package provides two classical variance estimators based on the
information matrix:

### Observed Information Matrix (`oim`) — Default

The default variance estimator in `mlmodels` is the **Observed
Information Matrix** (`type = "oim"`). It is used automatically unless
you specify otherwise.

``` r

data(mroz)
mroz$incthou <- mroz$faminc / 1000

# Fit a homoskedastic linear model
fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

# Summary calls vcov() which uses oim as default.
summary(fit)
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
#>  all: Chisq(5) = 972.786, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -29.4778     8.6369  -3.413 0.000643 ***
#>   value::age           1.2623     0.3996   3.159 0.001585 ** 
#>   value::I(age^2)     -0.0136     0.0046  -2.943 0.003250 ** 
#>   value::huswage       1.9566     0.0733  26.691  < 2e-16 ***
#>   value::educ          0.9636     0.1356   7.106 1.19e-12 ***
#>   value::unem         -0.2538     0.0957  -2.651 0.008028 ** 
#> Scale (log(sigma)):  
#>   scale::lnsigma       2.0853     0.0258  80.924  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 5291.36  BIC: 5323.72 
#> Residual standard error (sigma): 8.047

v_oim <- vcov(fit, type = "oim")
sqrt(diag(v_oim))
#> value::(Intercept)         value::age    value::I(age^2)     value::huswage 
#>        8.636925227        0.399627579        0.004609191        0.073304149 
#>        value::educ        value::unem     scale::lnsigma 
#>        0.135599149        0.095744322        0.025768404
```

Notice how the standard errors shown in
[`summary()`](https://rdrr.io/r/base/summary.html) are identical to
those obtained directly from `vcov(fit, type = "oim")`.
[`summary()`](https://rdrr.io/r/base/summary.html) displays the type of
variance used clearly, above the coefficient table.

The `oim` estimator is computationally efficient and statistically
attractive **when the model is correctly specified** (i.e., the outcome
follows the assumed distribution and both the value and scale, if
present, equations are appropriate).

### Outer Product of Gradients (`opg`)

An alternative classical estimator is the **Outer Product of Gradients**
(also known as BHHH).

When a maximum likelihood estimation model is well specified, the
negative Hessian matrix (which is the base for the `oim` variance), and
the outer product of the gradients matrix (which is the basis for the
`opg` variance) whould be equal. This is what is called the
**Information Matrix Identity**.

``` r

v_opg <- vcov(fit, type = "opg")
```

When a maximum likelihood model is correctly specified, the
**Information Matrix Identity** holds. This means that the negative
expected Hessian matrix (basis of the `oim`variance) should be equal (in
expectation) to the outer product of the gradients matrix (basis of the
`opg` variance).

In practice, you can compare the two estimators to get a sense of model
adequacy:

``` r

comp <- data.frame(
  oim = sqrt(diag(v_oim)),
  opg = sqrt(diag(v_opg))
)

comp
#>                            oim         opg
#> value::(Intercept) 8.636925227 9.281669188
#> value::age         0.399627579 0.431781719
#> value::I(age^2)    0.004609191 0.005002672
#> value::huswage     0.073304149 0.044304744
#> value::educ        0.135599149 0.116741374
#> value::unem        0.095744322 0.096801039
#> scale::lnsigma     0.025768404 0.014396434
```

If the model is well identified and correctly specified, the `oim` and
`opg` standard errors should be very similar. Significant differences
between them can be an indication of model misspecification or weak
identification.

This is one of the reasons `mlmodels` provides more robust variance
estimators (`robust`, `boot`, and `jack`) — they do not rely on the
Information Matrix Identity holding.

### The “Hessian may be singular” warning

When computing the `oim` variance, you may occasionally see the
following warning:

> `! OIM variance may be unreliable due to singularity in the Hessian.`

This warning indicates that the observed information matrix is close to
singular (nearly non-invertible).

**What this means**:  
It is a **strong diagnostic signal** of model misspecification or
numerical instability. Common causes include:

- Multicollinearity among predictors.
- Overparameterized models.
- Inappropriate specification of the value or scale equations.
- Poor identification of some parameters.

**What it does NOT mean**:

- It does **not** automatically mean that all coefficient estimates are
  meaningless.
- It does **not** imply that all misspecified models trigger this
  warning. Only that those that do are very likely misspecified.

#### Example: Heteroskedastic Logit Model

A common situation where this warning appears is when the scale
parameters are **counteracting** with the value parameters, creating
very different levels of curvature in the log-likelihood.

We illustrate it with a binary logit example:

``` r


fit_het <- ml_logit(inlf ~ age + I(age^2) + huswage + educ + unem,
                    scale = ~ repwage,
                    data = mroz)
#> ℹ Improving initial values by scaling (factor = 2).
#> ℹ Initial log-likelihood: -477.993
#> ℹ Final scaled log-likelihood: -426.248

summary(fit_het, vcov.type = "oim")
#> ! OIM variance may be unreliable due to singularity in the Hessian.
#> ℹ Relative singularity: "value::(Intercept)" has the largest eigenvalue (72495495908044.1), while "scale::repwage" has the smallest (15.666216).
#> ℹ Consider using `type = 'robust'` or revising the model specification.
#> 
#> Maximum Likelihood Model
#>  Type: Heteroskedastic Binary Logit 
#> ---------------------------------------
#> Call:
#> ml_logit(value = inlf ~ age + I(age^2) + huswage + educ + unem, 
#>     scale = ~repwage, data = mroz)
#> 
#> Log-Likelihood: -349.21 
#> 
#> Wald significance tests:
#>  all: Chisq(6) = 1194.371, Pr(>Chisq) = < 1e-8
#>  Mean: Chisq(5) = 6.844, Pr(>Chisq) = 0.2325
#>  Scale: Chisq(1) = 82.955, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (inlf):  
#>   value::(Intercept)  0.55597    0.23180   2.399   0.0165 *  
#>   value::age         -0.02582    0.01120  -2.306   0.0211 *  
#>   value::I(age^2)     0.00033    0.00014   2.327   0.0200 *  
#>   value::huswage     -0.00046    0.00030  -1.511   0.1308    
#>   value::educ        -0.00310    0.00248  -1.250   0.2115    
#>   value::unem         0.00153    0.00109   1.404   0.1604    
#> Scale (log(sigma)):  
#>   scale::repwage     -1.46321    0.16065  -9.108   <2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:753 (Successes: 428, Failures: 325)
#> Pseudo R-squared - Cor.Sq.: 0.467 McFadden: 0.3218
#> AIC: 712.42  BIC: 744.79 
#> 
#> Distribution of Std. Deviation (sigma):
#> ---------------------------------------
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0000  0.0053  1.0000  0.5603  1.0000  1.0000
```

In this case, [`summary()`](https://rdrr.io/r/base/summary.html)
triggers the warning. While the package still displays results, you
**should not trust** the wald tests, the standard errors, z-statistics,
or p-values shown under the `oim` variance when this warning is
triggered.

One of the recommendations shown in the warning message is to switch to
the `robust` variance estimator. We explore robust and other more
flexible variance options in the next section.

## Robust Variance Estimators

The use of **more robust** variance estimators becomes important when
the Information Matrix Identity may not hold. The singular Hessian
warning we saw in the previous section is one practical indication of
this problem. A more formal diagnostic is the **Information Matrix
Test**
([`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md)),
which we cover in the Model Diagnostics vignette.

### Robust (Sandwich / Huber-White) Estimator

The most commonly used robust estimator is the sandwich (or Huber-White)
estimator. It is the recommended default when you suspect the model may
be misspecified or when you see the singular Hessian warning.

``` r

  # summary() with the vcov.type argument
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

# or vcov() with the type argument
v_robust <- vcov(fit, type = "robust")
```

When you specify `vcov.type = "robust"` in
[`summary()`](https://rdrr.io/r/base/summary.html), the function uses
the robust variance not only for the standard errors, z-values, and
p-values in the coefficient table, but also for the Wald significance
tests shown at the top of the output.

You can also pass a pre-computed variance directly using the `vcov`
argument.

``` r

summary(fit, vcov = v_robust)
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

The displayed results will be identical because the variance matrix is
the same.

#### Clustering

If you suspect that it is the independence of the observations that is
violated, in that observations may be correlated within groups
(clusters), then the **cluster-robust** variance is appropriate:

``` r

# Suspect that income is correlated among education levels.
v_rob_clust <- vcov(fit, type = "robust", cl_var = "educ")
summary(fit, vcov = v_rob_clust)
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
#>  all: Chisq(5) = 1062.475, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust | Clusters: 13 (educ)
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept) -29.4778     6.2953  -4.683 2.83e-06 ***
#>   value::age           1.2623     0.3009   4.196 2.72e-05 ***
#>   value::I(age^2)     -0.0136     0.0033  -4.122 3.76e-05 ***
#>   value::huswage       1.9566     0.0884  22.132  < 2e-16 ***
#>   value::educ          0.9636     0.1770   5.444 5.22e-08 ***
#>   value::unem         -0.2538     0.0658  -3.857 0.000115 ***
#> Scale (log(sigma)):  
#>   scale::lnsigma       2.0853     0.1112  18.753  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations: 753 
#> Residual degrees of freedom: 747 
#> Multiple R-squared: 0.5637 Adjusted R-squared: 0.5608
#> AIC: 5291.36  BIC: 5323.72 
#> Residual standard error (sigma): 8.047
```

### Bootstrap and Jackknife Estimators

The package also provides two resampling-based variance estimators:

- `boot` - Bootstrap (with optional clustering).
- `jack`/`jacknife` - Jackknife (with optional clustering).

**Bootstrap** works by repeatedly resampling the estimation data **with
replacement** (default 999 times, but 500 is often sufficient),
re-fitting the model on each resampled dataset, and then computing the
variance from the resulting matrix of coefficient estimates.

**Jackknife** works by systematically leaving out one observation at a
time and re-fitting the model. The number of replications is therefore
equal to the sample size.

``` r

# Bootstrap with 500 repetitions
v_boot <- vcov(fit, type = "boot", repetitions = 500, seed = 123)
#> ℹ Bootstrap with 500 repetitions.
#>  0        10        20        30        40        50
#> ====================================================
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#> ====================================================
#> 
#> Bootstrapping finished - 100% of replications converged.
# Jackknife
v_jack <- vcov(fit, type = "jack")
#> ℹ Jackknife variance.
#>  0        10        20        30        40        50
#> ====================================================
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ..................................................
#>  ...
#> ====================================================
#> 
#> Jackknife finished - 100% of replications converged.
```

By default, [`vcov()`](https://rdrr.io/r/stats/vcov.html) shows a
progress table, with green dots where the iterations’ fits are
successful, or red crosses where they’re not, when `type` equals `boot`
or `jack`/`jacknife`. You can stop
[`vcov()`](https://rdrr.io/r/stats/vcov.html) from displaying it by
setting `progress = FALSE`.

When you set `vcov.type` to either `"boot"` or `"jack"` in
[`summary()`](https://rdrr.io/r/base/summary.html), or other functions
that also have a `vcov.type` argument, the progress display is not shown
by default (`progress = FALSE` by default).

You can also estimate clustered versions of these variances:

``` r

# Bootstrap with 500 repetitions
v_boot_clust <- vcov(fit, type = "boot", cl_var = "educ",
               repetitions = 500, seed = 123,
               progress = FALSE)
# Jackknife
v_jack_clust <- vcov(fit, type = "jack", cl_var = "educ",
               progress = FALSE)
```

#### Why use resampling methods?

Traditional robust estimators (like the sandwich estimator) still rely
on **asymptotic normality** for inference - that is, they assume that
the sampling distribution of the parameter estimates becomes
approximately normal in large samples.

Bootstrap and jackknife take a more direct approach: they **estimate the
sampling distribution** of the parameters empirically from the data
itself, without assuming normality. This makes them:

- More robust to violations of model assumptions.
- Better at capturing the actual variability in finite samples.
- Particularly useful when the model is complex or the asymptotic
  approximations are questionable.

#### Asymptotic Equivalence

- The regular bootstrap and jackknife are asymptotically equivalent to
  the robust (sandwich) estimator.
- The clustered versions are asymptotically equivalent to the
  cluster-robust estimator.

Here we compare the standard errors produced by each method.

``` r

# Robust Comparison
comp<- data.frame(
  robust = sqrt(diag(v_robust)),
  boot = sqrt(diag(v_boot)),
  jack = sqrt(diag(v_jack))
)
comp
#>                         robust       boot        jack
#> value::(Intercept) 8.411271366 8.31747787 8.507232115
#> value::age         0.382284787 0.38167461 0.386935880
#> value::I(age^2)    0.004391868 0.00440611 0.004446869
#> value::huswage     0.135952916 0.13824719 0.144060112
#> value::educ        0.165176906 0.16856717 0.167837527
#> value::unem        0.098676875 0.10157004 0.100966601
#> scale::lnsigma     0.055550631 0.05434783 0.057408605

# Clustered comparison
comp <- data.frame(
  robust = sqrt(diag(v_rob_clust)),
  boot = sqrt(diag(v_boot_clust)),
  jack = sqrt(diag(v_jack_clust))
)
comp
#>                        robust        boot        jack
#> value::(Intercept) 6.29528438 8.683355306 6.310413395
#> value::age         0.30086228 0.410357979 0.297872890
#> value::I(age^2)    0.00329108 0.004470342 0.003268108
#> value::huswage     0.08840541 0.110632230 0.093001512
#> value::educ        0.17701278 0.181024534 0.212512715
#> value::unem        0.06580827 0.088824945 0.065578436
#> scale::lnsigma     0.11119600 0.112927477 0.154391247
```

The benefit of the resampling methods is that they make fewer parametric
assumptions, which often leads to better finite-sample performance.

#### Important Feature of Our Implementation

When some bootstrap or jackknife replications fail to converge, our
[`vcov()`](https://rdrr.io/r/stats/vcov.html) methods compute the
variance using **only the successful replications**. This is
statistically more appropriate than including failed optimizations
(which is what
[`sandwich::vcovBS()`](https://sandwich.R-Forge.R-project.org/reference/vcovBS.html)
does).

## Practical Recommendations

Choosing the right variance estimator depends on your goals and concerns
about model specification. Here is a practical guide:

### Quick Decision Guide

| Situation | Recommended Variance | Why |
|----|----|----|
| Model is well-specified | `oim` | Most efficient |
| OIM singular Hessian warning | `robust` / `boot` | Robust to violations |
| Suspect heteroskedasticity or other violations | `robust` / `boot` | Robust to violations |
| Suspect correlation within groups | `robust` / `boot` + `cl_var` | Accounts for dependence |

### Best Practices

1.  **Start with `robust`** – In most applied work,
    `vcov.type = "robust"` is expected in journals, even if you’re
    modeling heteroskedasticty.
2.  **Pre-compute expensive variances** – This is mostly for `boot` and
    `jack`, although it also applies to the other ones. When you want to
    use one type of variance for your estimates, predictions, and
    marginal effects, it is very convenient to store the variance and
    feed it to the different functions.
3.  **Check model diagnostics** – Use
    [`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md)
    to assess whether classical estimators (`oim`, `opg`) are suitable.
    See the vignette on diagnostics for a discussion about this test and
    its methods.
4.  **When to prefer bootstrap/jackknife**
    - Small to moderate sample sizes.
    - Complex models.
    - Avoid strong parametric assumptions.
5.  **Reporting** – Always state which variance estimator, or estimators
    if you use several, you used in your paper or report.

## Concluding Remarks

The `mlmodels` package produces a comprehensive and coherent suite of
variance-covariance estimators through the unified
[`vcov()`](https://rdrr.io/r/stats/vcov.html) function. From the
classical information matrix estimation (`oim`, `opg`), to robust and
cluster-robust methods, allowing you to get these via parametric
assumptions (`robust`) or via more flexible resampling methods (`boot`
and `jack`).

Our design philosophy has been to give you both **power and guidance**:

- Clear informative warnings when OIM estimator may be unreliable.
- Robust default behavior and sensible options for advanced users.
- Consistent interfaces for switching between variance stypes is
  straightforward.
- Careful handling of edge cases, such as using only successful
  replications in bootstrap and jackknife procedures.

We built these tools because we observed that other existing packages
did not provide this level of detail and/or solid estimation methods.
Our goal was to give you the right standard errors for your estimations,
predictions, and marginal effects, by simply setting one or two
arguments (`type` and `cl_var`).

We hope his vignette helps you understand the strengths and appropriate
use of each estimator, so you can conduct your analyses with confidence
and transparency.

Happy modeling!
