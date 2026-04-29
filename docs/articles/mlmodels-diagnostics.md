# Diagnostic Tools in \`mlmodels\`

``` r

library(mlmodels)
```

## Introduction

No estimation package would be complete without providing the means to
test the validity of its estimates. Tests that practitioners need range
from misspecification, to significance of the parameters, as well as
model comparison, whether nested or not.

The `mlmodels` package provides a suite of testing functions, several of
which are general for all models, and a couple of which that are
designed for count data models. In this vignette we explore them.

## Information Matrix Test

With maximum likelihood estimators, the main question that we usually
ask ourselves is: **does the information matrix equality hold?** The
**information matrix test** – White (1982) – attempts to answer this
question. Even though it is a test of misspecification, from a
practitioner’s point of view, it is more a test of how to do inference
with your model. If the test rejects the equality, inference should be
done using robust methods.

The `mlmodels`
[`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md)
function offers four methods of the test:

- `quad` – **Analyitical** quadratic orthogonalized (Davidson and
  MacKinnon, 1992).
- `boot_quad` – **Analyitical** quadratic orthogonalized plus
  bootstrapped *p*-values.
- `opg` – **Analytical** Chesher/Lancester score version of the test
  (Chesher, 1983; Lancaster, 1984).
- `boot_opg` – **Analytical** Chesher/Lancaster score version of the
  test plus bootstrapped *p*-values.

Even though both analytical versions of the test are asymptotically
equivalent, each has shown practical limitations. The
**Chesher/Lancaster (OPG)** version is known for **extreme size
distortion**, often rejecting correctly specified models (Horowitz,
1994). The **orthogonalized quadratic version** is more robust to size
problems, but its construction can lead to a **loss of power**, failing
to detect actual misspecification.

Bootstrapping the tests – an approach pioneered for the IM test by
Horowitz (1994) and further championed by Davidson and MacKinnon (1999)
– bridges these gaps. Rather than relying on asymptotic critical values,
our implementation provides bootstrap *p*-values by estimating the test
statistic’s distribution directly from the sample data, which
significantly improves the reliability of the test in finite samples.

### Implementing the Test

Let’s consider a linear model, and how to implement the test.

``` r

data(mroz)
mroz$incthou <- mroz$faminc / 1000

# Fit a homoskedastic linear model
fit_lm <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

# Default is quad (quadratic orthognalized)
IMtest(fit_lm)
#> Information Matrix Test
#>  Method: Orthogonalized Quadratic Form 
#>  Model:  Homoskedastic Linear Model 
#> --------------------------------------------
#>  Chisq(28) = 122.732    Pr(>Chisq) = 0.0000
#> --------------------------------------------

# Let's check the Chesher/Lancaster
IMtest(fit_lm, method = "opg")
#> Information Matrix Test
#>  Method: Chesher/Lancaster OPG 
#>  Model:  Homoskedastic Linear Model 
#> --------------------------------------------
#>  Chisq(28) = 298.563    Pr(>Chisq) = 0.0000
#> --------------------------------------------
```

We see that both analytical versions of the test reject the null, and
state that the information matrix equality doesn’t hold. This is not
surprising, the outcome variable is skewed to the right, so that may be
what the tests are capturing. Let us check what happens when we add the
bootsrapped *p*-values.

``` r

# Bootstrapped quadratic (low repetitions for speed)
IMtest(fit_lm, method = "boot_quad", repetitions = 100)
#> Information Matrix Test
#>  Method: Orthogonalized Quadratic Form + Model-based bootstrap 
#>  Model:  Homoskedastic Linear Model 
#> --------------------------------------------
#>  Repetitions: Total 100 Successful 100 
#>  Chisq(28) = 122.732
#>  P(>Chisq): Analytical   = 0.0000 
#>             Bootstrapped = 0.8100
#> --------------------------------------------

# Boot opg (low repetitions for speed)
IMtest(fit_lm, method = "boot_opg", repetitions = 100)
#> Information Matrix Test
#>  Method: Chesher/Lancaster OPG + Model-based bootstrap 
#>  Model:  Homoskedastic Linear Model 
#> --------------------------------------------
#>  Repetitions: Total 100 Successful 100 
#>  Chisq(28) = 298.563
#>  P(>Chisq): Analytical   = 0.0000 
#>             Bootstrapped = 0.9600
#> --------------------------------------------
```

We do the bootstrapping with a low number of repetitions, for
illustration purposes. For actual research you should consider 500
repetitions or more, depending on the size of your dataset.

Both bootstrapped *p*-values oppose the analytical test, and tell us
that we cannot reject that the model is well specified. Which of the two
should we believe? You, actually, should take this as an indication that
there is some level of misspecification because both analytical versions
rejected the null. From a practical point of view, however, the
bootstrap *p*-values are telling you that you are probably good using
the model for predictions. Finally, the contrast between analytic and
bootstrapped *p*-values in both versions should not be ignored, and you
should proceed **with robust inference** to be safe.

### Illustrating the Size and Power Problems

Let us fit a Gamma model with the same specification, since the Gamma
model is specific for positive outcomes with long right tails, and test
it.

``` r

fit_gam <- ml_gamma(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

IMtest(fit_gam, method = "boot_quad", repetitions = 100)
#> Information Matrix Test
#>  Method: Orthogonalized Quadratic Form + Model-based bootstrap 
#>  Model:  Homoskedastic Gamma Model 
#> --------------------------------------------
#>  Repetitions: Total 100 Successful 100 
#>  Chisq(28) = 9.481
#>  P(>Chisq): Analytical   = 0.9996 
#>             Bootstrapped = 0.7800
#> --------------------------------------------

IMtest(fit_gam, method = "boot_opg", repetitions = 100)
#> Information Matrix Test
#>  Method: Chesher/Lancaster OPG + Model-based bootstrap 
#>  Model:  Homoskedastic Gamma Model 
#> --------------------------------------------
#>  Repetitions: Total 100 Successful 100 
#>  Chisq(28) = 81.513
#>  P(>Chisq): Analytical   = 0.0000 
#>             Bootstrapped = 0.9800
#> --------------------------------------------
```

We did the bootstrapped versions of the tests, because they allow us to
see the analytical result and the bootstrapped *p*-values. The
*p*-values of the analytical versions of the test, illustrate clearly
the weaknesses of each of the tests.

The model is probably well specified, but a *p*-value of 0.9996 is very
high, and it shows that the quadratic orthogonalized version may be
losing power with this model. The opposite *p*-value is on the
analytical Chesher/Lancaster version, showing that it may be strongly
rejecting a perhaps valid model (the size problem).

When we add the bootstrapped *p*-values, the picture becomes clearer.
They both fail to reject that the model is well specified. But look at
how much more reasonable the bootstrapped *p*-value for the quadratic
orthogonalized model is with respect to its analytical counterpart (0.78
versus 0.9996).

In this case you’re probably safe to proceed with `oim` standard errors,
but you may want to check the results with `robust` ones as well, to see
if there are major differences in significance, just in case.

``` r

# Show results with oim standard errors
summary(fit_gam)
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Gamma Model 
#> ---------------------------------------
#> Call:
#> ml_gamma(value = incthou ~ age + I(age^2) + huswage + educ + 
#>     unem, data = mroz)
#> 
#> Log-Likelihood: -2556.78 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 669.603, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept)  0.59821    0.37344   1.602 0.109174    
#>   value::age          0.06588    0.01725   3.818 0.000134 ***
#>   value::I(age^2)    -0.00071    0.00020  -3.567 0.000361 ***
#>   value::huswage      0.07594    0.00357  21.291  < 2e-16 ***
#>   value::educ         0.04362    0.00581   7.511 5.85e-14 ***
#>   value::unem        -0.01169    0.00420  -2.786 0.005337 ** 
#> Scale (log(nu)):  
#>   scale::lnnu         2.11298    0.05053  41.815  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:753 Deg. of freedom: 747
#> Pseudo R-squared - Cor.Sq.: 0.3446 McFadden: 0.09782
#> AIC: 5127.57  BIC: 5159.94 
#> Shape Param.: 8.27  - Coef.Var.: 0.35

summary(fit_gam, vcov.type = "robust")
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Gamma Model 
#> ---------------------------------------
#> Call:
#> ml_gamma(value = incthou ~ age + I(age^2) + huswage + educ + 
#>     unem, data = mroz)
#> 
#> Log-Likelihood: -2556.78 
#> 
#> Wald significance tests:
#>  all: Chisq(5) = 347.459, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Robust
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept)  0.59821    0.36392   1.644 0.100220    
#>   value::age          0.06588    0.01686   3.907 9.34e-05 ***
#>   value::I(age^2)    -0.00071    0.00019  -3.658 0.000254 ***
#>   value::huswage      0.07594    0.00528  14.389  < 2e-16 ***
#>   value::educ         0.04362    0.00706   6.182 6.33e-10 ***
#>   value::unem        -0.01169    0.00432  -2.709 0.006747 ** 
#> Scale (log(nu)):  
#>   scale::lnnu         2.11298    0.07231  29.223  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:753 Deg. of freedom: 747
#> Pseudo R-squared - Cor.Sq.: 0.3446 McFadden: 0.09782
#> AIC: 5127.57  BIC: 5159.94 
#> Shape Param.: 8.27  - Coef.Var.: 0.35
```

We see that the bigger differences are in the standard errors for
`huswage` and `educ`, but that they’re not enough to change the
conclusions about the significance of the parameters. For the rest of
the coefficients, the standard errors are very close across results. We
also see that Wald’s overall significance chi-squared has decreased, but
the significance of the test statistic hasn’t. This confirms our
interpretation of the information matrix tests on the model.

## Testing Nested Models

With maximum likelihood estimators, there are two usual ways to test
across nested models:

- **Likelihood-ratio** test – You estimate the unrestricted and the
  restricted (nested) model, and do the test.
- **Wald** test – You only estimate the unrestricted model, and the
  tests are tests restrictions of the coefficients.

As you can see, we present Wald test(s) in our estimation results, not
likelihood-ratio test(s), of significance. The reason is that the
**likelihood-ratio test is not robust** to the information matrix
equality not holding, whereas the **Wald test is**, when it uses the
variance that is robust to that equality not holding. Our
[`summary()`](https://rdrr.io/r/base/summary.html) function uses either
the type of variance you set in the `vcov.type` argument, which defaults
to `oim`, or the variance matrix you pass through the `vcov` argument,
to produce the standard errors of the coefficients and to do the Wald
test(s). That way it produces consistent and robust tests, as long a you
specify an asymptotically consistent type of variance for the model.

The limitation of the likelihood-ratio test is an indication that you
should use the insights provided by the information matrix test to guide
you in selecting what type of test to do with nested models.

### Likelihood-ratio Test

We fit a gamma model with more predictors, to have the unrestricted
model, and then test the nested models using `mlmodels`
[`lrtest()`](https://alfisankipan.github.io/mlmodels/reference/lrtest.md):

``` r

fit_gam_ur <- ml_gamma(incthou ~ age + I(age^2) + huswage + educ + unem + 
                         kidslt6 + kidsge6, 
                       data = mroz)

summary(fit_gam_ur)
#> 
#> Maximum Likelihood Model
#>  Type: Homoskedastic Gamma Model 
#> ---------------------------------------
#> Call:
#> ml_gamma(value = incthou ~ age + I(age^2) + huswage + educ + 
#>     unem + kidslt6 + kidsge6, data = mroz)
#> 
#> Log-Likelihood: -2554.40 
#> 
#> Wald significance tests:
#>  all: Chisq(7) = 677.722, Pr(>Chisq) = < 1e-8
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#>                       Estimate Std. Error z value Pr(>|z|)     
#> Value (incthou):  
#>   value::(Intercept)  0.87421    0.39479   2.214  0.02680 *  
#>   value::age          0.05400    0.01835   2.942  0.00326 ** 
#>   value::I(age^2)    -0.00059    0.00021  -2.778  0.00547 ** 
#>   value::huswage      0.07662    0.00357  21.453  < 2e-16 ***
#>   value::educ         0.04395    0.00582   7.547 4.47e-14 ***
#>   value::unem        -0.01158    0.00419  -2.763  0.00573 ** 
#>   value::kidslt6     -0.05910    0.02730  -2.165  0.03038 *  
#>   value::kidsge6      0.00382    0.01107   0.345  0.72995    
#> Scale (log(nu)):  
#>   scale::lnnu         2.11906    0.05054  41.930  < 2e-16 ***
#> ---------------------------------------
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Number of observations:753 Deg. of freedom: 745
#> Pseudo R-squared - Cor.Sq.: 0.3397 McFadden: 0.09866
#> AIC: 5126.80  BIC: 5168.42 
#> Shape Param.: 8.32  - Coef.Var.: 0.35

# The order you enter the models doesn't matter
lrtest(fit_gam_ur, fit_gam)
#> ℹ `object_2` is the restricted model (nested in `object_1`).
#> Likelihood Ratio Test
#> --------------------------------------------
#> Chisq(2) = 4.766    Pr(>Chisq) = 0.0923
#> --------------------------------------------
#> LogLik (restricted)   : -2556.784   (df = 9)
#> LogLik (unrestricted) : -2554.401   (df = 7)
```

In the results we see that the function determines which model is nested
in which, so the order that you enter them in the function doesn’t
matter. The test indicates that can reject the null at the 10% level, so
at least one of the two added variables seems to be significant. Looking
at the estimation result, it seems to be `kidslt6`.

### Wald Test

The
[`waldtest()`](https://alfisankipan.github.io/mlmodels/reference/waldtest.md)
function gives you two arguments to perform this same test:

- `coef_names` allows you to pass the vector with the names of the
  coefficients that you want to test.
- `indices` allows you to pass the vector with the order indices of the
  coefficients that you want to test.

Looking at the estimation results, we see that the names are
`value::kidslt6` and `value::kidsge6`, and they are in 7^(th) and 8^(th)
position, respectively, in the coefficients’ vector.

``` r

# Using indices
waldtest(fit_gam_ur, indices = 7:8)
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#> Restrictions:
#>   1: value::kidslt6 = 0
#>   2: value::kidsge6 = 0
#> Chisq(2) = 4.853    Pr(>Chisq) = 0.08836
#> --------------------------------------------

# Using coef_names
waldtest(fit_gam_ur, coef_names = c("value::kidslt6", "value::kidsge6"))
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#> Restrictions:
#>   1: value::kidslt6 = 0
#>   2: value::kidsge6 = 0
#> Chisq(2) = 4.853    Pr(>Chisq) = 0.08836
#> --------------------------------------------
```

We see that both methods return the same result, and that it is
consistent with that of the likelihood-ratio test.

The
[`waldtest()`](https://alfisankipan.github.io/mlmodels/reference/waldtest.md)
function also has a `rhs` argument, that defaults to 0, which allows you
to specify if you want to test that a coefficient equals something else
than 0. If you set it to a single number, it will use that number as the
right hand side of the equalities for all the coefficients you pass. If
you want to have different right hand side values for different
coefficients you pass a vector of the values, in the same order as the
vector of coefficients you pass in either `indices` or `coef_names`.

``` r

# Testing that age equals .05 and I(age^2) equals 0
waldtest(fit_gam_ur, coef_names = c("value::age", "value::I(age^2)"), rhs = c(.05, 0))
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#> Restrictions:
#>   1: value::age = 0.05
#>   2: value::I(age^2) = 0
#> Chisq(2) = 605.194    Pr(>Chisq) = < 1e-8
#> --------------------------------------------

# Testing age and educ both equal .05
waldtest(fit_gam_ur, coef_names = c("value::age", "value::educ"), rhs = .05)
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#> Restrictions:
#>   1: value::age = 0.05
#>   2: value::educ = 0.05
#> Chisq(2) = 1.131    Pr(>Chisq) = 0.5682
#> --------------------------------------------
```

This last test allows us to illustrate the use of yet another argument
that allows you to setup the restrictions to be tested: `rest_matrix`,
the actual **matrix of restrictions**.

This matrix becomes useful when you want to test linear combinations of
the coefficients. To do that you have to setup a vector of the values
multiplying each of the coefficients in each restriction. Let us
illustrate. The last test makes us want to test if
`value::age = value::educ`. This is equivalent to
`value::age - value::educ = 0`. Even though you don’t see them this is
the same as `1 * value::age - 1 * value::educ = 0`, and we could extend
this to include all other parameters (including `scale::lnnu`) adding
them but pre-multiplied by zero. So that means that the coefficients for
all other parameters are zero, for age it’s 1 and for educ it’s -1. Then
all you have to do is form the row vector in matrix form with the
coefficients for the parameters in place:

``` r

# One restriction, one row in the matrix with the vector.
R <- rbind(c(0, 1, 0, 0, -1, 0, 0, 0, 0))

# Do the test (rhs is 0, which is the default)
waldtest(fit_gam_ur, rest_matrix = R)
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#> Restrictions:
#>   1: value::age - value::educ = 0
#> Chisq(1) = 0.273    Pr(>Chisq) = 0.601
#> --------------------------------------------

# Test both that they're both equal and that age equals .05
R <- rbind(R,                             # age - value (from before)
           c(0, 1, 0, 0, 0, 0, 0, 0, 0))  # only age

# Don't forget the rhs are different
waldtest(fit_gam_ur, rest_matrix = R, rhs = c(0, .05))
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#> Restrictions:
#>   1: value::age - value::educ = 0
#>   2: value::age = 0.05
#> Chisq(2) = 1.131    Pr(>Chisq) = 0.5682
#> --------------------------------------------
```

We should notice that the Chi-squared statistic of the first of these
two last tests, since there is only one restriction tested, is the
squared of the *z* statistic you would get for that same test. The
*p*-value is identical, since it’s the same test.

Now we illustrate how to replicate the Wald test of joint significance
of the parameters in the value equation:

``` r

# Get the indices of the coefficients for the value equation and drop the intercept.
idx <- grep("value::", names(coef(fit_gam_ur)))[-1]

waldtest(fit_gam_ur, indices = idx)
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Original Information Matrix
#> ---------------------------------------
#> Restrictions:
#>   1: value::age = 0
#>   2: value::I(age^2) = 0
#>   3: value::huswage = 0
#>   4: value::educ = 0
#>   5: value::unem = 0
#>   6: value::kidslt6 = 0
#>   7: value::kidsge6 = 0
#> Chisq(7) = 677.722    Pr(>Chisq) = < 1e-8
#> --------------------------------------------
```

You can see that the Chi-squared and the *p*-value match those of the
overall test of significance in the results of the model. The reason is
that it is the same test, because the scale is constant.

This last example illustrates the convenience of having the prefixes in
the names of the coefficient, so you can extract them by equation.

Finally, we just show how to do the same test with robust standard
errors, to illustrate that you can do the Wald test with any
asymptotically correct variance for your models:

``` r

# robust (sandwich) variance
waldtest(fit_gam_ur, indices = idx, vcov.type = "robust")
#> 
#> Wald Test of Linear Restrictions
#> 
#> Variance type: Robust
#> ---------------------------------------
#> Restrictions:
#>   1: value::age = 0
#>   2: value::I(age^2) = 0
#>   3: value::huswage = 0
#>   4: value::educ = 0
#>   5: value::unem = 0
#>   6: value::kidslt6 = 0
#>   7: value::kidsge6 = 0
#> Chisq(7) = 345.310    Pr(>Chisq) = < 1e-8
#> --------------------------------------------
```

## Testing Non-nested Models

What if we want to compare two models that are not nested? Vuong (1989)
designed such test. We can use it to test different specifications that
are not nested, or even across models. The only issue is that Vuong’s
test is not robust to the information matrix equality not holding, so
you must be fairly certain that the models are well specified.

Since the analytical versions of the information test rejected that the
linear model was well specified, even though you can probably do
inference with robust standard errors with it, we cannot really use
Vuong’s test to compare it to the restricted Gamma model we fit. For
that reason, we fit a lognormal model, with the same predictors for the
value function as that of the linear model.

``` r

fit_ln <- ml_lm(log(incthou) ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

IMtest(fit_ln, method = "boot_quad", repetitions = 100)
#> Information Matrix Test
#>  Method: Orthogonalized Quadratic Form + Model-based bootstrap 
#>  Model:  Homoskedastic Lognormal Model 
#> --------------------------------------------
#>  Repetitions: Total 100 Successful 97 
#>  Chisq(28) = 26.141
#>  P(>Chisq): Analytical   = 0.5653 
#>             Bootstrapped = 0.6907
#> --------------------------------------------

vuongtest(fit_gam, fit_ln)
#> 
#> Vuong's (1989) Test
#> --------------------------------------------------
#>   Model 1: Homoskedastic Gamma Model 
#>   Model 2: Homoskedastic Lognormal Model 
#> --------------------------------------------------
#>   z-stat:  1.592
#>   p-value: 0.1114
#> --------------------------------------------------
#>  Inconclusive test: neither model is preferred.
```

We see that the quadratic orthogonalized version of the information test
and its bootstrapped *p*-value, cannot reject that the lognormal model
is well specified. Therefore, we perform Vuong’s test against the
restricted gamma model we estimated first.

The test is inconclusive and shows that neither model is preferred. If
you’re interested in a deeper comparison between the lognormal and Gamma
models, see the [Gamma versus
Lognormal](https://alfisankipan.github.io/mlmodels/articles/mlmodels-gamma-lognormal.md)
vignette.

## Tests for Count Data Models

The last couple of models that we have implemented in the `mlmodels`
package, are specific for count data models. The first one is a test of
overdispersion proposed by Cameron and Trivedi (1990). The null
hypothesis is that there is no overdispersion, and we test against the
two forms of dispersion that our negative binomial estimator models: NB1
and NB2.

``` r

data("docvis")
fit_poi <- ml_poisson(docvis ~ private + medicaid + age + I(age^2) + educyr + actlim + totchr,
                      data = docvis)

OVDtest(fit_poi)
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

We see that the test rejects the null hypothesis of equal dispersion
(Poisson) against either of the two versions of dispersion. You should
be encouraged to estimate a negative binomial model with any of the two
dispersion types: `mlmodels` offers you both in
[`ml_negbin()`](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md),
as well as modeling the actual dispersion parameter for **both**
dispersion types, through `scale`.

The second test that we offer for models of count data, is a
goodness-of-fit test proposed by Manjón and Martínez (2014). This test
allows you to select the bins you want to form, to compare the observed
proportion of observations whose counts fall in the bin, to the
predicted probabilities of the count falling in that bin, and performs
the test. Because we report the Pearson statistic for every bin, it not
only allows you to test if the model predicts the probabilites well or
badly, but also to identify with what particular bins the model has
trouble with.

``` r

# Default bins are individual values from 0 to 5.
GOFtest(fit_poi)
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
```

Clearly, the Poisson model has trouble with most of the bins, due to the
overdispersion, but where it really fails is in predicting 0 and 1
counts, as you can see from the Pearson statistic we include.

We don’t estimate the negative binomial models in this vignette and test
them, but we do that in the [Introduction to Count
Data](https://alfisankipan.github.io/mlmodels/articles/mlmodels-countintro.md)
vignette, that we encourage you to read if you want to learn more about
what you can do with `mlmodels` count data models.

## Concluding Remarks

The `mlmodels` package provides a comprehensive suite of diagnostic
tools designed to help you assess model adequacy and choose between
competing specifications:

- [`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md)
  – Tests for general model misspecification through the Information
  Matrix equality.
- [`lrtest()`](https://alfisankipan.github.io/mlmodels/reference/lrtest.md)
  – Compares nested models using the likelihood ratio test.
- [`waldtest()`](https://alfisankipan.github.io/mlmodels/reference/waldtest.md)
  — Tests linear hypotheses about the parameters, which also allow you
  to test nested models from the unrestricted one.
- [`vuongtest()`](https://alfisankipan.github.io/mlmodels/reference/vuongtest.md)
  – Compares non-nested models based on their relative fit.
- [`OVDtest()`](https://alfisankipan.github.io/mlmodels/reference/OVDtest.md)
  – Tests for overdispersion in Poisson models against both NB1 and NB2
  alternatives.
- [`GOFtest()`](https://alfisankipan.github.io/mlmodels/reference/GOFtest.md)
  – Evaluates how well a count model reproduces the observed frequency
  distribution of counts.

While this collection is not exhaustive of all possible diagnostic
procedures available in the literature, it is **flexible and powerful
enough** to address the majority of practical diagnostic needs when
working with the models included in `mlmodels`.

Used together with the rich set of variance-covariance estimators and
post-estimation tools (predictions and marginal effects), these
diagnostics give you a solid foundation for responsible and transparent
modeling.

Happy modeling!

## References

Cameron, A. C., & Trivedi, P. K. (1990). Regression-based tests for
overdispersion in the Poisson model. *Journal of Econometrics*, 46(3),
347-364. <https://doi.org/10.1016/0304-4076(90)90014-K>

Chesher, A. (1983). The information matrix test: Simplified calculation
via a score test interpretation. *Econometric Letters*, 13(1), 45–48.
<https://doi.org/10.1016/0165-1765(83)90009-5>

Davidson, R., & MacKinnon, J. G. (1992). A New Form of the Information
Matrix Test. *Econometrica*, 60(1), 145–157.
<https://doi.org/10.2307/2951680>

Davidson, R., & MacKinnon, J. G. (1999). The Size Distortion of
Bootstrap Tests. *Econometric Theory*, 15(3), 361-376.
<https://www.jstor.org/stable/3533339>

Horowitz, J. L. (1994). Bootstrap-based critical values for the
information matrix test. *Journal of Econometrics*, 61(2), 395–411.
<https://doi.org/10.1016/0304-4076(94)90092-2>

Lancaster, T. (1984). The covariance matrix of the information matrix
test. *Econometrica*, 52(4), 1051–1053.
<https://doi.org/10.2307/1911198>

Manjón, M., & Martínez, O. (2014). The chi-squared goodness-of-fit test
for count-data models. *The Stata Journal*, 14(4), 798–816.
<https://doi.org/10.1177/1536867X1401400406>

Vuong, Q. H. (1989). Likelihood Ratio Tests for Model Selection and
Non-Nested Hypotheses. *Econometrica*, 57(2), 307-333.
<https://doi.org/10.2307/1912557>

White, H. (1982). Maximum likelihood estimation of misspecified models.
*Econometrica*, 50(1), 1–25. <https://doi.org/10.2307/1912526>
