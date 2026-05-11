# mlmodels

**Maximum Likelihood Estimation with a Consistent S3 Interface**

[![CRAN status](https://www.r-pkg.org/badges/version/mlmodels)](https://cran.r-project.org/package=mlmodels)
[![CRAN monthly downloads](https://cranlogs.r-pkg.org/badges/mlmodels)](https://cran.r-project.org/package=mlmodels)
[![CRAN total downloads](https://cranlogs.r-pkg.org/badges/grand-total/mlmodels)](https://cran.r-project.org/package=mlmodels)

## Current CRAN Version: 0.1.2

The package is now available on CRAN (0.1.2). It includes:
- A full suite of maximum likelihood estimators (Gaussian, logit, probit, Poisson, negative binomial, gamma, and beta)
- Flexible modeling of the scale (variance, dispersion, shape, precision) parameter
- Multiple variance-covariance estimators
- Hypothesis tests (Wald, LR, Information Matrix, Vuong, overdispersion, and goodness-of-fit)

## Development Version (main branch)

We are currently working on:
- **Truncated models** (`ml_trunc_lm()`, `ml_trunc_poisson()`, etc.)
- Implementing Clarke's test for non-nested models.
- Implementing bootstrapping in Vuong's and Clarke's test.
- Improvements to estimation, and displaying statistics with weighted models.

## Documentation

- [Getting started](https://alfisankipan.github.io/mlmodels/articles/mlmodels-basics.html) — basic usage, model specification, and examples.
- [Reference manual](https://alfisankipan.github.io/mlmodels/reference/index.html) — complete function documentation.
- [NEWS](https://alfisankipan.github.io/mlmodels/news/index.html) — changelog.

## Installation

```r
# From CRAN (stable)
install.packages("mlmodels")

# Latest development version
remotes::install_github("alfisankipan/mlmodels")
```