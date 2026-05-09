# mlmodels

**Maximum Likelihood Estimation with a Consistent S3 Interface**

[![CRAN status](https://www.r-pkg.org/badges/version/mlmodels)](https://cran.r-project.org/package=mlmodels)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/mlmodels)](https://cran.r-project.org/package=mlmodels)

## Current CRAN Version: 0.1.2

The package is now available on CRAN (0.1.2). It includes:
- A full suite of maximum likelihood estimators (Gaussian, logit, probit, Poisson, negative binomial, gamma, beta, etc.)
- Flexible modeling of the scale/dispersion parameter
- Multiple variance-covariance estimators
- Hypothesis tests (Wald, LR, Information Matrix, Vuong)

## Development Version (main branch)

We are currently working on:
- **Truncated models** (`ml_trunc_lm()`, `ml_trunc_poisson()`, etc.)
- Further improvements to the test suite and documentation

---

## Installation

```r
# From CRAN (stable)
install.packages("mlmodels")

# Latest development version
remotes::install_github("alfisankipan/mlmodels")
```