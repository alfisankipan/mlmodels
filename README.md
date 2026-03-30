# mlmodels: Maximum Likelihood Models for R

<!-- badges: start -->
[![R-CMD-check](https://github.com/alfisankipan/mlmodels/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alfisankipan/mlmodels/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**mlmodels** provides a clean, consistent interface for maximum-likelihood estimation in R.
All models share the same `mlmodel` class and support:

- Flexible formulas (including heteroskedasticity in the scale equation)
- Rich `predict()` methods (response, mean, median, sigma, variance, etc.)
- Multiple variance-covariance estimators (OIM, OPG, robust, cluster-robust, bootstrap, cluster-bootstrap)
- Wald, LR, and Information-Matrix tests
- Full **marginaleffects** compatibility (`avg_slopes()`, `avg_predictions()`, etc.)

## Installation

```r
# Development version from GitHub
devtools::install_github("alfisankipan/mlmodels")
