# mlmodels: Maximum Likelihood Models and Tools for Estimation, Prediction, and Testing

Provides a collection of maximum likelihood estimators with a consistent
S3 interface. Supported models include Gaussian (linear and log-normal),
logit, probit, Poisson, negative binomial (NB1 and NB2), gamma, and beta
regression. A distinctive feature is flexible modeling of the scale
parameter (variance, dispersion, precision, or shape) alongside the
location/mean parameters. The package offers unified predict() methods,
multiple variance-covariance estimators (observed information, outer
product of gradients, robust/Huber-White, cluster-robust, bootstrap,
jackknife), and a full suite of hypothesis tests (Wald, likelihood
ratio, information matrix, Vuong, overdispersion, and goodness-of-fit).
It is fully compatible with 'marginaleffects' for post-estimation
analysis. Methods implemented include Cameron and Trivedi (1990)
[doi:10.1016/0304-4076(90)90014-K](https://doi.org/10.1016/0304-4076%2890%2990014-K)
, for Poisson overdispersion testing, Manjon and Martinez (2014)
[doi:10.1177/1536867X1401400406](https://doi.org/10.1177/1536867X1401400406)
, for goodness-of-fit testing of count data models, Vuong (1989)
[doi:10.2307/1912557](https://doi.org/10.2307/1912557) , for non-nested
likelihood ratio testing, and White (1982)
[doi:10.2307/1912526](https://doi.org/10.2307/1912526) , for information
matrix tests.

## See also

Useful links:

- <https://alfisankipan.github.io/mlmodels/>

- Report bugs at <https://github.com/alfisankipan/mlmodels/issues>

## Author

**Maintainer**: Alfonso Sanchez-Penalver <oneiros_spain@yahoo.com>
([ORCID](https://orcid.org/0000-0001-8491-4632))
