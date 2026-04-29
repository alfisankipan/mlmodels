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
ratio, information matrix). It is fully compatible with marginaleffects
for post-estimation analysis.

## Author

**Maintainer**: Alfonso Sanchez-Penalver <oneiros_spain@yahoo.com>
([ORCID](https://orcid.org/0000-0001-8491-4632))
