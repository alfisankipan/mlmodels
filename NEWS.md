# mlmodels 0.2.0

## New Features

* **Clarke test** (`clarketest()`): Full implementation of Clarke's (2007) 
  distribution-free test for non-nested models. Includes both binomial/sign 
  and signed-rank (Wilcoxon) variants, three penalty options 
  (`none`, `akaike`, `schwarz`), model-based bootstrap, and an informative 
  print method with expected values under the null and robustness notes.

* **Vuong test** (`vuongtest()`): Now supports the same three penalty options 
  (`none`, `akaike`, `schwarz`) for consistency with Clarke's test and 
  Vuong (1989). Bootstrap functionality has been added and aligned with the 
  rest of the test suite.

* **Wald test** (`waldtest()`): Major upgrade — new `constraints` argument 
  (replaces `coef_names`) that supports arbitrary linear combinations of 
  coefficients.

## Improvements

* Much better handling of weights across all tests and estimators. Statistics 
  are now properly scaled to ensure correct asymptotic behavior.
* Enhanced `logLik()`, `AIC()`, and `BIC()` methods to return effective or 
  scaled statistics depending on user request.
* Improved interaction with `maxLik` optimizer, leading to more reliable 
  convergence in many models.
* Enhanced print methods for hypothesis tests (clearer conclusions, reference 
  values under the null, bootstrap robustness warnings).
* Better print methods for summaries of weighted models, showing both 
  effective and scaled (to sample size) measures where appropriate.
* Updated documentation and examples for non-nested model comparison.

## Bug fixes & minor changes

* Fixed bootstrap storage (numeric vectors instead of character).
* Various internal cleanups and consistency improvements across the test suite.

# mlmodels 0.1.2

* Added return values in the documentation of exported functions that were missing them.
* Added references to implemented methods in the description.

# mlmodels 0.1.1

* Fixed weighted log-likelihood calculation in `ml_logit()` (both homoskedastic and heteroskedastic versions). 
  This bug previously caused incorrect log-likelihood, AIC, BIC, and convergence issues in weighted logit models.
* All other models were already handling weights correctly.

# mlmodels 0.1.0

* First public release - Initial CRAN submission
* Provides maximum likelihood estimation for Gaussian (linear and log-normal), logit, probit, Poisson, negative binomial (NB1 and NB2), gamma, and beta models.
* Consistent S3 interface with support for modeling scale parameters.
* Multiple variance-covariance estimators (OIM, OPG, robust, cluster-robust, bootstrap, jackknife).
* Full suite of post-estimation tools and hypothesis tests.
* Compatible with `marginaleffects` for marginal effects and predictions.
* Comprehensive vignettes covering main model families and diagnostics.