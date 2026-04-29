# Package index

## Model Fitting

Main model estimation functions

- [`ml_beta()`](https://alfisankipan.github.io/mlmodels/reference/ml_beta.md)
  : Fit Beta Model by Maximum Likelihood
- [`ml_gamma()`](https://alfisankipan.github.io/mlmodels/reference/ml_gamma.md)
  : Fit Gamma Model by Maximum Likelihood
- [`ml_lm()`](https://alfisankipan.github.io/mlmodels/reference/ml_lm.md)
  : Fit linear model by Maximum Likelihood
- [`ml_logit()`](https://alfisankipan.github.io/mlmodels/reference/ml_logit.md)
  : Fit Binary Logit Model by Maximum Likelihood
- [`ml_negbin()`](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md)
  : Fit negative binomial models by Maximum Likelihood
- [`ml_poisson()`](https://alfisankipan.github.io/mlmodels/reference/ml_poisson.md)
  : Fit Poisson model by Maximum Likelihood
- [`ml_probit()`](https://alfisankipan.github.io/mlmodels/reference/ml_probit.md)
  : Fit Binary Probit Model by Maximum Likelihood

## Post-Estimation

Prediction, summaries, and marginal effects

- [`predict(`*`<ml_beta>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_gamma>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_lm>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_logit>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_negbin>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_poisson>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_probit>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  : Predictions for mlmodel models
- [`summary(`*`<ml_beta>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_gamma>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_lm>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_logit>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_negbin>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_poisson>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_probit>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  : Summary for mlmodel objects
- [`coef(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/coef.mlmodel.md)
  : Extract Model Coefficients
- [`vcov(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md)
  : Variance-Covariance Matrix for mlmodel Objects
- [`logLik(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/logLik.mlmodel.md)
  [`logLik(`*`<summary.mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/logLik.mlmodel.md)
  : Extract Log-Likelihood from mlmodel objects
- [`fitted(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/fitted.mlmodel.md)
  [`fitted(`*`<values.mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/fitted.mlmodel.md)
  : Extract Fitted Values from mlmodel
- [`residuals(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/residuals.mlmodel.md)
  : Extract Model Residuals
- [`confint(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/confint.mlmodel.md)
  : Confidence Intervals for mlmodel Coefficients
- [`update(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/update.mlmodel.md)
  [`update(`*`<ml_poisson>`*`)`](https://alfisankipan.github.io/mlmodels/reference/update.mlmodel.md)
  : Update an mlmodel Call

## Hypothesis Testing

Model comparison and specification tests

- [`waldtest()`](https://alfisankipan.github.io/mlmodels/reference/waldtest.md)
  : Wald Test for Linear Restrictions
- [`lrtest()`](https://alfisankipan.github.io/mlmodels/reference/lrtest.md)
  : Likelihood Ratio Test for Nested mlmodel Objects
- [`vuongtest()`](https://alfisankipan.github.io/mlmodels/reference/vuongtest.md)
  : Vuong's Test for Non-Nested Models
- [`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md)
  : Information Matrix Test for Model Misspecification
- [`OVDtest()`](https://alfisankipan.github.io/mlmodels/reference/OVDtest.md)
  : Overdispersion Tests for Count Models
- [`GOFtest()`](https://alfisankipan.github.io/mlmodels/reference/GOFtest.md)
  : Goodness-of-Fit Test for Count Models

## Internal & Helper Functions

Not usually called directly by users

- [`AIC(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/AIC.mlmodel.md)
  [`AIC(`*`<summary.mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/AIC.mlmodel.md)
  : Extract AIC from mlmodel objects
- [`BIC(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/BIC.mlmodel.md)
  [`BIC(`*`<summary.mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/BIC.mlmodel.md)
  : Extract BIC from mlmodel objects
- [`GOFtest()`](https://alfisankipan.github.io/mlmodels/reference/GOFtest.md)
  : Goodness-of-Fit Test for Count Models
- [`IMtest()`](https://alfisankipan.github.io/mlmodels/reference/IMtest.md)
  : Information Matrix Test for Model Misspecification
- [`OVDtest()`](https://alfisankipan.github.io/mlmodels/reference/OVDtest.md)
  : Overdispersion Tests for Count Models
- [`coef(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/coef.mlmodel.md)
  : Extract Model Coefficients
- [`confint(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/confint.mlmodel.md)
  : Confidence Intervals for mlmodel Coefficients
- [`docvis`](https://alfisankipan.github.io/mlmodels/reference/docvis.md)
  : U.S. Medical Expenditure Panel Survey
- [`find_predictors(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/find_predictors.mlmodel.md)
  : Extract the predictors used in the model (for
  insight/marginaleffects compatibility)
- [`find_variables.mlmodel()`](https://alfisankipan.github.io/mlmodels/reference/find_variables.mlmodel.md)
  : Extract the variables used in the model (for insight/marginaleffects
  compatibility)
- [`fitted(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/fitted.mlmodel.md)
  [`fitted(`*`<values.mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/fitted.mlmodel.md)
  : Extract Fitted Values from mlmodel
- [`formula(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/formula.mlmodel.md)
  : Extract value formula from mlmodel objects
- [`get_data(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/get_data.mlmodel.md)
  [`get_modeldata.mlmodel()`](https://alfisankipan.github.io/mlmodels/reference/get_data.mlmodel.md)
  : Extract data used to fit the model (for insight/marginaleffects
  compatibility)
- [`gradientObs()`](https://alfisankipan.github.io/mlmodels/reference/gradientObs.md)
  : Gradient (Score) by Observation
- [`hessianObs()`](https://alfisankipan.github.io/mlmodels/reference/hessianObs.md)
  : Hessian by Observation
- [`logLik(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/logLik.mlmodel.md)
  [`logLik(`*`<summary.mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/logLik.mlmodel.md)
  : Extract Log-Likelihood from mlmodel objects
- [`loglikeObs()`](https://alfisankipan.github.io/mlmodels/reference/loglikeObs.md)
  : Log-Likelihood by Observation
- [`lrtest()`](https://alfisankipan.github.io/mlmodels/reference/lrtest.md)
  : Likelihood Ratio Test for Nested mlmodel Objects
- [`ml_beta()`](https://alfisankipan.github.io/mlmodels/reference/ml_beta.md)
  : Fit Beta Model by Maximum Likelihood
- [`ml_gamma()`](https://alfisankipan.github.io/mlmodels/reference/ml_gamma.md)
  : Fit Gamma Model by Maximum Likelihood
- [`ml_lm()`](https://alfisankipan.github.io/mlmodels/reference/ml_lm.md)
  : Fit linear model by Maximum Likelihood
- [`ml_logit()`](https://alfisankipan.github.io/mlmodels/reference/ml_logit.md)
  : Fit Binary Logit Model by Maximum Likelihood
- [`ml_negbin()`](https://alfisankipan.github.io/mlmodels/reference/ml_negbin.md)
  : Fit negative binomial models by Maximum Likelihood
- [`ml_poisson()`](https://alfisankipan.github.io/mlmodels/reference/ml_poisson.md)
  : Fit Poisson model by Maximum Likelihood
- [`ml_probit()`](https://alfisankipan.github.io/mlmodels/reference/ml_probit.md)
  : Fit Binary Probit Model by Maximum Likelihood
- [`mroz`](https://alfisankipan.github.io/mlmodels/reference/mroz.md) :
  University of Michigan Panel Study of Income Dynamics (PSID)
- [`nobs(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/nobs.mlmodel.md)
  : Extract the Number of Observations from an mlmodel
- [`` `%||%` ``](https://alfisankipan.github.io/mlmodels/reference/null-default.md)
  : Null default operator
- [`predict(`*`<ml_beta>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_gamma>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_lm>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_logit>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_negbin>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_poisson>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  [`predict(`*`<ml_probit>`*`)`](https://alfisankipan.github.io/mlmodels/reference/predict.mlmodel.md)
  : Predictions for mlmodel models
- [`pw401k`](https://alfisankipan.github.io/mlmodels/reference/pw401k.md)
  : 401(k) Participation Rates
- [`residuals(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/residuals.mlmodel.md)
  : Extract Model Residuals
- [`se()`](https://alfisankipan.github.io/mlmodels/reference/se.md) :
  Extract Standard Errors from mlmodel Objects
- [`smoke`](https://alfisankipan.github.io/mlmodels/reference/smoke.md)
  : 1979 National Health Interview Survey
- [`summary(`*`<ml_beta>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_gamma>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_lm>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_logit>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_negbin>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_poisson>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  [`summary(`*`<ml_probit>`*`)`](https://alfisankipan.github.io/mlmodels/reference/summary.mlmodel.md)
  : Summary for mlmodel objects
- [`update(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/update.mlmodel.md)
  [`update(`*`<ml_poisson>`*`)`](https://alfisankipan.github.io/mlmodels/reference/update.mlmodel.md)
  : Update an mlmodel Call
- [`vcov(`*`<mlmodel>`*`)`](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md)
  : Variance-Covariance Matrix for mlmodel Objects
- [`vuongtest()`](https://alfisankipan.github.io/mlmodels/reference/vuongtest.md)
  : Vuong's Test for Non-Nested Models
- [`waldtest()`](https://alfisankipan.github.io/mlmodels/reference/waldtest.md)
  : Wald Test for Linear Restrictions
- [`mlmodels`](https://alfisankipan.github.io/mlmodels/reference/mlmodels-package.md)
  [`mlmodels-package`](https://alfisankipan.github.io/mlmodels/reference/mlmodels-package.md)
  : mlmodels: Maximum Likelihood Models and Tools for Estimation,
  Prediction, and Testing
