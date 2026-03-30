#' mlmodels: Maximum Likelihood Models for R
#'
#' @description
#' A collection of maximum-likelihood estimators with a consistent S3/S4
#' interface. Currently includes Gaussian linear models (homoskedastic and
#' heteroskedastic) with full support for log-transformed outcomes and
#' parametric retransformation. Provides unified `predict()`, multiple
#' variance-covariance estimators (OIM, OPG, robust, cluster, bootstrap,
#' cluster-bootstrap), and a full suite of hypothesis tests. Fully compatible
#' with **marginaleffects** for average marginal effects and slopes.
#'
#' @details
#' All models inherit from the base class `mlmodel` and follow the same
#' workflow:
#'
#' - Fitting: `ml_lm()`, `ml_logit()`, …
#' - Prediction: `predict()` with many types (`response`, `mean`, `median`,
#'   `sigma`, `variance`, `link`, etc.)
#' - Inference: `vcov()`, `summary()`, `waldtest()`, `lrtest()`, `IMtest()`
#' - Bootstrap and cluster-robust standard errors
#' - Marginal effects via `marginaleffects::avg_slopes()` and friends
#'
#' @author Alfonso Sanchez-Penalver
#'
#' @references
#' Sanchez-Penalver, A. (2026). mlmodels: Maximum Likelihood Models for R.
#' R package version 0.1.0.
#'
#' @keywords package
#' @name mlmodels
#' @docType package
#'
#' @importFrom cli cli_abort cli_warn cli_alert_info
#' @importFrom hardhat mold forge default_formula_blueprint
#' @importFrom maxLik maxLik
#' @importFrom rlang arg_match enquo eval_tidy
#' @importFrom tibble tibble
#'
#' @importFrom stats coef
#' @importFrom stats complete.cases
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats cov2cor
#' @importFrom stats dnorm
#' @importFrom stats formula
#' @importFrom stats lm
#' @importFrom stats logLik
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @importFrom stats predict
#' @importFrom stats printCoefmat
#' @importFrom stats update
#' @importFrom stats update.formula
#' @importFrom stats var
#' @importFrom stats vcov
NULL
