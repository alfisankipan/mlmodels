# tests/testthat/test-summary.R

library(testthat)

test_that("summary.ml_lm returns expected structure", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)

  s <- summary(fit)

  expect_s3_class(s, "summary.ml_lm")
  expect_true("coefficients" %in% names(s))
  expect_true("nobs" %in% names(s))
  expect_true("df.residual" %in% names(s))
  expect_true("logLik" %in% names(s))
  expect_true("AIC" %in% names(s))
  expect_true("BIC" %in% names(s))
})

test_that("summary.ml_lm includes joint significance tests", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)

  s <- summary(fit)

  expect_type(s$significance, "list")
  expect_true("all" %in% names(s$significance))
  expect_s3_class(s$significance$all, "waldtest.mlmodel")
})

test_that("print.summary.ml_lm works without error", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)

  expect_output(print(summary(fit)), "Maximum Likelihood Model")
  expect_output(print(summary(fit)), "Joint significance tests")
})

test_that("summary.ml_lm handles heteroskedastic models", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, scale = ~ wt, data = mtcars)

  s <- summary(fit)

  expect_true(s$is_heteroskedastic)
  expect_equal(s$model_type, "Heteroskedastic Gaussian Linear Model")
})
