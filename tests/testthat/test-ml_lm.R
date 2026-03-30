# tests/testthat/test-ml_lm.R

library(testthat)

test_that("ml_lm fits homoskedastic model correctly", {
  data(mtcars)

  fit <- ml_lm(mpg ~ wt + hp + qsec, data = mtcars)

  expect_s3_class(fit, "ml_lm")
  expect_s3_class(fit, "mlmodel")
  expect_true(!is.null(fit$model$fitted.values))
  expect_true(!is.null(fit$model$residuals))
  expect_true(is.null(fit$model$scale_formula))
})

test_that("ml_lm fits heteroskedastic model correctly", {
  data(mtcars)

  fit <- ml_lm(value = mpg ~ wt + hp,
               scale = ~ wt,
               data = mtcars)

  expect_s3_class(fit, "ml_lm")
  expect_true(!is.null(fit$model$scale_formula))
})

test_that("ml_lm removes NAs correctly", {
  data(mtcars)
  mtcars[1:3, "mpg"] <- NA

  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)

  expect_equal(fit$model$n_used, nrow(mtcars) - 3)
  expect_equal(length(fit$model$fitted.values), fit$model$n_used)
})
