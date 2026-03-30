# tests/testthat/test-predict.R
library(testthat)

# Basic functionality
test_that("predict.ml_lm returns correct length for in-sample predictions", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp + qsec, data = mtcars)
  pred <- predict(fit)
  expect_equal(length(pred), nrow(mtcars))
  expect_type(pred, "double")
})

test_that("predict.ml_lm works with newdata", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)
  new_data <- mtcars[1:10, ]
  pred <- predict(fit, newdata = new_data)
  expect_equal(length(pred), 10)
})

# Different types
test_that("predict.ml_lm supports different prediction types", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, scale = ~ wt, data = mtcars)
  expect_silent(predict(fit, type = "response"))
  expect_silent(predict(fit, type = "mean"))
  expect_silent(predict(fit, type = "fitted"))
  expect_silent(predict(fit, type = "sigma"))
  expect_silent(predict(fit, type = "sd"))
  expect_silent(predict(fit, type = "variance"))
  expect_silent(predict(fit, type = "link"))
})

test_that("predict.ml_lm returns expected values for different types", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, scale = ~ wt, data = mtcars)
  mu <- predict(fit, type = "mean")
  sigma <- predict(fit, type = "sigma")
  var <- predict(fit, type = "variance")
  link <- predict(fit, type = "link")
  expect_equal(length(mu), nrow(mtcars))
  expect_equal(length(sigma), nrow(mtcars))
  expect_equal(length(var), nrow(mtcars))
  expect_equal(length(link), nrow(mtcars))
  expect_true(all(sigma > 0))
  expect_true(all(var > 0))
})

# Log-transformed outcome (new — critical for the retransformation logic)
test_that("predict.ml_lm correctly handles log(y) models and retransformation", {
  data(mtcars)
  fit_log <- ml_lm(log(mpg) ~ wt + hp, data = mtcars)   # pure log
  fit_shift <- ml_lm(log(mpg + 2) ~ wt, data = mtcars)  # with shift

  # point predictions
  expect_silent(mean_log <- predict(fit_log, type = "mean"))
  expect_silent(med_log <- predict(fit_log, type = "median"))
  expect_silent(var_y <- predict(fit_log, type = "variance_y"))

  expect_true(all(mean_log > 0))
  expect_true(all(med_log > 0))
  expect_true(all(var_y > 0))

  # shift warning
  expect_warning(
    predict(fit_shift, type = "mean"),
    regexp = "Outcome was transformed as log\\(y \\+ d\\)"
  )
})

# se.fit functionality
test_that("predict.ml_lm returns standard errors when se.fit = TRUE", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)
  pred <- predict(fit, se.fit = TRUE)
  expect_type(pred, "list")
  expect_named(pred, c("fit", "se.fit"))
  expect_equal(length(pred$fit), nrow(mtcars))
  expect_equal(length(pred$se.fit), nrow(mtcars))
  expect_true(all(pred$se.fit > 0, na.rm = TRUE))
})

test_that("predict.ml_lm se.fit works for heteroskedastic models and all types", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, scale = ~ wt, data = mtcars)

  # different types should give different SEs (sanity check on gradients)
  se_mean <- predict(fit, type = "mean", se.fit = TRUE)$se.fit
  se_sigma <- predict(fit, type = "sigma", se.fit = TRUE)$se.fit
  se_var <- predict(fit, type = "variance", se.fit = TRUE)$se.fit

  expect_true(all(se_mean > 0))
  expect_true(all(se_sigma > 0))
  expect_true(all(se_var > 0))
  expect_false(isTRUE(all.equal(se_mean, se_sigma)))  # gradients differ
})

# NA handling
test_that("predict.ml_lm handles NA predictions gracefully (in-sample)", {
  data(mtcars)
  mtcars[1:5, "mpg"] <- NA
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)
  pred <- predict(fit)
  expect_equal(sum(is.na(pred)), 5)
  expect_equal(length(pred), nrow(mtcars))
})

test_that("predict.ml_lm propagates NAs in newdata (no dropping message)", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)
  new_data <- mtcars[1:10, ]
  new_data[3:5, "wt"] <- NA
  pred <- expect_silent(predict(fit, newdata = new_data))  # no message
  expect_equal(sum(is.na(pred)), 3)
  expect_equal(length(pred), 10)
})

# Edge cases
test_that("predict.ml_lm returns NA for empty newdata", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)
  empty_data <- mtcars[0, ]
  pred <- predict(fit, newdata = empty_data)
  expect_equal(length(pred), 0)
})
