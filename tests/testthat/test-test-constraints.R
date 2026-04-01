# tests/testthat/test-constraints.R

library(testthat)

test_that("parse_constraints works with simple equality constraints", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, scale = ~ wt, data = mtcars)
  coef_names <- names(coef(fit))

  # Test simple constraints
  cons <- c("value::wt = 0",
            "scale::(Intercept) = -2.5",
            "value::hp + value::wt = 0")

  result <- .parse_constraints(cons, coef_names)

  expect_type(result, "list")
  expect_true("eqA" %in% names(result))
  expect_true("eqB" %in% names(result))

  # Check dimensions
  expect_equal(nrow(result$eqA), 3)
  expect_equal(ncol(result$eqA), length(coef_names))
})

test_that("parse_constraints handles scaled coefficients", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, scale = ~ wt, data = mtcars)
  coef_names <- names(coef(fit))

  cons <- c("value::hp * 2 = 0",
            "value::wt / 5 >= 0.1")

  result <- .parse_constraints(cons, coef_names)

  expect_true("eqA" %in% names(result))
  expect_true("ineqA" %in% names(result))
})

test_that("parse_constraints throws nice error on invalid input", {
  data(mtcars)
  fit <- ml_lm(mpg ~ wt + hp, data = mtcars)
  coef_names <- names(coef(fit))

  expect_error(
    .parse_constraints("value::wt = abc", coef_names),
    "Right-hand side must be a number"
  )

  expect_error(
    .parse_constraints("value::wt = value::hp", coef_names),
    "Right-hand side must be a number"
  )
})
