# =============================================================================
# Tests for ml_logit
# =============================================================================

library(testthat)

# Basic fitting
test_that("ml_logit fits homoskedastic and heteroskedastic models", {
  data(mtcars)
  fit_hom <- ml_logit(mpg > 20 ~ wt + hp, data = mtcars)
  fit_het <- ml_logit(mpg > 20 ~ wt + hp, scale = ~ wt, data = mtcars)
  
  expect_s3_class(fit_hom, "ml_logit")
  expect_s3_class(fit_het, "ml_logit")
  expect_true(!is.null(fit_het$model$scale))
  expect_true(is.null(fit_hom$model$scale))  # homoskedastic
})

# Summary and print
test_that("summary.ml_logit produces expected output", {
  data(mtcars)
  fit <- ml_logit(mpg > 20 ~ wt + hp, scale = ~ wt, data = mtcars)
  s <- summary(fit)
  expect_s3_class(s, "summary.ml_logit")
  expect_true(!is.null(s$r.squared$cor))
  expect_true(!is.null(s$AIC))
})

# Predict with all types
test_that("predict.ml_logit supports all types and se.fit", {
  data(mtcars)
  fit <- ml_logit(mpg > 20 ~ wt + hp, scale = ~ wt, data = mtcars)
  
  types <- c("response", "prob", "fitted", "prob0", "link", "odds", "xb", "sigma", "variance", "zd")
  
  for (t in types) {
    p <- predict(fit, type = t)
    expect_type(p, "double")
    expect_length(p, nrow(mtcars))
  }
  
  # se.fit
  p_se <- predict(fit, type = "response", se.fit = TRUE)
  expect_named(p_se, c("fit", "se.fit"))
})

# vcov types
test_that("vcov.ml_logit supports all variance types", {
  library(wooldridge)
  data("smoke")
  
  fit <- ml_logit(cigs > 0 ~ cigpric + income + age,
                   data = smoke)
  
  suppressWarnings({
    v <- vcov(fit, type = "oim")
    v_rob <- vcov(fit, type = "robust")
    v_boot <- vcov(fit, type = "boot", repetitions = 50, progress = FALSE)
    v_jack <- vcov(fit, type = "jack", progress = FALSE)
  })
  
  expect_true(is.matrix(v))
  expect_true(is.matrix(v_rob))
  expect_true(is.matrix(v_boot))
  expect_true(is.matrix(v_jack))
})

# Constraints
test_that("ml_logit supports constraints", {
  data(mtcars)
  fit <- ml_logit(mpg > 20 ~ wt + hp, scale = ~ wt,
                  constraints = c("value::wt = -0.5"),
                  start = c(0, -0.5, 0, 0),
                  data = mtcars)
  expect_true(!is.null(fit$model$constraints))
})

# NA handling and subset
test_that("ml_logit handles NAs and subset", {
  data(mtcars)
  mtcars$mpg[1:5] <- NA
  fit <- ml_logit(mpg > 20 ~ wt + hp, data = mtcars)
  expect_equal(length(predict(fit)), nrow(mtcars))
})

# Marginaleffects compatibility
test_that("ml_logit works with marginaleffects", {
  library(wooldridge)
  data("smoke")
  
  fit <- ml_logit(cigs > 0 ~ cigpric + income + age,
                   data = smoke)
  expect_silent(predictions(fit))
  expect_silent(avg_slopes(fit, variables = "income"))
})