# =============================================================================
# Tests for ml_probit
# =============================================================================

library(testthat)

# Basic fitting
test_that("ml_probit fits homoskedastic and heteroskedastic models", {
  library(wooldridge)
  data("smoke")
  
  fit_hom <- ml_probit(cigs > 0 ~ cigpric + income + age,
                   data = smoke)
  fit_het <- ml_probit(cigs > 0 ~ cigpric + income + age,
                       scale = ~ educ,
                       data = smoke)
  
  expect_s3_class(fit_hom, "ml_probit")
  expect_s3_class(fit_het, "ml_probit")
  expect_true(!is.null(fit_het$model$scale))
  expect_true(is.null(fit_hom$model$scale))
})

# Summary and print
test_that("summary.ml_probit produces expected output", {
  library(wooldridge)
  data("smoke")
  
  fit <- ml_probit(cigs > 0 ~ cigpric + income + age,
                   data = smoke)
  s <- summary(fit)
  expect_s3_class(s, "summary.ml_probit")
  expect_true(!is.null(s$r.squared$cor))
  expect_true(!is.null(s$r.squared$mczav))
  expect_true(!is.null(s$AIC))
})

# Predict with all types
test_that("predict.ml_probit supports all types and se.fit", {
  library(wooldridge)
  data("smoke")
  
  fit <- ml_probit(cigs > 0 ~ cigpric + income + age,
                   data = smoke)
  
  types <- c("response", "prob", "fitted", "prob0", "link", "odds", "xb", "sigma", "variance", "zd")
  
  for (t in types) {
    p <- predict(fit, type = t)
    expect_type(p, "double")
    expect_length(p, nrow(smoke))
  }
  
  p_se <- predict(fit, type = "response", se.fit = TRUE)
  expect_named(p_se, c("fit", "se.fit"))
})

# vcov types
test_that("vcov.ml_probit supports all variance types", {
  library(wooldridge)
  data("smoke")
  
  fit <- ml_probit(cigs > 0 ~ cigpric + income + age,
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
test_that("ml_probit supports constraints", {
  data(mtcars)
  fit <- ml_probit(mpg > 20 ~ wt + hp, scale = ~ wt,
                   constraints = c("value::wt = -0.5"),
                   start = c(0, -0.5, 0, 0),
                   data = mtcars)
  expect_true(!is.null(fit$model$constraints))
})

# NA handling and subset
test_that("ml_probit handles NAs and subset", {
  library("wooldridge")
  data("smoke")
  smoke$cigs[1:5] <- NA
  fit <- ml_probit(cigs > 0 ~ cigpric + income + age,
                   data = smoke)
  expect_equal(length(predict(fit)), nrow(smoke))
})

# Marginaleffects compatibility
test_that("ml_probit works with marginaleffects", {
  library(marginaleffects)
  library(wooldridge)
  data("smoke")
  
  fit <- ml_probit(cigs > 0 ~ cigpric + income + age,
                   data = smoke)
  expect_silent(predictions(fit))
  expect_silent(avg_slopes(fit, variables = "income"))
})