# Clear everything in the environment.
rm(list = ls(all.names = TRUE))
# Clear console
cat("\f")

library(devtools)
library(wooldridge)
library(sandwich)
library(microbenchmark)

load_all(".")

data("smoke")
data("mtcars", package = "datasets")
data("PetersenCL")

# fit <- ml_lm(mpg ~ wt + hp, data = mtcars, weights = wt)
# v <- vcov(fit, type = "boot", repetitions = 300, seed = 123)
# print(sqrt(diag(v)))
# v <- vcov(fit, type = "robust")
# print(sqrt(diag(v)))


# ols <- lm(y ~ x, data = PetersenCL)
ml <- ml_lm(y ~x, data = PetersenCL)

# set.seed(123)
# boot_sand <- vcovBS(ols, R = 300, cluster = ~ firm)


time_fast <- system.time({
  v_fast <- vcov(ml,
                 type = "boot", repetitions = 500,
                 cl_var = "firm",
                 seed = 123,
                 progress = FALSE)
})

time_general <- system.time({
  v_general <- vcov(ml,
                 type = "boot", repetitions = 500,
                 cl_var = "firm",
                 seed = 123,
                 progress = FALSE)
})


cat("Fast version (vcov_boot.ml_lm)  :", round(time_fast["elapsed"], 3), "seconds\n")
cat("General version (vcov_boot.mlmodel):", round(time_general["elapsed"], 3), "seconds\n")
cat("Speedup factor:", round(time_general["elapsed"] / time_fast["elapsed"], 2), "x faster\n")

# boot_ml <- vcov(ml, type = "boot", repetitions = 500, cl_var = "firm", seed = 123)
#
# print(sqrt(diag(boot_sand)))
# print(sqrt(diag(boot_ml)))


# print(summary(fit,
#               vcov.type = "boot",
#               repetitions = 50,
#               seed = 123))
#
# print(waldtest(fit, indices = 7:8,
#                vcov.type = "boot",
#                repetitions = 50,
#                seed = 123))
#
# predict <- predict(fit,
#                    vcov.type = "boot",
#                    repetitions = 50,
#                    se.fit = TRUE)

# Our bootstrap (excludes non-converged)
# boot_our <- vcov(fit, type = "boot",
#                  repetitions = 50,
#                  seed = 123,
#                  progress = FALSE)
#
#
# boot_robust <- vcov(fit, type = "robust")
#
# set.seed(123)

# suppressMessages(
#   suppressWarnings({
# # sandwich::vcovBS (includes whatever update() returns on failure)
# boot_sandwich <- sandwich::vcovBS(fit, R = 500)
# }))

# Or look at diagonals (standard errors)
# print(round(sqrt(diag(boot_our)), 5))
# print(round(sqrt(diag(boot_robust)), 5))

#
# # Assuming you have a fitted model
# fit <- ml_lm(value = mpg ~ wt + hp + factor(cyl), data = mtcars)
#
#
# # Try the default bootstrap covariance
# vcov_bs <- vcovBS(fit, R = 100)   # start with small R for testing
#
# # See if it works
# print(dim(vcov_bs))
# print(rownames(vcov_bs))
#
# library(sandwich)
#
# # Fit a heteroskedastic model
# fit_het <- ml_lm(
#   value = mpg ~ wt + hp + factor(cyl),
#   scale = ~ wt + hp,           # heteroskedastic scale
#   data = mtcars
# )
#
# # Try bootstrap covariance
# vcov_bs <- vcovBS(fit_het, R = 100)   # small R for speed
#
# # Check the result
# print(dim(vcov_bs))
# print(rownames(vcov_bs))
