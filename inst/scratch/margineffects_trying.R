# Clear everything in the environment.
rm(list = ls(all.names = TRUE))
# Clear console
cat("\f")

library(devtools)
library(marginaleffects)
library(wooldridge)

load_all(".")

data("mtcars")

library(marginaleffects)

# Fit a model
fit <- ml_lm(log(mpg) ~ wt + hp,
             scale = ~ hp,
             data = mtcars)

print(summary(fit))

# insight::get_data(fit)                  # should return the full mtcars data.frame
# dy <- marginaleffects::avg_slopes(fit, variables = "hp")        # no newdata argument needed anymore!
# pred <- predict(fit)

# # Original coefficients
# coef_original <- coef(fit)
#
# # Modify coefficients using set_coef
# fit_new <- set_coef(fit, coefs = coef(fit) * 1.1)
#
# # Check that the original model was NOT modified
# print(identical(coef(fit), coef_original))        # should be TRUE
#
# print(coef_original)
#
# # Check that the new model has the modified coefficients
# print(coef(fit_new))

# Try a marginaleffects function that uses set_coef internally
# slopes(fit, variables = "wt", newdata = mtcars)

# print(get_predict(fit, newdata = mtcars))

# print(avg_slopes(fit, variables = "hp", newdata = mtcars,
#                  type = "variance_y"))
#
# print(avg_slopes(fit, variables = "hp", newdata = mtcars,
#                  type = "fitted"))
#
# print(avg_slopes(fit, variables = "hp", newdata = mtcars,
#                  type = "response"))
#
# print(summary(fit))

