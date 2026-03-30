# Clear everything in the environment.
rm(list = ls(all.names = TRUE))
# Clear console
cat("\f")

library(devtools)
library(wooldridge)
library(devtools)

load_all(".")

data("smoke")

# Fit a heteroskedastic Gaussian model on cigs (we expect it to fail the IM test)
fit <- ml_lm(
  value = cigs ~ log(income) + age + I(age^2) + educ + restaurn,
  scale = ~ age + I(age^2) + restaurn,
  data = smoke
)

summary(fit)

# Run the improved OPG bootstrap IM test
im <- IMtest(fit,
             method = "boot_opg",
             repetitions = 10,     # feel free to lower to 200 if it's slow
             seed = 123)

print(im)

# Let's do it with logit
smoke$smokes <- smoke$cigs > 0

fit_log <- ml_logit(smokes ~ cigpric + income + age,
                    data = smoke)

im <- IMtest(fit_log,
             method = "quad")
print(im)

im <- IMtest(fit_log,
             method = "opg")
print(im)
