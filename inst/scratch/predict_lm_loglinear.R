# Clear everything in the environment.
rm(list = ls(all.names = TRUE))
# Clear console
cat("\f")


library(devtools)
library(wooldridge)

load_all(".")

data(mtcars)

# 1. Fit with logged outcome
fit_lm  <- lm(log(mpg) ~ wt + hp, data = mtcars)
fit_ml  <- ml_lm(log(mpg) ~ wt + hp,
                 scale = ~ hp,
                 data = mtcars)

# 2. Compare predictions on the original scale of mpg
newdata <- mtcars[1:5, c("wt", "hp")]

pred_lm_raw   <- predict(fit_lm,  newdata = newdata)           # log scale
pred_lm       <- exp(pred_lm_raw)                              # manual back-transform (median only)

pred_ml_mean  <- predict(fit_ml, newdata = newdata, type = "mean")      # correct E[mpg]
pred_ml_med   <- predict(fit_ml, newdata = newdata, type = "median")    # should match exp()
pred_ml_link  <- predict(fit_ml, newdata = newdata, type = "link")      # log-scale linear predictor

# 3. Print comparison
cat("lm() raw (log scale):     ", round(pred_lm_raw, 3), "\n")
cat("lm() manual exp():        ", round(pred_lm, 3), "\n")
cat("ml_lm() median:           ", round(pred_ml_med, 3), "\n")
cat("ml_lm() mean (E[mpg]):    ", round(pred_ml_mean, 3), "\n")
cat("ml_lm() link (log scale): ", round(pred_ml_link, 3), "\n")
