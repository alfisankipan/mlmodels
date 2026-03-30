# Clear everything in the environment.
rm(list = ls(all.names = TRUE))
# Clear console
cat("\f")

library(devtools)
library(wooldridge)

load_all(".")

data("smoke")

smoke$smokes <- smoke$cigs > 0

fit_log <- ml_logit(smokes ~ cigpric + income + age,
                    data = smoke)

# Fit with glm (standard binary logit)
fit_glm <- glm(smokes ~ cigpric + income + age,
               data = smoke,
               family = binomial(link = "logit"))

# Comparison
cat("=== Coefficient Comparison ===\n")
cat(" ml_logit     glm",
  sprintf("%.6f      %.6f",
          coef(fit_log),
          coef(fit_glm)),
  sep = "\n")

cat("\n=== Standard Error Comparison ===\n")
cat(" ml_logit     glm",
    sprintf("%.6f      %.6f",
            sqrt(diag(vcov(fit_log))),
            sqrt(diag(vcov(fit_glm)))),
    sep = "\n")

cat("\n=== Log-Likelihood ===\n")
cat("ml_logit :", logLik(fit_log), "\n")
cat("glm      :", logLik(fit_glm), "\n")

cat("\n=== Number of observations ===\n")
cat("ml_logit :", fit_log$model$n_used, "\n")
cat("glm      :", nobs(fit_glm), "\n")


fit_het  <- ml_logit(smokes ~ cigpric + income + age,
                     scale = ~ restaurn,
                     data = smoke)

cat("\n=== Heteroskedastic Results ===\n")
cat("        Estimate     Std.Error",
    sprintf("%s  %.4f   %.4f",
            names(coef(fit_het)),
            coef(fit_het),
            sqrt(diag(vcov(fit_het)))),
    sep = "\n")
