# Clear everything in the environment.
rm(list = ls(all.names = TRUE))
# Clear console
cat("\f")


library(devtools)
library(wooldridge)

load_all(".")


data("smoke")
data("mtcars", package = "datasets")

fit <- ml_lm(log(cigs) ~ educ,
             data = smoke)



model <- hardhat::mold(log(cigs) ~ educ,
                       data = smoke)

# fit <- ml_lm(cigs ~ educ, scale ~, data = smoke)

# model <- hardhat::mold(log(cigs) ~ income + educ,
#                                                data = smoke)

# fit <- ml_lm(log(cigs) ~ income + educ,
#              data = smoke)

# fit <- ml_lm(value = mpg ~ wt + hp,
#              scale = ~ wt,
#              data = mtcars)
#
# # Test update with a subset
# new_fit <- update(fit, data = mtcars[1:20, ])
#
# print(coef(new_fit))
#
#
# smoke$education <- 1 * (smoke$educ <= 8) +
#   2 * (smoke$educ > 8 & smoke$educ <= 12) +
#   3 * (smoke$educ > 12)
#
# educvector <- smoke$educ
#
# smoke$education <- factor(smoke$education,
#                            labels = c("Primary",
#                                       "Secondary",
#                                       "Other"))
#
# smoke$educ <- smoke$education
# smoke$education <- educvector

# model <- hardhat::mold(cigs ~ income + I((educ*education)^2),
#                        data = smoke)
#
# build_factor_mapping(list(value = model))

# fit<- ml_lm(cigs ~ income + education + poly(educ, 2),
#             scale = ~ restaurn,
#             data = smoke)

# fit<- ml_lm(cigs ~ income + cigpric + educ + I(cigpric / income),
#             scale = ~ restaurn,
#             data = smoke)
#
#
# print(summary(fit))

# print(summary(fit, vcov.type = "robust"))


# mtcars$engtrans <- 1 * (mtcars$vs == 0 & mtcars$am == 0) +
#   2 * (mtcars$vs == 0 & mtcars$am == 1) +
#   3 * (mtcars$vs == 1 & mtcars$am == 0) +
#   4 * (mtcars$vs == 1 & mtcars$am == 1)
#
# mtcars$engtrans2 <- mtcars$engtrans
# mtcars$engtrans3 <- mtcars$engtrans
# mtcars$engtrans <- factor(mtcars$engtrans,
#                           labels = c("Vaut", "Vstick", "Saut", "Sstick"))
#
#
# # This works, one column for enginetrans2.
# # model_mult <- hardhat::mold(mpg ~ I(engtrans^2) * (wt + hp) + I(engtrans2^2) * wt, data = mtcars)
# # Works too.
# # model_mult <- hardhat::mold(mpg ~ I(engtrans^2) * (wt + hp) + I(engtrans2^2) + I(engtrans3^2) * wt, data = mtcars)
# # This doesn't throw an error so it works too.
# # model_mult <- hardhat::mold(mpg ~ I(engtrans2^2) * (wt + hp), data = mtcars)
# # rm(model_mult)
# model_mult <- hardhat::mold(mpg ~ I((engtrans * engtrans2)^2) + wt + hp, data = mtcars)
# build_factor_map(model_mult)
# model_int <- hardhat::mold(mpg ~ (wt + hp) * engtrans, data = mtcars)
# model_int_2 <- hardhat::mold(mpg ~ engtrans / (wt + hp), data = mtcars)
# model_int_3 <- hardhat::mold(mpg ~ I(engtrans^2) * (wt + hp), data = mtcars)
# # model_int_4 <- hardhat::mold(mpg ~ I(log(engtrans)) * (wt + hp), data = mtcars)
# model_int_5 <- hardhat::mold(mpg ~ I(engtrans * wt) + hp, data = mtcars)
# model_poly <- hardhat::mold(mpg ~ poly(engtrans, 3) + wt + hp, data = mtcars)
# # model_cut <- hardhat::mold(mpg ~ cut(engtrans, 3) + wt + hp, data = mtcars)
# # model_bs <- hardhat::mold(mpg ~ bs(engtrans, 3) + wt + hp, data = mtcars)
# # model_ns <= hardhat::mold(mpg ~ ns(engtrans, 3) + wt + hp, data = mtcars)
#
# cat("int 1\n")
# build_factor_map(model_int) # should return nothing.
# cat("int 2\n")
# build_factor_map(model_int_2) # should return nothing.
# cat("int 3\n")
# build_factor_map(model_int_3) # should return the names with I(.
# cat("int 5\n")
# build_factor_map(model_int_5) # should return the names with I(.
# cat("poly\n")
# build_factor_map(model_poly) # should return the names with poly(.
#

