# Confidence Intervals for mlmodel Coefficients

Confidence Intervals for mlmodel Coefficients

## Usage

``` r
# S3 method for class 'mlmodel'
confint(
  object,
  parm,
  level = 0.95,
  vcov = NULL,
  vcov.type = "oim",
  cl_var = NULL,
  repetitions = 999,
  seed = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- object:

  An `mlmodel` object.

- parm:

  A specification of which parameters are to be given confidence
  intervals (names or numeric indices). If missing, all parameters are
  used.

- level:

  The confidence level required. Default is 0.95.

- vcov:

  Optional user-supplied variance-covariance matrix.

- vcov.type:

  Type of variance-covariance matrix to use. See
  [vcov.mlmodel](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md).

- cl_var:

  Clustering variable (name as string or vector).

- repetitions:

  Number of bootstrap replications when `vcov.type = "boot"`.

- seed:

  Random seed for bootstrap/jackknife.

- progress:

  Show progress bar? Default `FALSE`.

- ...:

  Further arguments passed to methods.

## Value

Matrix with the confidence intervals for the requested parameters.

## Details

Confidence intervals are constructed as \\\hat{\beta} \pm
z\_{1-\alpha/2} \times SE(\hat{\beta})\\, where the standard errors come
from the requested variance-covariance matrix.

The function supports all variance types available in
[vcov.mlmodel](https://alfisankipan.github.io/mlmodels/reference/vcov.mlmodel.md),
including robust, clustered, bootstrap, and jackknife estimators.

## Examples

``` r

data(mroz)
mroz$incthou <- mroz$faminc / 1000

fit <- ml_lm(incthou ~ age + I(age^2) + huswage + educ + unem, 
             data = mroz)

# Default 95% confidence intervals (using OIM)
confint(fit)
#>                           2.5 %        97.5 %
#> value::(Intercept) -46.40582832 -12.549703550
#> value::age           0.47905887   2.045570195
#> value::I(age^2)     -0.02259875  -0.004531053
#> value::huswage       1.81288004   2.100227030
#> value::educ          0.69783722   1.229376112
#> value::unem         -0.44146127  -0.066150422
#> scale::lnsigma       2.03477652   2.135786802

# 90% confidence intervals
confint(fit, level = 0.90)
#>                             5 %          95 %
#> value::(Intercept) -43.68424372 -15.271288147
#> value::age           0.60498566   1.919643406
#> value::I(age^2)     -0.02114634  -0.005983457
#> value::huswage       1.83597894   2.077128134
#> value::educ          0.74056591   1.186647416
#> value::unem         -0.41129124  -0.096320450
#> scale::lnsigma       2.04289641   2.127666911

# Confidence intervals for specific parameters
confint(fit, parm = c("value::educ", "value::huswage"))
#>                    2.5 %   97.5 %
#> value::huswage 1.8128800 2.100227
#> value::educ    0.6978372 1.229376
confint(fit, parm = 4:5)                     # by position
#>                    2.5 %   97.5 %
#> value::huswage 1.8128800 2.100227
#> value::educ    0.6978372 1.229376

# Using different variance types
confint(fit, vcov.type = "robust")
#>                          2.5 %        97.5 %
#> value::(Intercept) -45.9635549 -12.991976991
#> value::age           0.5130501   2.011578948
#> value::I(age^2)     -0.0221728  -0.004956998
#> value::huswage       1.6900907   2.223016356
#> value::educ          0.6398659   1.287347451
#> value::unem         -0.4472090  -0.060402725
#> scale::lnsigma       1.9764044   2.194158895

# Clustered confidence intervals
confint(fit, vcov.type = "robust", cl_var = "age")
#>                           2.5 %       97.5 %
#> value::(Intercept) -43.66514935 -15.29038252
#> value::age           0.62392349   1.90070558
#> value::I(age^2)     -0.02107893  -0.00605087
#> value::huswage       1.67335723   2.23974985
#> value::educ          0.68879938   1.23841395
#> value::unem         -0.41141057  -0.09620112
#> scale::lnsigma       1.96934455   2.20121877

# Using a pre-computed bootstrap variance matrix
v_boot <- vcov(fit, type = "boot", repetitions = 100, seed = 123)
#> ℹ Bootstrap with 100 repetitions.
#>  0        10        20        30        40        50
#> ====================================================
#>  ..................................................
#>  ..................................................
#> ====================================================
#> 
#> Bootstrapping finished - 100% of replications converged.
confint(fit, vcov = v_boot)
#>                           2.5 %        97.5 %
#> value::(Intercept) -44.00792493 -14.947606932
#> value::age           0.59339474   1.931234322
#> value::I(age^2)     -0.02133841  -0.005791393
#> value::huswage       1.66080920   2.252297880
#> value::educ          0.66735078   1.259862545
#> value::unem         -0.45997057  -0.047641119
#> scale::lnsigma       1.98213086   2.188432458
```
