# =============================================================================
# marginaleffects compatibility (optional)
# =============================================================================

#' Implementation of `marginaleffects` package's `get_predict()` function.
#' 
#' @param model An `mlmodel` estimation object (from one of the package's
#'    estimation functions)
#' @param newdata A data.frame with the values of the predictors you want to use.
#'    It defaults to `NULL`, and in such case, it will use the data used in the
#'    estimation for the predictions (in sample prediction).
#' @param type A string with the type of prediction you want. This will depend
#'    on the types available for the model you estimated. See
#'    [predict][mlmodels::predict.mlmodel].
#' @param ... Additional arguments sent to methods.
#'
#' @returns A data.frame with two variables: `rowid` with the row number, and
#'   `estimate` with the value of the prediction requested.
#' 
#' @export
get_predict.mlmodel <- function(model,
                                newdata = NULL,
                                type = "response",
                                ...)
{
  pred <- predict(model, newdata = newdata, type = type, ...)
  out <- data.frame(
    rowid = seq_len(length(pred)),
    estimate = pred
  )
}

#' Implementation of `marginaleffects` package's `get_vcov()` function.
#' 
#' @param model An `mlmodel` estimation object (from one of the package's
#'    estimation functions)
#' @param ... Additional arguments passed to [vcov][mlmodels::vcov.mlmodel]
#' 
#' @returns A matrix from the package's `vcov()` function. See
#'    `[vcov][mlmodels::vcov.mlmodel]` for more details.
#' 
#' @export
get_vcov <- function(model, ...) {
  UseMethod("get_vcov")
}

#' @rdname get_vcov
#' @export
get_vcov.mlmodel <- function(model, ...)
{
  vcov(model, ...)
}

#' Implementation of `marginaleffects` package's `get_coef()` function.
#' 
#' @param model An `mlmodel` estimation object (from one of the package's
#'    estimation functions)
#' @param ... Additional arguments sent to methods.
#' 
#' @returns A named vector with the coefficients in the model.
#' 
#' @export
get_coef.mlmodel <- function(model, ...)
{
  coef(model, ...)
}

#' Implementation of `marginaleffects` package's `set_coef()` function.
#' 
#' @param model An `mlmodel` estimation object (from one of the package's
#'    estimation functions).
#' @param coefs A vector with the coefficients you want to se in the model.
#' @param ... Additional arguments sent to methods.
#' 
#' @returns A new `mlmodel` object, with the coefficients set to those in
#'   `coefs`/
#' 
#' @export
set_coef.mlmodel <- function(model, coefs, ...)
{
  out <- model
  out$estimate <- coefs
  return(out)
}

#' Implementation of `insight` package's `get_data()` function.
#' 
#' @param x An `mlmodel` estimation object (from one of the package's
#'    estimation functions).
#' @param ... Additional arguments sent to methods.
#' 
#' @returns A data.frame with the data the user passed when estimating the
#'    model.
#' 
#' @importFrom insight get_data
#' @export
get_data.mlmodel <- function(x, ...) {
  # Data is stored inside the $model element.
  if (!is.null(x$model$data) && is.data.frame(x$model$data)) {
    return(x$model$data)
  }

  # Safety fallback (in case we ever store it at root level in future models)
  if (!is.null(x$data) && is.data.frame(x$data)) {
    return(x$data)
  }

  cli::cli_abort("Could not recover the original data from the mlmodel object.",
                 call = NULL)
}

#' @rdname get_data.mlmodel
#' @export
get_modeldata.mlmodel <- get_data.mlmodel

# Register the mlmodel class when marginaleffects is available
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("marginaleffects", quietly = TRUE)) {

    # Safely turn NULL into an empty character vector
    options("marginaleffects_model_classes" = "mlmodel")
  }
  invisible(NULL)
}
