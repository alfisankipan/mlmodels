# =============================================================================
# marginaleffects compatibility (optional)
# =============================================================================

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

#' @export
get_vcov <- function(model, ...) {
  UseMethod("get_vcov")
}

#' @export
get_vcov.mlmodel <- function(model, ...)
{
  vcov(model, ...)
}

#' @export
get_coef.mlmodel <- function(model, ...)
{
  coef(model, ...)
}

#' @export
set_coef.mlmodel <- function(model, coefs, ...)
{
  # This is the most important one for marginaleffects
  # It returns a new model object with the coefficients replaced
  out <- model
  out$estimate <- coefs
  return(out)
}

#' @importFrom insight get_data
#' @export
get_data.mlmodel <- function(model, ...) {
  # Data is stored inside the $model component (as we do in ml_lm())
  if (!is.null(model$model$data) && is.data.frame(model$model$data)) {
    return(model$model$data)
  }

  # Safety fallback (in case we ever store it at root level in future models)
  if (!is.null(model$data) && is.data.frame(model$data)) {
    return(model$data)
  }

  cli::cli_abort("Could not recover the original data from the mlmodel object.",
                 call = NULL)
}

#' @export
get_modeldata.mlmodel <- get_data.mlmodel


# Register the mlmodel class when marginaleffects is available
.onLoad <- function(libname, pkgname) {
  # Explicitly register our S3 method so it works even if marginaleffects masks the generic
  #registerS3method("get_vcov", "mlmodel", get_vcov.mlmodel, envir = asNamespace("mlmodels"))

  # if (exists("get_vcov", envir = asNamespace("marginaleffects"), inherits = FALSE)) {
  #   assign("get_vcov", get("get_vcov", envir = asNamespace("marginaleffects")),
  #          envir = asNamespace("mlmodels"))
  # }

  # Also register the other marginaleffects methods for safety
  # registerS3method("get_predict", "mlmodel", get_predict.mlmodel, envir = asNamespace("mlmodels"))
  # registerS3method("get_coef", "mlmodel", get_coef.mlmodel, envir = asNamespace("mlmodels"))
  # registerS3method("set_coef", "mlmodel", set_coef.mlmodel, envir = asNamespace("mlmodels"))
  
  if (requireNamespace("marginaleffects", quietly = TRUE)) {

    # Safely turn NULL into an empty character vector
    options("marginaleffects_model_classes" = "mlmodel")
  }
  invisible(NULL)
}
