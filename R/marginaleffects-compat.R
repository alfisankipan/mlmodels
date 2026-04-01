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

  cli::cli_abort("Could not recover the original data from the mlmodel object.")
}

#' @export
get_modeldata.mlmodel <- get_data.mlmodel


# Register the mlmodel class when marginaleffects is available
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("marginaleffects", quietly = TRUE)) {

    # Get current list of registered classes (may be NULL)
    current <- getOption("marginaleffects_model_classes")

    # Safely turn NULL into an empty character vector
    if (is.null(current))
      options(marginaleffects_model_classes = "mlmodel")
    else
      options(marginaleffects_model_classes = c(current, "mlmodel"))
  }
  invisible(NULL)
}


# usethis::use_git_remote(
#   name = "origin",
#   url  = "https://github.com/alfisankipan/mlmodels.git",
#   overwrite = TRUE
# )
