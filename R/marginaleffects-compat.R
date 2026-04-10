# =============================================================================
# marginaleffects compatibility layer
# This file must be processed early so generics are defined before methods
# =============================================================================

# --- S3 methods for mlmodel class ----------------------------------------
#' @export
get_predict.mlmodel <- function(model,
                                newdata = NULL,
                                type = "response",
                                vcov.type = "oim",
                                ...)
{
  pred <- predict(model, newdata = newdata, type = type, vcov.type = vcov.type, ...)
  data.frame(
    rowid = seq_len(length(pred)),
    estimate = pred
  )
}

#' @export
get_vcov.mlmodel <- function(model, vcov = NULL, ...)
{
  if(!is.null(vcov))
    return(vcov)
  vcov(model)
}

#' @export
get_coef.mlmodel <- function(model, ...)
{
  coef(model, ...)
}

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