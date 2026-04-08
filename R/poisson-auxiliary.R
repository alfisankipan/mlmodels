## HESSIANS ====================================================================
.ml_poisson_hessianObs <- function(object)
{
  if (!inherits(object, "ml_poisson"))
    cli::cli_abort("`object` must be a model of class 'ml_poisson' (from ml_poisson).",
                   call = NULL)
  
  # -- 1. Common elements for both hetero and homo -----------------------------
  b <- coef(object)
  y <- object$model$value$outcomes[[1]]
  x <- as.matrix(object$model$value$predictors)
  w <- object$model$weights %||% rep(1, length(y))   # default to 1 if NULL
  
  xb <- x %*% cbind(b)
  mu <- exp(xb)
  
  # Pre-allocate list to collect individual Hessians
  H_list <- vector("list", length(y))
  
  s <- as.vector(- mu * w)
  
  for (i in seq_len(nrow(x))) {
    xi  <- cbind(x[i, ])
    H_list[[i]] <- s[i] * tcrossprod(xi)
  }
  
  # Stack all individual Hessians
  do.call(rbind, H_list)
}



## ML EVALUATOR ==============================================================
#' @keywords internal
.ml_poisson_ll <- function(b, y, x, w = NULL)
{
  if(is.null(w))
    w <- rep(1, nrow(x))
  
  xb <- x %*% cbind(b)
  mu <- exp(xb)
  
  # Log-likelihood
  ll <- w * (y * xb - mu - lfactorial(y))
  
  # gradient
  g <- as.vector(w * (y - mu)) * x
  
  # Hessian
  s <- as.vector(- mu * w)
  H <- crossprod(x * s, x)
  
  # Attach gradient and Hessian as attributes
  attr(ll, "gradient") <- g
  attr(ll, "hessian")  <- H
  
  return(ll)
}