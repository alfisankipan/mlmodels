.ml_logit_ll <- function(b, y, x, z = NULL, w)
{
  if(is.null(w))
    w <- rep(1, nrow(x))
  if (is.null(z)) {
    # ── Homoskedastic binary logit ─────────────────────────────────────
    xb <- as.vector(x %*% cbind(b))
    py <- exp(xb) / (1 + exp(xb))

    # Weighted log-likelihood
    ll <- y * xb - log(1 + exp(xb))

    # Gradient
    g <- w * as.vector(y - py) * x

    # Hessian
    H <- matrix(0, nrow = length(b), ncol = length(b))
    for (i in seq_len(nrow(x))) {
      pyi <- py[i]
      xi  <- cbind(x[i, ])
      H   <- H - w[i] * pyi * (1 - pyi) * tcrossprod(xi)
    }

  } else {
    # ── Heteroskedastic binary logit ───────────────────────────────────
    k1   <- ncol(x)
    k    <- k1 + ncol(z)
    beta <- b[1:k1]
    delta <- b[(k1+1):k]

    zd <- as.vector(z %*% cbind(delta))
    xz <- x / exp(zd)
    xb <- as.vector(xz %*% cbind(beta))
    py <- exp(xb) / (1 + exp(xb))

    # Weighted log-likelihood
    ll <- w * (y * xb - log(1 + exp(xb)))

    # Gradient
    gb <- w * as.vector(y - py) * xz
    gd <- w * as.vector((py - y) * xb) * z
    g  <- cbind(gb, gd)

    # Hessian
    H <- matrix(0, nrow = k, ncol = k)
    for (i in seq_len(nrow(x))) {
      pyi <- py[i]
      xbi <- xb[i]
      xzi <- cbind(xz[i, ])
      zi  <- cbind(z[i, ])
      wi  <- w[i]

      hbb <- -wi * pyi * (1 - pyi) * tcrossprod(xzi)
      hbd <- wi * (pyi * (1 - pyi) * xbi - (y[i] - pyi)) * tcrossprod(xzi, zi)
      hdd <- wi * ((y[i] - pyi) * xbi - pyi * (1 - pyi) * xbi^2) * tcrossprod(zi)

      H <- H + rbind(cbind(hbb, hbd),
                     cbind(t(hbd), hdd))
    }
  }

  # Attach gradient and Hessian as attributes
  attr(ll, "gradient") <- g
  attr(ll, "hessian")  <- H

  return(ll)
}


# HESSIANS ---------------------------------------------------------------------
.ml_logit_hessianObs <- function(object)
{
  if (!inherits(object, "ml_logit"))
    cli::cli_abort("`object` must be a model of class 'ml_logit' (from ml_logit).",
                   call = NULL)

  # -- 1. Common elements for both hetero and homo -----------------------------
  b <- coef(object)
  y <- object$model$value$outcomes[[1]]
  x <- as.matrix(object$model$value$predictors)
  w <- object$model$weights %||% rep(1, length(y))   # default to 1 if NULL

  k1   <- ncol(x)
  beta <- b[1:k1]

  # Pre-allocate list to collect individual Hessians
  H_list <- vector("list", length(y))

  if (is.null(object$model$scale)) {
    # ── Homoskedastic case ─────────────────────────────────────
    xb <- as.vector(x %*% cbind(beta))
    py <- exp(xb) / (1 + exp(xb))

    for (i in seq_along(y)) {
      pyi <- py[i]
      xi  <- cbind(x[i, ])
      H_list[[i]] <- -w[i] * pyi * (1 - pyi) * tcrossprod(xi)
    }

  } else {
    # ── Heteroskedastic case ───────────────────────────────────
    z     <- as.matrix(object$model$scale$predictors)
    k     <- k1 + ncol(z)
    delta <- b[(k1+1):k]

    zg <- as.vector(z %*% cbind(delta))
    xz <- x / exp(zg)
    xb <- as.vector(xz %*% cbind(beta))
    py <- exp(xb) / (1 + exp(xb))

    for (i in seq_along(y)) {
      pyi <- py[i]
      xbi <- xb[i]
      xzi <- cbind(xz[i, ])
      zi  <- cbind(z[i, ])
      wi  <- w[i]

      hbb <- -wi * pyi * (1 - pyi) * tcrossprod(xzi)
      hbd <- wi * (pyi * (1 - pyi) * xbi - (y[i] - pyi)) * tcrossprod(xzi, zi)
      hdd <- wi * ((y[i] - pyi) * xbi - pyi * (1 - pyi) * xbi^2) * tcrossprod(zi)

      H_list[[i]] <- rbind(cbind(hbb, hbd),
                           cbind(t(hbd), hdd))
    }
  }

  # Stack all individual Hessians
  do.call(rbind, H_list)
}
