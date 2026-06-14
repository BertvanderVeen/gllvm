## ---------------------------------------------------------------------------
## HO-specific helpers for predictSR se.fit > 0
## ---------------------------------------------------------------------------

## Draw R parameter sets from the VA posterior of a gllvmHO object.
## Returns a list suitable for passing to perturb_gllvmHO().
simulate_params_gllvmHO <- function(object, R, seed = 42, level = 1, n = NULL) {
  set.seed(seed)

  if (is.null(object$Hess) || is.null(object$Hess$incl))
    stop("Hessian not available. Run se.gllvm(object) before predictSR with se.fit > 0.")

  p     <- ncol(object$y)
  n_fit <- nrow(object$y)
  if (is.null(n)) n <- n_fit

  par_hat <- object$TMBfn$env$last.par.best
  par_nms <- names(par_hat)

  num.RR  <- object$num.RR  %||% 0L
  num.lvc <- object$num.lv.c %||% 0L
  num.lv  <- object$num.lv  %||% 0L
  d       <- num.RR + num.lvc + num.lv
  d_active <- num.RR + num.lvc        # cols used by b_z / b_gamma
  Kz <- if (!is.null(object$lv.X)) ncol(object$lv.X) else 0L
  Kt <- if (!is.null(object$TR))   ncol(as.matrix(object$TR)) else 0L

  d_va_z <- ncol(.ho_diag(object$A))
  d_va_a <- ncol(.ho_diag(object$B))

  ## --- fixed params (structural) from Hessian MVN ---
  incl <- object$Hess$incl
  ffs  <- NULL
  if (sum(incl) > 0L) {
    Vf  <- vcov(object)
    ffs <- try(MASS::mvrnorm(R, object$TMBfn$par[incl], Vf), silent = TRUE)
    if (inherits(ffs, "try-error"))
      stop("Fixed-effects covariance matrix is not semi positive-definite.")
    colnames(ffs) <- names(par_hat)[incl]
  }

  ## --- VA draws: all independent normals ---
  u_sim   <- NULL
  a_sim   <- NULL
  bz_sim  <- NULL
  bg_sim  <- NULL
  r0r_sim <- NULL

  if (level > 0L && d_va_z > 0L && n == n_fit) {
    u_hat  <- matrix(par_hat[par_nms == "u"], nrow = n_fit, ncol = d_va_z)
    u_sd   <- sqrt(pmax(.ho_diag(object$A), 0))
    u_sim  <- matrix(rnorm(R * n_fit * d_va_z), R, n_fit * d_va_z)
    u_sim  <- sweep(u_sim, 2L, as.vector(u_sd), `*`) +
              matrix(as.vector(u_hat), R, n_fit * d_va_z, byrow = TRUE)
  }

  if (level > 0L && d_va_a > 0L) {
    a_hat  <- matrix(par_hat[par_nms == "a_lv_sp"], nrow = p, ncol = d_va_a)
    a_sd   <- sqrt(pmax(.ho_diag(object$B), 0))
    a_sim  <- matrix(rnorm(R * p * d_va_a), R, p * d_va_a)
    a_sim  <- sweep(a_sim, 2L, as.vector(a_sd), `*`) +
              matrix(as.vector(a_hat), R, p * d_va_a, byrow = TRUE)
  }

  ## b_z (LvXcoef unscaled) — only draw VA uncertainty when randomB is set
  if (Kz > 0L && d_active > 0L && !isFALSE(object$randomB)) {
    bz_hat  <- matrix(par_hat[par_nms == "b_z"], nrow = Kz, ncol = d)[, seq_len(d_active), drop = FALSE]
    Ab_z_raw <- par_hat[par_nms == "Ab_z"][seq_len(Kz * d_active)]
    ## C++ diagonal: AB_z(l)(q,q) = exp(Ab_z(q*d_c + l)), so row-major reshape
    bz_sd   <- sqrt(pmax(exp(matrix(Ab_z_raw, nrow = d_active, ncol = Kz)), 0))  # d_active x Kz
    bz_sd   <- t(bz_sd)   # Kz x d_active
    bz_sim  <- matrix(rnorm(R * Kz * d_active), R, Kz * d_active)
    bz_sim  <- sweep(bz_sim, 2L, as.vector(bz_sd), `*`) +
               matrix(as.vector(bz_hat), R, Kz * d_active, byrow = TRUE)
  }

  ## b_gamma (LoadTRcoef) — only draw VA uncertainty when randomT is set
  if (Kt > 0L && d_active > 0L && !isFALSE(object$randomT)) {
    bg_hat   <- matrix(par_hat[par_nms == "b_gamma"], nrow = Kt, ncol = d)[, seq_len(d_active), drop = FALSE]
    Ab_g_raw <- par_hat[par_nms == "Ab_gamma"][seq_len(Kt * d_active)]
    bg_sd    <- sqrt(pmax(exp(matrix(Ab_g_raw, nrow = d_active, ncol = Kt)), 0))  # d_active x Kt
    bg_sd    <- t(bg_sd)  # Kt x d_active
    bg_sim   <- matrix(rnorm(R * Kt * d_active), R, Kt * d_active)
    bg_sim   <- sweep(bg_sim, 2L, as.vector(bg_sd), `*`) +
               matrix(as.vector(bg_hat), R, Kt * d_active, byrow = TRUE)
  }

  ## r0r (random row effects)
  if (!is.null(object$params$row.params.random)) {
    r0r_hat <- par_hat[par_nms == "r0r"]
    lg_Ar   <- par_hat[par_nms == "lg_Ar"]
    r0r_sim <- matrix(rnorm(R * length(r0r_hat)), R, length(r0r_hat))
    r0r_sim <- sweep(r0r_sim, 2L, exp(lg_Ar), `*`) +
               matrix(r0r_hat, R, length(r0r_hat), byrow = TRUE)
  }

  list(ffs = ffs, incl = incl,
       u_sim = u_sim, a_sim = a_sim,
       bz_sim = bz_sim, bg_sim = bg_sim, r0r_sim = r0r_sim,
       d_va_z = d_va_z, d_va_a = d_va_a,
       d_active = d_active, Kz = Kz, Kt = Kt, n_fit = n_fit)
}


## Rebuild a gllvmHO object for simulation draw r, ready for predict().
perturb_gllvmHO <- function(object, params, r, skeleton = NULL, template = NULL) {
  if (is.null(skeleton)) skeleton <- object$TMBfn$env$parList()

  p       <- ncol(object$y)
  num.RR  <- object$num.RR  %||% 0L
  num.lvc <- object$num.lv.c %||% 0L
  num.lv  <- object$num.lv  %||% 0L
  d       <- num.RR + num.lvc + num.lv
  d_active <- params$d_active
  Kz       <- params$Kz
  Kt       <- params$Kt
  n_fit    <- params$n_fit
  d_va_z   <- params$d_va_z
  d_va_a   <- params$d_va_a

  Kz_eff <- if (!is.null(object$lv.X)) ncol(object$lv.X) else 0L
  rr_va_a <- if (Kz_eff > 0L) 0L else num.RR  # VA offset into a_lv_sp for non-RR dims

  newobj <- if (!is.null(template)) template else {
    tmp <- object
    tmp$params <- lapply(object$params, function(x) { if (is.numeric(x)) x[] <- NA; x })
    if (!is.null(tmp$lvs)) tmp$lvs[] <- NA
    tmp
  }

  sigma_new <- object$params$sigma.lv  # fallback if no fixed draw

  ## --- Fixed structural params ---
  if (!is.null(params$ffs) && sum(params$incl) > 0L) {
    newpars <- relist_gllvm(params$ffs[r, ], skeleton)

    newobj$params$beta0 <- newpars$b[1L, ]
    names(newobj$params$beta0) <- names(object$params$beta0)
    if (nrow(newpars$b) > 1L) {
      newobj$params$Xcoef <- t(newpars$b[-1L, , drop = FALSE])
      rownames(newobj$params$Xcoef) <- rownames(object$params$Xcoef)
      colnames(newobj$params$Xcoef) <- colnames(object$params$Xcoef)
    }

    ## sigma.lv from sigmaLV: sigma[1] = exp(slv[1]); sigma[k] = sigma[k-1]*exp(-exp(slv[k]))
    if (d > 0L && !is.null(newpars$sigmaLV)) {
      slv <- newpars$sigmaLV
      sigma_new <- numeric(d)
      sigma_new[1L] <- exp(slv[1L])
      if (d > 1L)
        for (k in 2L:d) sigma_new[k] <- sigma_new[k - 1L] * exp(-exp(slv[k]))
      names(sigma_new) <- names(object$params$sigma.lv)
      newobj$params$sigma.lv <- sigma_new
    }

    ## phi / dispersion (ZIP/ZINB: logistic, others: exp)
    if (any(object$family %in% c("ZIP", "ZINB")) && !is.null(newpars$lg_phi)) {
      lp0 <- newpars$lg_phi[object$disp.group]
      newobj$params$phi <- object$params$phi
      newobj$params$phi[object$family %in% c("ZIP", "ZINB")] <-
        (exp(lp0) / (1 + exp(lp0)))[object$family %in% c("ZIP", "ZINB")]
    } else if (!is.null(newpars$lg_phi)) {
      newobj$params$phi <- exp(newpars$lg_phi)
      names(newobj$params$phi) <- names(object$params$phi)
    }

    ## fixed row effects
    if (!is.null(newpars$r0f)) {
      newobj$params$row.params.fixed <- c(newpars$r0f)
      names(newobj$params$row.params.fixed) <- names(object$params$row.params.fixed)
    }

    ## LvXcoef (b_z sigma-scaled) when b_z is fixed (randomB=FALSE)
    if (Kz > 0L && d_active > 0L && isFALSE(object$randomB) && !is.null(newpars$b_z)) {
      bz_full <- matrix(c(newpars$b_z), nrow = Kz, ncol = d)
      lvx_cols_ho <- if (num.RR > 0L && num.lvc > 0L) {
        c(seq(num.RR + 1L, num.RR + num.lvc), seq_len(num.RR))
      } else seq_len(d_active)
      bz_act <- bz_full[, lvx_cols_ho, drop = FALSE]
      newobj$params$LvXcoef <- t(t(bz_act) * sigma_new[lvx_cols_ho])
      rownames(newobj$params$LvXcoef) <- rownames(object$params$LvXcoef)
      colnames(newobj$params$LvXcoef) <- colnames(object$params$LvXcoef)
    }

    ## LoadTRcoef (b_gamma) when fixed (randomT=FALSE)
    if (Kt > 0L && d_active > 0L && isFALSE(object$randomT) && !is.null(newpars$b_gamma)) {
      bg_full <- matrix(c(newpars$b_gamma), nrow = Kt, ncol = d)
      newobj$params$LoadTRcoef <- bg_full[, seq_len(d_active), drop = FALSE]
      rownames(newobj$params$LoadTRcoef) <- rownames(object$params$LoadTRcoef)
      colnames(newobj$params$LoadTRcoef) <- colnames(object$params$LoadTRcoef)
    }
  }

  ## --- VA: site scores (u) → lvs ---
  if (!is.null(params$u_sim) && d_va_z > 0L) {
    newobj$lvs <- matrix(params$u_sim[r, ], nrow = n_fit, ncol = d_va_z)
  }

  ## --- VA: species loadings (a_lv_sp) → theta ---
  if (!is.null(params$a_sim) && d_va_a > 0L) {
    a_r <- matrix(params$a_sim[r, ], nrow = p, ncol = d_va_a)
    ## Diagonal is stored on log scale for sign identification; back-transform here.
    for (ia in seq_len(d_va_a)) if (ia <= p) a_r[ia, ia] <- exp(a_r[ia, ia])

    ## Determine LoadTRcoef to use for RR dims (fixed or VA-drawn)
    bg_curr <- if (!is.null(params$bg_sim)) {
      matrix(params$bg_sim[r, ], nrow = Kt, ncol = d_active)
    } else if (!is.null(newobj$params$LoadTRcoef)) {
      newobj$params$LoadTRcoef
    } else {
      object$params$LoadTRcoef
    }

    TR_mat <- if (!is.null(object$TR) && Kt > 0L) as.matrix(object$TR) else NULL

    loadings_full <- matrix(0, p, d)
    for (k in seq_len(d)) {
      if (k <= num.RR) {
        if (!is.null(TR_mat) && !is.null(bg_curr)) {
          loadings_full[, k] <- TR_mat %*% bg_curr[, k]
        } else if (d_va_a > 0L) {
          loadings_full[, k] <- a_r[, k]
        }
      } else {
        ia <- rr_va_a + (k - num.RR)
        if (ia >= 1L && ia <= d_va_a) loadings_full[, k] <- a_r[, ia]
      }
    }

    ## Reorder to standard [lvc|RR|lv]
    ho_to_std_idx <- if (num.RR > 0L && num.lvc > 0L) {
      c(seq(num.RR + 1L, num.RR + num.lvc), seq_len(num.RR),
        if (num.lv > 0L) seq(num.RR + num.lvc + 1L, d) else integer(0))
    } else seq_len(d)

    newobj$params$theta <- loadings_full[, ho_to_std_idx, drop = FALSE]
    rownames(newobj$params$theta) <- rownames(object$params$theta)
    colnames(newobj$params$theta) <- colnames(object$params$theta)
  }

  ## --- VA: LvXcoef (b_z sigma-scaled) when randomB=TRUE ---
  if (!is.null(params$bz_sim) && Kz > 0L && d_active > 0L) {
    bz_r <- matrix(params$bz_sim[r, ], nrow = Kz, ncol = d_active)
    lvx_cols_ho <- if (num.RR > 0L && num.lvc > 0L) {
      c(seq(num.RR + 1L, num.RR + num.lvc), seq_len(num.RR))
    } else seq_len(d_active)
    newobj$params$LvXcoef <- t(t(bz_r) * sigma_new[lvx_cols_ho])
    rownames(newobj$params$LvXcoef) <- rownames(object$params$LvXcoef)
    colnames(newobj$params$LvXcoef) <- colnames(object$params$LvXcoef)
  }

  ## --- VA: LoadTRcoef (b_gamma) when randomT=TRUE ---
  if (!is.null(params$bg_sim) && Kt > 0L && d_active > 0L) {
    bg_r <- matrix(params$bg_sim[r, ], nrow = Kt, ncol = d_active)
    newobj$params$LoadTRcoef <- bg_r
    rownames(newobj$params$LoadTRcoef) <- rownames(object$params$LoadTRcoef)
    colnames(newobj$params$LoadTRcoef) <- colnames(object$params$LoadTRcoef)
  }

  ## --- VA: random row effects ---
  if (!is.null(params$r0r_sim)) {
    newobj$params$row.params.random <- params$r0r_sim[r, ]
    names(newobj$params$row.params.random) <- names(object$params$row.params.random)
  }

  newobj
}


#' @title Predict Method for gllvmHO Fits
#' @description Obtains predictions from a fitted hierarchical ordination model.
#'
#' @param object   A \code{gllvmHO} object.
#' @param newX     Optional data frame of environmental variables for new sites.
#' @param newTR    Ignored (trait effects are absorbed into loadings).
#' @param newLV    Optional matrix of latent variable scores for new sites,
#'   \eqn{n_\mathrm{new} \times (\mathrm{num.lv.c} + \mathrm{num.lv})} in
#'   standard gllvm order.
#' @param type     \code{"link"} (default) or \code{"response"}.
#' @param level    1 (default) uses VA posterior means for existing sites; 0 sets
#'   latent variables to zero (only covariate-driven RR/lvc part retained).
#' @param offset   Logical or matrix. \code{TRUE} (default) includes the
#'   training offsets; \code{FALSE} ignores them; a matrix of new offsets.
#' @param ...      Not used.
#'
#' @return An \eqn{n \times p} matrix of predictions.
#'
#' @method predict gllvmHO
#' @export
predict.gllvmHO <- function(object, newX = NULL, newTR = NULL, newLV = NULL,
                             type = "link", level = 1, offset = TRUE, ...) {
  p       <- ncol(object$y)
  n_orig  <- nrow(object$y)
  d       <- object$num.RR + object$num.lv.c + object$num.lv
  sigma   <- object$params$sigma.lv    # HO order [RR|lvc|lv], length d
  num.RR  <- object$num.RR
  num.lvc <- object$num.lv.c
  num.lv.unc <- object$num.lv          # unconstrained

  n <- if (!is.null(newX))  nrow(newX)
       else if (!is.null(newLV)) nrow(newLV)
       else n_orig

  ## --- beta0 ---------------------------------------------------------------
  eta <- matrix(object$params$beta0, n, p, byrow = TRUE)

  ## --- fixed X effects -----------------------------------------------------
  if (!is.null(object$X) && !is.null(object$params$Xcoef)) {
    if (is.null(newX)) {
      ## X.design includes the intercept column; Xcoef has no intercept row
      X.d <- object$X.design[, -1L, drop = FALSE]
    } else {
      ## For new sites: match columns to training design (no intercept)
      xcols <- colnames(object$X.design)[-1L]  # drop intercept name
      newX_mat <- as.matrix(as.data.frame(newX))
      if (all(xcols %in% colnames(newX_mat))) {
        X.d <- newX_mat[, xcols, drop = FALSE]
      } else {
        X.d <- newX_mat
      }
    }
    eta <- eta + X.d %*% t(object$params$Xcoef)
  }

  ## --- bilinear LV contribution --------------------------------------------
  .lv_x_for_new <- function(newdata) {
    lv_X_cols <- colnames(object$lv.X)
    if (is.data.frame(newdata) || is.matrix(newdata)) {
      mat <- as.matrix(newdata)
      if (!is.null(lv_X_cols) && all(lv_X_cols %in% colnames(mat)))
        return(mat[, lv_X_cols, drop = FALSE])
      return(mat)
    }
    return(NULL)
  }

  if (level == 0 && is.null(newLV)) {
    ## Marginal (z_i = 0): only the covariate-driven lvc/RR part
    if ((num.RR + num.lvc) > 0L && !is.null(object$params$LvXcoef)) {
      lv.X_use <- if (!is.null(newX)) .lv_x_for_new(newX) else object$lv.X.design
      if (!is.null(lv.X_use)) {
        theta_lvc_rr <- object$params$theta[, seq_len(num.lvc + num.RR), drop = FALSE]
        eta <- eta + lv.X_use %*% object$params$LvXcoef %*% t(theta_lvc_rr)
      }
    }

  } else if (level == 1) {
    if (!is.null(newLV)) {
      ## User-supplied VA scores: n_new × (num.lvc + num.lv.unc), standard order
      if (ncol(newLV) != num.lvc + num.lv.unc)
        stop("newLV must have ", num.lvc + num.lv.unc, " columns ",
             "(num.lv.c + num.lv).")
      ## Scale by sigma in standard [lvc|lv] order
      sigma_std <- sigma[c(seq(num.RR + 1L, num.RR + num.lvc),
                           seq(num.RR + num.lvc + 1L, d))]
      lvs_sc <- t(t(newLV) * sigma_std)
      ## Expand to full-d standard [lvc|RR|lv] order to match theta
      lvs_full_std <- cbind(
        if (num.lvc > 0L)    lvs_sc[, seq_len(num.lvc), drop = FALSE] else NULL,
        matrix(0.0, n, num.RR),
        if (num.lv.unc > 0L) lvs_sc[, num.lvc + seq_len(num.lv.unc), drop = FALSE] else NULL
      )
      eta <- eta + lvs_full_std %*% t(object$params$theta)
      ## Add lv.X covariate contribution for lvc/RR dims
      if ((num.lvc + num.RR) > 0L && !is.null(object$params$LvXcoef)) {
        lv.X_use <- if (!is.null(newX)) .lv_x_for_new(newX) else object$lv.X.design
        if (!is.null(lv.X_use)) {
          theta_lvc_rr <- object$params$theta[, seq_len(num.lvc + num.RR), drop = FALSE]
          eta <- eta + lv.X_use %*% object$params$LvXcoef %*% t(theta_lvc_rr)
        }
      }
    } else if (is.null(newX)) {
      ## Training data: use posterior means directly
      ## eta_bilinear = sigma*z (std order) %*% t(theta_std)
      eta <- eta + getLV(object, type = "conditional") %*% t(object$params$theta)
    } else {
      ## New sites without newLV
      if ((num.lv.unc + num.lvc) > 0L)
        stop("Level-1 predictions for new sites require 'newLV' ",
             "when the model has VA dims (num.lv or num.lv.c > 0). ",
             "Use level = 0 for marginal predictions.")
      ## Pure RR model: fully deterministic
      if (num.RR > 0L && !is.null(object$params$LvXcoef)) {
        lv.X_use <- .lv_x_for_new(newX)
        if (!is.null(lv.X_use)) {
          theta_lvc_rr <- object$params$theta[, seq_len(num.lvc + num.RR), drop = FALSE]
          eta <- eta + lv.X_use %*% object$params$LvXcoef %*% t(theta_lvc_rr)
        }
      }
    }
  }

  ## --- row effects ---------------------------------------------------------
  if (!is.null(object$params$row.params.fixed) && is.null(newX)) {
    xr_use <- if (!is.null(object$xr)) object$xr else object$TMBfn$env$data$xr
    if (!is.null(xr_use) && ncol(xr_use) > 0L && nrow(xr_use) == n) {
      r0 <- as.vector(xr_use %*% as.matrix(object$params$row.params.fixed))
      eta <- eta + matrix(r0, n, p)
    }
  }
  if (!is.null(object$params$row.params.random) && level > 0 && is.null(newX)) {
    dr0_use <- object$TMBfn$env$data$dr0
    if (!is.null(dr0_use) && ncol(dr0_use) > 0L && nrow(dr0_use) == n) {
      r0 <- as.vector(dr0_use %*% as.matrix(object$params$row.params.random))
      eta <- eta + matrix(r0, n, p)
    }
  }

  ## --- offset --------------------------------------------------------------
  if (!isFALSE(offset)) {
    off_mat <- if (is.matrix(offset)) {
      offset
    } else if (!is.null(object$offset) && isTRUE(offset)) {
      object$offset
    } else {
      NULL
    }
    if (!is.null(off_mat) && NROW(off_mat) == n)
      eta <- eta + off_mat
  }

  ## --- type transform ------------------------------------------------------
  if (type == "link") return(eta)

  fam  <- if (length(object$family) == p) object$family else rep(object$family[1L], p)
  lnk  <- if (!is.null(object$link) && length(object$link) == p) object$link else rep("log", p)
  out  <- eta

  for (j in seq_len(p)) {
    out[, j] <- switch(fam[j],
      poisson         = ,
      negative.binomial = ,
      gamma           = ,
      tweedie         = ,
      exponential     = exp(eta[, j]),
      ZIP             = ,
      ZINB            = exp(eta[, j]),
      gaussian        = eta[, j],
      binomial        = ,
      beta            = ,
      ZIB             = ,
      ZNIB            = switch(lnk[j],
                               logit   = plogis(eta[, j]),
                               probit  = pnorm(eta[, j]),
                               cloglog = 1 - exp(-exp(eta[, j])),
                               plogis(eta[, j])),
      eta[, j]  # default: identity
    )
  }
  out
}
