#' @title Extract latent variables from a hierarchical ordination model
#' @description
#' Extracts site scores from a \code{gllvmHO} object.  The three types mirror
#' \code{\link{getLV.gllvm}} but are computed from the HO parameterisation
#' (\eqn{z_i^\top \Sigma \gamma_j}):
#'
#' \describe{
#'   \item{\code{"conditional"}}{Sigma-scaled full site scores for all d dims
#'     (posterior means of \eqn{\sigma_k z_{ik}}), ordered in standard
#'     \[lvc | RR | lv\] column order.  Default when \code{num.lv.c > 0}.}
#'   \item{\code{"marginal"}}{Sigma-scaled deterministic (covariate-driven)
#'     site scores: \eqn{\sigma_k \cdot X_{lv,i} b_{z,k}} for lvc and RR dims.
#'     Requires canonical covariates.  Default when only \code{num.RR > 0}.}
#'   \item{\code{"residual"}}{Sigma-scaled VA residual site scores
#'     (\eqn{\sigma_k (u_{ik} - X_{lv,i} b_{z,k})} for lvc,
#'      \eqn{\sigma_k u_{ik}} for unconstrained lv).
#'     Default when only \code{num.lv > 0}.}
#' }
#'
#' @param object A \code{gllvmHO} object.
#' @param type  One of \code{"conditional"}, \code{"marginal"},
#'   \code{"residual"}.  \code{NULL} picks a sensible default.
#' @param ... Not used.
#'
#' @return A site-score matrix (\eqn{n \times k}) with appropriate column names.
#'
#' @method getLV gllvmHO
#' @export
#' @export getLV.gllvmHO
getLV.gllvmHO <- function(object, type = NULL, ...) {
  n      <- nrow(object$y)
  num.RR  <- as.integer(object$num.RR  %||% 0L)
  num.lvc <- as.integer(object$num.lv.c %||% 0L)
  num.lv  <- as.integer(object$num.lv  %||% 0L)  # unconstrained
  d       <- num.RR + num.lvc + num.lv
  sigma   <- object$params$sigma.lv
  lv_X    <- if (!is.null(object$lv.X.design)) as.matrix(object$lv.X.design) else NULL
  Kz      <- if (!is.null(lv_X)) ncol(lv_X) else 0L
  ## VA-index offset: when Kz==0, RR dims are also in VA space
  rr_va_z <- if (Kz > 0L) 0L else num.RR
  ## LvXcoef: sigma-scaled b_z, std [lvc|RR] order
  lvxc <- object$params$LvXcoef  # Kz × (num.lvc + num.RR)

  if (d == 0L) stop("No latent variables in model.")
  if (!is.null(type) && !type %in% c("residual", "conditional", "marginal"))
    stop("type must be one of: residual, conditional, marginal.")

  ## --- Default type ---------------------------------------------------------
  if (is.null(type)) {
    if (num.lvc == 0L && num.RR == 0L) {
      type <- "residual"
    } else if (num.lvc > 0L) {
      type <- "conditional"
    } else {
      type <- "marginal"
    }
  }

  ## --- Validate -------------------------------------------------------------
  if (type == "conditional" && num.lvc == 0L && num.RR == 0L && num.lv == 0L)
    stop("'conditional' scores require at least one latent variable.")
  if (type == "residual" && num.lvc == 0L && num.lv == 0L)
    stop("'residual' scores require num.lv.c > 0 or num.lv > 0.")
  if (type == "marginal" && num.lvc == 0L && num.RR == 0L)
    stop("'marginal' scores require num.RR > 0 or num.lv.c > 0.")
  if (type == "marginal" && Kz == 0L)
    stop("'marginal' scores require canonical covariates (lv.X / lv.formula).")

  ## --- Compute scores -------------------------------------------------------
  if (type == "conditional") {
    ## Full unscaled z in HO [RR|lvc|lv] order, then scale and reorder to std
    z_full <- matrix(0.0, n, d)

    if (num.RR > 0L) {
      if (Kz > 0L && !is.null(lvxc) && !is.null(lv_X)) {
        ## deterministic: z = lv_X * b_z  (unscale LvXcoef by sigma)
        rr_cols_std <- num.lvc + seq_len(num.RR)
        b_z_rr <- sweep(lvxc[, rr_cols_std, drop = FALSE],
                        2L, sigma[seq_len(num.RR)], `/`)
        z_full[, seq_len(num.RR)] <- lv_X %*% b_z_rr
      } else {
        ## VA: stored as first rr_va_z cols of object$lvs
        z_full[, seq_len(num.RR)] <- object$lvs[, seq_len(num.RR), drop = FALSE]
      }
    }

    if (num.lvc > 0L) {
      iz <- rr_va_z + seq_len(num.lvc)
      z_full[, num.RR + seq_len(num.lvc)] <- object$lvs[, iz, drop = FALSE]
    }

    if (num.lv > 0L) {
      iz <- rr_va_z + num.lvc + seq_len(num.lv)
      z_full[, num.RR + num.lvc + seq_len(num.lv)] <- object$lvs[, iz, drop = FALSE]
    }

    ## Scale by sigma (HO order), then reorder to standard [lvc|RR|lv]
    lvs_ho <- t(t(z_full) * sigma)
    if (num.RR > 0L && num.lvc > 0L) {
      idx <- c(seq(num.RR + 1L, num.RR + num.lvc),
               seq_len(num.RR),
               if (num.lv > 0L) seq(num.RR + num.lvc + 1L, d) else integer(0))
      lvs <- lvs_ho[, idx, drop = FALSE]
    } else {
      lvs <- lvs_ho
    }

    n_clv <- num.lvc + num.RR
    if (n_clv > 0L && num.lv > 0L)
      colnames(lvs) <- c(paste0("CLV", seq_len(n_clv)), paste0("LV", seq_len(num.lv)))
    else if (n_clv > 0L)
      colnames(lvs) <- paste0("CLV", seq_len(n_clv))
    else
      colnames(lvs) <- paste0("LV", seq_len(num.lv))

  } else if (type == "marginal") {
    ## LvXcoef = sigma * b_z, standard [lvc|RR] order — use directly
    lvs  <- lv_X %*% lvxc   # n x (num.lvc + num.RR), already sigma-scaled
    colnames(lvs) <- paste0("CLV", seq_len(ncol(lvs)))

  } else {  ## "residual"
    ## Sigma-scaled residual: u_hat - lv.X * b_z for lvc; raw u_hat for lv

    lvc_part <- if (num.lvc > 0L) {
      iz <- rr_va_z + seq_len(num.lvc)
      raw <- object$lvs[, iz, drop = FALSE]
      if (Kz > 0L && !is.null(lvxc) && !is.null(lv_X)) {
        b_z_lvc <- sweep(lvxc[, seq_len(num.lvc), drop = FALSE],
                         2L, sigma[num.RR + seq_len(num.lvc)], `/`)
        raw <- raw - lv_X %*% b_z_lvc
      }
      sweep(raw, 2L, sigma[num.RR + seq_len(num.lvc)], `*`)
    } else matrix(0.0, n, 0L)

    lv_part <- if (num.lv > 0L) {
      iz <- rr_va_z + num.lvc + seq_len(num.lv)
      sweep(object$lvs[, iz, drop = FALSE],
            2L, sigma[num.RR + num.lvc + seq_len(num.lv)], `*`)
    } else matrix(0.0, n, 0L)

    lvs <- cbind(lvc_part, lv_part)

    n_clv <- num.lvc
    if (n_clv > 0L && num.lv > 0L)
      colnames(lvs) <- c(paste0("CLV", seq_len(n_clv)), paste0("LV", seq_len(num.lv)))
    else if (n_clv > 0L)
      colnames(lvs) <- paste0("CLV", seq_len(n_clv))
    else
      colnames(lvs) <- paste0("LV", seq_len(num.lv))
  }

  rownames(lvs) <- rownames(object$y)
  lvs
}
