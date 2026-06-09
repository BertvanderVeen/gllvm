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
  d      <- ncol(object$lvs.full)   # total dims
  num.RR  <- as.integer(object$num.RR  %||% 0L)
  num.lvc <- as.integer(object$num.lv.c %||% 0L)
  num.lv  <- as.integer(object$num.lv  %||% 0L)  # unconstrained
  sigma   <- object$params$sigma.lv
  Kz <- if (!is.null(object$lv.X.design)) ncol(as.matrix(object$lv.X.design)) else 0L

  if ((num.RR + num.lvc + num.lv) == 0L)
    stop("No latent variables in model.")
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
  if (type == "conditional" && num.lvc == 0L)
    stop("'conditional' scores require num.lv.c > 0.")
  if (type == "residual" && num.lvc == 0L && num.lv == 0L)
    stop("'residual' scores require num.lv.c > 0 or num.lv > 0.")
  if (type == "marginal" && num.lvc == 0L && num.RR == 0L)
    stop("'marginal' scores require num.RR > 0 or num.lv.c > 0.")
  if (type == "marginal" && Kz == 0L)
    stop("'marginal' scores require canonical covariates (lv.X / lv.formula).")

  ## --- Compute scores -------------------------------------------------------
  if (type == "conditional") {
    ## Full sigma-scaled z (posterior means), all d dims, HO [RR|lvc|lv] order
    lvs <- t(t(object$lvs.full) * sigma)
    ## Reorder to standard [lvc | RR | lv]
    if (num.RR > 0L && num.lvc > 0L) {
      idx <- c(seq(num.RR + 1L, num.RR + num.lvc),
               seq_len(num.RR),
               if (num.lv > 0L) seq(num.RR + num.lvc + 1L, d) else integer(0))
      lvs <- lvs[, idx, drop = FALSE]
    }
    n_clv <- num.lvc + num.RR
    if (n_clv > 0L && num.lv > 0L)
      colnames(lvs) <- c(paste0("CLV", seq_len(n_clv)), paste0("LV", seq_len(num.lv)))
    else if (n_clv > 0L)
      colnames(lvs) <- paste0("CLV", seq_len(n_clv))
    else
      colnames(lvs) <- paste0("LV", seq_len(num.lv))

  } else if (type == "marginal") {
    ## Sigma-scaled covariate prior mean: sigma_k * lv_X * b_z[:,k]
    lv_X <- as.matrix(object$lv.X.design)
    b_z  <- object$params$b_z   # Kz x d, HO [RR|lvc|lv] order
    d_active <- num.RR + num.lvc
    lvs_det <- lv_X %*% b_z[, seq_len(d_active), drop = FALSE]  # n x d_active
    lvs_det <- t(t(lvs_det) * sigma[seq_len(d_active)])
    ## Reorder from HO [RR|lvc] to standard [lvc|RR]
    if (num.RR > 0L && num.lvc > 0L) {
      lvs <- lvs_det[, c(seq(num.RR + 1L, d_active), seq_len(num.RR)), drop = FALSE]
    } else {
      lvs <- lvs_det
    }
    colnames(lvs) <- paste0("CLV", seq_len(ncol(lvs)))

  } else {  ## "residual"
    ## Sigma-scaled VA residuals (posterior mean minus prior mean for lvc)
    ## object$lvs stores these in [lvc | lv] order, NOT sigma-scaled
    lvs_va <- object$lvs   # n x (num.lvc + num.lv)
    ## sigma for lvc and lv dims (in [lvc|lv] standard order)
    sigma_lvc <- if (num.lvc > 0L) sigma[seq(num.RR + 1L, num.RR + num.lvc)] else numeric(0)
    sigma_lv  <- if (num.lv  > 0L) sigma[seq(num.RR + num.lvc + 1L, d)]      else numeric(0)
    sigma_va  <- c(sigma_lvc, sigma_lv)
    lvs <- t(t(lvs_va) * sigma_va)
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
