#' @title Random coefficient plot for hierarchical ordination models
#' @description Caterpillar plot of per-species environmental slopes for
#'   \code{gllvmHO} objects with \code{randomB = "LV"} (random canonical
#'   coefficients). The effective slope for predictor \eqn{k} and species
#'   \eqn{j} is \eqn{\beta_{jk} = \mathbf{L}_{k\cdot}^\top \boldsymbol{\gamma}_j}
#'   where \eqn{\mathbf{L} = } \code{LvXcoef}. SE accounts for VA uncertainty
#'   in both \eqn{\boldsymbol{\gamma}_j} (species loadings) and
#'   \eqn{\mathbf{b}_z} (canonical coefficients) via the law of total variance
#'   (van der Veen et al. 2022, Appendix S1 eq 4).
#'
#' @param object a fitted \code{gllvmHO} object with \code{randomB = "LV"}.
#' @param y.label logical; if \code{TRUE} species names are printed on the y-axis.
#' @param which.Xcoef character or integer vector selecting predictors to plot.
#'   Defaults to all.
#' @param order logical; whether to order species by point estimate.
#' @param cex.ylab magnification for y-axis labels.
#' @param cex.xlab magnification for x-axis labels.
#' @param mfrow same as \code{par(mfrow)}.
#' @param mar same as \code{par(mar)}.
#' @param xlim.list list of x-axis limits, one per predictor.
#' @param ind.spp integer vector of species to include.
#' @param ... additional graphical arguments passed to \code{plot}.
#'
#' @return Invisibly returns a list with elements \code{coef} (p×Kz matrix of
#'   point estimates) and \code{se} (p×Kz matrix of standard errors).
#'
#' @seealso \code{\link{coefplot.gllvmHO}}
#' @method randomCoefplot gllvmHO
#' @export
randomCoefplot.gllvmHO <- function(object,
                                    y.label      = TRUE,
                                    which.Xcoef  = NULL,
                                    order        = TRUE,
                                    cex.ylab     = 0.5,
                                    cex.xlab     = 1.3,
                                    mfrow        = NULL,
                                    mar          = c(4, 6, 2, 1),
                                    xlim.list    = NULL,
                                    ind.spp      = NULL,
                                    ...) {
  if (!identical(object$randomB, "LV"))
    stop("randomCoefplot for gllvmHO requires randomB = 'LV'.")
  if (is.null(object$params$LvXcoef))
    stop("No canonical coefficient matrix (LvXcoef) found in the model.")
  if (is.null(object$Ab.lv))
    stop("Ab.lv (VA covariance of b_z) not found — refit with randomB = 'LV'.")

  p      <- ncol(object$y)
  nRR    <- object$num.RR
  nlvc   <- object$num.lv.c
  d_act  <- nRR + nlvc
  Kz     <- nrow(object$params$LvXcoef)

  if (is.null(ind.spp)) ind.spp <- seq_len(p)

  ## Point estimates: p × Kz matrix (beta = theta[,1:d_act] %*% t(LvXcoef))
  theta_act <- object$params$theta[, seq_len(d_act), drop = FALSE]
  coef_mat  <- theta_act %*% t(object$params$LvXcoef)   # p × Kz

  ## Standard errors: total variance over both b_z and theta_j
  se_mat <- .ho_randomRRse(object)   # p × Kz

  ## Predictor selection
  pred_names <- rownames(object$params$LvXcoef)
  if (is.null(pred_names)) pred_names <- paste0("cov", seq_len(Kz))
  if (is.null(which.Xcoef)) which.Xcoef <- seq_len(Kz)
  if (is.character(which.Xcoef)) which.Xcoef <- match(which.Xcoef, pred_names)

  coef_mat <- coef_mat[ind.spp, which.Xcoef, drop = FALSE]
  se_mat   <- se_mat[ind.spp, which.Xcoef, drop = FALSE]
  cnames   <- pred_names[which.Xcoef]
  spp_names <- if (!is.null(colnames(object$y))) colnames(object$y)[ind.spp] else
               paste0("sp", ind.spp)
  m <- length(ind.spp)
  k <- length(which.Xcoef)

  if (is.null(mfrow) && k > 1) mfrow <- c(1, k)
  if (!is.null(mfrow)) par(mfrow = mfrow, mar = mar) else par(mar = mar)

  for (i in seq_len(k)) {
    Xc    <- coef_mat[, i]
    lower <- Xc - 1.96 * se_mat[, i]
    upper <- Xc + 1.96 * se_mat[, i]

    if (order) {
      ord <- order(Xc)
      Xc  <- Xc[ord]; lower <- lower[ord]; upper <- upper[ord]
      sn  <- spp_names[ord]
    } else {
      sn <- spp_names
    }

    col.seq <- rep("black", m)
    col.seq[lower < 0 & upper > 0] <- "grey"
    At.y <- seq_len(m)

    xlim <- if (!is.null(xlim.list[[i]])) xlim.list[[i]] else
            c(min(lower, na.rm = TRUE), max(upper, na.rm = TRUE))
    plot(x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq,
         xlab = cnames[i], xlim = xlim, pch = "x", cex.lab = cex.xlab, ...)
    segments(x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq)
    abline(v = 0, lty = 1)
    if (y.label) axis(2, at = At.y, labels = sn, las = 1, cex.axis = cex.ylab)
  }

  invisible(list(coef = coef_mat, se = se_mat))
}

## Full variance for beta_jk = sum_l f_{kl} * g_l, f = sigma*b_z (random), g = theta_j (random),
## independent under factored q (van der Veen et al. 2022 Appendix S1 eq 4):
##
##   Var(beta_jk)
##     = sum_{l,l'} L[k,l]*L[k,l'] * B_out[j,ia(l),ia(l')]   [theta uncertainty]
##     + sum_l sig_l^2 * Ab.lv[ablv(l),k,k] * B_diag[j,ia(l)] [var_f × var_g]
##     + sum_l sig_l^2 * Ab.lv[ablv(l),k,k] * theta[j,l]^2    [var_f × E[g]^2]
##
## Indices:
##   ia(l)    — B_out column for standard dim l (0 = deterministic theta, skip var_g term)
##   ablv(l)  — Ab.lv first-dim index for standard dim l (HO order: [RR|lvc])
#' @keywords internal
.ho_randomRRse <- function(object) {
  p    <- ncol(object$y)
  nRR  <- object$num.RR
  nlvc <- object$num.lv.c
  d_act <- nRR + nlvc
  Kz   <- nrow(object$params$LvXcoef)
  Kt   <- if (!is.null(object$TR)) ncol(as.matrix(object$TR)) else 0L
  rr_va_a <- if (Kt > 0L) 0L else nRR

  ## Standard column l → B_out VA column ia (0 = deterministic)
  ia_for_std <- c(
    rr_va_a + seq_len(nlvc),
    if (nRR > 0L) {
      if (Kt == 0L) seq_len(nRR) else rep(0L, nRR)
    } else integer(0)
  )

  ## Standard column l → Ab.lv slice index (HO order: [RR|lvc])
  ## l=1..nlvc (lvc): ablv = nRR + l
  ## l=nlvc+1..nlvc+nRR (RR): ablv = l - nlvc
  ablv_for_std <- c(nRR + seq_len(nlvc),
                    if (nRR > 0L) seq_len(nRR) else integer(0))

  ## sigma.lv in standard [lvc|RR|lv] order; first d_act entries cover the active dims
  sigma_std <- object$params$sigma.lv[seq_len(d_act)]
  sig2      <- sigma_std^2

  L      <- object$params$LvXcoef    # Kz × d_act, L[k,l] = sigma_l * b_z_hat[k,l]
  theta  <- object$params$theta[, seq_len(d_act), drop = FALSE]  # p × d_act

  Ab.lv   <- object$Ab.lv   # (d_c, Kz, Kz) where d_c = nRR + nlvc
  B_out   <- object$B
  is_full <- length(dim(B_out)) == 3L
  B_diag  <- .ho_diag(B_out)          # p × d_va_a

  ## Precompute Ab.lv diagonal per (l, k): var_f_kl = sig_l^2 * Ab.lv[ablv(l), k, k]
  ## Shape: d_act × Kz — row l, col k
  var_f <- matrix(0, d_act, Kz)
  for (l in seq_len(d_act)) {
    al <- ablv_for_std[l]
    if (al >= 1L && al <= dim(Ab.lv)[1L])
      var_f[l, ] <- sig2[l] * Ab.lv[al, , ][ cbind(seq_len(Kz), seq_len(Kz)) ]
  }
  # var_f[l,k] = Var_q(sigma_l * b_z[k,l])

  se_mat <- matrix(0, p, Kz)
  rownames(se_mat) <- if (!is.null(colnames(object$y))) colnames(object$y) else paste0("sp", seq_len(p))
  colnames(se_mat) <- if (!is.null(rownames(L))) rownames(L) else paste0("cov", seq_len(Kz))

  ## Valid (non-deterministic theta) dims
  valid  <- which(ia_for_std > 0L)
  ia_val <- ia_for_std[valid]
  L_val  <- L[, valid, drop = FALSE]   # Kz × n_valid

  for (k in seq_len(Kz)) {
    lk <- L_val[k, ]   # n_valid elements

    ## Term 1: theta uncertainty (same as .ho_RRse for the valid dims)
    if (is_full) {
      var_theta <- vapply(seq_len(p), function(j) {
        Bj <- B_out[j, ia_val, ia_val, drop = FALSE]
        as.numeric(t(lk) %*% matrix(Bj, length(ia_val), length(ia_val)) %*% lk)
      }, numeric(1))
    } else {
      var_theta <- as.vector(B_diag[, ia_val, drop = FALSE] %*% lk^2)
    }

    ## Term 2: var_f × var_g (only for valid ia dims)
    var_fvg <- numeric(p)
    for (v in seq_along(valid)) {
      l  <- valid[v]
      ia <- ia_val[v]
      var_fvg <- var_fvg + var_f[l, k] * B_diag[, ia]
    }

    ## Term 3: var_f × E[g]^2 (all dims — even deterministic theta from traits)
    var_feg2 <- numeric(p)
    for (l in seq_len(d_act)) {
      var_feg2 <- var_feg2 + var_f[l, k] * theta[, l]^2
    }

    se_mat[, k] <- sqrt(pmax(0, var_theta + var_fvg + var_feg2))
  }
  se_mat
}
