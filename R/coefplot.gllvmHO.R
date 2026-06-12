#' @title Coefficient plot for hierarchical ordination models
#' @description Caterpillar plot of per-species environmental slopes for
#'   \code{gllvmHO} objects with \code{randomB = FALSE} (fixed canonical
#'   coefficients). The effective species slope for predictor \eqn{k} and
#'   species \eqn{j} is \eqn{\beta_{jk} = \mathbf{L}_{k\cdot}^\top \boldsymbol{\gamma}_j},
#'   where \eqn{\mathbf{L} = } \code{LvXcoef} and \eqn{\boldsymbol{\gamma}_j} are
#'   the species loadings.  SE comes from the VA posterior covariance of
#'   \eqn{\boldsymbol{\gamma}_j} (van der Veen et al. 2022, Appendix S1 eq 4
#'   simplified for fixed B).
#'
#' @param object a fitted \code{gllvmHO} object with \code{randomB = FALSE}.
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
#' @seealso \code{\link{randomCoefplot.gllvmHO}}
#' @method coefplot gllvmHO
#' @export
coefplot.gllvmHO <- function(object,
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
  if (!isFALSE(object$randomB))
    stop("coefplot for gllvmHO requires randomB = FALSE.")
  if (is.null(object$params$LvXcoef))
    stop("No canonical coefficient matrix (LvXcoef) found in the model.")

  p      <- ncol(object$y)
  nRR    <- object$num.RR
  nlvc   <- object$num.lv.c
  d_act  <- nRR + nlvc
  Kz     <- nrow(object$params$LvXcoef)   # number of env predictors
  Kt     <- if (!is.null(object$TR)) ncol(as.matrix(object$TR)) else 0L

  if (is.null(ind.spp)) ind.spp <- seq_len(p)

  ## Point estimates: p × Kz matrix (beta = theta[,1:d_act] %*% t(LvXcoef))
  theta_act <- object$params$theta[, seq_len(d_act), drop = FALSE]
  coef_mat  <- theta_act %*% t(object$params$LvXcoef)   # p × Kz

  ## Standard errors using VA posterior covariance of theta_j
  se_mat <- .ho_RRse(object)   # p × Kz

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

## Internal SE helper: Var_q(beta_{jk}) = L[k,]^T B_j L[k,] (fixed B = b_z)
## where B_j = VA posterior covariance of theta_j (from object$B = A_lv_out)
#' @keywords internal
.ho_RRse <- function(object) {
  p    <- ncol(object$y)
  nRR  <- object$num.RR
  nlvc <- object$num.lv.c
  d_act <- nRR + nlvc
  Kz   <- nrow(object$params$LvXcoef)
  Kt   <- if (!is.null(object$TR)) ncol(as.matrix(object$TR)) else 0L
  rr_va_a <- if (Kt > 0L) 0L else nRR

  ## Map standard [lvc|RR] column l → VA column ia in B_out
  ## lvc dims (l=1..nlvc): ia = rr_va_a + l
  ## RR  dims (l=nlvc+1..nlvc+nRR, Kt=0): ia = l - nlvc; Kt>0: ia = 0 (deterministic)
  ia_for_std <- c(
    rr_va_a + seq_len(nlvc),
    if (nRR > 0L) {
      if (Kt == 0L) seq_len(nRR) else rep(0L, nRR)
    } else integer(0)
  )

  valid   <- which(ia_for_std > 0L)
  ia_val  <- ia_for_std[valid]
  L       <- object$params$LvXcoef[, valid, drop = FALSE]  # Kz × n_valid
  B_out   <- object$B                                        # p × d_va_a or p×d_va_a×d_va_a
  is_full <- length(dim(B_out)) == 3L
  B_diag  <- .ho_diag(B_out)                                # p × d_va_a

  se_mat <- matrix(0, p, Kz)
  rownames(se_mat) <- if (!is.null(colnames(object$y))) colnames(object$y) else paste0("sp", seq_len(p))
  colnames(se_mat) <- if (!is.null(rownames(object$params$LvXcoef))) rownames(object$params$LvXcoef) else paste0("cov", seq_len(Kz))

  for (k in seq_len(Kz)) {
    lk <- L[k, ]  # n_valid-vector
    if (is_full) {
      for (j in seq_len(p)) {
        Bj <- B_out[j, ia_val, ia_val, drop = FALSE]
        Bj <- matrix(Bj, length(ia_val), length(ia_val))
        se_mat[j, k] <- sqrt(pmax(0, as.numeric(t(lk) %*% Bj %*% lk)))
      }
    } else {
      ## Diagonal B: SE^2 = sum_l L[k,l]^2 * B_diag[j, ia(l)]
      se_mat[, k] <- sqrt(pmax(0, as.vector(B_diag[, ia_val, drop = FALSE] %*% lk^2)))
    }
  }
  se_mat
}
