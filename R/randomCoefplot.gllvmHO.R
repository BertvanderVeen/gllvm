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

## Full variance for beta_jk = sum_l f_l * g_l, f_l = sigma_l*b_z[k,l] (random),
## g_l = gamma_j[l] (random), independent under factored q.
## Var(f_l * g_l) = E[f_l]^2 Var(g_l) + E[g_l]^2 Var(f_l) + Var(f_l)Var(g_l).
## Since f and g dims are VA-independent, summing gives:
##
##   Var(beta_jk) = L_k^T Cov(gamma_j) L_k          [Term 1: full quadratic in species covariance]
##                + theta_j^T C_k theta_j            [Term 3: full quadratic in b_z covariance]
##                + sum_l C_k[l,l] * Var(gamma_j[l]) [Term 2: diagonal because g dims independent]
##
## The full Isserlis formula requires the joint CMSEP of (b_z, a_lv_sp) from
## CMSEPf_HO_bilinear() — the cross-covariance Cov(b_z[k,l], a_lv_sp[j,l'])
## cannot be recovered from separate marginal CMSEPs.  Falls back to
## VA-only independence formula (no CMSEP correction) when Hessian is absent.
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

  ## CMSEP-corrected gamma (species loading) variances: p × d_total (HO order SDs)
  ## getPredictErr adds the Hessian Schur-complement correction to B_out when a
  ## Hessian is available.  Fall back to raw B_out diagonal if Hessian absent.
  d_total <- nRR + nlvc + (object$num.lv %||% 0L)
  d_va_a  <- ncol(.ho_diag(B_out))
  if (!is.null(object$Hess)) {
    pe     <- getPredictErr(object, CMSEP = TRUE, cov = FALSE)
    ## pe$loadings is p × d_total SDs in HO order; we need the d_va_a VA-active columns.
    ## VA-active species columns sit at the END of the d_total HO-order block.
    ho_idx_a <- seq(d_total - d_va_a + 1L, d_total)
    B_diag   <- pe$loadings[, ho_idx_a, drop = FALSE]^2   # p × d_va_a variances
  } else {
    B_diag  <- .ho_diag(B_out)   # p × d_va_a (raw VA variances)
  }

  ## Joint CMSEP for (b_z, a_lv_sp): cannot be computed separately because the
  ## cross-covariance Cov_CMSEP(b_z[k,l], a_lv_sp[j,l']) only appears when both
  ## are included jointly in C_va = H[{b_z,a_lv_sp}, fixed].
  jcmsep <- if (!is.null(object$Hess))
    tryCatch(CMSEPf_HO_bilinear(object), error = function(e) NULL)
  else NULL
  have_joint <- !is.null(jcmsep)

  d_c_ab <- dim(Ab.lv)[1L]

  ## Which standard dims have a valid Ab.lv entry (b_z VA covariance exists)
  valid_l  <- ablv_for_std >= 1L & ablv_for_std <= d_c_ab
  bz_dims  <- which(valid_l)   # indices in standard order with b_z

  ## Valid (non-deterministic theta) dims — where a_lv_sp VA covariance exists
  valid  <- which(ia_for_std > 0L)
  ia_val <- ia_for_std[valid]
  L_val  <- L[, valid, drop = FALSE]   # Kz × n_valid

  ## Position of valid dims within bz_dims (for the Isserlis tr terms)
  ## valid ⊆ 1..d_act = bz_dims in all typical cases
  valid_in_bz <- match(valid, bz_dims)

  se_mat <- matrix(0, p, Kz)
  rownames(se_mat) <- if (!is.null(colnames(object$y))) colnames(object$y) else paste0("sp", seq_len(p))
  colnames(se_mat) <- if (!is.null(rownames(L))) rownames(L) else paste0("cov", seq_len(Kz))

  for (k in seq_len(Kz)) {
    lk <- L_val[k, ]   # LvXcoef[k, valid] = sigma*b_z_hat[k, valid]

    if (have_joint) {
      n_bz <- jcmsep$n_bz
      jmat <- jcmsep$joint

      ## Indices of b_z[k, bz_dims] in the joint CMSEP (b_z block, col-major)
      bz_idx_k <- (ablv_for_std[bz_dims] - 1L) * Kz + k   # length(bz_dims)

      ## sigma-scaled b_z sub-block (constant across species j):
      ## Sff[l,l'] = sigma_l * sigma_l' * Cov_CMSEP(b_z[k,l], b_z[k,l'])
      sigma_f <- sigma_std[bz_dims]
      Sff     <- outer(sigma_f, sigma_f) * jmat[bz_idx_k, bz_idx_k, drop = FALSE]

      ## Per-species quantities that depend on species j
      var_jk <- vapply(seq_len(p), function(j) {
        ## Indices of a_lv_sp[j, ia_val] in the joint CMSEP (asp block, col-major)
        asp_j_idx <- n_bz + (ia_val - 1L) * p + j   # length(valid)

        ## a_lv_sp sub-block (NOT sigma-scaled; theta is already unscaled a_lv_sp)
        Sgg <- jmat[asp_j_idx, asp_j_idx, drop = FALSE]

        ## Cross-block: Sfg[l,l'] = sigma_l * Cov_CMSEP(b_z[k,l], a_lv_sp[j,ia(l')])
        ## (sigma_l on the b_z side only; a_lv_sp enters unscaled via lk)
        Sfg <- outer(sigma_f, rep(1, length(valid))) *
               jmat[bz_idx_k, asp_j_idx, drop = FALSE]

        ## Gradients of beta_jk = sum_l LvXcoef[k,l] * theta[j,l]:
        ##   d/d(b_z[k,l])      = sigma_l * theta[j,l]  (stored in sigma_f * theta[j,l])
        ##   d/d(a_lv_sp[j,ia]) = LvXcoef[k,l] = lk
        grad_f <- theta[j, bz_dims]   # unscaled a_lv_sp (sigma absorbed in Sff)
        grad_g <- lk                  # sigma-scaled b_z means

        ## Isserlis formula for Var(f^T g) with jointly Gaussian (f, g):
        ##   Delta terms: grad_f^T Sff grad_f  +  grad_g^T Sgg grad_g  +  2 grad_f^T Sfg grad_g
        ##   Trace terms: sum(Sff_valid * Sgg)  +  sum(Sfg_valid * t(Sfg_valid))
        ## where "valid" restricts to dims where BOTH f and g have VA variance
        t3  <- as.numeric(t(grad_f) %*% Sff %*% grad_f)
        t1  <- as.numeric(t(grad_g) %*% Sgg %*% grad_g)
        tc  <- 2 * as.numeric(t(grad_f) %*% Sfg %*% grad_g)
        ## Isserlis trace terms — only over dims where g has VA variance (valid dims)
        Sff_v   <- Sff[valid_in_bz, valid_in_bz, drop = FALSE]
        Sfg_v   <- Sfg[valid_in_bz,  , drop = FALSE]
        t2a <- sum(Sff_v * Sgg)
        t2b <- sum(Sfg_v * t(Sfg_v))
        t1 + t3 + tc + t2a + t2b
      }, numeric(1))

    } else {
      ## Fallback: VA-only, no CMSEP correction, dims independent.
      ## Var(beta_jk) = sum_l [E[f_l]^2 Var(g_l) + E[g_l]^2 Var(f_l) + Var(f_l)Var(g_l)]
      var_f_k <- numeric(d_act)   # Var(sigma_l * b_z[k,l]) from Ab.lv diagonal
      for (l in bz_dims) {
        al <- ablv_for_std[l]
        var_f_k[l] <- sig2[l] * matrix(Ab.lv[al, , ], Kz, Kz)[k, k]
      }
      var_jk <- vapply(seq_len(p), function(j) {
        t1  <- sum(lk^2 * B_diag[j, ia_val])                       # E[f]^2 Var(g)
        t3  <- sum(theta[j, bz_dims]^2 * var_f_k[bz_dims])         # E[g]^2 Var(f)
        t2  <- sum(var_f_k[valid] * B_diag[j, ia_val])             # Var(f)Var(g)
        t1 + t2 + t3
      }, numeric(1))
    }

    se_mat[, k] <- sqrt(pmax(0, var_jk))
  }
  se_mat
}
