#' Ordination plot for hierarchical ordination models
#'
#' Plots site scores and (optionally) species scores from a \code{gllvmHO}
#' model.  \code{alpha} distributes the singular values \eqn{\Sigma} between
#' sites (\eqn{\sigma^\alpha}) and species (\eqn{\sigma^{1-\alpha}}), so the
#' inner product recovers \eqn{z_i^\top \Sigma a_j}.  No post-hoc rotation is
#' applied to constrained axes; unconstrained axes can optionally be rotated to
#' their principal direction via SVD.
#'
#' @param object   A \code{gllvmHO} object.
#' @param biplot   Logical; add species scores?  Default \code{FALSE}.
#' @param ind.spp  Number of species to show (default: all).
#' @param alpha    In [0,1]. Distributes \eqn{\Sigma}: sites carry
#'   \eqn{\sigma^\alpha}, species carry \eqn{\sigma^{1-\alpha}}.
#'   \code{alpha = 0.5} (default) gives a symmetric biplot.
#' @param main     Plot title.
#' @param which.lvs Integer vector of length 2 selecting axes to plot within
#'   the type-selected dimension set.
#' @param type     Which ordination type to plot.
#'   \code{"conditional"} (default): all \eqn{d} dimensions (VA posterior means
#'   for all dims; deterministic \eqn{X b_z} for RR dims when both X and TR are
#'   present).
#'   \code{"residual"}: only VA (random) dimensions — \code{num.lv.c} and
#'   \code{num.lv} dims; excludes deterministic RR dims.
#'   \code{"marginal"}: only the deterministic RR dimensions (\eqn{X b_z}).
#'   Only meaningful when \code{num.RR > 0} and X is present.
#' @param rotate   Logical (default \code{TRUE}).  If \code{TRUE}, rotate the
#'   selected site scores via SVD before plotting (as in standard ordination).
#'   For constrained or deterministic axes the rotation is suppressed.
#' @param predict.region  \code{TRUE}/\code{"sites"}: ellipses for sites.
#'   \code{"species"}: ellipses for species (implies biplot).
#' @param CMSEP    Logical (default \code{TRUE}).  Use CMSEPf-corrected variances
#'   for ellipses when a Hessian is available.
#' @param level    Ellipse coverage level.  Default 0.95.
#' @param jitter   Logical; jitter site positions?
#' @param jitter.amount  Half-range of jitter.
#' @param s.colors  Colour(s) for site labels/points.
#' @param s.cex     Character expansion for site labels.
#' @param symbols   Use points instead of text for sites?
#' @param cex.spp   Character expansion for species labels.
#' @param spp.colors  Colour for species labels.
#' @param spp.arrows  Plot species scores outside the site-score range as
#'   arrows?  Default \code{FALSE}.
#' @param spp.arrows.lty  Line type for out-of-range species arrows.
#' @param arrow.spp.scale Positive scalar; scale for species arrows.
#' @param lwd.ellips  Line width for site ellipses.
#' @param col.ellips  Colour(s) for site ellipses.
#' @param lty.ellips  Line type for site ellipses.
#' @param col.spp.ellips  Colour for species ellipses.
#' @param lty.spp.ellips  Line type for species ellipses.
#' @param arrow.scale  Positive scalar; arrows are scaled to this fraction of
#'   the half-axis range.  Default 0.8.
#' @param arrow.ci   Logical (default \code{TRUE}).  If \code{TRUE} and standard
#'   errors are available for \eqn{b_z}/\eqn{b_\gamma}, arrows whose 95\% CI
#'   excludes zero are drawn in the full colour; others are drawn in a lighter
#'   shade.
#' @param arrow.lty  Line type for arrows.  Default \code{"solid"}.
#' @param cex.env  Character expansion for covariate and trait arrow labels.
#'   Default 0.7.
#' @param lab.dist  Gap between arrowhead and label as a fraction of arrow
#'   length.  Default 0.05.
#' @param col.arrow.bz    Colour for canonical covariate (\eqn{b_z}) arrows.
#'   Default \code{"red"}.
#' @param col.arrow.bgamma  Colour for trait (\eqn{b_\gamma}) arrows.
#'   Default \code{"darkblue"}.
#' @param ...  Further graphical parameters passed to the plot function.
#'
#' @return Invisibly, a list with \code{sites} (n x 2) and \code{species}
#'   (p x 2) plot coordinates.
#'
#' @method ordiplot gllvmHO
#' @export
ordiplot.gllvmHO <- function(object,
                              biplot         = FALSE,
                              ind.spp        = NULL,
                              alpha          = 0.5,
                              main           = NULL,
                              which.lvs      = c(1, 2),
                              type           = NULL,
                              predict.region = FALSE,
                              CMSEP          = TRUE,
                              level          = 0.95,
                              jitter         = FALSE,
                              jitter.amount  = 0.2,
                              s.colors       = 1,
                              s.cex          = 1.2,
                              symbols        = FALSE,
                              cex.spp        = 0.7,
                              spp.colors     = "blue",
                              spp.arrows     = FALSE,
                              spp.arrows.lty = "dashed",
                              arrow.spp.scale = 0.8,
                              lwd.ellips     = 0.5,
                              col.ellips     = 4,
                              lty.ellips     = 1,
                              col.spp.ellips = "blue",
                              lty.spp.ellips = 2,
                              arrow.scale    = 0.8,
                              arrow.ci       = TRUE,
                              arrow.lty      = "solid",
                              cex.env        = 0.7,
                              lab.dist       = 0.05,
                              col.arrow.bz     = "red",
                              col.arrow.bgamma = "darkblue",
                              ...) {

  n <- nrow(object$y)
  p <- ncol(object$y)
  d <- as.integer(object$num.RR %||% 0L) +
       as.integer(object$num.lv.c %||% 0L) +
       as.integer(object$num.lv %||% 0L)
  if (d < 1L) stop("No latent variables in model.")

  num.RR  <- if (!is.null(object$num.RR))  as.integer(object$num.RR)  else 0L
  num.lvc <- if (!is.null(object$num.lv.c)) as.integer(object$num.lv.c) else 0L
  num.lv  <- object$num.lv   # unconstrained LV dims

  Kz <- if (!is.null(object$lv.X)) ncol(as.matrix(object$lv.X)) else 0L
  Kt <- if (!is.null(object$TR))   ncol(as.matrix(object$TR))   else 0L

  rr_det  <- num.RR > 0L && Kz > 0L
  lvc_det <- num.lvc > 0L && Kz > 0L   # concurrent dims with covariate structure

  sigma <- object$params$sigma.lv   # length d, ordered desc

  ## --- Determine plot type --------------------------------------------------
  if (is.null(type)) {
    if ((rr_det || lvc_det) && (num.lvc + num.lv) > 0L && !(rr_det || lvc_det)) {
      type <- "conditional"
    } else if (rr_det && (num.lvc + num.lv) == 0L) {
      type <- "marginal"
    } else {
      type <- "conditional"
    }
  }
  type <- match.arg(type, c("conditional", "residual", "marginal"))

  if (type == "marginal" && !(rr_det || lvc_det))
    stop("type = 'marginal' requires X covariates (lv.X) with RR or lvc dims.")
  if (type == "residual" && (num.lvc + num.lv) == 0L)
    stop("type = 'residual' requires at least one VA (lv.c or lv) dimension.")

  ## --- Reconstruct full unscaled z in HO [RR|lvc|lv] order on the fly -------
  lv_X    <- if (!is.null(object$lv.X.design)) as.matrix(object$lv.X.design) else NULL
  rr_va_z <- if (Kz > 0L) 0L else num.RR
  lvxc    <- object$params$LvXcoef
  loadings_ho <- getLoadings(object)   # p x d, HO [RR|lvc|lv] order

  .ho_z_full <- function() {
    ## Returns n x d unscaled z in HO order from stored u_hat (object$lvs)
    z <- matrix(0.0, n, d)
    if (num.RR > 0L) {
      if (Kz > 0L && !is.null(lvxc) && !is.null(lv_X)) {
        b_z_rr <- sweep(lvxc[, num.lvc + seq_len(num.RR), drop = FALSE],
                        2L, sigma[seq_len(num.RR)], `/`)
        z[, seq_len(num.RR)] <- lv_X %*% b_z_rr
      } else {
        z[, seq_len(num.RR)] <- object$lvs[, seq_len(num.RR), drop = FALSE]
      }
    }
    if (num.lvc > 0L)
      z[, num.RR + seq_len(num.lvc)] <-
        object$lvs[, rr_va_z + seq_len(num.lvc), drop = FALSE]
    if (num.lv > 0L)
      z[, num.RR + num.lvc + seq_len(num.lv)] <-
        object$lvs[, rr_va_z + num.lvc + seq_len(num.lv), drop = FALSE]
    z
  }

  ## --- Build site scores and loadings for selected type ---------------------
  if (type == "conditional") {
    lv_sel    <- .ho_z_full()
    theta_sel <- loadings_ho
    sigma_sel <- sigma
    col_idx   <- seq_len(d)

  } else if (type == "marginal") {
    ## Covariate-driven scores: z = lv.X * b_z (deterministic part only)
    z_full <- .ho_z_full()
    rr_scores  <- if (num.RR > 0L) z_full[, seq_len(num.RR), drop = FALSE] else NULL
    lvc_scores <- if (num.lvc > 0L) {
      ## lv.X * b_z = full z - residual u_hat
      lvc_full   <- z_full[, num.RR + seq_len(num.lvc), drop = FALSE]
      ## subtract the VA residual (u_hat - lv.X*b_z) to recover lv.X*b_z
      iz_lvc     <- rr_va_z + seq_len(num.lvc)
      resid_part <- if (Kz > 0L && !is.null(lvxc) && !is.null(lv_X)) {
        b_z_lvc <- sweep(lvxc[, seq_len(num.lvc), drop = FALSE],
                         2L, sigma[num.RR + seq_len(num.lvc)], `/`)
        object$lvs[, iz_lvc, drop = FALSE] - lv_X %*% b_z_lvc
      } else {
        matrix(0.0, n, num.lvc)
      }
      lvc_full - resid_part
    } else NULL

    lv_sel    <- cbind(rr_scores, lvc_scores)
    col_ho    <- c(if (num.RR > 0L)  seq_len(num.RR) else NULL,
                   if (num.lvc > 0L) num.RR + seq_len(num.lvc) else NULL)
    theta_sel <- loadings_ho[, col_ho, drop = FALSE]
    sigma_sel <- sigma[col_ho]
    col_idx   <- col_ho

  } else {
    ## type == "residual": VA random part (lvc residual + lv)
    ## Use getLV to get sigma-scaled residuals, then unscale for alpha distribution
    sigma_va <- sigma[c(if (num.lvc > 0L) num.RR + seq_len(num.lvc) else NULL,
                        if (num.lv  > 0L) num.RR + num.lvc + seq_len(num.lv) else NULL)]
    lv_res   <- getLV(object, type = "residual")  # sigma-scaled, std [lvc|lv]
    ## unscale so the alpha sweep below applies cleanly
    lv_sel   <- sweep(lv_res, 2L, sigma_va, `/`)
    col_ho   <- c(if (num.lvc > 0L) num.RR + seq_len(num.lvc) else NULL,
                  if (num.lv  > 0L) num.RR + num.lvc + seq_len(num.lv) else NULL)
    theta_sel <- loadings_ho[, col_ho, drop = FALSE]
    sigma_sel <- sigma[col_ho]
    col_idx   <- col_ho
  }

  d_sel <- ncol(lv_sel)
  if (d_sel < 1L) stop("No ordination dimensions selected for type = '", type, "'.")

  ## --- which.lvs ------------------------------------------------------------
  if (length(which.lvs) == 1L) which.lvs <- c(which.lvs, which.lvs)
  if (max(which.lvs) > d_sel)
    stop("which.lvs index exceeds the number of selected dimensions (", d_sel, ").")
  xl <- which.lvs[1L]; yl <- which.lvs[2L]

  ## --- Alpha-scale sites and species ----------------------------------------
  lv    <- sweep(lv_sel,    2L, sigma_sel^alpha,       `*`)
  theta <- sweep(theta_sel, 2L, sigma_sel^(1 - alpha), `*`)

  spp_names <- colnames(object$y)
  if (is.null(spp_names))
    spp_names <- rownames(object$params$theta) %||% paste0("sp", seq_len(p))
  if (!is.null(ind.spp)) ind.spp <- min(p, ind.spp) else ind.spp <- p

  if (length(spp.colors) == 1L) spp.colors <- rep(spp.colors, p)
  if (length(cex.spp)    == 1L) cex.spp    <- rep(cex.spp,    p)

  ## --- Visual norm balance --------------------------------------------------
  col_norms_lv <- pmax(sqrt(colSums(lv^2)),                                    1e-10)
  col_norms_th <- pmax(sqrt(colSums(theta[seq_len(ind.spp), , drop=FALSE]^2)), 1e-10)
  bothnorms    <- sqrt(col_norms_lv * col_norms_th)

  scale_lv <- sqrt(bothnorms) / col_norms_lv
  scale_th <- sqrt(bothnorms) / col_norms_th

  lv_sc <- t(t(lv) * scale_lv)
  th_sc <- t(t(theta) * scale_th)

  ## Transform multipliers for ellipses
  ell_lv <- (sigma_sel^alpha       * scale_lv)^2
  ell_th <- (sigma_sel^(1 - alpha) * scale_th)^2

  ## --- Ellipse data ---------------------------------------------------------
  need_ellips_sites   <- !isFALSE(predict.region) &&
                         (isTRUE(predict.region) || identical(predict.region, "sites"))
  need_ellips_species <- identical(predict.region, "species")
  show_spp            <- biplot || need_ellips_species

  if (need_ellips_sites || need_ellips_species) {
    if (isFALSE(object$sd))
      warning("No standard errors in model; prediction ellipses may be inaccurate.")
    pe          <- getPredictErr(object, CMSEP = CMSEP, cov = FALSE)
    Pvar_sites_full   <- pe$lvs^2
    Pvar_species_full <- pe$loadings^2

    .sel_pvar <- function(mat, cix) {
      nc_mat <- ncol(mat)
      out <- matrix(0, nrow(mat), length(cix))
      for (ki in seq_along(cix)) {
        k <- cix[ki]
        if (k <= nc_mat) out[, ki] <- mat[, k] else out[, ki] <- Inf
      }
      out
    }
    Pvar_sites   <- .sel_pvar(Pvar_sites_full,   col_idx)
    Pvar_species <- .sel_pvar(Pvar_species_full, col_idx)
  }

  ## --- Plot frame -----------------------------------------------------------
  gr_par_list <- list(...)
  plotfun <- if (("ylim" %in% names(gr_par_list)) || ("xlim" %in% names(gr_par_list)))
    plot else MASS::eqscplot

  if (is.null(main)) {
    ## Title reflects the ordination type:
    ##   num.lv.c > 0                           → Hierarchical ordination
    ##   num.RR > 0, both lv.X and TR present   → Double constrained ordination
    ##   num.RR > 0, only lv.X or only TR       → Constrained ordination
    has_lvc <- (object$num.lv.c %||% 0L) > 0L
    has_rr  <- (object$num.RR   %||% 0L) > 0L
    double  <- has_rr && !has_lvc &&
               !is.null(object$lv.X.design) && !is.null(object$TR)
    ord_label <- if (has_lvc)   "Hierarchical ordination"
                 else if (double) "Double constrained ordination"
                 else             "Constrained ordination"
    main <- if (d_sel == 1L) ord_label else
      paste0(ord_label, " (type = '", type, "')")
  }

  xlab_base <- switch(type,
    marginal  = "Constrained axis",
    residual  = "Residual LV",
    "LV"
  )
  site_labs <- rownames(object$y) %||% rownames(object$lvs) %||% as.character(seq_len(n))

  if (d_sel == 1L) {
    plotfun(seq_len(n), lv_sc[, 1L],
            ylab = paste0(xlab_base, " 1"), xlab = "Row index",
            main = main, type = "n", ...)
    if (symbols) {
      points(seq_len(n), lv_sc[, 1L], col = s.colors, cex = s.cex)
    } else {
      text(seq_len(n), lv_sc[, 1L], labels = site_labs, cex = s.cex, col = s.colors)
    }
    return(invisible(list(sites   = cbind(seq_len(n), lv_sc[, 1L]),
                          species = cbind(seq_len(p), th_sc[, 1L]))))
  }

  ## Axis limits: include species when biplot is active so loadings don't fall off.
  ## b_gamma arrows are always drawn (like b_z) so include their scale-contribution
  ## in the range even when biplot=FALSE.
  has_bgamma_arrows <- !is.null(object$params$LoadTRcoef) && type != "residual"
  range_pts <- lv_sc[, c(xl, yl), drop = FALSE]
  if (show_spp)
    range_pts <- rbind(range_pts, th_sc[seq_len(ind.spp), c(xl, yl), drop = FALSE])
  plotfun(range_pts, type = "n",
          xlab = paste(xlab_base, xl), ylab = paste(xlab_base, yl),
          main = main, ...)

  ## --- Site prediction ellipses ---------------------------------------------
  if (need_ellips_sites) {
    if (length(col.ellips) != n) col.ellips <- rep(col.ellips, length.out = n)
    rad <- sqrt(qchisq(level, df = 2L))
    for (i in seq_len(n)) {
      covm <- diag(ell_lv[c(xl, yl)] * Pvar_sites[i, c(xl, yl)], nrow = 2L)
      ellipse(lv_sc[i, c(xl, yl)], covM = covm, rad = rad,
              col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
    }
  }

  ## --- Site labels / symbols ------------------------------------------------
  if (jitter) {
    lv_sc[, xl] <- lv_sc[, xl] + runif(n, -jitter.amount, jitter.amount)
    lv_sc[, yl] <- lv_sc[, yl] + runif(n, -jitter.amount, jitter.amount)
  }
  if (symbols) {
    points(lv_sc[, c(xl, yl)], col = s.colors, cex = s.cex)
  } else {
    text(lv_sc[, c(xl, yl)], labels = site_labs, cex = s.cex,
         col = if (length(s.colors) == 1L) s.colors else s.colors[seq_len(n)])
  }

  ## --- Species scores -------------------------------------------------------
  if (show_spp) {
    if (need_ellips_species) {
      rad_sp <- sqrt(qchisq(level, df = 2L))
      for (j in seq_len(ind.spp)) {
        covm_sp <- diag(ell_th[c(xl, yl)] * Pvar_species[j, c(xl, yl)], nrow = 2L)
        ellipse(th_sc[j, c(xl, yl)], covM = covm_sp, rad = rad_sp,
                col = col.spp.ellips, lwd = lwd.ellips, lty = lty.spp.ellips)
      }
    }

    if (spp.arrows) {
      lv_rng   <- range(lv_sc[, c(xl, yl)])
      in_range <- apply(th_sc[seq_len(ind.spp), c(xl, yl), drop=FALSE], 1,
                        function(r) all(r >= lv_rng[1L] & r <= lv_rng[2L]))
      marg <- par("usr")
      origin_sp <- c(mean(marg[1:2]), mean(marg[3:4]))
      half_x <- diff(marg[1:2]) / 2; half_y <- diff(marg[3:4]) / 2
      for (j in seq_len(ind.spp)) {
        if (in_range[j]) {
          text(th_sc[j, c(xl, yl)], labels = spp_names[j],
               cex = cex.spp[j], col = spp.colors[j])
        } else {
          e  <- th_sc[j, c(xl, yl)]
          sc <- min(half_x, half_y) * arrow.spp.scale / sqrt(sum(e^2))
          arrows(origin_sp[1L], origin_sp[2L],
                 e[1L] * sc + origin_sp[1L], e[2L] * sc + origin_sp[2L],
                 col = spp.colors[j], lty = spp.arrows.lty, length = 0.1)
          text(e[1L] * sc * (1 + lab.dist) + origin_sp[1L],
               e[2L] * sc * (1 + lab.dist) + origin_sp[2L],
               labels = spp_names[j], cex = cex.spp[j], col = spp.colors[j])
        }
      }
    } else {
      text(th_sc[seq_len(ind.spp), c(xl, yl)],
           labels = spp_names[seq_len(ind.spp)],
           cex    = cex.spp[seq_len(ind.spp)],
           col    = spp.colors[seq_len(ind.spp)])
    }
  }

  ## --- Helper: draw arrows from origin -------------------------------------
  .draw_arrows <- function(coef_mat, col, scale, cex_lab, lab_dist,
                            lty = "solid", sig = NULL) {
    if (is.null(coef_mat) || nrow(coef_mat) == 0L) return(invisible(NULL))
    marg   <- par("usr")
    origin <- c(mean(marg[1:2]), mean(marg[3:4]))
    half_x <- diff(marg[1:2]) / 2
    half_y <- diff(marg[3:4]) / 2
    max_len <- max(sqrt(rowSums(coef_mat^2)), na.rm = TRUE)
    if (max_len < .Machine$double.eps) return(invisible(NULL))
    ends <- coef_mat / max_len * min(half_x, half_y) * scale
    nms  <- rownames(coef_mat)
    col_vec <- if (!is.null(sig))
      ifelse(sig, col, adjustcolor(col, alpha.f = 0.4))
    else
      rep(col, nrow(ends))
    for (k in seq_len(nrow(ends))) {
      ex <- ends[k, 1L]; ey <- ends[k, 2L]
      arrows(origin[1L], origin[2L],
             ex + origin[1L], ey + origin[2L],
             col = col_vec[k], lty = lty, length = 0.1)
      text(ex * (1 + lab_dist) + origin[1L],
           ey * (1 + lab_dist) + origin[2L],
           labels = if (!is.null(nms)) nms[k] else k,
           col = col_vec[k], cex = cex_lab)
    }
  }

  ## --- b_z arrows (canonical covariates) -----------------------------------
  if (!is.null(object$params$LvXcoef) && type != "residual") {
    ## Reconstruct unscaled b_z in HO [RR|lvc|lv] order from LvXcoef
    ## LvXcoef is Kz x d_c, standard [lvc|RR] order, sigma-scaled
    d_c  <- num.RR + num.lvc
    lvxc <- object$params$LvXcoef
    b_z_full <- matrix(0, nrow(lvxc), d, dimnames = list(rownames(lvxc), NULL))
    if (num.lvc > 0L)
      b_z_full[, num.RR + seq_len(num.lvc)] <-
        sweep(lvxc[, seq_len(num.lvc), drop = FALSE],
              2L, sigma[num.RR + seq_len(num.lvc)], `/`)
    if (num.RR > 0L)
      b_z_full[, seq_len(num.RR)] <-
        sweep(lvxc[, num.lvc + seq_len(num.RR), drop = FALSE],
              2L, sigma[seq_len(num.RR)], `/`)
    if (ncol(b_z_full) >= max(col_idx)) {
      bz_sel <- b_z_full[, col_idx, drop = FALSE]
      bz_sel <- sweep(bz_sel, 2L, sigma_sel^alpha, `*`)
      if (!is.null(object$lv.X)) {
        sds <- apply(as.matrix(object$lv.X), 2L, sd)
        sds[sds < .Machine$double.eps] <- 1
        bz_sel <- bz_sel / sds
      }
      rownames(bz_sel) <- if (!is.null(object$lv.X) &&
                               !is.null(colnames(object$lv.X)))
        colnames(object$lv.X)
      else
        rownames(b_z_full) %||% paste0("cov", seq_len(nrow(bz_sel)))
      bz_plot <- bz_sel[, c(xl, yl), drop = FALSE]

      sig_bz <- NULL
      if (arrow.ci && !isFALSE(object$sd) &&
          !is.null(object$sd) && !is.null(object$sd$b_z)) {
        se_bz <- object$sd$b_z
        if (ncol(se_bz) >= max(col_idx)) {
          bz_raw  <- b_z_full[, col_idx, drop = FALSE][, c(xl, yl), drop = FALSE]
          se_plot <- se_bz[,   col_idx, drop = FALSE][, c(xl, yl), drop = FALSE]
          sig_bz  <- apply(abs(bz_raw) > 1.96 * abs(se_plot), 1, any)
        }
      }
      .draw_arrows(bz_plot, col.arrow.bz, arrow.scale, cex.env, lab.dist,
                   lty = arrow.lty, sig = sig_bz)
    }
  }

  ## --- b_gamma arrows (traits) — shown by default, same as b_z arrows ------
  if (!is.null(object$params$LoadTRcoef) && type != "residual") {
    b_gamma_full <- object$params$LoadTRcoef   # Kt x d (HO order), same as LoadTRcoef
    if (ncol(b_gamma_full) >= max(col_idx)) {
      bg_sel <- b_gamma_full[, col_idx, drop = FALSE]
      bg_sel <- sweep(bg_sel, 2L, sigma_sel^(1 - alpha), `*`)
      if (!is.null(object$TR)) {
        sds <- apply(as.matrix(object$TR), 2L, sd)
        sds[sds < .Machine$double.eps] <- 1
        bg_sel <- bg_sel / sds
      }
      rownames(bg_sel) <- if (!is.null(object$TR) &&
                               !is.null(colnames(object$TR)))
        colnames(object$TR)
      else
        rownames(b_gamma_full) %||% paste0("trait", seq_len(nrow(bg_sel)))
      bg_plot <- bg_sel[, c(xl, yl), drop = FALSE]

      sig_bg <- NULL
      if (arrow.ci && !isFALSE(object$sd) &&
          !is.null(object$sd) && !is.null(object$sd$LoadTRcoef)) {
        se_bg <- object$sd$LoadTRcoef
        if (ncol(se_bg) >= max(col_idx)) {
          bg_raw  <- b_gamma_full[, col_idx, drop = FALSE][, c(xl, yl), drop = FALSE]
          se_plot <- se_bg[,        col_idx, drop = FALSE][, c(xl, yl), drop = FALSE]
          sig_bg  <- apply(abs(bg_raw) > 1.96 * abs(se_plot), 1, any)
        }
      }
      .draw_arrows(bg_plot, col.arrow.bgamma, arrow.scale, cex.env, lab.dist,
                   lty = arrow.lty, sig = sig_bg)
    }
  }

  invisible(list(sites   = lv_sc[, c(xl, yl)],
                 species = th_sc[, c(xl, yl)]))
}

## Infix NULL-coalesce: x %||% y  → x if not NULL, else y
`%||%` <- function(x, y) if (!is.null(x)) x else y
