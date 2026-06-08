#' Ordination plot for hierarchical ordination models
#'
#' Plots site scores and (optionally) species scores from a \code{gllvmHO}
#' model.  Uncertainty ellipses for sites use the VA posterior variances
#' \code{object$A}; ellipses for species use \code{sigma^2 * object$A_lv},
#' which is the variance of the biplot-space quantity \eqn{\Sigma \gamma_j}.
#'
#' @param object   A \code{gllvmHO} object from \code{\link{gllvm}} with
#'   \code{random.loadings = TRUE}.
#' @param biplot   Logical; add species scores to the plot?  Defaults to
#'   \code{FALSE}.
#' @param ind.spp  Number of species to show (default: all).
#' @param alpha    Biplot balance parameter in [0,1].  \code{alpha = 0.5}
#'   (default) balances site and species spread.
#' @param main     Plot title.
#' @param which.lvs Integer vector of length 2 selecting which latent variable
#'   axes to plot.  Default \code{c(1, 2)}.
#' @param predict.region  If \code{TRUE} or \code{"sites"}, draw VA-based
#'   uncertainty ellipses for sites.  If \code{"species"}, draw uncertainty
#'   ellipses for species loadings (implies \code{biplot = TRUE}).
#'   Default \code{FALSE}.
#' @param level    Coverage level for ellipses.  Default 0.95.
#' @param jitter   Logical; add jitter to site positions?
#' @param jitter.amount  Half-range of jitter.
#' @param s.colors   Colour(s) for site labels/points.
#' @param s.cex      Character expansion for site labels.
#' @param symbols    Logical; use points instead of labels for sites?
#' @param cex.spp    Character expansion for species labels.
#' @param spp.colors  Colour(s) for species labels.
#' @param lwd.ellips  Line width for site ellipses.
#' @param col.ellips  Colour(s) for site ellipses.
#' @param lty.ellips  Line type for site ellipses.
#' @param col.spp.ellips  Colour for species ellipses.
#' @param lty.spp.ellips  Line type for species ellipses.
#' @param rotate  Logical; apply post-hoc SVD rotation for interpretability?
#' @param ...  Further graphical parameters passed to \code{\link[MASS]{eqscplot}}.
#'
#' @return Invisibly, a list with \code{sites} (n x 2) and \code{species}
#'   (p x 2) plot coordinates.
#'
#' @method ordiplot gllvmHO
#' @export
ordiplot.gllvmHO <- function(object,
                              biplot        = FALSE,
                              ind.spp       = NULL,
                              alpha         = 0.5,
                              main          = NULL,
                              which.lvs     = c(1, 2),
                              predict.region = FALSE,
                              level         = 0.95,
                              jitter        = FALSE,
                              jitter.amount = 0.2,
                              s.colors      = 1,
                              s.cex         = 1.2,
                              symbols       = FALSE,
                              cex.spp       = 0.7,
                              spp.colors    = "blue",
                              lwd.ellips    = 0.5,
                              col.ellips    = 4,
                              lty.ellips    = 1,
                              col.spp.ellips = "blue",
                              lty.spp.ellips = 2,
                              rotate        = TRUE,
                              ...) {

  n <- nrow(object$y)
  p <- ncol(object$y)
  d <- object$num.lv
  if (d < 1) stop("No latent variables in model.")

  ## --- Sigma scaling --------------------------------------------------------
  ## params$theta stores unscaled loadings a_j (E[gamma_j]).
  ## The biplot quantity is Sigma * a_j, where sigma = object$params$sigma.lv.
  ## All species-side geometry must use this scaled version.
  sigma <- object$params$sigma.lv                       # length-d
  theta <- sweep(object$params$theta, 2, sigma, `*`)    # p x d  (Sigma * a_j)

  spp_names <- rownames(object$params$theta)
  if (is.null(spp_names))
    spp_names <- if (!is.null(colnames(object$y))) colnames(object$y) else
      paste0("sp", seq_len(p))
  rownames(theta) <- spp_names

  lv <- object$lvs   # n x d  (E[z_i])

  ## --- Axis selection -------------------------------------------------------
  if (!is.null(ind.spp)) ind.spp <- min(p, ind.spp) else ind.spp <- p
  if (length(which.lvs) == 1) which.lvs <- c(which.lvs, which.lvs)
  xl <- which.lvs[1]; yl <- which.lvs[2]
  if (xl > d || yl > d)
    stop("which.lvs contains an index exceeding the number of latent variables.")

  ## --- Post-hoc SVD rotation ------------------------------------------------
  ## Rotate in the full d-dimensional space so all axes are consistent.
  rot    <- if (rotate && d > 1) svd(lv)$v else diag(d)
  lv_rot <- lv    %*% rot    # n x d
  th_rot <- theta %*% rot    # p x d

  ## --- Biplot alpha-scaling -------------------------------------------------
  ## For each plotted axis k, scale so that lv_sc and th_sc have comparable
  ## spread.  Inner product lv_sc[i,k] * th_sc[j,k] approximates the
  ## ordination contribution z_ik * sigma_k * a_jk.
  col_norms_lv <- pmax(sqrt(colSums(lv_rot^2)), 1e-10)
  col_norms_th <- pmax(sqrt(colSums(th_rot[seq_len(ind.spp), , drop = FALSE]^2)), 1e-10)
  bothnorms    <- sqrt(col_norms_lv * col_norms_th)     # length-d

  lv_sc <- t(t(lv_rot) * bothnorms^alpha       / col_norms_lv)
  th_sc <- t(t(th_rot) * bothnorms^(1 - alpha) / col_norms_th)

  ## --- Plot frame -----------------------------------------------------------
  if (is.null(main)) main <- "Hierarchical ordination"
  MASS::eqscplot(lv_sc[, c(xl, yl)], type = "n",
                 xlab = paste("LV", xl), ylab = paste("LV", yl),
                 main = main, ...)

  ## --- Site VA ellipses -----------------------------------------------------
  ## Uncertainty for site i comes from the VA posterior: q(z_i) ~ N(mu_i, diag(A[i,])).
  ## In plot space the covariance is B %*% diag(A[i,]) %*% t(B), where B maps
  ## the d-dimensional LV space to the two plotted axes after scaling + rotation.
  if (!isFALSE(predict.region) &&
      (isTRUE(predict.region) || identical(predict.region, "sites"))) {
    rad    <- sqrt(qchisq(level, df = 2))
    sscale <- bothnorms^alpha / col_norms_lv             # per-axis scale
    B      <- (diag(sscale) %*% rot)[c(xl, yl), , drop = FALSE]   # 2 x d
    for (i in seq_len(n)) {
      covm <- B %*% diag(object$A[i, ], d) %*% t(B)
      ellipse(lv_sc[i, c(xl, yl)], covM = covm, rad = rad,
              col = col.ellips[min(i, length(col.ellips))],
              lwd = lwd.ellips, lty = lty.ellips)
    }
  }

  ## --- Site labels / symbols ------------------------------------------------
  if (jitter) {
    lv_sc[, xl] <- lv_sc[, xl] + runif(n, -jitter.amount, jitter.amount)
    lv_sc[, yl] <- lv_sc[, yl] + runif(n, -jitter.amount, jitter.amount)
  }
  site_labs <- if (!is.null(rownames(lv))) rownames(lv) else as.character(seq_len(n))
  if (symbols) {
    points(lv_sc[, c(xl, yl)], col = s.colors, cex = s.cex)
  } else {
    text(lv_sc[, c(xl, yl)], labels = site_labs, cex = s.cex,
         col = s.colors[if (length(s.colors) == 1) 1 else seq_len(n)])
  }

  ## --- Species scores -------------------------------------------------------
  ## Shown when biplot = TRUE or predict.region = "species".
  show_spp <- biplot || identical(predict.region, "species")
  if (show_spp) {

    ## Species VA ellipses.
    ## q(gamma_j) ~ N(a_j, diag(A_lv[j,])).
    ## In biplot space the relevant quantity is Sigma * gamma_j; its variance
    ## is diag(sigma^2 * A_lv[j,]).  The 2-d plot covariance is
    ##   Bsp %*% diag(sigma^2 * A_lv[j,]) %*% t(Bsp).
    if (identical(predict.region, "species")) {
      rad_sp <- sqrt(qchisq(level, df = 2))
      tscale <- bothnorms^(1 - alpha) / col_norms_th    # per-axis scale
      Bsp    <- (diag(tscale) %*% rot)[c(xl, yl), , drop = FALSE]   # 2 x d
      for (j in seq_len(ind.spp)) {
        covm_sp <- Bsp %*% diag(sigma^2 * object$A_lv[j, ], d) %*% t(Bsp)
        ellipse(th_sc[j, c(xl, yl)], covM = covm_sp, rad = rad_sp,
                col = col.spp.ellips, lwd = lwd.ellips, lty = lty.spp.ellips)
      }
    }

    text(th_sc[seq_len(ind.spp), c(xl, yl)],
         labels = spp_names[seq_len(ind.spp)],
         cex = cex.spp, col = spp.colors)
  }

  invisible(list(sites   = lv_sc[, c(xl, yl)],
                 species = th_sc[, c(xl, yl)]))
}
