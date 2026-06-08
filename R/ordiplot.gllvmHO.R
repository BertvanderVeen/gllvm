#' Ordination plot for hierarchical ordination models
#'
#' Plots site scores and (optionally) species scores from a \code{gllvmHO}
#' model.  \code{alpha} distributes the singular values \eqn{\Sigma} between
#' sites (\eqn{\sigma^\alpha}) and species (\eqn{\sigma^{1-\alpha}}), so the
#' inner product recovers \eqn{z_i^\top \Sigma a_j}.  No post-hoc rotation is
#' applied: the ordering constraint on \eqn{\sigma} already identifies the axes.
#'
#' @param object   A \code{gllvmHO} object.
#' @param biplot   Logical; add species scores?  Default \code{FALSE}.
#' @param ind.spp  Number of species to show (default: all).
#' @param alpha    In [0,1]. Distributes \eqn{\Sigma}: sites carry
#'   \eqn{\sigma^\alpha}, species carry \eqn{\sigma^{1-\alpha}}.
#'   \code{alpha = 0.5} (default) gives a symmetric biplot.
#' @param main     Plot title.
#' @param which.lvs Integer vector of length 2 selecting axes to plot.
#' @param predict.region  \code{TRUE}/\code{"sites"}: VA ellipses for sites.
#'   \code{"species"}: VA ellipses for species (implies biplot).
#' @param level    Ellipse coverage level.  Default 0.95.
#' @param jitter   Logical; jitter site positions?
#' @param jitter.amount  Half-range of jitter.
#' @param s.colors  Colour(s) for site labels/points.
#' @param s.cex     Character expansion for site labels.
#' @param symbols   Use points instead of text for sites?
#' @param cex.spp   Character expansion for species labels.
#' @param spp.colors  Colour for species labels.
#' @param lwd.ellips  Line width for site ellipses.
#' @param col.ellips  Colour(s) for site ellipses.
#' @param lty.ellips  Line type for site ellipses.
#' @param col.spp.ellips  Colour for species ellipses.
#' @param lty.spp.ellips  Line type for species ellipses.
#' @param ...  Further graphical parameters passed to \code{\link[MASS]{eqscplot}}.
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
                              predict.region = FALSE,
                              level          = 0.95,
                              jitter         = FALSE,
                              jitter.amount  = 0.2,
                              s.colors       = 1,
                              s.cex          = 1.2,
                              symbols        = FALSE,
                              cex.spp        = 0.7,
                              spp.colors     = "blue",
                              lwd.ellips     = 0.5,
                              col.ellips     = 4,
                              lty.ellips     = 1,
                              col.spp.ellips = "blue",
                              lty.spp.ellips = 2,
                              ...) {

  n <- nrow(object$y)
  p <- ncol(object$y)
  d <- object$num.lv
  if (d < 1) stop("No latent variables in model.")

  sigma <- object$params$sigma.lv   # length-d, ordered descending

  ## Distribute Sigma via alpha:
  ##   sites    <- z_i  * sigma^alpha          (n x d)
  ##   species  <- a_j  * sigma^(1-alpha)      (p x d)
  ## Inner product: (sigma^alpha * z_ik)(sigma^(1-alpha) * a_jk) = sigma_k * z_ik * a_jk  ✓
  ## No rotation: the ordering constraint on sigma identifies the axes.
  lv    <- sweep(object$lvs,          2L, sigma^alpha,       `*`)
  theta <- sweep(object$params$theta, 2L, sigma^(1 - alpha), `*`)

  spp_names <- rownames(object$params$theta)
  if (is.null(spp_names))
    spp_names <- if (!is.null(colnames(object$y))) colnames(object$y) else
      paste0("sp", seq_len(p))

  if (!is.null(ind.spp)) ind.spp <- min(p, ind.spp) else ind.spp <- p
  if (length(which.lvs) == 1L) which.lvs <- c(which.lvs, which.lvs)
  xl <- which.lvs[1L]; yl <- which.lvs[2L]
  if (xl > d || yl > d)
    stop("which.lvs index exceeds the number of latent variables.")

  ## Visual norm balance: scale each axis so sites and species have equal
  ## visual spread in the plot.  Sigma is already distributed; this is purely
  ## cosmetic and uses a fixed 0.5/0.5 split.
  col_norms_lv <- pmax(sqrt(colSums(lv^2)),                              1e-10)
  col_norms_th <- pmax(sqrt(colSums(theta[seq_len(ind.spp),,drop=FALSE]^2)), 1e-10)
  bothnorms    <- sqrt(col_norms_lv * col_norms_th)   # geometric mean, length-d

  scale_lv <- sqrt(bothnorms) / col_norms_lv   # visual scale per axis, sites
  scale_th <- sqrt(bothnorms) / col_norms_th   # visual scale per axis, species

  lv_sc <- t(t(lv) * scale_lv)
  th_sc <- t(t(theta) * scale_th)

  ## For ellipses: the complete linear transform from the VA mean to the plot
  ## coordinate is, per axis k:
  ##   z_ik   -> lv_sc[i,k]  via factor  sigma[k]^alpha     * scale_lv[k]
  ##   a_jk   -> th_sc[j,k]  via factor  sigma[k]^(1-alpha) * scale_th[k]
  ## VA variances A[i,k] and A_lv[j,k] are in the original z / a space, so
  ## the covariance in plot space is that factor squared times the VA variance.
  ## No rotation => covariance matrices are diagonal (axes are independent).
  ell_lv <- (sigma^alpha       * scale_lv)^2   # variance multipliers, sites
  ell_th <- (sigma^(1 - alpha) * scale_th)^2   # variance multipliers, species

  ## --- Plot frame -----------------------------------------------------------
  if (is.null(main)) main <- "Hierarchical ordination"
  MASS::eqscplot(lv_sc[, c(xl, yl)], type = "n",
                 xlab = paste("LV", xl), ylab = paste("LV", yl),
                 main = main, ...)

  ## --- Site VA ellipses -----------------------------------------------------
  if (!isFALSE(predict.region) &&
      (isTRUE(predict.region) || identical(predict.region, "sites"))) {
    rad <- sqrt(qchisq(level, df = 2))
    for (i in seq_len(n)) {
      covm <- diag(ell_lv[c(xl, yl)] * object$A[i, c(xl, yl)], nrow = 2)
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
  site_labs <- if (!is.null(rownames(object$lvs))) rownames(object$lvs) else
    as.character(seq_len(n))
  if (symbols) {
    points(lv_sc[, c(xl, yl)], col = s.colors, cex = s.cex)
  } else {
    text(lv_sc[, c(xl, yl)], labels = site_labs, cex = s.cex,
         col = s.colors[if (length(s.colors) == 1L) 1L else seq_len(n)])
  }

  ## --- Species scores -------------------------------------------------------
  show_spp <- biplot || identical(predict.region, "species")
  if (show_spp) {

    ## Species VA ellipses.
    ## Var(th_sc[j,k]) = ell_th[k] * A_lv[j,k]
    if (identical(predict.region, "species")) {
      rad_sp <- sqrt(qchisq(level, df = 2))
      for (j in seq_len(ind.spp)) {
        covm_sp <- diag(ell_th[c(xl, yl)] * object$A_lv[j, c(xl, yl)], nrow = 2)
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
