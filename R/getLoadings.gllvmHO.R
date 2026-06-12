#' @title Extract species loadings from a hierarchical ordination model
#' @description
#' Returns the p x d matrix of unscaled species scores (gamma) in HO
#' \[RR | lvc | lv\] column order, reconstructed from \code{params$theta}
#' (which is stored in standard \[lvc | RR | lv\] order).
#'
#' @param object A \code{gllvmHO} object.
#' @param ... Not used.
#'
#' @return A \eqn{p \times d} matrix with unscaled species scores in HO order.
#'
#' @method getLoadings gllvmHO
#' @export
#' @export getLoadings.gllvmHO
getLoadings.gllvmHO <- function(object, ...) {
  num.RR  <- as.integer(object$num.RR  %||% 0L)
  num.lvc <- as.integer(object$num.lv.c %||% 0L)
  num.lv  <- as.integer(object$num.lv  %||% 0L)

  theta_std <- object$params$theta  # p x d, std [lvc|RR|lv] order

  ## Reorder to HO [RR|lvc|lv]
  if (num.RR > 0L && num.lvc > 0L) {
    ho_idx <- c(num.lvc + seq_len(num.RR),
                seq_len(num.lvc),
                if (num.lv > 0L) num.lvc + num.RR + seq_len(num.lv) else integer(0))
    theta_std[, ho_idx, drop = FALSE]
  } else {
    theta_std
  }
}
