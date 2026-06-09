##############################################################################
## Hierarchical Ordination VA (HO_VA) fitting function
## Uses gllvm_HO TMB template built as a side-car DLL during R CMD INSTALL.
## Model: eta_ij = beta0_j + x_i^T beta_j + z_i^T Sigma gamma_j
##   z_i ~ N(0,I),  gamma_j ~ N(0,I)  (both random)
##############################################################################

#' @keywords internal
.ensure_ho_dll <- function() {
  if (is.element("gllvm_HO", names(getLoadedDLLs()))) return(invisible(NULL))
  ext <- .Platform$dynlib.ext
  so <- system.file("libs", paste0("gllvm_HO", ext), package = "gllvm")
  if (!nzchar(so))
    so <- system.file("src",  paste0("gllvm_HO", ext), package = "gllvm")
  if (!nzchar(so))
    stop("gllvm_HO shared library not found; please reinstall the gllvm package.")
  dyn.load(so)
  invisible(NULL)
}

# Family name → integer code (must match C++ enum in gllvm_HO.cpp)
.fam_code <- c(
  poisson      = 0L, negative.binomial = 1L, binomial    = 2L,
  gaussian     = 3L, gamma             = 4L, tweedie     = 5L,
  ZIP          = 6L, ordinal           = 7L, exponential = 8L,
  beta         = 9L, betaH             = 10L, ZINB        = 11L,
  ordered.beta = 12L, ZIB              = 13L, ZNIB        = 14L,
  beta.binomial = 15L
)

##############################################################################
#' Fit a Hierarchical Ordination model with Variational Approximation
#'
#' @param y  n x p response matrix.
#' @param X  Optional n x Kx covariate matrix (without intercept). An intercept
#'   is always included.
#' @param family  Character family name (length 1 or p). Currently tested:
#'   \code{"poisson"}, \code{"gaussian"}, \code{"negative.binomial"},
#'   \code{"binomial"}, \code{"gamma"}, \code{"tweedie"}.
#' @param num.lv  Ordination dimension d (default 2).
#' @param offset  Optional n x p offset matrix.
#' @param Ntrials  For binomial: n x p integer matrix of trial counts.
#' @param row.eff  Logical. If \code{TRUE}, include i.i.d. random row effects
#'   (variance estimated from the data).
#' @param maxit  Maximum iterations for nlminb.
#' @param reltol  Relative tolerance for convergence.
#' @param diag.iter  Number of diagonal-A inner iterations (0 = skip).
#' @param start.params  Optional list with \code{$lvs} (n x d) and
#'   \code{$loadings} (p x d) to override SVD-based starts.
#' @param trace  Print optimisation trace.
#'
#' @return An object of class \code{c("gllvmHO", "gllvm")}.
#'
#' @keywords internal
gllvm.HO.TMB <- function(
    y, X = NULL, lv.X = NULL, TR = NULL,
    family = "poisson", num.lv = 2, num.lv.c = 0L, num.RR = 0L,
    offset = NULL,
    Ntrials = matrix(1L), row.eff = FALSE, studyDesign = NULL,
    Lambda.struc = "unstructured", zeta.struc = "common",
    n.init = 1L, n.init.max = 10L,
    maxit = 2000, reltol = 1e-8,
    diag.iter = 1, start.params = NULL, trace = FALSE, call. = NULL,
    starting.val = "res", jitter.var = 0,
    randomB = "LV", randomT = "LV",
    csb_z = matrix(0L, 0L, 2L), csb_gamma = matrix(0L, 0L, 2L)
) {
  .ensure_ho_dll()

  ## ---- validate randomB / randomT -----------------------------------------
  if (!identical(randomB, "LV"))
    stop("randomB = '", randomB, "' is not yet supported for HO models. ",
         "Use randomB = 'LV'. Fixed canonical coefficients (randomB = FALSE) ",
         "are planned but not yet implemented.")
  if (!identical(randomT, "LV"))
    stop("randomT = '", randomT, "' is not yet supported for HO models. ",
         "Use randomT = 'LV'.")

  ## ---- dimensions & family setup ------------------------------------------
  n      <- nrow(y)
  p      <- ncol(y)
  d      <- as.integer(num.lv)      # total ordination dimension (num.RR + num.lv.c + num_lv_unc)
  num.RR  <- as.integer(num.RR)
  num.lv.c <- as.integer(num.lv.c)
  if (length(family) == 1L) family <- rep(family, p)
  fam_int <- .fam_code[family]
  if (any(is.na(fam_int)))
    stop("Unknown family: ", paste(family[is.na(fam_int)], collapse = ", "))

  ## extra link flags (0 = canonical/logit for binary, etc.)
  extra <- rep(0, p)

  ## ---- LV covariate / trait dimensions ------------------------------------
  ## These are needed early to compute VA array dimensions d_va_z / d_va_a.
  Kz <- if (!is.null(lv.X)) ncol(as.matrix(lv.X)) else 0L
  Kt <- if (!is.null(TR))   ncol(as.matrix(TR))   else 0L

  ## RR dims contribute to VA only when the corresponding predictor is absent
  rr_va_z <- if (Kz > 0L) 0L else num.RR
  rr_va_a <- if (Kt > 0L) 0L else num.RR
  d_va_z  <- rr_va_z + num.lv.c + (d - num.RR - num.lv.c)   # = rr_va_z + (d - num.RR)
  d_va_a  <- rr_va_a + num.lv.c + (d - num.RR - num.lv.c)   # = rr_va_a + (d - num.RR)

  ## Active constrained columns of b_z / b_gamma (lv dims never use covariates)
  d_active <- num.RR + num.lv.c
  d_c <- if (Kz == 0L) 0L else min(Kz, d_active)
  d_t <- if (Kt == 0L) 0L else min(Kt, d_active)

  ## ---- fixed-effect design matrix -----------------------------------------
  if (is.null(X)) {
    Xmat <- matrix(1, n, 1)
  } else {
    Xmat <- cbind(1, as.matrix(X))   # ensure matrix; prepend intercept
  }
  Kx <- ncol(Xmat)

  ## ---- offset, Ntrials -----------------------------------------------------
  if (is.null(offset)) offset_mat <- matrix(0, n, p) else offset_mat <- offset
  if (length(Ntrials) == 1L) Ntrials <- matrix(as.integer(Ntrials), n, p)

  ## ---- Starting values ------------------------------------------------
  link_y <- .ho_link_residuals(y, fam_int, Xmat)

  starting.val <- match.arg(starting.val, c("res", "zero", "random"))

  if (starting.val == "res") {
    sv <- .ho_svd_starts(link_y$residuals, d_va_z, d_va_a, link_y$beta0,
                         d_total = d)
  } else if (starting.val == "random") {
    sv <- .ho_svd_starts(link_y$residuals, d_va_z, d_va_a, link_y$beta0,
                         d_total = d)
    if (d_va_z > 0L)
      sv$u  <- sv$u  + matrix(rnorm(n * d_va_z, sd = sqrt(jitter.var + 0.5)),
                               n, d_va_z)
    if (d_va_a > 0L)
      sv$a_sp <- sv$a_sp + matrix(rnorm(p * d_va_a, sd = sqrt(jitter.var + 0.5)),
                                   p, d_va_a)
  } else {
    ## "zero"
    sv <- list(
      u      = matrix(0, n, d_va_z),
      a_sp   = matrix(0, p, d_va_a),
      beta0  = link_y$beta0,
      sigmaLV = rep(0, d)
    )
  }

  ## Jitter on top of "res" starting values
  if (starting.val == "res" && jitter.var > 0) {
    if (d_va_z > 0L)
      sv$u <- sv$u + matrix(rnorm(n * d_va_z, sd = sqrt(jitter.var)), n, d_va_z)
    if (d_va_a > 0L)
      sv$a_sp <- sv$a_sp + matrix(rnorm(p * d_va_a, sd = sqrt(jitter.var)),
                                    p, d_va_a)
  }

  if (!is.null(start.params)) {
    ## start.params$lvs is n×d (full); take VA-z columns
    if (!is.null(start.params$lvs)) {
      sp_lvs <- as.matrix(start.params$lvs)
      ## VA-z columns from the full d-column start:
      ## if Kz>0: skip first num.RR cols (det), keep rest → d_va_z cols
      ## if Kz==0: keep all d cols → but d_va_z < d only when num.RR>0 && Kz>0
      if (ncol(sp_lvs) == d_va_z) sv$u <- sp_lvs
      else if (ncol(sp_lvs) >= d_va_z)
        sv$u <- sp_lvs[, seq_len(d_va_z), drop = FALSE]
    }
    if (!is.null(start.params$loadings)) {
      sp_ld <- as.matrix(start.params$loadings)
      if (ncol(sp_ld) == d_va_a) sv$a_sp <- sp_ld
      else if (ncol(sp_ld) >= d_va_a)
        sv$a_sp <- sp_ld[, seq_len(d_va_a), drop = FALSE]
    }
  }

  ## ---- dispersion starting values -----------------------------------------
  lg_phi    <- rep(0, p)     # log(phi) for NB/Gaussian/Gamma
  lg_phiZINB<- rep(0, p)
  ePower    <- 0             # Tweedie power (logit-scale)

  ## ---- ordinal cutpoint (zeta) starting values ----------------------------
  ordinal_cols <- which(fam_int == 7L)
  has_ordinal  <- length(ordinal_cols) > 0L
  zetastruc_int <- if (zeta.struc == "species") 1L else 0L

  if (has_ordinal && min(y[, ordinal_cols, drop = FALSE], na.rm = TRUE) == 0L)
    y[, ordinal_cols] <- y[, ordinal_cols] + 1L

  if (!has_ordinal) {
    zeta <- numeric(0)
  } else if (zetastruc_int == 0L) {
    ## Common cutpoints: K = max(y over ordinal cols) - 1 cutpoints,
    ## first anchored at 0; K-1 raw (log-diff) parameters
    ymax_ord <- max(y[, ordinal_cols, drop = FALSE], na.rm = TRUE)
    K_cuts   <- as.integer(ymax_ord) - 1L      # number of cutpoints
    n_raw    <- K_cuts - 1L                    # raw parameters needed
    zeta     <- if (n_raw > 0L) rep(0, n_raw) else numeric(0)
  } else {
    ## Species-specific cutpoints: for each ordinal species j, Kj-1 raw params
    zeta <- unlist(lapply(ordinal_cols, function(j) {
      ymaxj <- max(y[, j], na.rm = TRUE)
      Kj    <- as.integer(ymaxj) - 1L
      n_rawj <- Kj - 1L
      if (n_rawj > 0L) rep(0, n_rawj) else numeric(0)
    }))
  }

  ## ---- row-effect setup ---------------------------------------------------
  ## row.eff = FALSE/"none"    : no row effects
  ## row.eff = "fixed"         : one fixed intercept per site (xr = I_n)
  ## row.eff = "random"/TRUE   : i.i.d. random row effects (xrr = I_n)
  ## row.eff = ~group          : fixed effect per group level
  ## row.eff = ~(1|group)      : grouped random intercepts (one VA mean per group)
  ## Mixed formula (fixed + random bars) also handled.
  ## studyDesign: data.frame supplying grouping variables for formula row effects.

  ## Canonicalise non-formula strings
  if (!inherits(row.eff, "formula")) {
    if (isFALSE(row.eff) || identical(row.eff, "none")) {
      row.eff <- FALSE
    } else if (identical(row.eff, "fixed")) {
      if (is.null(studyDesign))
        studyDesign <- data.frame(sample = factor(seq_len(n)))
      row.eff <- ~sample
    } else if (isTRUE(row.eff) || identical(row.eff, "random")) {
      if (is.null(studyDesign))
        studyDesign <- data.frame(sample = factor(seq_len(n)))
      row.eff <- ~(1|sample)
    } else {
      warning("Unrecognised row.eff value; row effects ignored.")
      row.eff <- FALSE
    }
  }

  xr          <- matrix(0, 0, 0)   # fixed row-effect design (n x Kr)
  r0f         <- matrix(0, 0, 1)   # fixed row-effect coefficients (Kr x 1)
  xrr         <- matrix(0, n, 0)   # random row-effect design (n x G)
  r0r         <- matrix(0, 0, 1)   # random row-effect VA means (G x 1)
  lg_Ar       <- numeric(0)        # random row-effect log-Chol diagonals (length G)
  log_sigma_r <- numeric(0)
  random_flag <- 0L

  if (inherits(row.eff, "formula")) {
    if (is.null(studyDesign))
      stop("studyDesign must be provided for formula-based row effects.")

    ## --- Fixed part (terms without | bars) ---
    fixed_form <- nobars1_(row.eff)
    if (!is.null(fixed_form) && length(all.vars(fixed_form)) > 0) {
      xr  <- model.matrix(fixed_form, data = studyDesign)
      r0f <- matrix(0, ncol(xr), 1)
    }

    ## --- Random part (terms with | bars) ---
    if (anyBars(row.eff)) {
      bar.f <- findbars1(row.eff)
      ## Build one indicator matrix per bar term, column-bind them
      xrr_list <- lapply(bar.f, function(b) {
        grp_var <- deparse(b[[3]])
        if (!grp_var %in% colnames(studyDesign))
          stop("Grouping variable '", grp_var, "' not found in studyDesign.")
        grp <- factor(studyDesign[[grp_var]])
        ## Dummy indicator: n x G (drop intercept = TRUE gives G columns)
        model.matrix(~ grp - 1)
      })
      xrr <- do.call(cbind, xrr_list)
      G   <- ncol(xrr)
      r0r <- matrix(0, G, 1)
      lg_Ar       <- rep(log(sqrt(0.1)), G)
      log_sigma_r <- log(0.3)
      random_flag <- 1L
    }
  }

  ## ---- VA covariance structure --------------------------------------------
  va_struct_int <- if (identical(Lambda.struc, "unstructured")) 1L else 0L
  tri_z <- d_va_z * (d_va_z + 1L) / 2L   # lower-triangle elements for site VA
  tri_a <- d_va_a * (d_va_a + 1L) / 2L   # lower-triangle elements for species VA

  ## Au(idx*n+i) for diagonal:     idx = k (0..d_va_z-1), value = log sqrt(A_i(k,k))
  ## Au(idx*n+i) for unstructured: idx = r*(r+1)/2+c (0..tri_z-1)
  if (va_struct_int == 0L) {
    Au_init    <- rep(log(sqrt(0.1)), d_va_z * n)
    Au_sp_init <- rep(log(sqrt(0.1)), d_va_a * p)
  } else {
    au_z_unit <- numeric(tri_z)
    for (r in seq_len(d_va_z)) au_z_unit[r*(r+1L)/2L] <- log(sqrt(0.1))
    au_a_unit <- numeric(tri_a)
    for (r in seq_len(d_va_a)) au_a_unit[r*(r+1L)/2L] <- log(sqrt(0.1))
    Au_init    <- rep(au_z_unit, each = n)   # length tri_z*n
    Au_sp_init <- rep(au_a_unit, each = p)   # length tri_a*p
  }

  ## ---- data list ----------------------------------------------------------
  data.list <- list(
    y         = as.matrix(y),
    x         = Xmat,
    xr        = if (nrow(xr) > 0) xr else matrix(0, 0, 0),
    xrr       = xrr,
    offset    = offset_mat,
    Ntrials   = Ntrials,
    family    = fam_int,
    extra     = extra,
    num_lv    = as.integer(d),
    num_RR    = num.RR,
    num_lvc   = num.lv.c,
    method    = 0L,              # VA
    zetastruc = zetastruc_int,
    p_betaH   = 0L,
    random    = random_flag,
    va_struct = va_struct_int,
    ## Canonical covariates and traits
    lv_X_env  = if (!is.null(lv.X)) as.matrix(lv.X) else matrix(0, n, 0),
    TR        = if (!is.null(TR))   as.matrix(TR)    else matrix(0, p, 0),
    ## randomB / randomT and correlation structures
    randomB   = 1L,                  # "LV" is the only supported option
    randomT   = 1L,
    csb_z     = csb_z,
    csb_gamma = csb_gamma
  )

  ## Ensure csb matrices are integer
  csb_z     <- matrix(as.integer(csb_z),     ncol = 2L)
  csb_gamma <- matrix(as.integer(csb_gamma), ncol = 2L)

  ## ---- parameter list -----------------------------------------------------
  param.list <- list(
    b         = rbind(sv$beta0, matrix(0, Kx - 1, p)),  # Kx x p
    u         = sv$u,          # n x d_va_z
    Au        = Au_init,       # n*d_va_z (diag) or n*tri_z (unstructured)
    a_lv_sp   = sv$a_sp,       # p x d_va_a
    Au_sp     = Au_sp_init,    # p*d_va_a (diag) or p*tri_a (unstructured)
    sigmaLV   = sv$sigmaLV,    # d
    lg_phi    = lg_phi,
    lg_phiZINB= lg_phiZINB,
    zeta      = zeta,
    ePower    = ePower,
    r0f       = r0f,
    r0r       = r0r,
    lg_Ar     = lg_Ar,
    log_sigma = log_sigma_r,
    b_z       = {
      ## When all dims are deterministic (d_va_z==0) and both lv.X and TR are present,
      ## b_z=0 & b_gamma=0 is a saddle point — initialise from SVD of residuals.
      bz0 <- matrix(0, Kz, d)
      if (Kz > 0L && Kt > 0L && d_va_z == 0L && d_va_a == 0L && d_c > 0L) {
        svd_r <- tryCatch(svd(link_y$residuals, nu = d_c, nv = d_c),
                          error = function(e) NULL)
        if (!is.null(svd_r) && length(svd_r$d) >= d_c) {
          lv_Xm <- as.matrix(lv.X)
          for (k in seq_len(d_c)) {
            u_k <- svd_r$u[, k] * sqrt(svd_r$d[k])
            bz0[, k] <- solve(crossprod(lv_Xm) + diag(Kz), t(lv_Xm) %*% u_k)
          }
        }
      }
      bz0
    },
    b_gamma   = {
      bg0 <- matrix(0, Kt, d)
      if (Kz > 0L && Kt > 0L && d_va_z == 0L && d_va_a == 0L && d_t > 0L) {
        svd_r <- tryCatch(svd(link_y$residuals, nu = d_t, nv = d_t),
                          error = function(e) NULL)
        if (!is.null(svd_r) && length(svd_r$d) >= d_t) {
          TRm <- as.matrix(TR)
          for (k in seq_len(d_t)) {
            v_k <- svd_r$v[, k] * sqrt(svd_r$d[k])
            bg0[, k] <- solve(crossprod(TRm) + diag(Kt), t(TRm) %*% v_k)
          }
        }
      }
      bg0
    },
    ## randomB VA parameters (mapped out when Kz=0 or randomB=FALSE)
    Ab_z          = {
      n_pairs_z <- nrow(csb_z)
      rep(log(sqrt(0.1)), max(Kz * d_c + n_pairs_z * d_c, 1L))
    },
    log_sigma_bz  = {
      n_pairs_z <- nrow(csb_z)
      ## All entries initialised at 0 (sigma_bz = 1); RR entries are mapped out below
      rep(0, max(d_c + n_pairs_z, 1L))
    },
    ## randomT VA parameters
    Ab_gamma       = {
      n_pairs_t <- nrow(csb_gamma)
      rep(log(sqrt(0.1)), max(Kt * d_t + n_pairs_t * d_t, 1L))
    },
    log_sigma_bgamma = {
      n_pairs_t <- nrow(csb_gamma)
      ## All entries initialised at 0; RR entries are mapped out below
      rep(0, max(d_t + n_pairs_t, 1L))
    }
  )

  ## ---- map: fix params that are unused ------------------------------------
  map.list <- list()
  if (Kz == 0L) {
    map.list$b_z <- factor(rep(NA, length(param.list$b_z)))
  } else if (d_c < d) {
    ## Padding: map out columns d_c+1..d of b_z (column-major: rows vary fastest)
    bz_map <- seq_len(Kz * d)
    bz_map[seq(d_c * Kz + 1L, d * Kz)] <- NA_integer_
    map.list$b_z <- factor(bz_map)
  }
  if (Kt == 0L) {
    map.list$b_gamma <- factor(rep(NA, length(param.list$b_gamma)))
  } else if (d_t < d) {
    bg_map <- seq_len(Kt * d)
    bg_map[seq(d_t * Kt + 1L, d * Kt)] <- NA_integer_
    map.list$b_gamma <- factor(bg_map)
  }
  ## Map out random-slope VA params when unused
  if (Kz == 0L || d_c == 0L) {
    map.list$Ab_z         <- factor(rep(NA, length(param.list$Ab_z)))
    map.list$log_sigma_bz <- factor(rep(NA, length(param.list$log_sigma_bz)))
  } else {
    ## For RR dims, sigma_bz is confounded with Sigma; fix to 1 via map
    n_rr_bz <- min(num.RR, d_c)
    if (n_rr_bz > 0L) {
      lsigbz_map <- seq_len(length(param.list$log_sigma_bz))
      lsigbz_map[seq_len(n_rr_bz)] <- NA_integer_
      map.list$log_sigma_bz <- factor(lsigbz_map)
    }
  }
  if (Kt == 0L || d_t == 0L) {
    map.list$Ab_gamma         <- factor(rep(NA, length(param.list$Ab_gamma)))
    map.list$log_sigma_bgamma <- factor(rep(NA, length(param.list$log_sigma_bgamma)))
  } else {
    n_rr_bt <- min(num.RR, d_t)
    if (n_rr_bt > 0L) {
      lsigbt_map <- seq_len(length(param.list$log_sigma_bgamma))
      lsigbt_map[seq_len(n_rr_bt)] <- NA_integer_
      map.list$log_sigma_bgamma <- factor(lsigbt_map)
    }
  }
  if (!has_ordinal)
    map.list$zeta   <- factor(rep(NA, length(zeta)))
  map.list$lg_phiZINB <- factor(rep(NA, length(lg_phiZINB)))
  map.list$ePower     <- factor(NA)
  # Fix dispersion for Poisson (no phi needed)
  poisson_cols <- which(fam_int == 0L)
  if (length(poisson_cols) > 0) {
    phi_map <- 1:p
    phi_map[poisson_cols] <- NA
    map.list$lg_phi <- factor(phi_map)
  }
  # Fix row-effect params that are unused
  has_fixed_re  <- nrow(xr) > 0 && ncol(xr) > 0
  has_random_re <- random_flag == 1L
  if (!has_fixed_re)  map.list$r0f       <- factor(rep(NA, length(r0f)))
  if (!has_random_re) {
    map.list$r0r       <- factor(rep(NA, length(r0r)))
    map.list$lg_Ar     <- factor(rep(NA, length(lg_Ar)))
    map.list$log_sigma <- factor(rep(NA, length(log_sigma_r)))
  }

  ## ---- MakeADFun ----------------------------------------------------------
  obj <- TMB::MakeADFun(
    data       = data.list,
    parameters = param.list,
    map        = map.list,
    DLL        = "gllvm_HO",
    silent     = !trace
  )

  ## ---- Diagonal-A inner iterations ----------------------------------------
  ## Rebuild full param list from obj$env$last.par.best (includes fixed params),
  ## then fix everything except Au/Au_sp using correct full-length factors.
  for (iter in seq_len(max(diag.iter, 0L))) {
    full_par <- obj$env$last.par.best

    # Reconstruct named param list from the full (all-param) vector
    cur_pl <- param.list
    fp_idx <- 1L
    for (nm in names(param.list)) {
      sz <- prod(dim(as.array(param.list[[nm]])))
      if (sz == 0L) next
      chunk <- full_par[fp_idx:(fp_idx + sz - 1L)]
      cur_pl[[nm]] <- if (is.matrix(param.list[[nm]]))
        matrix(chunk, nrow = nrow(param.list[[nm]])) else chunk
      fp_idx <- fp_idx + sz
    }

    # Fix everything except Au, Au_sp, Ab_z, Ab_gamma
    map_diag <- map.list
    for (nm in setdiff(names(param.list), c("Au", "Au_sp", "Ab_z", "Ab_gamma"))) {
      sz <- prod(dim(as.array(param.list[[nm]])))
      if (sz == 0L) next
      map_diag[[nm]] <- factor(rep(NA_integer_, sz))
    }

    obj_diag <- TMB::MakeADFun(
      data       = data.list,
      parameters = cur_pl,
      map        = map_diag,
      DLL        = "gllvm_HO",
      silent     = TRUE
    )
    .res_diag <- tryCatch(
      suppressWarnings(
        nlminb(obj_diag$par, obj_diag$fn, obj_diag$gr,
               control = list(rel.tol = 1e-6, iter.max = 50, eval.max = 200))
      ),
      error = function(e) NULL
    )
    if (!is.null(.res_diag)) {
      for (nm in c("Au", "Au_sp")) {
        to_idx   <- names(obj$par)      == nm
        from_idx <- names(obj_diag$par) == nm
        if (any(to_idx) && any(from_idx))
          obj$par[to_idx] <- obj_diag$par[from_idx]
      }
    }
  }

  ## ---- Full optimisation --------------------------------------------------
  opt <- try(
    nlminb(obj$par, obj$fn, obj$gr,
           control = list(rel.tol = reltol, iter.max = maxit,
                          eval.max = maxit * 5, trace = as.integer(trace))),
    silent = TRUE
  )
  if (inherits(opt, "try-error")) {
    warning("nlminb failed, trying optim BFGS")
    opt <- optim(obj$par, obj$fn, obj$gr, method = "BFGS",
                 control = list(reltol = reltol, maxit = maxit, trace = trace))
    opt$convergence <- opt$convergence == 0
  }

  ## ---- Extract parameters -------------------------------------------------
  par_hat <- obj$env$last.par.best
  names(par_hat) <- names(obj$env$par)
  par_lst <- obj$env$parList()   # full list with map applied (padded entries = initial value)

  b_hat      <- matrix(par_hat[names(par_hat) == "b"],    nrow = Kx, ncol = p)
  u_hat      <- matrix(par_hat[names(par_hat) == "u"],    nrow = n,  ncol = d_va_z)
  Au_hat     <- par_hat[names(par_hat) == "Au"]
  alv_hat    <- matrix(par_hat[names(par_hat) == "a_lv_sp"], nrow = p, ncol = d_va_a)
  Au_sp_hat  <- par_hat[names(par_hat) == "Au_sp"]
  sigLV_hat  <- par_hat[names(par_hat) == "sigmaLV"]
  zeta_hat   <- par_hat[names(par_hat) == "zeta"]
  lgphi_hat  <- par_hat[names(par_hat) == "lg_phi"]

  ## Recover ordered sigma from cumulative-sum parameterisation
  ## sigma(k) = sum_{l=k}^{d-1} exp(sigmaLV(l))
  sigma_hat       <- numeric(d)
  sigma_hat[d]    <- exp(sigLV_hat[d])
  if (d > 1) for (k in (d-1):1) sigma_hat[k] <- sigma_hat[k+1] + exp(sigLV_hat[k])

  ## Site VA covariances (d_va_z dimensions)
  if (va_struct_int == 0L) {
    Ai_diag <- matrix(exp(2 * Au_hat), nrow = n, ncol = d_va_z)    # n x d_va_z
    Aj_diag <- matrix(exp(2 * Au_sp_hat), nrow = p, ncol = d_va_a) # p x d_va_a
    A_out    <- Ai_diag
    A_lv_out <- Aj_diag
  } else {
    .build_cov_array <- function(Au_vec, m, dv, tv) {
      arr <- array(0, dim = c(m, dv, dv))
      for (unit in seq_len(m)) {
        L <- matrix(0, dv, dv)
        for (r in seq_len(dv)) {
          for (cc in seq_len(r)) {
            idx <- r*(r-1L)/2L + cc
            val <- Au_vec[(idx - 1L)*m + unit]
            L[r, cc] <- if (r == cc) exp(val) else val
          }
        }
        arr[unit,,] <- L %*% t(L)
      }
      arr
    }
    A_out    <- .build_cov_array(Au_hat,    n, d_va_z, tri_z)  # n x d_va_z x d_va_z
    A_lv_out <- .build_cov_array(Au_sp_hat, p, d_va_a, tri_a)  # p x d_va_a x d_va_a
    Ai_diag <- matrix(0, n, d_va_z)
    Aj_diag <- matrix(0, p, d_va_a)
    for (k in seq_len(d_va_z)) Ai_diag[,k] <- A_out[,k,k]
    for (k in seq_len(d_va_a)) Aj_diag[,k] <- A_lv_out[,k,k]
  }

  ## ---- Parameter extraction helpers -----------------------------------------
  b_z_hat     <- if (Kz > 0L) matrix(par_lst$b_z,     Kz, d) else NULL
  b_gamma_hat <- if (Kt > 0L) matrix(par_lst$b_gamma,  Kt, d) else NULL
  lv_X_mat    <- if (!is.null(lv.X)) as.matrix(lv.X) else NULL
  TR_mat      <- if (!is.null(TR))   as.matrix(TR)   else NULL
  num.lv_unc  <- d - num.RR - num.lv.c

  ## ---- lvs.full: n × d full padded z, HO [RR|lvc|lv] order -----------------
  ## RR dims: deterministic z = lv_X * b_z if Kz>0, else VA u_hat
  ## lvc/lv:  z = u_hat[:, rr_va_z + (k - num.RR)]
  lvs_full <- matrix(0, n, d)
  for (k in seq_len(d)) {
    if (k <= num.RR) {
      if (Kz > 0L && !is.null(b_z_hat) && !is.null(lv_X_mat)) {
        lvs_full[, k] <- lv_X_mat %*% b_z_hat[, k]
      } else {
        lvs_full[, k] <- u_hat[, k]
      }
    } else {
      iz <- rr_va_z + (k - num.RR)
      lvs_full[, k] <- u_hat[, iz]
    }
  }

  ## ---- loadings (p × d unscaled gamma, HO [RR|lvc|lv] order) ---------------
  loadings_full <- matrix(0, p, d)
  for (k in seq_len(d)) {
    if (k <= num.RR) {
      if (Kt > 0L && !is.null(b_gamma_hat) && !is.null(TR_mat)) {
        loadings_full[, k] <- TR_mat %*% b_gamma_hat[, k]
      } else if (d_va_a > 0L) {
        loadings_full[, k] <- alv_hat[, k]
      }
    } else {
      ia <- rr_va_a + (k - num.RR)
      if (ia >= 1L && ia <= d_va_a) loadings_full[, k] <- alv_hat[, ia]
    }
  }

  ## ---- lvs: residual VA scores, n × (num.lv.c + num.lv_unc) -----------------
  ## Matches standard gllvm convention: lvc cols store (u_hat - lv_X * b_z),
  ## lv cols store raw u_hat.  predict.gllvm then adds back lv_X * LvXcoef.
  d_va_lvc_lv <- num.lv.c + num.lv_unc
  if (d_va_lvc_lv == 0L || d_va_z == 0L) {
    lvs_std <- matrix(0, n, 0)
  } else {
    lvs_raw <- if (rr_va_z == 0L) {
      u_hat[, seq_len(d_va_z), drop = FALSE]
    } else {
      u_hat[, seq(rr_va_z + 1L, d_va_z), drop = FALSE]
    }
    if (num.lv.c > 0L && Kz > 0L && !is.null(b_z_hat) && !is.null(lv_X_mat)) {
      for (k in seq_len(num.lv.c)) {
        k_ho <- num.RR + k
        lvs_raw[, k] <- lvs_raw[, k] - lv_X_mat %*% b_z_hat[, k_ho, drop = FALSE]
      }
    }
    lvs_std <- lvs_raw
  }

  ## ---- params$theta: unscaled gamma, standard [lvc|RR|lv] order -------------
  if (num.RR > 0L && num.lv.c > 0L) {
    ho_to_std_idx <- c(seq(num.RR + 1L, num.RR + num.lv.c),
                       seq_len(num.RR),
                       if (num.lv_unc > 0L) seq(num.RR + num.lv.c + 1L, d) else integer(0))
  } else {
    ho_to_std_idx <- seq_len(d)
  }
  theta_std <- loadings_full[, ho_to_std_idx, drop = FALSE]

  ## ---- params$LvXcoef: sigma-scaled b_z, standard [lvc|RR] order ------------
  ## sigma scaling here means predict.gllvm computes the correct
  ## sigma_k * lv_X * b_z_k * gamma_jk without further scaling.
  if (Kz > 0L && (num.RR + num.lv.c) > 0L) {
    d_active <- num.RR + num.lv.c
    lvx_cols_ho <- if (num.RR > 0L && num.lv.c > 0L) {
      c(seq(num.RR + 1L, num.RR + num.lv.c), seq_len(num.RR))
    } else {
      seq_len(d_active)
    }
    LvXcoef_mat <- t(t(b_z_hat[, lvx_cols_ho, drop = FALSE]) * sigma_hat[lvx_cols_ho])
  } else {
    LvXcoef_mat <- NULL
  }

  ## ---- Build output -------------------------------------------------------
  out <- list(
    y            = y,
    X            = if (is.null(X)) NULL else X,
    X.design     = Xmat,
    family       = family,
    num.lv       = num.lv_unc,   # unconstrained LV count (not total d)
    method       = "VA",
    TMB          = TRUE,
    random.loadings = TRUE,
    ## lvs: residual VA scores n × (num.lv.c + num.lv_unc), standard gllvm format
    lvs          = lvs_std,
    ## lvs.full: full n × d unscaled z, HO [RR|lvc|lv] order (for ordiplot)
    lvs.full     = lvs_full,
    ## loadings: full p × d unscaled gamma, HO [RR|lvc|lv] order (for ordiplot)
    loadings     = loadings_full,
    params       = list(
      beta0    = b_hat[1, ],
      Xcoef    = if (Kx > 1) t(b_hat[-1, , drop = FALSE]) else NULL,
      sigma.lv = sigma_hat,
      ## theta: unscaled gamma in standard [lvc|RR|lv] order (used by predict.gllvm)
      theta    = theta_std,
      ## LvXcoef: sigma-scaled b_z in standard [lvc|RR] order (used by getLV/predict)
      LvXcoef  = LvXcoef_mat,
      zeta     = if (has_ordinal && length(zeta_hat) > 0L) zeta_hat else NULL,
      phi      = if (length(lgphi_hat) > 0L) exp(lgphi_hat) else NULL,
      inv.phi  = if (length(lgphi_hat) > 0L) 1 / exp(lgphi_hat) else NULL,
      b_z      = if (Kz > 0L) matrix(par_lst$b_z, Kz, d) else NULL,
      b_gamma  = if (Kt > 0L) matrix(par_lst$b_gamma, Kt, d) else NULL,
      ## sigma.bz / sigma.bgamma: only the lvc entries (RR entries are fixed to 1)
      sigma.bz = {
        n_rr_bz  <- min(num.RR, d_c)
        n_lvc_bz <- d_c - n_rr_bz
        if (Kz > 0L && n_lvc_bz > 0L)
          exp(par_hat[names(par_hat) == "log_sigma_bz"][seq(n_rr_bz + 1L, d_c)])
        else NULL
      },
      sigma.bgamma = {
        n_rr_bt  <- min(num.RR, d_t)
        n_lvc_bt <- d_t - n_rr_bt
        if (Kt > 0L && n_lvc_bt > 0L)
          exp(par_hat[names(par_hat) == "log_sigma_bgamma"][seq(n_rr_bt + 1L, d_t)])
        else NULL
      }
    ),
    ## Variational covariances (VA dims only: d_va_z / d_va_a)
    Lambda.struc = Lambda.struc,
    A            = A_out,
    A_lv         = A_lv_out,
    A_diag       = Ai_diag,
    A_lv_diag    = Aj_diag,
    ## Hessian / SE info — computed by gllvm() post-fitting
    Hess         = NULL,
    sd           = FALSE,
    ## Fields expected by generic S3 methods
    num.lv.c     = num.lv.c,
    num.RR       = num.RR,
    num.lvcor    = 0L,
    randomB      = randomB,
    randomT      = randomT,
    quadratic    = FALSE,
    lv.X         = lv.X,
    lv.X.design  = lv.X,
    TR           = TR,
    col.eff      = list(col.eff = FALSE),
    row.eff      = if (isFALSE(row.eff)) FALSE else row.eff,
    zeta.struc   = zeta.struc,
    ## Optimisation info
    logL         = -opt$objective,
    convergence  = if (is.numeric(opt$convergence)) opt$convergence == 0 else opt$convergence,
    TMBfn        = obj,
    optim.method = "nlminb",
    Ntrials      = Ntrials,
    offset       = offset_mat,
    call         = if (!is.null(call.)) call. else match.call()
  )
  if (!is.null(opt$par)) out$TMBfn$par <- opt$par

  class(out) <- c("gllvmHO", "gllvm")
  out
}

## Helper: compute link-scale residuals from column-wise GLM intercepts
#' @keywords internal
.ho_link_residuals <- function(y, fam_int, Xmat) {
  n <- nrow(y); p <- ncol(y)
  beta0 <- numeric(p)
  resid <- matrix(0, n, p)

  for (j in seq_len(p)) {
    yj <- y[, j]
    ok <- !is.na(yj)
    fj <- fam_int[j]

    # Fit species-wise GLM intercept
    if (fj == 0L || fj == 6L) {            # Poisson / ZIP
      mu_j  <- mean(yj[ok]) + 0.5
      beta0[j] <- log(mu_j)
      resid[ok, j] <- log(pmax(yj[ok], 0.5)) - beta0[j]
    } else if (fj == 3L) {                 # Gaussian
      beta0[j] <- mean(yj[ok])
      resid[ok, j] <- yj[ok] - beta0[j]
    } else if (fj == 1L || fj == 11L) {    # NB / ZINB
      mu_j <- mean(yj[ok]) + 0.5
      beta0[j] <- log(mu_j)
      resid[ok, j] <- log(pmax(yj[ok], 0.5)) - beta0[j]
    } else if (fj == 2L) {                 # Binomial (logit)
      pj <- mean(yj[ok]) + 1e-3
      pj <- min(pj, 1 - 1e-3)
      beta0[j] <- log(pj / (1 - pj))
      resid[ok, j] <- log((pmax(yj[ok], 1e-3)) / (1 - pmin(yj[ok], 1 - 1e-3))) - beta0[j]
    } else if (fj == 4L) {                 # Gamma (log link)
      mu_j <- mean(yj[ok]) + 0.5
      beta0[j] <- log(mu_j)
      resid[ok, j] <- log(pmax(yj[ok], 0.5)) - beta0[j]
    } else {
      beta0[j] <- 0
      resid[ok, j] <- yj[ok] - mean(yj[ok])
    }
    # Centre residuals
    resid[, j] <- resid[, j] - mean(resid[ok, j])
  }

  list(residuals = resid, beta0 = beta0)
}

## Helper: truncated SVD → initial u (n × d_va_z), a_sp (p × d_va_a), sigmaLV (d_total)
#' @keywords internal
.ho_svd_starts <- function(R, d_va_z, d_va_a, beta0, d_total = max(d_va_z, d_va_a)) {
  d_svd <- max(d_va_z, d_va_a, d_total)
  sv <- tryCatch(
    svd(R, nu = d_svd, nv = d_svd),
    error = function(e) {
      message("SVD failed, using random starts: ", conditionMessage(e))
      list(u = matrix(rnorm(nrow(R) * d_svd, sd = 0.1), nrow(R), d_svd),
           v = matrix(rnorm(ncol(R) * d_svd, sd = 0.1), ncol(R), d_svd),
           d = rep(1, d_svd))
    }
  )

  d_trunc <- min(d_svd, length(sv$d))
  svals   <- sv$d[seq_len(d_trunc)]
  U_full  <- sv$u[, seq_len(d_trunc), drop = FALSE]
  V_full  <- sv$v[, seq_len(d_trunc), drop = FALSE]

  if (d_trunc < d_svd) {
    U_full <- cbind(U_full, matrix(rnorm(nrow(R) * (d_svd - d_trunc), sd = 0.01), nrow(R)))
    V_full <- cbind(V_full, matrix(rnorm(ncol(R) * (d_svd - d_trunc), sd = 0.01), ncol(R)))
    svals  <- c(svals, rep(max(svals) * 0.1, d_svd - d_trunc))
  }

  U <- U_full[, seq_len(d_va_z), drop = FALSE]
  V <- V_full[, seq_len(d_va_a), drop = FALSE]

  svals_scaled <- pmin(svals / sqrt(nrow(R)), 1.0)
  svals_scaled <- pmax(svals_scaled, 1e-3)
  # Pad to d_total if d_svd > d_total (shouldn't happen) or truncate
  svals_d <- svals_scaled[seq_len(d_total)]

  sigmaLV <- numeric(d_total)
  sigmaLV[d_total] <- log(svals_d[d_total])
  if (d_total > 1)
    for (k in (d_total - 1):1)
      sigmaLV[k] <- log(max(svals_d[k] - svals_d[k + 1], 1e-6))

  list(u = U, a_sp = V, sigmaLV = sigmaLV, beta0 = beta0)
}

##############################################################################
## summary.gllvmHO
##############################################################################

#' Summary method for Hierarchical Ordination models
#' @export
#' @keywords internal
summary.gllvmHO <- function(object, ...) {
  cat("Hierarchical Ordination model (VA)\n")
  cat("Call:", deparse(object$call), "\n\n")
  cat("Family:", paste(unique(object$family), collapse = ", "), "\n")
  n <- nrow(object$lvs)
  p <- nrow(object$loadings)
  d <- object$num.lv
  cat("n =", n, "  p =", p, "  d =", d, "\n\n")

  ## Ordination scale
  sig <- object$params$sigma.lv
  cat("Ordination scale (sigma_1 >= ... >= sigma_d):\n")
  names(sig) <- paste0("dim", seq_along(sig))
  print(round(sig, 4))
  cat("\n")


  ## Fixed-effect intercepts
  cat("Species intercepts (beta0):\n")
  b0 <- object$params$beta0
  if (!is.null(object$sd) && !isFALSE(object$sd) && !is.null(object$sd$beta0)) {
    df_b0 <- data.frame(
      Estimate = round(b0, 4),
      Std.Err  = round(object$sd$beta0, 4)
    )
    rownames(df_b0) <- if (!is.null(colnames(object$y))) colnames(object$y) else
      paste0("sp", seq_len(p))
    print(df_b0)
  } else {
    print(round(b0, 4))
  }
  cat("\n")

  ## b_z
  if (!is.null(object$params$b_z)) {
    cat("Canonical covariate coefficients (b_z):\n")
    bz <- object$params$b_z
    rownames(bz) <- if (!is.null(colnames(object$lv.X))) colnames(object$lv.X) else
      paste0("cov", seq_len(nrow(bz)))
    colnames(bz) <- paste0("dim", seq_len(ncol(bz)))
    print(round(bz, 4))
    cat("\n")
  }

  ## b_gamma
  if (!is.null(object$params$b_gamma)) {
    cat("Trait coefficients (b_gamma):\n")
    bg <- object$params$b_gamma
    rownames(bg) <- if (!is.null(colnames(object$TR))) colnames(object$TR) else
      paste0("trait", seq_len(nrow(bg)))
    colnames(bg) <- paste0("dim", seq_len(ncol(bg)))
    print(round(bg, 4))
    cat("\n")
  }

  cat("Log-likelihood:", round(object$logL, 4), "\n")
  cat("Converged:", isTRUE(object$convergence), "\n")
  invisible(object)
}

## Helper: rebuild named parameter list from a flat vector + template list
#' @keywords internal
.par_to_list <- function(flat, template) {
  out  <- template
  lens <- sapply(template, length)
  idx  <- 1L
  for (nm in names(template)) {
    sz <- prod(dim(template[[nm]]))
    if (sz == 0L) next
    chunk <- flat[idx:(idx + sz - 1L)]
    if (is.matrix(template[[nm]]))
      out[[nm]] <- matrix(chunk, nrow = nrow(template[[nm]]))
    else
      out[[nm]] <- chunk
    idx <- idx + sz
  }
  out
}
