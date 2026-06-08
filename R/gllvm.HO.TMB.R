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
    y, X = NULL, family = "poisson", num.lv = 2, offset = NULL,
    Ntrials = matrix(1L), row.eff = FALSE, maxit = 2000, reltol = 1e-8,
    diag.iter = 1, start.params = NULL, trace = FALSE
) {
  .ensure_ho_dll()

  ## ---- dimensions & family setup ------------------------------------------
  n  <- nrow(y)
  p  <- ncol(y)
  d  <- num.lv
  if (length(family) == 1L) family <- rep(family, p)
  fam_int <- .fam_code[family]
  if (any(is.na(fam_int)))
    stop("Unknown family: ", paste(family[is.na(fam_int)], collapse = ", "))

  ## extra link flags (0 = canonical/logit for binary, etc.)
  extra <- rep(0, p)

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

  ## ---- SVD starting values ------------------------------------------------
  link_y <- .ho_link_residuals(y, fam_int, Xmat)
  sv     <- .ho_svd_starts(link_y$residuals, d, link_y$beta0)

  if (!is.null(start.params)) {
    if (!is.null(start.params$lvs))
      sv$u      <- as.matrix(start.params$lvs)
    if (!is.null(start.params$loadings))
      sv$a_sp   <- as.matrix(start.params$loadings)
  }

  ## ---- dispersion starting values -----------------------------------------
  lg_phi    <- rep(0, p)     # log(phi) for NB/Gaussian/Gamma
  lg_phiZINB<- rep(0, p)
  zeta      <- numeric(0)
  ePower    <- 0             # Tweedie power (logit-scale)

  ## ---- row-effect setup ---------------------------------------------------
  xr     <- matrix(0, 0, 0)
  r0f    <- matrix(0, 0, 1)
  r0r    <- matrix(0, 0, 1)
  lg_Ar  <- numeric(0)
  log_sigma_r <- numeric(0)
  random_flag <- 0L

  if (isTRUE(row.eff)) {
    r0r     <- matrix(0, n, 1)
    lg_Ar   <- rep(log(sqrt(0.1)), n)   # small initial row-effect variance
    log_sigma_r <- log(0.3)
    random_flag <- 1L
  }

  ## ---- VA covariance starting values (log-Chol diagonals) ----------------
  ## Au(k*n+i)   = log(sqrt(A_i(k,k)))  → A_i(k,k) = exp(2*Au(k*n+i))
  ## Au_sp(k*p+j) = log(sqrt(A_j(k,k))) → A_j(k,k) = exp(2*Au_sp(k*p+j))
  Au_init    <- rep(log(sqrt(0.1)), d * n)   # A_i(k,k) = 0.1 (small start keeps ck > 0)
  Au_sp_init <- rep(log(sqrt(0.1)), d * p)   # A_j(k,k) = 0.1

  ## ---- data list ----------------------------------------------------------
  data.list <- list(
    y         = as.matrix(y),
    x         = Xmat,
    xr        = xr,
    offset    = offset_mat,
    Ntrials   = Ntrials,
    family    = fam_int,
    extra     = extra,
    num_lv    = as.integer(d),
    method    = 0L,              # VA
    zetastruc = 0L,
    p_betaH   = 0L,
    random    = random_flag
  )

  ## ---- parameter list -----------------------------------------------------
  param.list <- list(
    b         = rbind(sv$beta0, matrix(0, Kx - 1, p)),  # Kx x p
    u         = sv$u,          # n x d
    Au        = Au_init,       # n*d
    a_lv_sp   = sv$a_sp,       # p x d
    Au_sp     = Au_sp_init,    # p*d
    sigmaLV   = sv$sigmaLV,    # d
    lg_phi    = lg_phi,
    lg_phiZINB= lg_phiZINB,
    zeta      = zeta,
    ePower    = ePower,
    r0f       = r0f,
    r0r       = r0r,
    lg_Ar     = lg_Ar,
    log_sigma = log_sigma_r
  )

  ## ---- map: fix params that are unused ------------------------------------
  map.list <- list()
  map.list$zeta       <- factor(rep(NA, length(zeta)))
  map.list$lg_phiZINB <- factor(rep(NA, length(lg_phiZINB)))
  map.list$ePower     <- factor(NA)
  # Fix dispersion for Poisson (no phi needed)
  poisson_cols <- which(fam_int == 0L)
  if (length(poisson_cols) > 0) {
    phi_map <- 1:p
    phi_map[poisson_cols] <- NA
    map.list$lg_phi <- factor(phi_map)
  }
  # Fix row-effect params if not used
  if (!isTRUE(row.eff)) {
    map.list$r0f       <- factor(rep(NA, length(r0f)))
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

    # Fix everything except Au and Au_sp (use full lengths from param.list)
    map_diag <- map.list
    for (nm in setdiff(names(param.list), c("Au", "Au_sp"))) {
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
    try(nlminb(obj_diag$par, obj_diag$fn, obj_diag$gr,
               control = list(rel.tol = 1e-6, iter.max = 50, eval.max = 200)),
        silent = TRUE)

    # Transfer updated Au and Au_sp back to main obj$par (by name)
    for (nm in c("Au", "Au_sp")) {
      to_idx   <- names(obj$par)      == nm
      from_idx <- names(obj_diag$par) == nm
      if (any(to_idx) && any(from_idx))
        obj$par[to_idx] <- obj_diag$par[from_idx]
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

  b_hat      <- matrix(par_hat[names(par_hat) == "b"],    nrow = Kx, ncol = p)
  u_hat      <- matrix(par_hat[names(par_hat) == "u"],    nrow = n,  ncol = d)
  Au_hat     <- par_hat[names(par_hat) == "Au"]
  alv_hat    <- matrix(par_hat[names(par_hat) == "a_lv_sp"], nrow = p, ncol = d)
  Au_sp_hat  <- par_hat[names(par_hat) == "Au_sp"]
  sigLV_hat  <- par_hat[names(par_hat) == "sigmaLV"]

  ## Recover ordered sigma from cumulative-sum parameterisation
  ## sigma(k) = sum_{l=k}^{d-1} exp(sigmaLV(l))
  sigma_hat       <- numeric(d)
  sigma_hat[d]    <- exp(sigLV_hat[d])
  if (d > 1) for (k in (d-1):1) sigma_hat[k] <- sigma_hat[k+1] + exp(sigLV_hat[k])

  ## Site VA covariances (diagonal, stored as n x d)
  Ai_diag <- matrix(exp(2 * Au_hat), nrow = n, ncol = d)   # A_i(k,k)
  Aj_diag <- matrix(exp(2 * Au_sp_hat), nrow = p, ncol = d) # A_j(k,k)

  ## ---- Build output -------------------------------------------------------
  out <- list(
    y            = y,
    X            = if (is.null(X)) NULL else X,
    X.design     = Xmat,
    family       = family,
    num.lv       = d,
    method       = "VA",
    random.loadings = TRUE,
    ## Site scores: VA means (a_i)
    lvs          = u_hat,
    ## Species loadings: VA means (a_j), unscaled
    loadings     = alv_hat,
    ## Fitted params
    ## theta is UNSCALED (= alv_hat); ordiplot/getResidualCov apply sigma.lv separately
    params       = list(
      beta0    = b_hat[1, ],
      Xcoef    = if (Kx > 1) t(b_hat[-1, , drop = FALSE]) else NULL,
      sigma.lv = sigma_hat,
      theta    = alv_hat           # unscaled species loadings (p x d)
    ),
    ## Variational covariance diagonals (stored as matrices, not 3D arrays)
    A            = Ai_diag,        # n x d  diag entries of A_i
    A_lv         = Aj_diag,        # p x d  diag entries of A_j
    ## Fields expected by generic S3 methods inherited from gllvm
    num.lv.c     = 0L,
    num.RR       = 0L,
    num.lvcor    = 0L,
    randomB      = FALSE,
    quadratic    = FALSE,
    lv.X         = NULL,
    lv.X.design  = NULL,
    col.eff      = list(col.eff = FALSE),
    row.eff      = if (isTRUE(row.eff)) "random" else FALSE,
    sd           = TRUE,           # non-FALSE: getPredictErr.gllvmHO provides errors
    ## Optimisation info
    logL         = -opt$objective,
    convergence  = if (is.numeric(opt$convergence)) opt$convergence == 0 else opt$convergence,
    TMBfn        = obj,
    optim.method = "nlminb",
    Ntrials      = Ntrials,
    offset       = offset_mat,
    call         = match.call()
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

## Helper: truncated SVD → initial u, a_sp, sigmaLV
#' @keywords internal
.ho_svd_starts <- function(R, d, beta0) {
  sv <- tryCatch(
    svd(R, nu = d, nv = d),
    error = function(e) {
      message("SVD failed, using random starts: ", conditionMessage(e))
      list(u = matrix(rnorm(nrow(R) * d, sd = 0.1), nrow(R), d),
           v = matrix(rnorm(ncol(R) * d, sd = 0.1), ncol(R), d),
           d = rep(1, d))
    }
  )

  d_trunc <- min(d, length(sv$d))
  svals   <- sv$d[seq_len(d_trunc)]
  U       <- sv$u[, seq_len(d_trunc), drop = FALSE]
  V       <- sv$v[, seq_len(d_trunc), drop = FALSE]

  # Pad if fewer singular values returned than d
  if (d_trunc < d) {
    U <- cbind(U, matrix(rnorm((nrow(R)) * (d - d_trunc), sd = 0.01), nrow(R), d - d_trunc))
    V <- cbind(V, matrix(rnorm((ncol(R)) * (d - d_trunc), sd = 0.01), ncol(R), d - d_trunc))
    svals <- c(svals, rep(max(svals) * 0.1, d - d_trunc))
  }

  # Scale down singular values so ck = 1 - sigma^2 * A_i * A_j > 0 at init.
  # With Au_init -> A_diag = 0.1, need sigma < 1/sqrt(0.01) = 10; cap at 1.
  svals_scaled <- pmin(svals / sqrt(nrow(R)), 1.0)
  svals_scaled <- pmax(svals_scaled, 1e-3)   # avoid log(0)

  # Cumulative-sum parameterisation: sigma(k) = sum_{l=k}^{d} exp(sigmaLV(l))
  # Invert: delta(d) = sigma(d);  delta(k) = sigma(k) - sigma(k+1)  for k < d
  # sigmaLV(k) = log(delta(k))
  sigmaLV <- numeric(d)
  sigmaLV[d] <- log(svals_scaled[d])
  if (d > 1) {
    for (k in (d - 1):1) {
      sigmaLV[k] <- log(max(svals_scaled[k] - svals_scaled[k + 1], 1e-6))
    }
  }

  list(u = U, a_sp = V, sigmaLV = sigmaLV, beta0 = beta0)
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
