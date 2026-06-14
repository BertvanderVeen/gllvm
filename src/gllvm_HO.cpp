#include <TMB.hpp>
#include <cmath>
#include "distrib.h"
#include "init.h"
#include "utils.h"
#include <R_ext/Error.h>

// Family codes — must match gllvm.cpp
enum Family : int {
  POISSON = 0, NEG_BINOMIAL = 1, BINOMIAL = 2, GAUSSIAN = 3, GAMMA = 4,
  TWEEDIE = 5, ZIP = 6, ORDINAL = 7, EXPONENTIAL = 8, BETA = 9,
  BETA_HURDLE = 10, ZINB = 11, ORDERED_BETA = 12, ZIB = 13, ZNIB = 14,
  BETA_BINOMIAL = 15
};

//-----------------------------------------------------------------------
// GLLVM Hierarchical Ordination (HO_VA)
//
// CORE INTEGRATION PATTERN: E[Z · Sigma · Gamma]
//   For every dim type (RR, lvc, lv), the ELBO contains the bilinear
//   integral E_{q(z_i)q(gamma_j)}[z_ik * sigma_k * gamma_jk].  The VA
//   approximation factorizes: q(z_i) and q(gamma_j) are Gaussian.  The
//   only thing that changes across dim types is how the VA means/variances
//   are set up:
//
//   num_lv (unconstrained):
//     z_i   ~ q(z_i)     = N(u_i,   A_i),   prior N(0, I)
//     gamma_j ~ q(gamma_j) = N(a_j,   B_j),   prior N(0, I)
//
//   num_lv.c (concurrent, lvc):
//     z_i   ~ q(z_i)     = N(u_i,   A_i),   prior N(lv_X_i * b_z,   I)
//     gamma_j ~ q(gamma_j) = N(a_j,   B_j),   prior N(TR_j   * b_gamma, I)
//     b_z and b_gamma are themselves VA parameters with own Gaussian posteriors
//     q(b_z[:,l]) = N(b_z_hat[:,l], AB_z(l)*AB_z(l)^T)
//
//   num_RR (reduced-rank, fully constrained):
//     z_i   = lv_X_i * b_z  (deterministic given b_z — no VA residual)
//     gamma_j = TR_j   * b_gamma (or VA a_j when no TR)
//     sigma_bz is not separately identified from sigma; fix log_sigma_bz to 0 for RR dims.
//
// CROSS-COVARIANCE CORRECTION (lvc dims only, when b_z/b_gamma are random VA params):
//   The factorized VA q(z_i) q(b_z) means the z_i KL must integrate over q(b_z):
//     E_{q(b_z)}[ ||z_i - lv_X_i*b_z||^2 ] = ||u_i - lv_X_i*b_z_hat||^2
//                                            + lv_X_i · Var_q(b_z) · lv_X_i^T
//   The second term is the cross-covariance correction.  Without it, the ELBO is
//   monotonically increasing in sigma_bz → sigma_bz explodes to infinity.
//   With it, the optimal sigma_bz* = (||b_z_hat||^2 / ||lv_X||_F^2)^{1/4} is finite.
//   See ms.pdf §3.2 for the full ELBO derivation.
//
// LV BLOCKS (ordered):  [RR | lvc | lv]
//   k = 0..num_RR-1               : RR dims
//   k = num_RR..num_RR+num_lvc-1  : concurrent (lvc) dims
//   k = num_RR+num_lvc..d-1       : unconstrained (lv) dims
//
// VA dimension indices:
//   z_i VA  (d_va_z dims):  lvc + lv dims always VA;
//                            RR dims VA iff Kz==0 (no X covariates)
//   gamma_j VA (d_va_a dims): lvc + lv always VA;
//                              RR dims VA iff Kt==0 (no TR traits)
//
// b_z/b_gamma active columns: 0..d_c-1 where d_c = min(Kz, num_RR+num_lvc)
//   lv dims (k >= num_RR+num_lvc) never use covariates/traits
//-----------------------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  // ===== DATA =====
  DATA_MATRIX(y);            // n x p responses
  DATA_MATRIX(x);            // n x Kx fixed covariates (first col = intercept)
  DATA_MATRIX(xr);           // fixed row-effect design matrix
  DATA_SPARSE_MATRIX(dr0);   // random row-effect design matrix (sparse)
  DATA_MATRIX(offset);       // n x p offsets
  DATA_IMATRIX(Ntrials);     // n x p binomial trials
  DATA_IVECTOR(family);      // family code per species, length p
  DATA_VECTOR(extra);        // link/variant flag per species, length p
  DATA_INTEGER(num_lv);      // total ordination dimension d = num_RR + num_lvc + num_lv_unc
  DATA_INTEGER(num_RR);      // number of RR (reduced-rank) dims
  DATA_INTEGER(num_lvc);     // number of concurrent (constrained) dims
  DATA_INTEGER(method);      // must be 0 (VA)
  DATA_INTEGER(zetastruc);   // ordinal cutoff structure: 0=common, 1=species-specific
  DATA_INTEGER(p_betaH);     // number of betaH columns (beta-hurdle)
  DATA_IVECTOR(random);      // (0)=row RE, (1)=TMBtrait Br, (2)=b_lv slopes, (3)=formula Br
  DATA_INTEGER(va_struct);   // 0 = diagonal VA covariance; 1 = unstructured (full Cholesky)
  DATA_MATRIX(lv_X_env);     // n x Kz  (0-column matrix when not used)
  DATA_MATRIX(TR);           // p x Kt  (0-column matrix when not used)
  DATA_INTEGER(randomB);     // 1 = "LV" prior on b_z
  DATA_INTEGER(randomT);     // 1 = "LV" prior on b_gamma
  DATA_IMATRIX(csb_z);       // (n_pairs x 2) 1-based predictor-pair indices for b_z
  DATA_IMATRIX(csb_gamma);   // same for b_gamma
  DATA_MATRIX(xb);           // n x spnr  formula design matrix (0-col when unused)
  // Species-specific random effects (shared with species_effects.h)
  DATA_IMATRIX(cs);          // correlation index pairs for species random slopes
  DATA_STRUCT(colMatBlocksI, gllvmutils::dclist); // phylogenetic correlation blocks
  DATA_IMATRIX(nncolMat);    // nearest-neighbour matrix for sparse phylo approx
  DATA_VECTOR(Abranks);      // ranks for species VA structure
  DATA_INTEGER(Abstruc);     // 0=diag/blockdiag, 1=MN, 2=CL2, 3=CL1, 4=CL2, 5=unstructured
  // Row random effects (shared with row_effects.h)
  DATA_IMATRIX(trmsize);     // 2-row: row1=term LHS size, row2=group count
  DATA_IMATRIX(csR);         // correlation pairs for row effects
  DATA_IVECTOR(cstruc);      // correlation structure codes
  DATA_STRUCT(proptoMats, gllvmutils::nesteddclist); // proportion matrices for row RE
  DATA_STRUCT(dc, gllvmutils::dclist); // site coordinates for exp-decay row RE

  // ===== PARAMETERS =====
  PARAMETER_MATRIX(b);       // Kx x p fixed effects
  PARAMETER_MATRIX(Br);      // spnr x p  species-specific random slopes for xb
  PARAMETER_MATRIX(B);       // community-level effects / RE means for formula
  PARAMETER_VECTOR(sigmaB);  // log-SDs (+ phylo range pars) for Br
  PARAMETER_VECTOR(sigmaij); // off-diagonal covariance terms for species slopes
  PARAMETER_VECTOR(Abb);     // VA covariance entries for Br

  // Site VA parameters — covers d_va_z dims (computed below)
  PARAMETER_MATRIX(u);       // n x d_va_z
  PARAMETER_VECTOR(Au);      // n*d_va_z (diag) or n*tri_z (unstructured)

  // Species VA parameters — covers d_va_a dims
  PARAMETER_MATRIX(a_lv_sp); // p x d_va_a
  PARAMETER_VECTOR(Au_sp);   // p*d_va_a (diag) or p*tri_a (unstructured)

  PARAMETER_MATRIX(b_z);     // Kz x d  (cols >= d_c and >= num_RR+num_lvc mapped to 0)
  PARAMETER_MATRIX(b_gamma); // Kt x d  (same padding)
  PARAMETER_VECTOR(Ab_z);
  PARAMETER_VECTOR(log_sigma_bz);
  PARAMETER_VECTOR(Ab_gamma);
  PARAMETER_VECTOR(log_sigma_bgamma);

  // sigma(0) = exp(sigmaLV(0)); sigma(k) = sigma(k-1)*exp(-exp(sigmaLV(k))), k>=1
  PARAMETER_VECTOR(sigmaLV); // length d

  PARAMETER_VECTOR(lg_phi);
  PARAMETER_VECTOR(lg_phiZINB);
  PARAMETER_VECTOR(zeta);
  PARAMETER(ePower);

  PARAMETER_MATRIX(r0f);
  PARAMETER_MATRIX(r0r);
  PARAMETER_VECTOR(lg_Ar);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER_VECTOR(sigmaijr); // correlation parameters for row random effects

  // ===== DERIVED SCALARS =====
  int n   = y.rows();
  int p   = y.cols();
  int truep = p - p_betaH;
  int d   = num_lv;             // total dims
  int num_lv_unc = d - num_RR - num_lvc;

  int Kz  = lv_X_env.cols();
  int Kt  = TR.cols();

  // RR dims in the VA-z / VA-a arrays:
  //   if Kz>0, RR z is deterministic → rr_va_z = 0
  //   if Kz==0, RR z is random       → rr_va_z = num_RR
  int rr_va_z = (Kz > 0) ? 0 : num_RR;
  int rr_va_a = (Kt > 0) ? 0 : num_RR;
  int d_va_z  = rr_va_z + num_lvc + num_lv_unc;
  int d_va_a  = rr_va_a + num_lvc + num_lv_unc;

  // Active b_z/b_gamma columns: lv dims never use covariates/traits
  int d_active = num_RR + num_lvc;
  int d_c = (Kz > 0 && Kz < d_active) ? Kz : (Kz > 0 ? d_active : 0);
  int d_t = (Kt > 0 && Kt < d_active) ? Kt : (Kt > 0 ? d_active : 0);

  // Sign identification: for each VA-a dim ia, a_lv_sp(ia, ia) is stored on the
  // log scale so that exp(a_lv_sp(ia, ia)) > 0 always, breaking the per-dimension
  // sign symmetry (z_ik, gamma_jk) -> (-z_ik, -gamma_jk).
  matrix<Type> a_lv_sp_id = a_lv_sp;
  for (int ia = 0; ia < d_va_a && ia < p; ia++)
    a_lv_sp_id(ia, ia) = exp(a_lv_sp(ia, ia));

  vector<Type> iphi = exp(lg_phi);

  parallel_accumulator<Type> nll(this);

  // ===== SIGMA construction =====
  // sigma(0) = exp(sigmaLV(0)); sigma(k) = sigma(k-1)*exp(-exp(sigmaLV(k))), k>=1
  vector<Type> sigma(d);
  sigma(0) = exp(sigmaLV(0));
  for (int k = 1; k < d; k++) sigma(k) = sigma(k-1) * exp(-exp(sigmaLV(k)));
  REPORT(sigma);
  vector<Type> sigma2(d);
  for (int k = 0; k < d; k++) sigma2(k) = sigma(k) * sigma(k);
  REPORT(sigma2);

  // Log-barrier promoting distinct sigma values: log(1 - ratio(k)) -> -inf as ratio -> 1
  for (int k = 1; k < d; k++)
    nll -= log(Type(1) - exp(-exp(sigmaLV(k))));

  // ===== VARIATIONAL COVARIANCES =====
  // Ai_diag: n x d_va_z  (diagonal variances for site VA)
  // Aj_diag: p x d_va_a  (diagonal variances for species VA)
  int tri_z = d_va_z*(d_va_z+1)/2;
  int tri_a = d_va_a*(d_va_a+1)/2;

  matrix<Type> Ai_diag(n, d_va_z);
  matrix<Type> Aj_diag(p, d_va_a);

  if (va_struct == 0) {
    for (int iz = 0; iz < d_va_z; iz++)
      for (int i = 0; i < n; i++) Ai_diag(i,iz) = exp(Type(2)*Au(iz*n+i));
    for (int ia = 0; ia < d_va_a; ia++)
      for (int j = 0; j < p; j++) Aj_diag(j,ia) = exp(Type(2)*Au_sp(ia*p+j));
  } else {
    for (int i = 0; i < n; i++) {
      for (int r = 0; r < d_va_z; r++) {
        Type acc = Type(0);
        for (int c = 0; c <= r; c++) {
          int idx = r*(r+1)/2+c;
          Type lval = Au(idx*n+i);
          if (c == r) lval = exp(lval);
          acc += lval*lval;
        }
        Ai_diag(i,r) = acc;
      }
    }
    for (int j = 0; j < p; j++) {
      for (int r = 0; r < d_va_a; r++) {
        Type acc = Type(0);
        for (int c = 0; c <= r; c++) {
          int idx = r*(r+1)/2+c;
          Type lval = Au_sp(idx*p+j);
          if (c == r) lval = exp(lval);
          acc += lval*lval;
        }
        Aj_diag(j,r) = acc;
      }
    }
  }

  // ===== INDEX MAPPING =====
  // va_z_idx(k): full dim k → VA-z index (-1 if z is deterministic for this dim)
  // va_a_idx(k): full dim k → VA-a index (-1 if gamma is deterministic)
  //
  // Layout: VA-z = [RR (if Kz==0) | lvc | lv]
  //         VA-a = [RR (if Kt==0) | lvc | lv]
  std::vector<int> va_z_idx(d), va_a_idx(d);
  for (int k = 0; k < d; k++) {
    if (k < num_RR) {
      va_z_idx[k] = (rr_va_z > 0) ? k : -1;
      va_a_idx[k] = (rr_va_a > 0) ? k : -1;
    } else {
      va_z_idx[k] = rr_va_z + (k - num_RR);
      va_a_idx[k] = rr_va_a + (k - num_RR);
    }
  }

  // Reverse mapping: VA-z index → full dim k  (for KL prior means)
  // VA-z[0..rr_va_z-1] → full k = 0..rr_va_z-1  (RR dims, only when Kz==0)
  // VA-z[rr_va_z..rr_va_z+num_lvc-1] → full k = num_RR..num_RR+num_lvc-1 (lvc)
  // VA-z[rr_va_z+num_lvc..] → full k = num_RR+num_lvc..d-1 (lv)
  std::vector<int> va_z_to_full(d_va_z), va_a_to_full(d_va_a);
  {
    int iz = 0;
    if (rr_va_z > 0) for (int k = 0; k < num_RR; k++) va_z_to_full[iz++] = k;
    for (int k = num_RR; k < d; k++) va_z_to_full[iz++] = k;
  }
  {
    int ia = 0;
    if (rr_va_a > 0) for (int k = 0; k < num_RR; k++) va_a_to_full[ia++] = k;
    for (int k = num_RR; k < d; k++) va_a_to_full[ia++] = k;
  }

  // ===== VA CHOLESKY FACTORS FOR RANDOM b_z AND b_gamma =====
  // Build AB_z(l) and AB_g(l) — Cholesky factors of the VA covariances for b_z[:,l]
  // and b_gamma[:,l] respectively.  These are needed in two places:
  //   (1) cross-covariance corrections in the z_i / gamma_j KL loops below
  //   (2) the b_z / b_gamma prior KL terms
  // Building them once here avoids duplication.
  //
  // Ab_z layout: diagonal entries at q*d_c+l; off-diag at Kz*d_c + k*d_c + l.
  // AB_z(l)(q,q) = exp(Ab_z(q*d_c+l));  off-diag entries stored raw (no exp).

  int dc_lvc = d_c - std::min(d_c, num_RR);  // active lvc cols in b_z
  int dt_lvc = d_t - std::min(d_t, num_RR);  // active lvc cols in b_gamma

  bool has_bz = (randomB == 1 && Kz > 0 && d_c > 0);
  bool has_bg = (randomT == 1 && Kt > 0 && d_t > 0);

  // Use size max(.,1) to avoid 0-size Eigen arrays when the parameter is absent.
  vector<matrix<Type>> AB_z(d_c > 0 ? d_c : 1);
  if (has_bz) {
    for (int l = 0; l < d_c; l++) { AB_z(l).resize(Kz, Kz); AB_z(l).setZero(); }
    for (int l = 0; l < d_c; l++)
      for (int q = 0; q < Kz; q++)
        AB_z(l)(q, q) = exp(Ab_z(q * d_c + l));
    if (Ab_z.size() > Kz * d_c) {
      int k = 0;
      for (int c = 0; c < Kz; c++) for (int r = c + 1; r < Kz; r++) {
        for (int l = 0; l < d_c; l++) AB_z(l)(r, c) = Ab_z(Kz * d_c + k * d_c + l);
        k++;
      }
    }
  }

  vector<matrix<Type>> AB_g(d_t > 0 ? d_t : 1);
  if (has_bg) {
    for (int l = 0; l < d_t; l++) { AB_g(l).resize(Kt, Kt); AB_g(l).setZero(); }
    for (int l = 0; l < d_t; l++)
      for (int q = 0; q < Kt; q++)
        AB_g(l)(q, q) = exp(Ab_gamma(q * d_t + l));
    if (Ab_gamma.size() > Kt * d_t) {
      int k = 0;
      for (int c = 0; c < Kt; c++) for (int r = c + 1; r < Kt; r++) {
        for (int l = 0; l < d_t; l++) AB_g(l)(r, c) = Ab_gamma(Kt * d_t + k * d_t + l);
        k++;
      }
    }
  }

  // bz_var_mat(Kz x Kz) = sum_{lvc dims l} AB_z(l) * AB_z(l)^T
  //   = total VA variance of b_z across all lvc dimensions, per covariate pair (q,q')
  // Cross-covariance correction for site i: lv_X_i^T * bz_var_mat * lv_X_i
  // (diagonal case: bz_var_mat is diagonal, entry q = sum_l exp(2*Ab_z(q*d_c+l)))
  matrix<Type> bz_var_mat(Kz > 0 ? Kz : 1, Kz > 0 ? Kz : 1); bz_var_mat.setZero();
  if (has_bz && dc_lvc > 0) {
    for (int iz = 0; iz < dc_lvc; iz++) {
      int l = num_RR + iz;
      bz_var_mat += AB_z(l) * AB_z(l).transpose();
    }
  }

  matrix<Type> bg_var_mat(Kt > 0 ? Kt : 1, Kt > 0 ? Kt : 1); bg_var_mat.setZero();
  if (has_bg && dt_lvc > 0) {
    for (int ia = 0; ia < dt_lvc; ia++) {
      int l = num_RR + ia;
      bg_var_mat += AB_g(l) * AB_g(l).transpose();
    }
  }

  // ===== KL(q(z_i) || N(mu_z_i, I)) =====
  // Only over d_va_z dims (det-z dims have no KL contribution).
  // Prior mean mu_z[iz]:
  //   iz < rr_va_z (RR, Kz==0): mean = 0
  //   rr_va_z <= iz < rr_va_z+num_lvc (lvc, k = num_RR + iz - rr_va_z):
  //     mean = lv_X_i * b_z_hat[:,k]   if k < d_c, else 0
  //   iz >= rr_va_z+num_lvc (lv dims): mean = 0
  //
  // dev2 = ||u_i - mu_z||^2 + bz_cross  (cross-covariance correction for lvc dims)
  for (int i = 0; i < n; i++) {
    vector<Type> mu_z(d_va_z); mu_z.setZero();
    if (Kz > 0) {
      // rr_va_z == 0; VA-z[0..num_lvc-1] are lvc dims
      int dc_lvc = d_c - std::min(d_c, num_RR); // active cols within lvc block
      for (int iz = 0; iz < dc_lvc; iz++) {
        int k = num_RR + iz;
        for (int l = 0; l < Kz; l++) mu_z(iz) += lv_X_env(i,l) * b_z(l,k);
      }
    }

    Type logdet_half = Type(0);
    Type trace_Ai   = Type(0);
    if (va_struct == 0) {
      for (int iz = 0; iz < d_va_z; iz++) {
        logdet_half += Au(iz*n+i);
        trace_Ai    += Ai_diag(i,iz);
      }
    } else {
      for (int r = 0; r < d_va_z; r++) logdet_half += Au((r*(r+1)/2+r)*n+i);
      for (int r = 0; r < d_va_z; r++) {
        for (int c = 0; c <= r; c++) {
          int idx = r*(r+1)/2+c;
          Type lval = Au(idx*n+i);
          if (c == r) lval = exp(lval);
          trace_Ai += lval*lval;
        }
      }
    }
    Type dev2 = Type(0);
    for (int iz = 0; iz < d_va_z; iz++) {
      Type dv = u(i,iz) - mu_z(iz);
      dev2 += dv*dv;
    }
    // Cross-covariance correction: accounts for q(b_z) uncertainty in the lvc prior mean.
    // E_{q(b_z)}[||z_i - lv_X_i*b_z||^2] = dev2_above + lv_X_i^T * bz_var_mat * lv_X_i
    // Without this term ELBO in sigma_bz is -0.5*||b_z||^2/sigma_bz^2 (monotone in sigma_bz).
    Type bz_cross = Type(0);
    if (has_bz && dc_lvc > 0) {
      matrix<Type> xi(Kz, 1);
      for (int q = 0; q < Kz; q++) xi(q, 0) = lv_X_env(i, q);
      bz_cross = (xi.transpose() * bz_var_mat * xi)(0, 0);
    }
    nll -= logdet_half - Type(0.5)*(trace_Ai + dev2 + bz_cross);
  }
  nll -= Type(0.5)*n*d_va_z;

  // ===== KL(q(gamma_j) || N(mu_g_j, I)) =====
  for (int j = 0; j < p; j++) {
    vector<Type> mu_g(d_va_a); mu_g.setZero();
    if (Kt > 0) {
      int dt_lvc = d_t - std::min(d_t, num_RR);
      for (int ia = 0; ia < dt_lvc; ia++) {
        int k = num_RR + ia;
        for (int l = 0; l < Kt; l++) mu_g(ia) += TR(j,l) * b_gamma(l,k);
      }
    }

    Type logdet_half = Type(0);
    Type trace_Aj   = Type(0);
    if (va_struct == 0) {
      for (int ia = 0; ia < d_va_a; ia++) {
        logdet_half += Au_sp(ia*p+j);
        trace_Aj    += Aj_diag(j,ia);
      }
    } else {
      for (int r = 0; r < d_va_a; r++) logdet_half += Au_sp((r*(r+1)/2+r)*p+j);
      for (int r = 0; r < d_va_a; r++) {
        for (int c = 0; c <= r; c++) {
          int idx = r*(r+1)/2+c;
          Type lval = Au_sp(idx*p+j);
          if (c == r) lval = exp(lval);
          trace_Aj += lval*lval;
        }
      }
    }
    Type dev2 = Type(0);
    for (int ia = 0; ia < d_va_a; ia++) {
      Type dv = a_lv_sp_id(j,ia) - mu_g(ia);
      dev2 += dv*dv;
    }
    // Cross-covariance correction: analogous to bz_cross for gamma_j with traits.
    // E_{q(b_gamma)}[||gamma_j - TR_j*b_gamma||^2] = dev2 + TR_j^T * bg_var_mat * TR_j
    Type bg_cross = Type(0);
    if (has_bg && dt_lvc > 0) {
      matrix<Type> xj(Kt, 1);
      for (int q = 0; q < Kt; q++) xj(q, 0) = TR(j, q);
      bg_cross = (xj.transpose() * bg_var_mat * xj)(0, 0);
    }
    nll -= logdet_half - Type(0.5)*(trace_Aj + dev2 + bg_cross);
  }
  nll -= Type(0.5)*p*d_va_a;

  // ===== KL FOR RANDOM b_z (randomB == 1, "LV") =====
  // AB_z already built above; used here for the prior KL on b_z[:,l].
  // KL(q(b_z[:,l]) || N(0, sigma_bz_l^2 Sigma_corr)) per active column l = 0..d_c-1.
  // sigma_bz for RR dims (l < num_RR) is mapped out (fixed = 1) in R; only lvc dims are free.
  if (has_bz) {
    if (csb_z.cols() < 2) {
      for (int l = 0; l < d_c; l++) {
        Type lsig = log_sigma_bz(l), sig2inv = exp(Type(-2)*lsig);
        nll -= AB_z(l).diagonal().array().log().sum()
             - Type(0.5)*sig2inv*AB_z(l).rowwise().squaredNorm().sum()
             - Type(0.5)*sig2inv*b_z.col(l).squaredNorm()
             + Type(0.5)*(Type(Kz) - Type(Kz)*Type(2)*lsig);
      }
    } else {
      int n_pairs_z = (Kz*(Kz-1))/2;
      vector<Type> corsb_z(n_pairs_z); corsb_z.setZero();
      for (int i = 0; i < csb_z.rows(); i++) {
        int idx = (csb_z(i,0)-1)*(csb_z(i,0)-2)/2 + csb_z(i,1)-1;
        corsb_z(idx) = log_sigma_bz(d_c + i);
      }
      matrix<Type> Sigmab_z_L  = gllvmutils::constructL(corsb_z);
      matrix<Type> Ikz         = matrix<Type>::Identity(Kz, Kz);
      matrix<Type> Sigmab_z_LI = Sigmab_z_L.template triangularView<Eigen::Lower>().solve(Ikz);
      matrix<Type> Sigmab_z_I  = Sigmab_z_LI.transpose() * Sigmab_z_LI;
      Type logdet_corr = Sigmab_z_L.diagonal().array().log().sum();
      for (int l = 0; l < d_c; l++) {
        Type lsig = log_sigma_bz(l), sig2inv = exp(Type(-2)*lsig);
        matrix<Type> SigI = sig2inv * Sigmab_z_I;
        nll -= AB_z(l).diagonal().array().log().sum()
             - Type(0.5)*(SigI*AB_z(l)*AB_z(l).transpose()).trace()
             - Type(0.5)*(b_z.col(l).transpose()*SigI*b_z.col(l)).value()
             + Type(0.5)*Type(Kz) - Type(Kz)*lsig - logdet_corr;
      }
    }
  }

  // ===== KL FOR RANDOM b_gamma (randomT == 1, "LV") =====
  // AB_g already built above; used here for the prior KL on b_gamma[:,l].
  if (has_bg) {
    if (csb_gamma.cols() < 2) {
      for (int l = 0; l < d_t; l++) {
        Type lsig = log_sigma_bgamma(l), sig2inv = exp(Type(-2)*lsig);
        nll -= AB_g(l).diagonal().array().log().sum()
             - Type(0.5)*sig2inv*AB_g(l).rowwise().squaredNorm().sum()
             - Type(0.5)*sig2inv*b_gamma.col(l).squaredNorm()
             + Type(0.5)*(Type(Kt) - Type(Kt)*Type(2)*lsig);
      }
    } else {
      int n_pairs_t = (Kt*(Kt-1))/2;
      vector<Type> corsb_g(n_pairs_t); corsb_g.setZero();
      for (int i = 0; i < csb_gamma.rows(); i++) {
        int idx = (csb_gamma(i,0)-1)*(csb_gamma(i,0)-2)/2 + csb_gamma(i,1)-1;
        corsb_g(idx) = log_sigma_bgamma(d_t + i);
      }
      matrix<Type> Sigmab_g_L  = gllvmutils::constructL(corsb_g);
      matrix<Type> Ikt         = matrix<Type>::Identity(Kt, Kt);
      matrix<Type> Sigmab_g_LI = Sigmab_g_L.template triangularView<Eigen::Lower>().solve(Ikt);
      matrix<Type> Sigmab_g_I  = Sigmab_g_LI.transpose() * Sigmab_g_LI;
      Type logdet_corr_g = Sigmab_g_L.diagonal().array().log().sum();
      for (int l = 0; l < d_t; l++) {
        Type lsig = log_sigma_bgamma(l), sig2inv = exp(Type(-2)*lsig);
        matrix<Type> SigI = sig2inv * Sigmab_g_I;
        nll -= AB_g(l).diagonal().array().log().sum()
             - Type(0.5)*(SigI*AB_g(l)*AB_g(l).transpose()).trace()
             - Type(0.5)*(b_gamma.col(l).transpose()*SigI*b_gamma.col(l)).value()
             + Type(0.5)*Type(Kt) - Type(Kt)*lsig - logdet_corr_g;
      }
    }
  }

  // ===== ETA and CQ =====
  matrix<Type> eta(n,p); eta.setZero();
  matrix<Type> cQ(n,p);  cQ.setZero();
  matrix<Type> mu(n,p);  mu.setZero();

  eta += x * b;
  if (offset.rows() == n) eta += offset;
  if (xr.cols() > 0 && r0f.rows() > 0) eta += (xr * r0f).replicate(1, p);
  int spnr = xb.cols();

  // -------------------------------------------------------------------
  // Ordination contribution to eta and cQ: E[Z · Sigma · Gamma]
  //
  // We always integrate the bilinear form z_ik * sigma_k * gamma_jk over the
  // joint VA posterior q(z_i) q(gamma_j).  The VA means and variances are:
  //
  //   RR (k < num_RR):
  //     Kz > 0, randomB=0: ai = lv_X_i*b_z[:,k] (det)   ui = 0
  //     Kz > 0, randomB=1: ai = lv_X_i*b_z_hat[:,k]     ui = x_i^T·AB_z(k)·AB_z(k)^T·x_i
  //     Kz = 0:            ai = u_i[iz]          (VA)    ui = A_i[iz]
  //     same logic for aj/vj with Kt/randomT
  //
  //   lvc (k in [num_RR, num_RR+num_lvc)):
  //     ai = u_i[iz]  (VA full z mean)   ui = A_i[iz] + x_i^T·AB_z(k)·AB_z(k)^T·x_i  (if randomB=1, k<d_c)
  //     aj = a_j[ia]  (VA full γ mean)   vj = B_j[ia] + t_j^T·AB_g(k)·AB_g(k)^T·t_j  (if randomT=1, k<d_t)
  //
  //   lv (k >= num_RR+num_lvc):
  //     ai = u_i[iz]  (VA, prior mean = 0)              ui = A_i[iz]
  //     aj = a_j[ia]  (VA, prior mean = 0)              vj = B_j[ia]
  //
  // eta(i,j) += sum_k sigma_k * ai_k * aj_k             (VA mean of bilinear form)
  // cQ(i,j):
  //   log-link:     sum_k cQ_k  (5-term MGF formula, ms.pdf eq 5, T1-T5 in loop below)
  //   non-log-link: 0.5 * Var_q(z_i^T Sigma gamma_j)
  //               = 0.5 * sum_k sigma_k^2 * (ui_k*vj_k + ui_k*aj_k^2 + vj_k*ai_k^2)
  // -------------------------------------------------------------------

  for (int i = 0; i < n; i++) {

    // Precompute deterministic z contributions for RR dims (only if Kz > 0)
    std::vector<Type> z_rr_det(num_RR, Type(0));
    if (Kz > 0) {
      for (int k = 0; k < num_RR; k++)
        for (int l = 0; l < Kz; l++) z_rr_det[k] += lv_X_env(i,l) * b_z(l,k);
    }

    // Precompute x_i as column vector for b_z variance calculations.
    matrix<Type> xi_bz(has_bz ? Kz : 1, 1); xi_bz.setZero();
    if (has_bz) for (int q = 0; q < Kz; q++) xi_bz(q, 0) = lv_X_env(i, q);

    // ---- Per-site b_z variance contributions ----
    // For RR dims: z_ik = x_i^T b_z[:,k] (no VA residual).  Var_q(z_ik) = x_i^T Var_q(b_z[:,k]) x_i.
    // For lvc dims: q(z_i)=N(u_i,A_eps_i); u_i is the full VA mean (= b_z_hat^T x_i + residual mean).
    //   Var_q(z_ik) = A_eps_i[k] + x_i^T Var_q(b_z[:,k]) x_i  (ms.Rmd eq lvcvar, randomB=1)
    //   Var_q(z_ik) = A_eps_i[k]                               (randomB=0; TODO Stiefel manifold)
    // Var_q(b_z[:,k]) = AB_z(k) AB_z(k)^T.
    int n_rr_bz = has_bz ? std::min(num_RR, d_c) : 0;
    std::vector<Type> rr_z_var_i(num_RR, Type(0));
    for (int k = 0; k < n_rr_bz; k++) {
      matrix<Type> Bk = AB_z(k) * AB_z(k).transpose();
      rr_z_var_i[k] = (xi_bz.transpose() * Bk * xi_bz)(0, 0);
    }
    std::vector<Type> lvc_z_var_i(dc_lvc, Type(0));
    if (has_bz) {
      for (int iz_lvc = 0; iz_lvc < dc_lvc; iz_lvc++) {
        int k = num_RR + iz_lvc;
        matrix<Type> Bk = AB_z(k) * AB_z(k).transpose();
        lvc_z_var_i[iz_lvc] = (xi_bz.transpose() * Bk * xi_bz)(0, 0);
      }
    }

    // ---- Full marginal moments of z_i: a_z[k]=E_q[z_ik], u_z[k]=Var_q[z_ik] ----
    // (ms.Rmd §"Posterior moments with covariates and traits", eqs lvcmean/lvcvar)
    vector<Type> a_z(d); a_z.setZero();
    vector<Type> u_z(d); u_z.setZero();
    for (int k = 0; k < d; k++) {
      int iz = va_z_idx[k];
      if (iz >= 0) {
        // lvc or lv (or RR when Kz=0): u(i,iz) is the full VA mean of z_ik;
        // Ai_diag(i,iz) is the residual (epsilon) variance.
        a_z(k) = u(i, iz);
        u_z(k) = Ai_diag(i, iz);
        // lvc dim k with active b_z, randomB=1: add x_i^T Var_q(b_z[:,k]) x_i
        if (has_bz && k >= num_RR && k < num_RR + dc_lvc)
          u_z(k) += lvc_z_var_i[k - num_RR];
      } else {
        // RR dim (Kz>0): no VA residual; z_ik = x_i^T b_z[:,k]
        a_z(k) = z_rr_det[k];
        u_z(k) = (k < n_rr_bz) ? rr_z_var_i[k] : Type(0);
      }
    }

    // ---- Build augmented Ai_full = Var_q(z_i) (d×d) and Mi = Σ Ai_full Σ ----
    // Diagonal VA (va_struct=0): Ai_full = diag(u_z) — u_z already holds the full augmented
    //   marginal diagonal (residual VA + b_z contribution), so no further work needed.
    // Unstructured VA (va_struct=1): build from the Cholesky Li_va, then augment at RR and
    //   lvc positions with rr_z_var_i and lvc_z_var_i (ms.Rmd eq lvcvar).
    // The two cases share Mi and the j-loop cQ formula: diagonal is just a special case of
    // the unstructured formula where Ai_full and Aj happen to be diagonal.
    matrix<Type> Ai_full(d,d); Ai_full.setZero();
    if (va_struct == 0) {
      for (int k = 0; k < d; k++) Ai_full(k, k) = u_z(k);
    } else {
      matrix<Type> Li_va(d_va_z, d_va_z); Li_va.setZero();
      for (int r = 0; r < d_va_z; r++) {
        for (int c = 0; c <= r; c++) {
          int idx = r*(r+1)/2+c;
          Type lval = Au(idx*n+i);
          Li_va(r,c) = (c==r) ? exp(lval) : lval;
        }
      }
      matrix<Type> Ai_va = Li_va * Li_va.transpose();
      for (int r = 0; r < d_va_z; r++) for (int c = 0; c < d_va_z; c++)
        Ai_full(va_z_to_full[r], va_z_to_full[c]) = Ai_va(r,c);
      for (int k = 0; k < n_rr_bz; k++)
        Ai_full(k, k) += rr_z_var_i[k];
      for (int iz_lvc = 0; iz_lvc < dc_lvc; iz_lvc++)
        if (has_bz) Ai_full(num_RR + iz_lvc, num_RR + iz_lvc) += lvc_z_var_i[iz_lvc];
    }
    matrix<Type> Mi(d,d);
    for (int r = 0; r < d; r++) for (int c = 0; c < d; c++)
      Mi(r,c) = sigma(r) * Ai_full(r,c) * sigma(c);

    for (int j = 0; j < p; j++) {

      // Precompute deterministic gamma for RR dims (only if Kt > 0)
      std::vector<Type> g_rr_det(num_RR, Type(0));
      if (Kt > 0) {
        for (int k = 0; k < num_RR; k++)
          for (int l = 0; l < Kt; l++) g_rr_det[k] += TR(j,l) * b_gamma(l,k);
      }

      // Precompute t_j as column vector for b_gamma variance calculations.
      matrix<Type> xj_bg(has_bg ? Kt : 1, 1); xj_bg.setZero();
      if (has_bg) for (int q = 0; q < Kt; q++) xj_bg(q, 0) = TR(j, q);

      // RR b_gamma variance: Var_q(gamma_jk) = t_j^T · AB_g(k)·AB_g(k)^T · t_j
      int n_rr_bg = has_bg ? std::min(num_RR, d_t) : 0;
      std::vector<Type> rr_g_var_j(num_RR, Type(0));
      for (int k = 0; k < n_rr_bg; k++) {
        matrix<Type> Bk = AB_g(k) * AB_g(k).transpose();
        rr_g_var_j[k] = (xj_bg.transpose() * Bk * xj_bg)(0, 0);
      }

      // ---- Per-species b_gamma variance contributions ----
      // Var_q(gamma_jk) for lvc dim k = Aj_diag_eps(j,ia) [residual] + t_j^T Var_q(b_gamma[:,k]) t_j
      // Var_q(b_gamma[:,k]) = AB_g(k) AB_g(k)^T.  Only the lvc-active block (k=num_RR..num_RR+dt_lvc-1).
      std::vector<Type> lvc_g_var_j(dt_lvc, Type(0));
      if (has_bg) {
        for (int ia_lvc = 0; ia_lvc < dt_lvc; ia_lvc++) {
          int k = num_RR + ia_lvc;
          matrix<Type> Bk = AB_g(k) * AB_g(k).transpose();
          lvc_g_var_j[ia_lvc] = (xj_bg.transpose() * Bk * xj_bg)(0, 0);
        }
      }

      // ---- Full marginal moments of gamma_j: a_g[k]=E_q[gamma_jk], v_g[k]=Var_q[gamma_jk] ----
      // (ms.Rmd §"Posterior moments with covariates and traits", eqs lvcmean/lvcvar, species side)
      vector<Type> a_g(d); a_g.setZero();
      vector<Type> v_g(d); v_g.setZero();
      for (int k = 0; k < d; k++) {
        int ia = va_a_idx[k];
        if (ia >= 0) {
          a_g(k) = a_lv_sp_id(j, ia);
          v_g(k) = Aj_diag(j, ia);
          // lvc dim k with active b_gamma, randomT=1: add t_j^T Var_q(b_gamma[:,k]) t_j
          if (has_bg && k >= num_RR && k < num_RR + dt_lvc)
            v_g(k) += lvc_g_var_j[k - num_RR];
        } else {
          // RR dim (Kt>0): no VA residual; gamma_jk = t_j^T b_gamma[:,k]
          a_g(k) = g_rr_det[k];
          v_g(k) = (k < n_rr_bg) ? rr_g_var_j[k] : Type(0);
        }
      }

      // ---- eta mean: E_q[z_i^T Sigma gamma_j] = a_z^T Sigma a_g  (ms.Rmd eq etamean) ----
      for (int k = 0; k < d; k++) eta(i,j) += sigma(k) * a_z(k) * a_g(k);

      int fam = family(j);
      bool log_link = (fam == POISSON || fam == NEG_BINOMIAL || fam == GAMMA ||
                       fam == TWEEDIE || fam == ZIP || fam == EXPONENTIAL ||
                       fam == ZINB    || fam == ZNIB);

      Type cq = Type(0);

      // ---- Build augmented Aj = Var_q(gamma_j) (d×d) and Mj = Σ Aj Σ ----
      // Diagonal VA (va_struct=0): Aj = diag(v_g) — v_g already has the full augmented variance.
      // Unstructured VA (va_struct=1): build from Lj_va, then augment at RR and lvc positions.
      matrix<Type> Aj(d,d); Aj.setZero();
      if (va_struct == 0) {
        for (int k = 0; k < d; k++) Aj(k, k) = v_g(k);
      } else {
        matrix<Type> Lj_va(d_va_a, d_va_a); Lj_va.setZero();
        for (int r = 0; r < d_va_a; r++) {
          for (int c = 0; c <= r; c++) {
            int idx = r*(r+1)/2+c;
            Type lval = Au_sp(idx*p+j);
            Lj_va(r,c) = (c==r) ? exp(lval) : lval;
          }
        }
        matrix<Type> Aj_va = Lj_va * Lj_va.transpose();
        for (int r = 0; r < d_va_a; r++) for (int c = 0; c < d_va_a; c++)
          Aj(va_a_to_full[r], va_a_to_full[c]) = Aj_va(r,c);
        for (int k = 0; k < n_rr_bg; k++)
          Aj(k, k) += rr_g_var_j[k];
        for (int ia_lvc = 0; ia_lvc < dt_lvc; ia_lvc++)
          if (has_bg) Aj(num_RR + ia_lvc, num_RR + ia_lvc) += lvc_g_var_j[ia_lvc];
      }

      // Mj = Σ Aj Σ;  va_m = a_z;  nu_m = a_g
      matrix<Type> Mj(d,d);
      for (int r = 0; r < d; r++) for (int c = 0; c < d; c++)
        Mj(r,c) = sigma(r) * Aj(r,c) * sigma(c);
      matrix<Type> va_m(d,1); for (int k=0; k<d; k++) va_m(k,0) = a_z(k);
      matrix<Type> nu_m(d,1); for (int k=0; k<d; k++) nu_m(k,0) = a_g(k);

      // ---- cQ: unified bilinear form formula for log-link and non-log-link ----
      // Uses Ai_full and Aj (both augmented with full marginal variances).
      // Diagonal VA is a special case (Ai_full = diag(u_z), Aj = diag(v_g)); the matrix
      // formula gives the same result as the per-dim T1-T5 scalar loop in that case.
      // (ms.Rmd eq solution for log-link; eq varbilinear for non-log-link)
      if (log_link) {
        // C_ij = I − Σ A_j Σ A_i;  cQ = -0.5 log|C| + 0.5 a_i^T A_i (...)
        // The integral converges iff C is positive definite (det(C) > 0).
        // Guard: clamp det(C) away from zero/negative so log() and C^{-1} never receive
        // degenerate input.  When det(C) <= 0 the cQ penalty grows large, driving the
        // optimizer back into the feasible region.
        matrix<Type> C = -(Mj * Ai_full);
        for (int k = 0; k < d; k++) C(k,k) += Type(1);
        Type detC      = C.determinant();
        Type detC_safe = CppAD::CondExpGt(detC, Type(1e-8), detC, Type(1e-8));
        Type logdetC   = log(detC_safe);
        // Regularize C diagonally when near-singular so inverse is always finite.
        matrix<Type> C_reg = C;
        Type reg = CppAD::CondExpLt(detC, Type(1e-4), Type(1e-4) - detC + Type(1e-4), Type(0));
        for (int k = 0; k < d; k++) C_reg(k,k) += reg;
        matrix<Type> Cinv = C_reg.inverse();

        matrix<Type> nu_sigma(d,1);
        for (int k = 0; k < d; k++) nu_sigma(k,0) = sigma(k) * nu_m(k,0);

        matrix<Type> Mj_va      = Mj * va_m;           // Σ A_j Σ a_z
        matrix<Type> Cinv_nu    = Cinv * nu_sigma;      // C^{-1} Σ a_g
        matrix<Type> Cinv_Mj_va = Cinv * Mj_va;        // C^{-1} Σ A_j Σ a_z
        matrix<Type> Ai_nu_sig  = Ai_full * nu_sigma;   // A_i Σ a_g

        Type t2_t4 = Type(0.5) * (va_m.transpose()     * Cinv_Mj_va)(0,0);  // 0.5 a_z^T C^{-1} Mj a_z
        Type t1_t5 = Type(0.5) * (Ai_nu_sig.transpose() * Cinv_nu   )(0,0);  // 0.5 (A_i Σ a_g)^T C^{-1} Σ a_g
        Type t3    =              (Ai_nu_sig.transpose() * Cinv_Mj_va)(0,0);  // (A_i Σ a_g)^T C^{-1} Mj a_z

        cq = -Type(0.5)*logdetC + t1_t5 + t2_t4 + t3;
      } else {
        // cQ = 0.5 Var_q(z_i^T Σ gamma_j) = 0.5 [ tr(Mi Aj) + a_g^T Mi a_g + a_z^T Mj a_z ]
        // (ms.Rmd eq varbilinear)
        Type trMAj = Type(0);
        for (int r = 0; r < d; r++) for (int c = 0; c < d; c++) trMAj += Mi(r,c)*Aj(c,r);
        matrix<Type> Mnu_m = Mi * nu_m;
        Type nu_Mnu   = (nu_m.transpose() * Mnu_m)(0,0);
        Type va_Mj_va = (va_m.transpose() * Mj * va_m)(0,0);
        cq = Type(0.5) * (trMAj + nu_Mnu + va_Mj_va);
      }

      cQ(i,j) = cq;
    }
  }

  // ===== SPECIES-SPECIFIC RANDOM EFFECTS (shared header) =====
  #include "species_effects.h"

  // ===== RANDOM ROW EFFECTS (shared header) =====
  #include "row_effects.h"

  // ===== VA FAMILY LIKELIHOODS =====
  if ((method < 1) || (method > 1)) {
    #include "family_va_nll.h"
  }

  return nll;
}
