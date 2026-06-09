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
// LV blocks (ordered):  [RR | lvc | lv]
//   k = 0..num_RR-1     : RR dims
//   k = num_RR..num_RR+num_lvc-1 : concurrent (lvc) dims
//   k = num_RR+num_lvc..d-1      : unconstrained (lv) dims
//
// VA parameters:
//   z_i VA  (d_va_z dims):  lvc + lv dims always VA;
//                            RR dims VA iff Kz==0 (no X)
//   gamma_j VA (d_va_a dims): lvc + lv always VA;
//                              RR dims VA iff Kt==0 (no TR)
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
  DATA_MATRIX(xrr);          // random row-effect design matrix
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
  DATA_INTEGER(random);      // bit 0 = random row effects
  DATA_INTEGER(va_struct);   // 0 = diagonal VA covariance; 1 = unstructured (full Cholesky)
  DATA_MATRIX(lv_X_env);     // n x Kz  (0-column matrix when not used)
  DATA_MATRIX(TR);           // p x Kt  (0-column matrix when not used)
  DATA_INTEGER(randomB);     // 1 = "LV" prior on b_z
  DATA_INTEGER(randomT);     // 1 = "LV" prior on b_gamma
  DATA_IMATRIX(csb_z);       // (n_pairs x 2) 1-based predictor-pair indices for b_z
  DATA_IMATRIX(csb_gamma);   // same for b_gamma

  // ===== PARAMETERS =====
  PARAMETER_MATRIX(b);       // Kx x p fixed effects

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

  // Ordination scale: sigma(k) = sum_{l=k}^{d-1} exp(sigmaLV(l))
  PARAMETER_VECTOR(sigmaLV); // length d

  PARAMETER_VECTOR(lg_phi);
  PARAMETER_VECTOR(lg_phiZINB);
  PARAMETER_VECTOR(zeta);
  PARAMETER(ePower);

  PARAMETER_MATRIX(r0f);
  PARAMETER_MATRIX(r0r);
  PARAMETER_VECTOR(lg_Ar);
  PARAMETER_VECTOR(log_sigma);

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

  vector<Type> iphi = exp(lg_phi);

  parallel_accumulator<Type> nll(this);

  // ===== SIGMA construction =====
  vector<Type> sigma(d);
  sigma(d-1) = exp(sigmaLV(d-1));
  for (int k = d-2; k >= 0; k--) sigma(k) = sigma(k+1) + exp(sigmaLV(k));
  vector<Type> sigma2(d);
  for (int k = 0; k < d; k++) sigma2(k) = sigma(k) * sigma(k);

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

  // ===== KL(q(z_i) || N(mu_z_i, I)) =====
  // Only over d_va_z dims (det-z dims have no KL contribution).
  // Prior mean mu_z[iz]:
  //   iz < rr_va_z (RR-TR-only, Kz==0): mean = 0
  //   rr_va_z <= iz < rr_va_z+num_lvc (lvc dims, k=num_RR+iz-rr_va_z):
  //     mean = lv_X_env_i * b_z[:,k] if k < d_c, else 0
  //   iz >= rr_va_z+num_lvc (lv dims): mean = 0
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
    nll -= logdet_half - Type(0.5)*(trace_Ai + dev2);
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
      Type dv = a_lv_sp(j,ia) - mu_g(ia);
      dev2 += dv*dv;
    }
    nll -= logdet_half - Type(0.5)*(trace_Aj + dev2);
  }
  nll -= Type(0.5)*p*d_va_a;

  // ===== KL FOR RANDOM b_z (randomB == 1, "LV") =====
  if (randomB == 1 && Kz > 0 && d_c > 0) {
    vector<matrix<Type>> AB_z(d_c);
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
  if (randomT == 1 && Kt > 0 && d_t > 0) {
    vector<matrix<Type>> AB_g(d_t);
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

  // -------------------------------------------------------------------
  // Ordination eta mean + cQ
  //
  // For each full dim k, the z_ik mean ("ai") and variance ("ui") and
  // gamma_jk mean ("aj") and variance ("vj") depend on which block k is in:
  //
  //   RR block (k < num_RR):
  //     z_ik: det  (ai = X*b_z[:,k], ui=0) if Kz>0; VA (ai=u[iz], ui=A_i[iz]) if Kz==0
  //     gamma_jk: det (aj = TR*b_g[:,k], vj=0) if Kt>0; VA (aj=a[ia], vj=A_j[ia]) if Kt==0
  //   lvc block (k in [num_RR, num_RR+num_lvc)):
  //     always VA for both z and gamma
  //   lv block (k >= num_RR+num_lvc):
  //     always VA for both z and gamma
  // -------------------------------------------------------------------

  for (int i = 0; i < n; i++) {

    // Precompute deterministic z contributions for RR dims (only if Kz > 0)
    std::vector<Type> z_rr_det(num_RR, Type(0));
    if (Kz > 0) {
      for (int k = 0; k < num_RR; k++)
        for (int l = 0; l < Kz; l++) z_rr_det[k] += lv_X_env(i,l) * b_z(l,k);
    }

    // Unstructured: build full d×d A_i (padded with zeros at det-z dims) and Mi
    matrix<Type> Ai_full(d,d); Ai_full.setZero();
    matrix<Type> Mi(d,d);      Mi.setZero();
    if (va_struct == 1) {
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
      for (int r = 0; r < d; r++) for (int c = 0; c < d; c++)
        Mi(r,c) = sigma(r) * Ai_full(r,c) * sigma(c);
    }

    for (int j = 0; j < p; j++) {

      // Precompute deterministic gamma for RR dims (only if Kt > 0)
      std::vector<Type> g_rr_det(num_RR, Type(0));
      if (Kt > 0) {
        for (int k = 0; k < num_RR; k++)
          for (int l = 0; l < Kt; l++) g_rr_det[k] += TR(j,l) * b_gamma(l,k);
      }

      // --- eta mean ---
      Type cross = Type(0);
      for (int k = 0; k < d; k++) {
        int iz = va_z_idx[k], ia = va_a_idx[k];
        Type ai_k = (iz >= 0) ? u(i,iz) : z_rr_det[k];
        Type aj_k = (ia >= 0) ? a_lv_sp(j,ia) : g_rr_det[k];
        cross += sigma(k) * ai_k * aj_k;
      }
      eta(i,j) += cross;

      int fam = family(j);
      bool log_link = (fam == POISSON || fam == NEG_BINOMIAL || fam == GAMMA ||
                       fam == TWEEDIE || fam == ZIP || fam == EXPONENTIAL ||
                       fam == ZINB    || fam == ZNIB);

      Type cq = Type(0);

      if (va_struct == 0) {
        // ---- diagonal VA ----
        // For dims with deterministic z (iz==-1): ui_k=0, ck=1 → no correction from those dims.
        // For dims with deterministic gamma (ia==-1): vj_k=0, ck=1 → same.
        if (log_link) {
          for (int k = 0; k < d; k++) {
            int iz = va_z_idx[k], ia = va_a_idx[k];
            Type sk  = sigma(k), sk2 = sigma2(k);
            Type ui_k = (iz >= 0) ? Ai_diag(i,iz) : Type(0);
            Type vj_k = (ia >= 0) ? Aj_diag(j,ia) : Type(0);
            Type ai_k = (iz >= 0) ? u(i,iz) : z_rr_det[k];
            Type aj_k = (ia >= 0) ? a_lv_sp(j,ia) : g_rr_det[k];
            Type ck  = Type(1) - sk2 * ui_k * vj_k;
            Type ck_safe = CppAD::CondExpGt(ck, Type(1e-6), ck, Type(1e-6));
            cq -= Type(0.5) * log(ck_safe);
            cq += Type(0.5) * sk2 * ui_k * aj_k * aj_k;
            cq += sk * sk2 * ui_k * vj_k * ai_k * aj_k / ck_safe;
            cq += Type(0.5) * sk * sk2 * vj_k * vj_k * ai_k * ai_k / ck_safe;
            cq += Type(0.5) * sk * sk2 * ui_k * ui_k * aj_k * aj_k / ck_safe;
          }
        } else {
          for (int k = 0; k < d; k++) {
            int iz = va_z_idx[k], ia = va_a_idx[k];
            Type sk2 = sigma2(k);
            Type ui_k = (iz >= 0) ? Ai_diag(i,iz) : Type(0);
            Type vj_k = (ia >= 0) ? Aj_diag(j,ia) : Type(0);
            Type ai_k = (iz >= 0) ? u(i,iz) : z_rr_det[k];
            Type aj_k = (ia >= 0) ? a_lv_sp(j,ia) : g_rr_det[k];
            cq += sk2 * (ui_k*vj_k + ui_k*aj_k*aj_k + vj_k*ai_k*ai_k);
          }
          cq *= Type(0.5);
        }
      } else {
        // ---- unstructured VA ----
        // Build full d×d A_j padded with zeros at det-a dims
        matrix<Type> Lj_va(d_va_a, d_va_a); Lj_va.setZero();
        for (int r = 0; r < d_va_a; r++) {
          for (int c = 0; c <= r; c++) {
            int idx = r*(r+1)/2+c;
            Type lval = Au_sp(idx*p+j);
            Lj_va(r,c) = (c==r) ? exp(lval) : lval;
          }
        }
        matrix<Type> Aj_va = Lj_va * Lj_va.transpose();
        matrix<Type> Aj(d,d); Aj.setZero();
        for (int r = 0; r < d_va_a; r++) for (int c = 0; c < d_va_a; c++)
          Aj(va_a_to_full[r], va_a_to_full[c]) = Aj_va(r,c);

        // Full d-vectors for means
        matrix<Type> va_m(d,1); va_m.setZero();
        matrix<Type> nu_m(d,1); nu_m.setZero();
        for (int k = 0; k < d; k++) {
          int iz = va_z_idx[k], ia = va_a_idx[k];
          va_m(k,0) = (iz >= 0) ? u(i,iz)         : z_rr_det[k];
          nu_m(k,0) = (ia >= 0) ? a_lv_sp(j,ia)   : g_rr_det[k];
        }
        matrix<Type> v_m(d,1);
        for (int k = 0; k < d; k++) v_m(k,0) = sigma(k) * va_m(k,0);
        matrix<Type> Mnu_m = Mi * nu_m;

        if (log_link) {
          matrix<Type> C = -(Mi * Aj);
          for (int k = 0; k < d; k++) C(k,k) += Type(1);
          Type logdetC = log(C.determinant() + Type(1e-15));
          matrix<Type> rhs_m  = v_m + Mnu_m;
          matrix<Type> w_m    = C.inverse() * (Aj * rhs_m);
          Type nu_Mnu = (nu_m.transpose() * Mnu_m)(0,0);
          Type rhs_w  = (rhs_m.transpose() * w_m)(0,0);
          cq = -Type(0.5)*logdetC + Type(0.5)*nu_Mnu + Type(0.5)*rhs_w;
        } else {
          matrix<Type> Mj(d,d);
          for (int r = 0; r < d; r++) for (int c = 0; c < d; c++)
            Mj(r,c) = sigma(r) * Aj(r,c) * sigma(c);
          Type trMAj = Type(0);
          for (int r = 0; r < d; r++) for (int c = 0; c < d; c++) trMAj += Mi(r,c)*Aj(c,r);
          Type nu_Mnu   = (nu_m.transpose() * Mnu_m)(0,0);
          Type va_Mj_va = (va_m.transpose() * Mj * va_m)(0,0);
          cq = Type(0.5) * (trMAj + nu_Mnu + va_Mj_va);
        }
      }

      cQ(i,j) = cq;
    }
  }

  // ===== RANDOM ROW EFFECTS =====
  if ((random & 1) > 0 && xrr.cols() > 0 && r0r.rows() > 0) {
    int G = r0r.rows();
    Type sigma_r  = exp(log_sigma(0));
    Type sigma_r2 = sigma_r * sigma_r;
    matrix<Type> r_contrib = xrr * r0r;
    for (int i = 0; i < n; i++) eta.row(i).array() += r_contrib(i, 0);
    for (int g = 0; g < G; g++) {
      Type Au_r = lg_Ar(g), Ar = exp(Type(2)*Au_r);
      nll -= Au_r - Type(0.5)*(Ar/sigma_r2 + r0r(g,0)*r0r(g,0)/sigma_r2);
    }
    nll -= Type(0.5)*G*(Type(1) - Type(2)*log_sigma(0));
  }

  // ===== VA FAMILY LIKELIHOODS =====
  if ((method < 1) || (method > 1)) {
    #include "family_va_nll.h"
  }

  return nll;
}
