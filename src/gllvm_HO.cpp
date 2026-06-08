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
// Model: eta_ij = beta0_j + x_i^T beta_j + z_i^T Sigma gamma_j + row_i
//   z_i ~ N(0, I_d),  gamma_j ~ N(0, I_d)
//   q(z_i)   = N(a_i,   A_i)  — site VA mean/cov (diagonal A_i)
//   q(gamma_j) = N(a_j, A_j)  — species VA mean/cov (diagonal A_j)
//   Sigma = diag(sigma),  sigma ordered via log-sigmoid parameterisation
//-----------------------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  // ===== DATA =====
  DATA_MATRIX(y);            // n x p responses
  DATA_MATRIX(x);            // n x Kx fixed covariates (first col = intercept)
  DATA_MATRIX(xr);           // row effect design matrix (n x Kr, or 0 rows if none)
  DATA_MATRIX(offset);       // n x p offsets
  DATA_IMATRIX(Ntrials);     // n x p binomial trials
  DATA_IVECTOR(family);      // family code per species, length p
  DATA_VECTOR(extra);        // link/variant flag per species, length p
  DATA_INTEGER(num_lv);      // ordination dimension d
  DATA_INTEGER(method);      // must be 0 (VA)
  DATA_INTEGER(zetastruc);   // ordinal cutoff structure: 0=common, 1=species-specific
  DATA_INTEGER(p_betaH);     // number of betaH columns (beta-hurdle)
  DATA_INTEGER(random);      // bit 0 = random row effects; bit 1 = future extensions

  // ===== PARAMETERS =====
  // Fixed effects: b is Kx x p (each column is a species' coefficient vector)
  PARAMETER_MATRIX(b);

  // Site VA parameters
  PARAMETER_MATRIX(u);        // n x d  — a_i means
  PARAMETER_VECTOR(Au);       // n*d    — log-Cholesky diagonals for A_i
                              //   Au(k*n+i) = log(L_i(k,k)), A_i(k,k) = exp(2*Au(k*n+i))

  // Species VA parameters
  PARAMETER_MATRIX(a_lv_sp);  // p x d  — a_j means
  PARAMETER_VECTOR(Au_sp);    // p*d    — log-Cholesky diagonals for A_j
                              //   Au_sp(k*p+j) = log(L_j(k,k)), A_j(k,k) = exp(2*Au_sp(k*p+j))

  // Ordination scale: Sigma = diag(sigma)
  // sigma(0) = exp(sigmaLV(0))
  // sigma(k) = sigma(k-1) * invlogit(sigmaLV(k)),  k >= 1  (enforces ordering)
  PARAMETER_VECTOR(sigmaLV);  // length d

  // Dispersion / zero-inflation
  PARAMETER_VECTOR(lg_phi);
  PARAMETER_VECTOR(lg_phiZINB);
  PARAMETER_VECTOR(zeta);     // ordinal cutoff raw params
  PARAMETER(ePower);          // Tweedie power (logit-scale)

  // Row effects
  PARAMETER_MATRIX(r0f);      // fixed row effects (xr.rows() x Kr)
  PARAMETER_MATRIX(r0r);      // random row effects (n x 1 or 0 rows if not used)
  PARAMETER_VECTOR(lg_Ar);    // log-Cholesky diagonals for VA row effect covariance
  PARAMETER_VECTOR(log_sigma); // log(SD) for random row effects

  // ===== DERIVED SCALARS =====
  int n = y.rows();
  int p = y.cols();
  int truep = p - p_betaH;
  int d = num_lv;

  vector<Type> iphi = exp(lg_phi);

  parallel_accumulator<Type> nll(this);

  // ===== SIGMA construction =====
  vector<Type> sigma(d);
  sigma(0) = exp(sigmaLV(0));
  for (int k = 1; k < d; k++) {
    sigma(k) = sigma(k-1) * invlogit(sigmaLV(k));
  }
  vector<Type> sigma2(d);
  for (int k = 0; k < d; k++) sigma2(k) = sigma(k) * sigma(k);

  // ===== DIAGONAL VARIATIONAL COVARIANCES =====
  // Ai_diag(i,k) = A_i(k,k) = exp(2*Au(k*n+i))
  // Aj_diag(j,k) = A_j(k,k) = exp(2*Au_sp(k*p+j))
  matrix<Type> Ai_diag(n, d);
  matrix<Type> Aj_diag(p, d);
  for (int k = 0; k < d; k++) {
    for (int i = 0; i < n; i++) Ai_diag(i,k) = exp(Type(2)*Au(k*n+i));
    for (int j = 0; j < p; j++) Aj_diag(j,k) = exp(Type(2)*Au_sp(k*p+j));
  }

  // ===== KL DIVERGENCES =====
  // KL(q(z_i)||N(0,I)) = 1/2 [tr(A_i) + a_i^T a_i - d - log|A_i|]
  // = 1/2 [sum_k Ai_diag(i,k) + ||u_i||^2 - d - 2*sum_k Au(k*n+i)]
  // In TMB convention (matches gllvm.cpp): nll -= logdet_half - 1/2*(tr + ||a||^2)
  for (int i = 0; i < n; i++) {
    Type logdet_half = Type(0);
    Type trace_Ai   = Type(0);
    for (int k = 0; k < d; k++) {
      logdet_half += Au(k*n+i);
      trace_Ai    += Ai_diag(i,k);
    }
    nll -= logdet_half - Type(0.5)*(trace_Ai + u.row(i).squaredNorm());
  }
  nll -= Type(0.5)*n*d;

  // KL(q(gamma_j)||N(0,I))
  for (int j = 0; j < p; j++) {
    Type logdet_half = Type(0);
    Type trace_Aj   = Type(0);
    for (int k = 0; k < d; k++) {
      logdet_half += Au_sp(k*p+j);
      trace_Aj    += Aj_diag(j,k);
    }
    nll -= logdet_half - Type(0.5)*(trace_Aj + a_lv_sp.row(j).squaredNorm());
  }
  nll -= Type(0.5)*p*d;

  // ===== ETA and CQ =====
  matrix<Type> eta(n,p);
  eta.setZero();
  matrix<Type> cQ(n,p);
  cQ.setZero();
  matrix<Type> mu(n,p);   // scratch for probit/cloglog families
  mu.setZero();

  // Fixed effects: eta += x * b  (x[:,0] = 1 gives intercept via b[0,:])
  eta += x * b;

  // Offset
  if (offset.rows() == n) eta += offset;

  // Fixed row effects
  if (xr.rows() == n) eta += (xr * r0f).replicate(1, p);

  // Ordination contribution
  // eta(i,j)  += sum_k sigma(k) * u(i,k) * a_lv_sp(j,k)  =  a_i^T Sigma a_j
  // cQ(i,j)   = depends on family link (see below)
  //
  // For log-link families: cQ = log E_q[exp(eta)] - eta   (eq 5, ms.pdf)
  //   With diagonal A_i, A_j, Sigma:
  //   c_k  = 1 - sigma2(k) * Ai_diag(i,k) * Aj_diag(j,k)
  //   cQ  = sum_k [ -1/2 log(c_k)
  //                + 1/2 sigma2(k)*Ai_diag(i,k)*a_j(k)^2
  //                + sigma(k)*sigma2(k)*Ai_diag(i,k)*Aj_diag(j,k)*a_i(k)*a_j(k)/c_k
  //                + 1/2 sigma(k)*sigma2(k)*Aj_diag(j,k)^2*a_i(k)^2/c_k
  //                + 1/2 sigma(k)*sigma2(k)*Ai_diag(i,k)^2*a_j(k)^2/c_k ]
  //
  // For identity/logit-link families: cQ = 1/2 Var_q[eta_ij]
  //   = 1/2 sum_k sigma2(k)*(Ai_diag(i,k)*Aj_diag(j,k)
  //                          + Ai_diag(i,k)*a_j(k)^2 + Aj_diag(j,k)*a_i(k)^2)

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      // Ordination mean contribution to eta
      Type cross = Type(0);
      for (int k = 0; k < d; k++) {
        cross += sigma(k) * u(i,k) * a_lv_sp(j,k);
      }
      eta(i,j) += cross;

      // Determine link type from family code
      int fam = family(j);
      bool log_link = (fam == POISSON || fam == NEG_BINOMIAL || fam == GAMMA ||
                       fam == TWEEDIE || fam == ZIP || fam == EXPONENTIAL ||
                       fam == ZINB    || fam == ZNIB);

      Type cq = Type(0);
      if (log_link) {
        // Full log-MGF correction (eq 5, ms.pdf), diagonal case
        for (int k = 0; k < d; k++) {
          Type sk  = sigma(k);
          Type sk2 = sigma2(k);
          Type ui  = Ai_diag(i,k);
          Type vj  = Aj_diag(j,k);
          Type ai  = u(i,k);
          Type aj  = a_lv_sp(j,k);
          Type ck  = Type(1) - sk2 * ui * vj;
          cq -= Type(0.5) * log(ck);
          cq += Type(0.5) * sk2 * ui * aj * aj;
          cq += sk * sk2 * ui * vj * ai * aj / ck;
          cq += Type(0.5) * sk * sk2 * vj * vj * ai * ai / ck;
          cq += Type(0.5) * sk * sk2 * ui * ui * aj * aj / ck;
        }
      } else {
        // Half-variance correction for non-log-link families
        for (int k = 0; k < d; k++) {
          Type sk2 = sigma2(k);
          Type ui  = Ai_diag(i,k);
          Type vj  = Aj_diag(j,k);
          Type ai  = u(i,k);
          Type aj  = a_lv_sp(j,k);
          cq += sk2 * (ui*vj + ui*aj*aj + vj*ai*ai);
        }
        cq *= Type(0.5);
      }
      cQ(i,j) = cq;
    }
  }

  // ===== RANDOM ROW EFFECTS KL =====
  // Simple independent N(0, sigma_r^2) prior; VA covariance diagonal via lg_Ar
  if ((random & 1) > 0) {
    // eta contribution already handled via dr0 in family_va_nll.h for standard gllvm.
    // For HO: add random row effects directly (xr as identity mapping or supplied matrix)
    // r0r has the VA means; lg_Ar has their log-Cholesky diagonals.
    // KL(q(r_i)||N(0,sigma_r^2)) = sum_i [ Au_r(i) - 1/2*(exp(2*Au_r(i))/sigma_r^2 + r0r(i)^2/sigma_r^2) ] + n/2*(1 - log(sigma_r^2))
    if (r0r.rows() == n) {
      Type sigma_r = exp(log_sigma(0));
      Type sigma_r2 = sigma_r * sigma_r;
      for (int i = 0; i < n; i++) {
        Type Au_r = lg_Ar(i);
        Type Ar   = exp(Type(2)*Au_r);
        nll -= Au_r - Type(0.5)*(Ar/sigma_r2 + r0r(i,0)*r0r(i,0)/sigma_r2);
        eta.row(i).array() += r0r(i,0);
      }
      nll -= Type(0.5)*n*(Type(1) - Type(2)*log_sigma(0));
    }
  }

  // ===== VA FAMILY LIKELIHOODS =====
  // family_va_nll.h uses: n, p, truep, method, zetastruc, family, extra, iphi,
  //   lg_phi, lg_phiZINB, zeta, ePower, y, eta, cQ, mu, Ntrials, nll
  // It declares internally: idx, has12
  if ((method < 1) || (method > 1)) {
    #include "family_va_nll.h"
  } // method end

  return nll;
}
