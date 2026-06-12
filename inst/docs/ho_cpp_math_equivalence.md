# Hierarchical Ordination (HO) — C++ to Math Mapping

**Source files:**
- C++: `src/gllvm_HO.cpp`
- Math: `~/Dropbox/HO_VA/ms.Rmd`

---

## 0. Notation Glossary

| Math symbol | C++ variable | Location |
|---|---|---|
| d | `d = num_lv` | line 145 |
| num_RR | `num_RR` | DATA_INTEGER, line 78 |
| num_lvc | `num_lvc` | DATA_INTEGER, line 79 |
| num_lv (pure) | `num_lv_unc = d - num_RR - num_lvc` | line 146 |
| K_z | `Kz = lv_X_env.cols()` | line 148 |
| K_t | `Kt = TR.cols()` | line 149 |
| σ_k | `sigma(k)` | lines 170–175: `sigma(0)=exp(sigmaLV(0))`, `sigma(k)=sigma(k-1)*exp(-exp(sigmaLV(k)))` |
| Σ | `sigma` vector; applied element-wise as `sigma(k)*...` | |
| E_q[z_{ik}] | `a_z(k)` | line 568 |
| Var_q[z_{ik}] | `u_z(k)` | line 569 |
| E_q[γ_{jk}] | `a_g(k)` | line 653 |
| Var_q[γ_{jk}] | `v_g(k)` | line 654 |
| Var_q(z_i) as d×d | `Ai_full` | lines 594–613 |
| Var_q(γ_j) as d×d | `Aj` | lines 683–702 |
| Σ A_i Σ | `Mi` | lines 614–616 |
| Σ A_j Σ | `Mj` | lines 705–707 |
| C_{ij} = I − Σ A_j Σ A_i | `C = -(Mj * Ai_full); C(k,k)+=1` | lines 722–723 |
| VA Cholesky of q(b_z[:,l]) | `AB_z(l)` (Kz×Kz lower-triangular) | lines 273–286 |
| VA Cholesky of q(b_γ[:,l]) | `AB_g(l)` (Kt×Kt lower-triangular) | lines 288–301 |
| Var_q(b_z[:,l]) | `AB_z(l) * AB_z(l).transpose()` | line 554 |
| Σ_{lvc} Var_q(b_z) | `bz_var_mat` (Kz×Kz) | lines 307–313 |
| Σ_{lvc} Var_q(b_γ) | `bg_var_mat` (Kt×Kt) | lines 315–320 |
| x_i^T Var_q(b_z[:,k]) x_i | `lvc_z_var_i[k-num_RR]` | lines 558–564 |
| t_j^T Var_q(b_γ[:,k]) t_j | `lvc_g_var_j[k-num_RR]` | lines 642–649 |

---

## 1. Dimension-type layout and VA index mapping

The full d-dimensional block is ordered `[RR | lvc | lv]`:

```
k = 0 .. num_RR-1               : RR dims
k = num_RR .. num_RR+num_lvc-1  : concurrent (lvc) dims
k = num_RR+num_lvc .. d-1       : pure latent (lv) dims
```

The VA arrays `u` (sites) and `a_lv_sp` (species) cover only **d_va_z** and **d_va_a** columns
respectively, which exclude deterministic RR dims:

```
rr_va_z = (Kz > 0) ? 0 : num_RR   // RR z is deterministic iff Kz > 0
rr_va_a = (Kt > 0) ? 0 : num_RR
d_va_z  = rr_va_z + num_lvc + num_lv_unc
d_va_a  = rr_va_a + num_lvc + num_lv_unc
```

`va_z_idx[k]` and `va_a_idx[k]` return the VA column index (or −1 if deterministic).

Active b_z / b_γ columns:

```
d_active = num_RR + num_lvc   // lv dims never use covariates/traits
d_c      = min(Kz, d_active)  // active b_z columns
d_t      = min(Kt, d_active)  // active b_γ columns
dc_lvc   = d_c - min(d_c, num_RR)
dt_lvc   = d_t - min(d_t, num_RR)
```

---

## 2. Sigma construction

ms.Rmd: Σ = diag(σ_1,...,σ_d), σ_1 > σ_2 > ... > σ_d > 0.

C++ (lines 170–175):
```
sigma(0) = exp(sigmaLV(0))
sigma(k) = sigma(k-1) * exp(-exp(sigmaLV(k)))   for k >= 1
```
The multiplicative shrinkage `exp(-exp(sigmaLV(k)))` ∈ (0,1) guarantees strict ordering.
A log-barrier (lines 179–180) penalises the ratio approaching 1.

**Note — sigma_bz identifiability for RR dims.** The R interface pins `log_sigma_bz[l] = 0`
for `l < num_RR`, so sigma_bz = 1 for all RR dims. Without this pin, sigma_k × sigma_bz_k
is unidentified. Same logic for sigma_bgamma at RR dims.

---

## 3. VA covariance matrices

### 3.1 Diagonal VA (va_struct = 0)

`Ai_diag(i, iz) = exp(2 * Au(iz*n + i))` (line 193).
`Aj_diag(j, ia) = exp(2 * Au_sp(ia*p + j))` (line 195).

### 3.2 Unstructured VA (va_struct = 1)

Cholesky factor `Li_va` (d_va_z × d_va_z) assembled from `Au` with `exp()` on the diagonal.
`Ai_va = Li_va * Li_va^T`. Diagonal entries used for `Ai_diag`. Same for `Lj_va`.

---

## 4. Per-case breakdown

### Case 1: Pure lv (num_RR = 0, num_lvc = 0, Kz = 0, Kt = 0)

**a_z[k], u_z[k]** (lines 568–585):
- `va_z_idx[k] = k` for all k
- `a_z(k) = u(i, k)` — VA mean
- `u_z(k) = Ai_diag(i, k)` — residual VA variance (no b_z correction)

Math: a_z(k) = a_i[k]; u_z(k) = (A_i)_{kk} (ms.Rmd §"Variational likelihood").

**a_g[k], v_g[k]** (lines 653–668): symmetric, from `a_lv_sp` and `Aj_diag`.

**Ai_full**: `diag(u_z)` (diagonal VA) or full Cholesky-based d×d (unstructured VA).

**η mean** (line 671): `Σ_k σ_k * a_z(k) * a_g(k)`  
Math: ē_{ij} = β_{0j} + a_z^T Σ a_g (ms.Rmd prose definition, "bar_eta_ij").

**cQ — log-link** (lines 716–745):
```cpp
C   = I - Mj * Ai_full       // C_{ij} = I - Σ A_j Σ A_i
cq  = -0.5*log|C| + t1_t5 + t2_t4 + t3
```
where `t2_t4 = 0.5 * a_z^T C^{-1} Mj a_z`, `t1_t5 = 0.5 * (A_i Σ a_g)^T C^{-1} Σ a_g`,
`t3 = (A_i Σ a_g)^T C^{-1} Mj a_z`.  
Math: ms.Rmd eq. `\label{solutionfin}`.

**cQ — non-log-link** (lines 746–755):
```cpp
cq = 0.5 * (tr(Mi*Aj) + a_g^T Mi a_g + a_z^T Mj a_z)
```
Math: ms.Rmd eq. `\label{varbilinear}` — `= 0.5 Var_q(z_i^T Σ γ_j)`, exact result.

**KL terms**: standard Gaussian KL over d_va_z / d_va_a dims (no b_z/b_γ corrections).

---

### Case 2: lvc X-only (num_lvc > 0, Kz > 0, Kt = 0, randomB = 1)

**u_z[k]** (lines 558–584):
- `u_z(k) = Ai_diag(i, iz) + lvc_z_var_i[k-num_RR]`

Math: ms.Rmd eq. `\label{lvcvar}`:  
`Var(z_i) = A_i + Var(B_z^T x_i)` where `lvc_z_var_i[k] = x_i^T AB_z(k) AB_z(k)^T x_i`
(eq. `\label{varbzLV}` with single active column).

**KL(q(z_i) || N(μ_z, I))** (lines 332–376):
- `μ_z(iz) = lv_X_i * b_z_hat[:,k]` (lines 337–340)
- `bz_cross = x_i^T bz_var_mat x_i` (line 373)

Math (ms.pdf §3.2 / C++ header comment lines 41–48):
```
E_{q(b_z)}[||z_i - lv_X_i b_z||^2] = ||u_i - lv_X_i b_z_hat||^2 + x_i^T Var_q(b_z) x_i
```

**KL(q(b_z[:,l]) || N(0, σ_{bz}^2 I))** (lines 425–458):
```
-log|diag(AB_z(l))| + 0.5/σ^2 * (||AB_z(l)||_F^2 + ||b_z_hat[:,l]||^2) - 0.5*(Kz - Kz*2*log_σ)
```
Math: KL(N(μ, V) || N(0, σ²I)) = 0.5[tr(V)/σ² + μ^Tμ/σ² − Kz + Kz log(σ²) − log|V|].

---

### Case 3: lvc TR-only (num_lvc > 0, Kt > 0, Kz = 0, randomT = 1)

Mirror of Case 2 on the species side.

**v_g[k]** = `Aj_diag(j, ia) + lvc_g_var_j[k-num_RR]`  
Math: ms.Rmd eq. `\label{lvcvar}` applied to γ_j.

**KL(q(γ_j) || N(μ_g, I))** (lines 379–423): `bg_cross = t_j^T bg_var_mat t_j`.

**KL(q(b_γ[:,l]))** (lines 461–492): symmetric to b_z KL.

---

### Case 4: lvc X+TR (num_lvc > 0, Kz > 0, Kt > 0, randomB = 1, randomT = 1)

Both bz_cross and bg_cross active. Ai_full and Aj both augmented with their respective
b_z / b_γ variance corrections. All four KL terms operate simultaneously.

---

### Case 5: RR with fixed b_z (num_RR > 0, Kz > 0, randomB = 0)

**a_z[k] for k < num_RR** (lines 580–584):
- `va_z_idx[k] = -1` (deterministic branch)
- `a_z(k) = z_rr_det[k] = Σ_l lv_X_env(i,l) * b_z(l,k)`
- `u_z(k) = 0` (fixed b_z → no VA uncertainty)

Math: Var_q(z_{ik}) = 0 since ε_i = 0 for RR and b_z fixed.

**Ai_full**: zero at RR positions (u_z = 0). C = I − Mj Ai_full has identity blocks at those positions.

**Scale TODO**: The column norms of b_z[:,k] and σ_k are not separately identified. The proper
fix — orthonormal constraints via the **Stiefel manifold** — is marked TODO at line 549 of
`gllvm_HO.cpp`: `// randomB=0; TODO Stiefel manifold`.

---

### Case 6: RR with random b_z (num_RR > 0, Kz > 0, randomB = 1)

**u_z[k] for k < num_RR** (lines 547–556):
- `u_z(k) = rr_z_var_i[k] = x_i^T AB_z(k) AB_z(k)^T x_i` (line 555)

Math: Var_q(z_{ik}) = x_i^T Var_q(b_z[:,k]) x_i (eq. `\label{varbzLV}` with A_i = 0;
no ε_i residual for RR dims).

**Ai_full**: `Ai_full(k,k) = rr_z_var_i[k]` at RR positions (lines 609–610).

**KL(q(b_z[:,l]))**: active for l = 0..d_c−1, including RR columns. For l < num_RR, sigma_bz
is pinned to 1 by the R-side map (identifiability constraint, see §2 note).

---

### Case 7: RR with X+TR, randomB = 1, randomT = 1

**v_g[k] for k < num_RR** (lines 631–637):
- `v_g(k) = rr_g_var_j[k] = t_j^T AB_g(k) AB_g(k)^T t_j`

Aj at RR positions: `Aj(k,k) += rr_g_var_j[k]` (lines 698–699).

Both bz and bg KL active. Sigma_bgamma for RR cols pinned to 1 (same as sigma_bz).

---

### Case 8: Mixed (num_RR > 0, num_lvc > 0, num_lv_unc > 0)

All cases are active simultaneously. The `va_z_idx[k]` / `va_a_idx[k]` dispatch handles each
dim type within a single loop over k = 0..d−1.

| k range | Type | va_z_idx | a_z source | u_z source |
|---|---|---|---|---|
| 0..num_RR-1, Kz>0 | RR | −1 | z_rr_det[k] | rr_z_var_i[k] (0 if randomB=0) |
| 0..num_RR-1, Kz=0 | RR | k | u(i,k) | Ai_diag(i,k) |
| num_RR..num_RR+num_lvc−1 | lvc | rr_va_z+(k−num_RR) | u(i,iz) | Ai_diag(i,iz)+lvc_z_var_i (if has_bz) |
| num_RR+num_lvc..d−1 | lv | rr_va_z+(k−num_RR) | u(i,iz) | Ai_diag(i,iz) |

Same table for γ side with `va_a_idx`, `a_g`, `v_g`, `a_lv_sp`, `Aj_diag`, `lvc_g_var_j`.

**Ai_full construction** (lines 594–613): diagonal entries directly from u_z (which already
incorporates all per-dim-type corrections); unstructured VA augments the Cholesky block with
rr_z_var_i (lines 609–610) and lvc_z_var_i (lines 611–612) on the diagonal.

---

## 5. Numerical guard: CppAD::CondExpGt on det(C)

**Location**: lines 724–731

```cpp
Type detC      = C.determinant();
Type detC_safe = CppAD::CondExpGt(detC, Type(1e-8), detC, Type(1e-8));
Type logdetC   = log(detC_safe);
Type reg       = CppAD::CondExpLt(detC, Type(1e-4), Type(1e-4)-detC+Type(1e-4), Type(0));
for (int k = 0; k < d; k++) C_reg(k,k) += reg;
matrix<Type> Cinv = C_reg.inverse();
```

**Purpose**: the log-link formula (ms.Rmd eq. `\label{solutionfin}`) requires C_{ij} ≻ 0.
This is satisfied at convergence, but can fail during optimisation.

**Two guards**:
1. `CppAD::CondExpGt`: clamps `detC_safe ≥ 1e-8` before `log()`.
2. Diagonal regularisation: when `detC < 1e-4`, adds a positive amount to the diagonal before
   inversion, preventing a numerically singular system.

`CppAD::CondExpGt(a, b, t, f)` is the AD-compatible conditional: `(a > b) ? t : f`, smooth in
the tape so gradients remain valid even when the condition triggers.

Math: ms.Rmd eq. (solutionfin) comment (line 63): "this form is not numerically stable as the
inverse for C_{ij} only exists if all eigenvalues of Σ A_j Σ A_i are below one."

---

## 6. Identifiability constraints and TODOs

### 6.1 Sigma ordering constraint

Strict ordering enforced via log-barrier (lines 179–180). Identifies ordination axes by
variation explained (most important dimension first).

### 6.2 sigma_bz pinned to 1 for RR dims (randomB = 1)

C++ comment line 39: "sigma_bz is not separately identified from sigma; fix log_sigma_bz to 0
for RR dims." Enforced in the R interface (not in C++). Same for sigma_bgamma (randomT = 1).

### 6.3 TODO: Stiefel manifold for randomB = 0 lvc dims

Line 549 of `gllvm_HO.cpp`:
```cpp
//   Var_q(z_ik) = A_eps_i[k]    (randomB=0; TODO Stiefel manifold)
```

When `randomB = 0` and `num_lvc > 0`, b_z columns are fixed parameters and the product
`σ_k * b_z[:,k]` has a rotational ambiguity. Proper fix: orthonormal columns via thin QR /
Stiefel manifold constraint. Not yet implemented.

### 6.4 etamean label

C++ comment at line 670 references `ms.Rmd eq etamean`. This label does not yet exist in
ms.Rmd; the η̄_{ij} = β_{0j} + a_z^T Σ a_g formula is given in prose (line ~91 of ms.Rmd)
without a `\label{}`. Should be added as `\label{etamean}`.

---

## 7. VA covariance parameterisation for q(b_z) and q(b_γ)

The VA covariance for `b_z[:,l]` is stored as lower-triangular Cholesky `AB_z(l)` (Kz×Kz):

```
diagonal of AB_z(l): Ab_z(q * d_c + l)  for q=0..Kz-1   (stored as log)
off-diagonal (r>c):  Ab_z(Kz*d_c + k*d_c + l)             (stored raw)
```

`Var_q(b_z[:,l]) = AB_z(l) * AB_z(l)^T`

`bz_var_mat` (Kz×Kz) accumulates over lvc-active l (lines 308–313):
```cpp
bz_var_mat = Σ_{l in dc_lvc} AB_z(l) * AB_z(l)^T
```
Used in `bz_cross = x_i^T bz_var_mat x_i` (line 373), matching ms.Rmd §"CROSS-COVARIANCE
CORRECTION" (C++ header comment, lines 41–48).

Note: `bz_var_mat` accumulates only over dc_lvc lvc columns. The RR columns contribute
to `rr_z_var_i` (used in Ai_full) but not to `bz_cross` (RR dims have no VA z residual).

---

## 8. ELBO structure summary

The full negative ELBO in `nll` consists of:

1. **Log-barrier for sigma ordering** (lines 179–180)
2. **KL(q(z_i) || N(μ_z, I))** for i=1..n over d_va_z dims (lines 332–376)  
   — includes bz_cross for lvc dims when randomB=1
3. **KL(q(γ_j) || N(μ_g, I))** for j=1..p over d_va_a dims (lines 379–423)  
   — includes bg_cross for lvc dims when randomT=1
4. **Constants** `0.5*n*d_va_z` and `0.5*p*d_va_a` (lines 377, 423)
5. **KL(q(b_z[:,l]) || N(0, σ_{bz}^2 I))** for l=0..d_c−1 when has_bz (lines 425–458)
6. **KL(q(b_γ[:,l]) || N(0, σ_{bγ}^2 I))** for l=0..d_t−1 when has_bg (lines 461–492)
7. **η and cQ construction** (lines 496–759)
8. **Species-specific random effects** (`species_effects.h`)
9. **Row random effects** (`row_effects.h`)
10. **Family VA likelihoods** (`family_va_nll.h`) — consumes η and cQ

---

## 9. Equation-label cross-reference table

| Computation | ms.Rmd label | C++ lines |
|---|---|---|
| Numerically stable cQ (log-link) | `\label{solutionfin}` | 716–745 |
| Probit VA likelihood | `\label{probitVA}` | family_va_nll.h |
| Var(η_{ij}) — non-log-link | `\label{varbilinear}` | 746–755 |
| E(z_i) with covariates | `\label{lvcmean}` | 337–340, 568–585 |
| Var(z_i) with covariates | `\label{lvcvar}` | 548–584 |
| Var(B_z^T x_i) — row prior | `\label{varbzLV}` | 553–563 |
| η mean (prose, no label yet) | should be `\label{etamean}` | 671 |
| Cross-covariance correction | ms.pdf §3.2 / C++ header lines 41–48 | 366–375, 413–420 |
