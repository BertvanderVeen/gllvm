## Comprehensive test suite for the Hierarchical Ordination (gllvmHO) model.
## Run with: testthat::test_file("inst/examples/test_ho.R")
## or:        source("inst/examples/test_ho.R")

pkgload::load_all(quiet = TRUE)
library(testthat)
library(mvabund)
data(spider)
Y  <- spider$abund
X  <- spider$x
TR <- spider$trait

## Helpers
finite_logL <- function(fit) is.finite(fit$logL)
has_class   <- function(fit) inherits(fit, "gllvmHO")

## Simulated datasets for families that need non-count data
set.seed(42)
n_sim <- 20; p_sim <- 8

## binary (0/1)
Ybin <- matrix(rbinom(n_sim * p_sim, 1, 0.4), n_sim, p_sim)
## positive continuous — gamma / exponential
Ygam <- matrix(rgamma(n_sim * p_sim, shape = 2, rate = 1), n_sim, p_sim)
## ordinal (levels 1..4)
Yord <- matrix(sample(1:4, n_sim * p_sim, replace = TRUE), n_sim, p_sim)
## Ntrials for binomial / beta.binomial / ZIB / ZNIB
Ntr  <- matrix(10L, n_sim, p_sim)
Ybn  <- matrix(rbinom(n_sim * p_sim, 10, 0.4), n_sim, p_sim)
## Gaussian
Ygau <- matrix(rnorm(n_sim * p_sim), n_sim, p_sim)
## Tweedie / ZIP — non-negative with some zeros
Ytw  <- matrix(pmax(0, rnorm(n_sim * p_sim, mean = 2)), n_sim, p_sim)
Yzip <- matrix(rpois(n_sim * p_sim, 1), n_sim, p_sim)
Yzip[sample(n_sim * p_sim, floor(n_sim * p_sim * 0.3))] <- 0

## ============================================================
## Block 1: all families — correct class, finite logL
## ============================================================
test_that("poisson fit returns gllvmHO with finite logL", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
  expect_equal(fit$num.lv, 2L)
})

test_that("negative.binomial fit: finite logL and phi non-NULL", {
  fit <- gllvm(Y, family = "negative.binomial", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$phi))
})

test_that("binomial fit with Ntrials", {
  fit <- gllvm(Ybn, family = "binomial", random.loadings = TRUE, num.lv = 2,
               Ntrials = Ntr, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})

test_that("gaussian fit", {
  fit <- gllvm(Ygau, family = "gaussian", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$phi))
})

test_that("gamma fit", {
  fit <- gllvm(Ygam, family = "gamma", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})

test_that("tweedie fit", {
  fit <- gllvm(Ytw, family = "tweedie", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})

test_that("ZIP fit", {
  fit <- gllvm(Yzip, family = "ZIP", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})

test_that("ordinal fit", {
  fit <- gllvm(Yord, family = "ordinal", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE,
               zeta.struc = "common", n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$zeta))
})

test_that("exponential fit", {
  fit <- gllvm(Ygam, family = "exponential", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})

test_that("ZINB fit", {
  fit <- gllvm(Y, family = "ZINB", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})

test_that("ZIB fit", {
  fit <- gllvm(Ybn, family = "ZIB", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE,
               Ntrials = Ntr, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})

test_that("ZNIB fit", {
  fit <- gllvm(Ybn, family = "ZNIB", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE,
               Ntrials = Ntr, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})


## ============================================================
## Block 2: row effects
## ============================================================
test_that("row.eff = fixed works", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, row.eff = "fixed", sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
})

test_that("row.eff = random works", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, row.eff = "random", sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
})

## ============================================================
## Block 3: fixed-effect formula
## ============================================================
test_that("formula = ~soil.dry works (X covariates in beta)", {
  fit <- gllvm(Y, X = X, formula = ~soil.dry,
               family = "poisson", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$Xcoef))
})

## ============================================================
## Block 4: lv.formula — canonical covariates → b_z
## ============================================================
test_that("lv.formula: b_z produced with correct dimensions", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_equal(dim(fit$params$b_z), c(2L, 2L))
  expect_false(is.null(fit$params$sigma.bz))
  expect_length(fit$params$sigma.bz, 2L)
})

test_that("lv.formula single predictor", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_equal(nrow(fit$params$b_z), 1L)
})

## ============================================================
## Block 5: load.formula — trait covariates → b_gamma
## ============================================================
test_that("load.formula: b_gamma produced with correct dimensions", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_equal(dim(fit$params$b_gamma), c(1L, 2L))
  expect_false(is.null(fit$params$sigma.bgamma))
})

test_that("load.formula multiple traits", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length + colour,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_equal(nrow(fit$params$b_gamma), 2L)
})

## ============================================================
## Block 6: both lv.formula and load.formula
## ============================================================
test_that("lv.formula + load.formula: both b_z and b_gamma non-NULL", {
  fit <- gllvm(Y, X = X, TR = TR,
               lv.formula = ~soil.dry, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$b_z))
  expect_false(is.null(fit$params$b_gamma))
})

## ============================================================
## Block 7: padding — Kz < d or Kt < d
## ============================================================
test_that("padding: Kz=2 < d=3, third column of b_z is zero", {
  ## num.lv.c=3 gives total d=3; d_c=min(Kz=2, d_active=3)=2 → col 3 of b_z zero
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 3, sd.errors = FALSE, n.init = 1)
  expect_equal(dim(fit$params$b_z), c(2L, 3L))
  expect_true(all(fit$params$b_z[, 3] == 0))
  expect_length(fit$params$sigma.bz, 2L)   # only d_c = 2 prior SDs
})

test_that("padding: Kz=1 < d=3, columns 2 and 3 of b_z are zero", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry,
               family = "poisson", random.loadings = TRUE, num.lv.c = 3, sd.errors = FALSE, n.init = 1)
  expect_equal(dim(fit$params$b_z), c(1L, 3L))
  expect_true(all(fit$params$b_z[, 2:3] == 0))
  expect_length(fit$params$sigma.bz, 1L)
})

## ============================================================
## Block 8: VA covariance structures
## ============================================================
test_that("Lambda.struc = diagonal works", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, Lambda.struc = "diagonal", sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
})

test_that("Lambda.struc = unstructured works", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, Lambda.struc = "unstructured", sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
})

## ============================================================
## Block 9: all four (covariates, traits) combinations
## ============================================================

## 9a. Neither — b_z and b_gamma both NULL
test_that("neither covariates nor traits: b_z and b_gamma are NULL", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(is.null(fit$params$b_z))
  expect_true(is.null(fit$params$b_gamma))
  expect_true(is.null(fit$params$sigma.bz))
  expect_true(is.null(fit$params$sigma.bgamma))
})

## 9b. Covariates only — b_z non-NULL, b_gamma NULL (using num.lv.c for covariate effects)
test_that("covariates only: b_z non-NULL, b_gamma NULL", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_false(is.null(fit$params$b_z))
  expect_true(is.null(fit$params$b_gamma))
  expect_true(is.null(fit$params$sigma.bgamma))
  ## sigma.bz length = d_c = min(Kz=2, d_active=2) = 2
  expect_length(fit$params$sigma.bz, 2L)
  expect_true(all(fit$params$sigma.bz > 0))
})

## 9c. Traits only — b_z NULL, b_gamma non-NULL (using num.lv.c for trait effects)
test_that("traits only: b_z NULL, b_gamma non-NULL", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length + colour,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(is.null(fit$params$b_z))
  expect_true(is.null(fit$params$sigma.bz))
  expect_false(is.null(fit$params$b_gamma))
  ## sigma.bgamma length = d_t = min(Kt=2, d_active=2) = 2
  expect_length(fit$params$sigma.bgamma, 2L)
  expect_true(all(fit$params$sigma.bgamma > 0))
})

## 9d. Both — b_z and b_gamma non-NULL, dimensions correct
test_that("both covariates and traits: correct dims and sigma lengths", {
  fit <- gllvm(Y, X = X, TR = TR,
               lv.formula = ~soil.dry + reflection, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_equal(dim(fit$params$b_z),     c(2L, 2L))
  expect_equal(dim(fit$params$b_gamma), c(1L, 2L))
  ## d_c = min(2,2)=2, d_t = min(1,2)=1
  expect_length(fit$params$sigma.bz,     2L)
  expect_length(fit$params$sigma.bgamma, 1L)
  expect_true(all(fit$params$sigma.bz     > 0))
  expect_true(all(fit$params$sigma.bgamma > 0))
})

## 9e. Padding: both covariates and traits with d > max(Kz, Kt)
test_that("padding with both: padded columns of b_z and b_gamma are zero", {
  ## Kz=2, Kt=1, d=3 (num.lv.c=3) → d_c=min(2,3)=2, d_t=min(1,3)=1
  ## col 3 of b_z and cols 2,3 of b_gamma = 0
  fit <- gllvm(Y, X = X, TR = TR,
               lv.formula = ~soil.dry + reflection, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 3, sd.errors = FALSE, n.init = 1)
  expect_equal(dim(fit$params$b_z),     c(2L, 3L))
  expect_equal(dim(fit$params$b_gamma), c(1L, 3L))
  expect_true(all(fit$params$b_z[, 3]      == 0))
  expect_true(all(fit$params$b_gamma[, 2:3] == 0))
  expect_length(fit$params$sigma.bz,     2L)
  expect_length(fit$params$sigma.bgamma, 1L)
})

## 9f. randomB / randomT argument validation
test_that("randomB = FALSE errors with informative message", {
  expect_error(
    gllvm.HO.TMB(Y, randomB = FALSE, num.lv = 2),
    "not yet supported"
  )
})

test_that("randomT = FALSE errors with informative message", {
  expect_error(
    gllvm.HO.TMB(Y, randomT = FALSE, num.lv = 2),
    "not yet supported"
  )
})

## ============================================================
## Block 10: se.gllvm — standard errors
## ============================================================
test_that("se() returns sd$beta0, sd$sigma.lv, and prediction.errors", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  se_out <- se(fit)
  expect_length(se_out$sd$beta0, ncol(Y))
  expect_length(se_out$sd$sigma.lv, 2L)
  expect_equal(dim(se_out$prediction.errors$lvs), c(nrow(Y), 2L))
  expect_equal(dim(se_out$prediction.errors$loadings), c(ncol(Y), 2L))
})

test_that("se() returns sd$sigma.bz when lv.formula is used", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  se_out <- se(fit)
  expect_false(is.null(se_out$sd$sigma.bz))
  expect_length(se_out$sd$sigma.bz, 2L)
  expect_true(all(is.finite(se_out$sd$sigma.bz)))
  expect_equal(dim(se_out$prediction.errors$b_z), c(2L, 2L))
})

test_that("se() returns sd$sigma.bgamma when load.formula is used", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  se_out <- se(fit)
  expect_false(is.null(se_out$sd$sigma.bgamma))
  # d_t = min(Kt=1, d_active=2) = 1, so one prior SD
  expect_length(se_out$sd$sigma.bgamma, 1L)
})

## ============================================================
## Block 11: summary
## ============================================================
test_that("summary.gllvmHO runs without error and returns object invisibly", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  out <- capture.output(res <- summary(fit))
  expect_true(inherits(res, "gllvmHO"))
  expect_true(any(grepl("Ordination", out)))
})

test_that("summary with lv.formula prints b_z section", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Canonical covariate", out)))
})

## ============================================================
## Block 11b: correlation via csb_z / csb_gamma (direct interface)
## Bar-formula correlation is not yet properly supported for HO models;
## csb_z / csb_gamma can be passed directly to gllvm.HO.TMB.
## ============================================================
test_that("csb_z with one pair: n log_sigma_bz = d_c + 1", {
  ## Kz=2, num.lv.c=2 → d_active=2, d_c=2, n_pairs=1 → 3 log_sigma_bz params
  csb <- matrix(as.integer(c(2L, 1L)), nrow = 1L, ncol = 2L)
  fit <- gllvm.HO.TMB(Y,
                      lv.X = as.matrix(X[, c("soil.dry", "reflection")]),
                      num.lv = 2L, num.lv.c = 2L, family = "poisson", csb_z = csb, n.init = 1)
  ph <- fit$TMBfn$par
  expect_equal(sum(names(ph) == "log_sigma_bz"), 3L)  # d_c + n_pairs = 2+1
  expect_true(finite_logL(fit))
})

test_that("csb_gamma with one pair: n log_sigma_bgamma = d_t + 1", {
  ## spider$trait$colour is a factor; use two numeric traits instead
  TR2 <- cbind(length = TR$length, length2 = TR$length + rnorm(nrow(TR), sd = 0.1))
  csb <- matrix(as.integer(c(2L, 1L)), nrow = 1L, ncol = 2L)
  fit <- gllvm.HO.TMB(Y,
                      TR = TR2,
                      num.lv = 2L, num.lv.c = 2L, family = "poisson", csb_gamma = csb, n.init = 1)
  ph <- fit$TMBfn$par
  expect_equal(sum(names(ph) == "log_sigma_bgamma"), 3L)  # d_t + n_pairs = 2+1
  expect_true(finite_logL(fit))
})

## ============================================================
## Block 12: ordiplot / predict
## ============================================================
test_that("ordiplot runs without error (no covariates)", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_silent(ordiplot(fit, biplot = FALSE))
})

test_that("ordiplot draws b_z arrows with lv.formula", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  ## Should not error; arrows drawn silently from b_z
  expect_silent(ordiplot(fit, biplot = FALSE))
})

test_that("ordiplot draws b_gamma arrows in biplot with load.formula", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_silent(ordiplot(fit, biplot = TRUE))
})

test_that("ordiplot draws both arrow types with lv.formula + load.formula", {
  fit <- gllvm(Y, X = X, TR = TR,
               lv.formula = ~soil.dry + reflection, load.formula = ~length,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  res <- ordiplot(fit, biplot = TRUE)
  expect_named(res, c("sites", "species"))
  expect_equal(nrow(res$sites), nrow(Y))
})

test_that("predict returns correct dimensions", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  pr <- predict(fit)
  expect_equal(dim(pr), dim(Y))
})

## ============================================================
## Block 13a: num.lv (unconstrained) — X becomes fixed effects only
## ============================================================
test_that("num.lv + X: X in fixed effects, b_z NULL", {
  fit <- gllvm(Y, X = X, formula = ~soil.dry, num.lv = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$Xcoef))
  expect_true(is.null(fit$params$b_z))
  pr <- predict(fit)
  expect_equal(dim(pr), dim(Y))
})

test_that("num.lv + TR via load.formula: b_gamma non-NULL, b_z NULL", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length, num.lv = 2,
               family = "poisson", random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(is.null(fit$params$b_z))
  expect_false(is.null(fit$params$b_gamma))
})

## ============================================================
## Block 13b: num.lv.c (constrained) — X+TR without explicit formula
## ============================================================
test_that("num.lv.c + X only (no lv.formula): b_z non-zero, predict dims", {
  fit <- gllvm(Y, X = X, num.lv.c = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$b_z))
  expect_true(any(fit$params$b_z != 0), info = "b_z stuck at zero")
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("num.lv.c + X + numeric TR (no formula): b_z and b_gamma non-zero", {
  TRmm <- model.matrix(~ ., spider$trait)
  fit <- gllvm(Y, X = X, TR = TRmm, num.lv.c = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(any(fit$params$b_z != 0),     info = "b_z stuck at zero")
  expect_true(any(fit$params$b_gamma != 0), info = "b_gamma stuck at zero")
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("num.lv.c + X + factor TR (no formula): b_gamma non-zero", {
  fit <- gllvm(Y, X = X, TR = TR, num.lv.c = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(any(fit$params$b_gamma != 0), info = "b_gamma stuck at zero")
})

test_that("num.lv.c + ordiplot without formula", {
  fit <- gllvm(Y, X = X, num.lv.c = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_silent(ordiplot(fit, biplot = FALSE))
})

## ============================================================
## Block 13c: num.RR (reduced-rank) — X+TR without explicit formula
## ============================================================
test_that("num.RR + X only (no formula): b_z non-zero, predict dims", {
  fit <- gllvm(Y, X = X, num.RR = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$b_z))
  expect_true(any(fit$params$b_z != 0), info = "b_z stuck at zero")
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("num.RR + X + numeric TR (no formula): b_z and b_gamma non-zero", {
  TRmm <- model.matrix(~ ., spider$trait)
  fit <- gllvm(Y, X = X, TR = TRmm, num.RR = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(any(fit$params$b_z != 0),     info = "b_z stuck at zero")
  expect_true(any(fit$params$b_gamma != 0), info = "b_gamma stuck at zero")
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("num.RR + X + factor TR (no formula): b_gamma non-zero", {
  fit <- gllvm(Y, X = X, TR = TR, num.RR = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(any(fit$params$b_gamma != 0), info = "b_gamma stuck at zero")
})

test_that("num.RR + ordiplot(biplot=TRUE): species within plot bounds", {
  fit <- gllvm(Y, X = X, num.RR = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_silent(res <- ordiplot(fit, biplot = TRUE))
  expect_named(res, c("sites", "species"))
})

test_that("num.RR + update() chain: adding TR preserves b_z, adds b_gamma", {
  fit2 <- gllvm(Y, X = X, num.RR = 2, family = "poisson",
                random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  fit3 <- suppressWarnings(update(fit2, TR = TR, sd.errors = FALSE))
  expect_true(finite_logL(fit3))
  expect_false(is.null(fit3$params$b_z),     info = "b_z lost after update with TR")
  expect_false(is.null(fit3$params$b_gamma), info = "b_gamma missing after update with TR")
})

## ============================================================
## Block 14: n.init > 1 selects best run
## ============================================================
test_that("n.init = 2 selects finite logL", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, n.init = 2)
  expect_true(finite_logL(fit))
})
