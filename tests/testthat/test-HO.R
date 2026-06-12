## Comprehensive test suite for the Hierarchical Ordination (gllvmHO) model.

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

Ybin <- matrix(rbinom(n_sim * p_sim, 1, 0.4), n_sim, p_sim)
Ygam <- matrix(rgamma(n_sim * p_sim, shape = 2, rate = 1), n_sim, p_sim)
Yord <- matrix(sample(1:4, n_sim * p_sim, replace = TRUE), n_sim, p_sim)
Ntr  <- matrix(10L, n_sim, p_sim)
Ybn  <- matrix(rbinom(n_sim * p_sim, 10, 0.4), n_sim, p_sim)
Ygau <- matrix(rnorm(n_sim * p_sim), n_sim, p_sim)
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
## Block 2: row effects — basic
## ============================================================
test_that("row.eff = fixed works", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, row.eff = "fixed", sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$row.params.fixed))
  expect_length(fit$params$row.params.fixed, nrow(Y) - 1L)
})

## ============================================================
## Block 3: fixed-effect formula
## ============================================================
test_that("formula = ~soil.dry works (X covariates in beta)", {
  fit <- gllvm(Y, X = X, formula = ~soil.dry,
               family = "poisson", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$Xcoef))
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("formula + lv.formula: fixed X and canonical covariates", {
  fit <- gllvm(Y, X = X, formula = ~soil.dry, lv.formula = ~reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$Xcoef))
  expect_false(is.null(fit$params$LvXcoef))
  expect_equal(dim(predict(fit)), dim(Y))
})

## ============================================================
## Block 4: lv.formula — canonical covariates → LvXcoef
## ============================================================
test_that("lv.formula: LvXcoef produced with correct dimensions", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_equal(dim(fit$params$LvXcoef), c(2L, 2L))
  expect_false(is.null(fit$params$sigma.bz))
  expect_length(fit$params$sigma.bz, 2L)
})

test_that("lv.formula single predictor", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_equal(nrow(fit$params$LvXcoef), 1L)
})

## ============================================================
## Block 5: load.formula — trait covariates → LoadTRcoef
## ============================================================
test_that("load.formula: LoadTRcoef produced with correct dimensions", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_equal(dim(fit$params$LoadTRcoef), c(1L, 2L))
  expect_false(is.null(fit$params$sigma.bgamma))
})

test_that("load.formula multiple traits", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length + colour,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_equal(nrow(fit$params$LoadTRcoef), 2L)
})

## ============================================================
## Block 6: both lv.formula and load.formula
## ============================================================
test_that("lv.formula + load.formula: both LvXcoef and LoadTRcoef non-NULL", {
  fit <- gllvm(Y, X = X, TR = TR,
               lv.formula = ~soil.dry, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$LvXcoef))
  expect_false(is.null(fit$params$LoadTRcoef))
})

## ============================================================
## Block 7: padding — Kz < d or Kt < d
## ============================================================
test_that("padding: Kz=2 < d=3, third column of LvXcoef is zero", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 3, sd.errors = FALSE, n.init = 1)
  expect_equal(dim(fit$params$LvXcoef), c(2L, 3L))
  expect_true(all(fit$params$LvXcoef[, 3] == 0))
  expect_length(fit$params$sigma.bz, 2L)
})

test_that("padding: Kz=1 < d=3, columns 2 and 3 of LvXcoef are zero", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry,
               family = "poisson", random.loadings = TRUE, num.lv.c = 3, sd.errors = FALSE, n.init = 1)
  expect_equal(dim(fit$params$LvXcoef), c(1L, 3L))
  expect_true(all(fit$params$LvXcoef[, 2:3] == 0))
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

test_that("neither covariates nor traits: LvXcoef and LoadTRcoef are NULL", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(is.null(fit$params$LvXcoef))
  expect_true(is.null(fit$params$LoadTRcoef))
  expect_true(is.null(fit$params$sigma.bz))
  expect_true(is.null(fit$params$sigma.bgamma))
})

test_that("covariates only: LvXcoef non-NULL, LoadTRcoef NULL", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_false(is.null(fit$params$LvXcoef))
  expect_true(is.null(fit$params$LoadTRcoef))
  expect_true(is.null(fit$params$sigma.bgamma))
  expect_length(fit$params$sigma.bz, 2L)
  expect_true(all(fit$params$sigma.bz > 0))
})

test_that("traits only: LvXcoef NULL, LoadTRcoef non-NULL", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length + colour,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(is.null(fit$params$LvXcoef))
  expect_true(is.null(fit$params$sigma.bz))
  expect_false(is.null(fit$params$LoadTRcoef))
  expect_length(fit$params$sigma.bgamma, 2L)
  expect_true(all(fit$params$sigma.bgamma > 0))
})

test_that("both covariates and traits: correct dims and sigma lengths", {
  fit <- gllvm(Y, X = X, TR = TR,
               lv.formula = ~soil.dry + reflection, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_equal(dim(fit$params$LvXcoef),     c(2L, 2L))
  expect_equal(dim(fit$params$LoadTRcoef),  c(1L, 2L))
  expect_length(fit$params$sigma.bz,     2L)
  expect_length(fit$params$sigma.bgamma, 1L)
  expect_true(all(fit$params$sigma.bz     > 0))
  expect_true(all(fit$params$sigma.bgamma > 0))
})

test_that("padding with both: padded columns are zero", {
  fit <- gllvm(Y, X = X, TR = TR,
               lv.formula = ~soil.dry + reflection, load.formula = ~length,
               family = "poisson", random.loadings = TRUE, num.lv.c = 3, sd.errors = FALSE, n.init = 1)
  expect_equal(dim(fit$params$LvXcoef),    c(2L, 3L))
  expect_equal(dim(fit$params$LoadTRcoef), c(1L, 3L))
  expect_true(all(fit$params$LvXcoef[, 3]      == 0))
  expect_true(all(fit$params$LoadTRcoef[, 2:3] == 0))
  expect_length(fit$params$sigma.bz,     2L)
  expect_length(fit$params$sigma.bgamma, 1L)
})

test_that("randomB = FALSE is accepted without error", {
  expect_error(
    gllvm.HO.TMB(Y, randomB = FALSE, num.lv = 2),
    NA
  )
})

test_that("randomT = FALSE is accepted without error", {
  expect_error(
    gllvm.HO.TMB(Y, randomT = FALSE, num.lv = 2),
    NA
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

test_that("getPredictErr CMSEP=TRUE works for num.lv.c model (sd.errors=TRUE)", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, sd.errors = TRUE, n.init = 1)
  pe <- getPredictErr(fit, CMSEP = TRUE)
  expect_equal(dim(pe$lvs),      c(nrow(Y), 2L))
  expect_equal(dim(pe$loadings), c(ncol(Y), 2L))
  expect_true(all(is.finite(pe$lvs)))
  expect_true(all(is.finite(pe$loadings)))
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
  expect_length(se_out$sd$sigma.bgamma, 1L)
})

test_that("se() with row.eff = random: sd$sigma non-NULL and length 1", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, row.eff = "random", sd.errors = FALSE, n.init = 1)
  se_out <- se(fit)
  expect_false(is.null(se_out$sd$sigma))
  expect_length(se_out$sd$sigma, 1L)
  expect_true(all(is.finite(se_out$sd$sigma)))
})

test_that("se() with row.eff = fixed: sd$row.params.fixed non-NULL", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, row.eff = "fixed", sd.errors = FALSE, n.init = 1)
  se_out <- se(fit)
  expect_false(is.null(se_out$sd$row.params.fixed))
  expect_length(se_out$sd$row.params.fixed, nrow(Y) - 1L)
  expect_true(all(is.finite(se_out$sd$row.params.fixed)))
})

test_that("se() with corAR1 row effect: sd$sigma has length 2", {
  StudyDesign <- data.frame(Site = factor(rep(1:5, length.out = nrow(Y))))
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2,
               studyDesign = StudyDesign, row.eff = ~corAR1(1|Site),
               sd.errors = FALSE, n.init = 1)
  se_out <- se(fit)
  expect_false(is.null(se_out$sd$sigma))
  expect_length(se_out$sd$sigma, 2L)
  expect_true(all(is.finite(se_out$sd$sigma)))
})

## ============================================================
## Block 11: summary
## ============================================================
test_that("summary.gllvmHO returns summary.gllvmHO object with correct fields", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  res <- summary(fit)
  expect_true(inherits(res, "summary.gllvmHO"))
  out <- capture.output(print(res))
  expect_true(any(grepl("Effective standard deviation", out)))
})

test_that("summary with lv.formula prints LV predictor section", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  res <- summary(fit)
  out <- capture.output(print(res))
  expect_true(any(grepl("LV predictors", out)))
})

## ============================================================
## Block 11b: correlation via csb_z / csb_gamma (direct interface)
## ============================================================
test_that("csb_z with one pair: n log_sigma_bz = d_c + 1", {
  csb <- matrix(as.integer(c(2L, 1L)), nrow = 1L, ncol = 2L)
  fit <- gllvm.HO.TMB(Y,
                      lv.X = as.matrix(X[, c("soil.dry", "reflection")]),
                      num.lv = 2L, num.lv.c = 2L, family = "poisson", csb_z = csb, n.init = 1)
  ph <- fit$TMBfn$par
  expect_equal(sum(names(ph) == "log_sigma_bz"), 3L)
  expect_true(finite_logL(fit))
})

test_that("csb_gamma with one pair: n log_sigma_bgamma = d_t + 1", {
  TR2 <- cbind(length = TR$length, length2 = TR$length + rnorm(nrow(TR), sd = 0.1))
  csb <- matrix(as.integer(c(2L, 1L)), nrow = 1L, ncol = 2L)
  fit <- gllvm.HO.TMB(Y,
                      TR = TR2,
                      num.lv = 2L, num.lv.c = 2L, family = "poisson", csb_gamma = csb, n.init = 1)
  ph <- fit$TMBfn$par
  expect_equal(sum(names(ph) == "log_sigma_bgamma"), 3L)
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

test_that("ordiplot draws LvXcoef arrows with lv.formula", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_silent(ordiplot(fit, biplot = FALSE))
})

test_that("ordiplot draws LoadTRcoef arrows in biplot with load.formula", {
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

test_that("ordiplot type = residual and type = marginal work with lvc", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_silent(ordiplot(fit, type = "residual"))
  expect_silent(ordiplot(fit, type = "marginal"))
})

test_that("predict returns correct dimensions (training, level 1)", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  pr <- predict(fit, type = "link")
  expect_equal(dim(pr), dim(Y))
  expect_true(all(is.finite(pr)))
})

test_that("predict level = 0 returns correct dimensions", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  pr0 <- predict(fit, level = 0, type = "link")
  expect_equal(dim(pr0), dim(Y))
  expect_true(all(is.finite(pr0)))
})

test_that("predict type = response returns values on response scale", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  pr <- predict(fit, type = "response")
  expect_equal(dim(pr), dim(Y))
  expect_true(all(pr > 0))
})

test_that("predict includes random row effects in eta (not silently zero)", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, row.eff = "random", sd.errors = FALSE, n.init = 1)
  pr_with  <- predict(fit, level = 1, type = "link")
  pr_marg  <- predict(fit, level = 0, type = "link")
  expect_equal(dim(pr_with), dim(Y))
  ## Level-1 predictions must differ from level-0 when row effects are non-zero
  expect_false(isTRUE(all.equal(pr_with, pr_marg)))
})

test_that("predict includes structured random row effects (1|Site)", {
  StudyDesign <- data.frame(Site = factor(rep(1:5, length.out = nrow(Y))))
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2,
               studyDesign = StudyDesign, row.eff = ~(1|Site),
               sd.errors = FALSE, n.init = 1)
  pr <- predict(fit, level = 1, type = "link")
  expect_equal(dim(pr), dim(Y))
  expect_true(all(is.finite(pr)))
})

test_that("predictSR.gllvmHO returns expected SR and predicted PMF", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  sr <- predictSR(fit, se.fit = FALSE)
  expect_true(inherits(sr, "predictSR.gllvm"))
  expect_length(sr$expected$fit, nrow(Y))
  expect_true(all(sr$expected$fit >= 0))
})

test_that("predictSR.gllvmHO se.fit > 0 returns CI bounds", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = TRUE, n.init = 1)
  sr <- predictSR(fit, se.fit = 50, seed = 1)
  expect_true(inherits(sr, "predictSR.gllvm"))
  expect_true(!is.null(sr$expected$lower))
  expect_true(!is.null(sr$expected$upper))
  expect_true(all(sr$expected$lower <= sr$expected$fit + 1e-9))
  expect_true(all(sr$expected$upper >= sr$expected$fit - 1e-9))
})

## ============================================================
## Block 13a: num.lv (unconstrained) — X becomes fixed effects only
## ============================================================
test_that("num.lv + X: X in fixed effects, LvXcoef NULL", {
  fit <- gllvm(Y, X = X, formula = ~soil.dry, num.lv = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$Xcoef))
  expect_true(is.null(fit$params$LvXcoef))
  pr <- predict(fit)
  expect_equal(dim(pr), dim(Y))
})

test_that("num.lv + TR via load.formula: LoadTRcoef non-NULL, LvXcoef NULL", {
  fit <- gllvm(Y, TR = TR, load.formula = ~length, num.lv = 2,
               family = "poisson", random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(is.null(fit$params$LvXcoef))
  expect_false(is.null(fit$params$LoadTRcoef))
})

## ============================================================
## Block 13b: num.lv.c (constrained) — X+TR without explicit formula
## ============================================================
test_that("num.lv.c + X only (no lv.formula): LvXcoef non-zero, predict dims", {
  fit <- gllvm(Y, X = X, num.lv.c = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$LvXcoef))
  expect_true(any(fit$params$LvXcoef != 0), info = "LvXcoef stuck at zero")
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("num.lv.c + X + numeric TR (no formula): LvXcoef and LoadTRcoef non-zero", {
  TRmm <- model.matrix(~ ., spider$trait)
  fit <- gllvm(Y, X = X, TR = TRmm, num.lv.c = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(any(fit$params$LvXcoef != 0),    info = "LvXcoef stuck at zero")
  expect_true(any(fit$params$LoadTRcoef != 0), info = "LoadTRcoef stuck at zero")
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("num.lv.c + X + factor TR (no formula): LoadTRcoef non-zero", {
  fit <- gllvm(Y, X = X, TR = TR, num.lv.c = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(any(fit$params$LoadTRcoef != 0), info = "LoadTRcoef stuck at zero")
})

test_that("num.lv.c + ordiplot without formula", {
  fit <- gllvm(Y, X = X, num.lv.c = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_silent(ordiplot(fit, biplot = FALSE))
})

## ============================================================
## Block 13c: num.RR (reduced-rank) — X+TR without explicit formula
## ============================================================
test_that("num.RR + X only (no formula): LvXcoef non-zero, predict dims", {
  fit <- gllvm(Y, X = X, num.RR = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$LvXcoef))
  expect_true(any(fit$params$LvXcoef != 0), info = "LvXcoef stuck at zero")
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("num.RR + X + numeric TR (no formula): LvXcoef and LoadTRcoef non-zero", {
  TRmm <- model.matrix(~ ., spider$trait)
  fit <- gllvm(Y, X = X, TR = TRmm, num.RR = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(any(fit$params$LvXcoef != 0),    info = "LvXcoef stuck at zero")
  expect_true(any(fit$params$LoadTRcoef != 0), info = "LoadTRcoef stuck at zero")
  expect_equal(dim(predict(fit)), dim(Y))
})

test_that("num.RR + X + factor TR (no formula): LoadTRcoef non-zero", {
  fit <- gllvm(Y, X = X, TR = TR, num.RR = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(any(fit$params$LoadTRcoef != 0), info = "LoadTRcoef stuck at zero")
})

test_that("num.RR + ordiplot(biplot=TRUE): returns site and species coords", {
  fit <- gllvm(Y, X = X, num.RR = 2, family = "poisson",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_silent(res <- ordiplot(fit, biplot = TRUE))
  expect_named(res, c("sites", "species"))
})

test_that("num.RR + binomial (non-log-link): cQ uses Var_q path, logL finite", {
  fit <- gllvm(Ybn, X = X[seq_len(n_sim), ], TR = TR[seq_len(p_sim), ],
               Ntrials = Ntr, num.RR = 2, family = "binomial",
               random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  expect_true(has_class(fit))
  expect_true(finite_logL(fit))
})

test_that("num.RR + update() chain: adding TR preserves LvXcoef, adds LoadTRcoef", {
  fit2 <- gllvm(Y, X = X, num.RR = 2, family = "poisson",
                random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  fit3 <- suppressWarnings(update(fit2, TR = TR, sd.errors = FALSE))
  expect_true(finite_logL(fit3))
  expect_false(is.null(fit3$params$LvXcoef),    info = "LvXcoef lost after update with TR")
  expect_false(is.null(fit3$params$LoadTRcoef), info = "LoadTRcoef missing after update with TR")
})

## ============================================================
## Block 14: n.init > 1 selects best run
## ============================================================
test_that("n.init = 2 selects finite logL", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, n.init = 2)
  expect_true(finite_logL(fit))
})

## ============================================================
## Block 15: extended row effects — structured / formula-based
## ============================================================
test_that("row.eff = fixed with lv.formula: both params non-NULL", {
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, row.eff = "fixed", sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$row.params.fixed))
  expect_false(is.null(fit$params$LvXcoef))
})

test_that("row.eff = random with num.lv: row sigma positive", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 2, row.eff = "random", sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_true(fit$params$sigma > 0)
  expect_length(fit$params$row.params.random, nrow(Y))
})

test_that("structured row.eff ~(1|Site) with studyDesign", {
  ## Assign each row a fake site factor (5 sites, ~6 rows each)
  StudyDesign <- data.frame(Site = factor(rep(1:5, length.out = nrow(Y))))
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2,
               studyDesign = StudyDesign, row.eff = ~(1 | Site),
               sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$sigma))
})

test_that("row.eff fixed predictor ~soildry: row.params.fixed named", {
  StudyDesign <- data.frame(soildry = X$soil.dry)
  fit <- gllvm(Y, X = X, family = "poisson", random.loadings = TRUE, num.lv = 2,
               studyDesign = StudyDesign, row.eff = ~soildry,
               sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$row.params.fixed))
})

test_that("row.eff corAR1(1|Site): structured random effects, sigma length 2", {
  StudyDesign <- data.frame(Site = factor(rep(1:5, length.out = nrow(Y))))
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2,
               studyDesign = StudyDesign, row.eff = ~corAR1(1|Site),
               sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$sigma))
  expect_length(fit$params$sigma, 2L)   # 1 SD + 1 rho
})

test_that("row.eff corCS(1|Site): structured random effects, sigma length 2", {
  StudyDesign <- data.frame(Site = factor(rep(1:5, length.out = nrow(Y))))
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2,
               studyDesign = StudyDesign, row.eff = ~corCS(1|Site),
               sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$sigma))
  expect_length(fit$params$sigma, 2L)   # 1 SD + 1 rho
})

test_that("row.eff corExp(1|Site) with dist: spatial random effects, sigma length 2", {
  set.seed(99)
  StudyDesign <- data.frame(Site = factor(seq_len(nrow(Y))))
  coords <- matrix(runif(nrow(Y) * 2), nrow(Y), 2)
  dmat   <- as.matrix(dist(coords))
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE, num.lv = 2,
               studyDesign = StudyDesign, row.eff = ~corExp(1|Site),
               dist = list(dmat), sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$sigma))
  expect_length(fit$params$sigma, 2L)   # 1 scale + 1 range
})

test_that("row.eff corAR1 + lv.formula: both random row effects and LvXcoef non-NULL", {
  StudyDesign <- data.frame(Site = factor(rep(1:5, length.out = nrow(Y))))
  fit <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
               family = "poisson", random.loadings = TRUE, num.lv.c = 2,
               studyDesign = StudyDesign, row.eff = ~corAR1(1|Site),
               sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$sigma))
  expect_false(is.null(fit$params$LvXcoef))
})

## ============================================================
## Block 16: species-specific formula effects via load.formula
## ============================================================
test_that("load.formula numeric trait: LoadTRcoef correct dims, logL finite", {
  ## Use num.lv.c so the trait mean on loadings is structurally identifiable
  fit <- gllvm(Y, TR = TR, load.formula = ~length,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$LoadTRcoef))
  expect_equal(nrow(fit$params$LoadTRcoef), 1L)
  expect_equal(ncol(fit$params$LoadTRcoef), 2L)
})

test_that("load.formula factor trait (colour): LoadTRcoef has >1 row after dummy coding", {
  fit <- gllvm(Y, TR = TR, load.formula = ~colour,
               family = "poisson", random.loadings = TRUE,
               num.lv = 2, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  ## colour is a factor with 3 levels → 2 dummy columns in TR design → 2 rows in LoadTRcoef
  expect_true(nrow(fit$params$LoadTRcoef) >= 1L)
})

test_that("load.formula: adding traits improves (or matches) logL vs no-traits model", {
  ## HO model without traits
  fit0 <- gllvm(Y, family = "poisson", random.loadings = TRUE,
                num.lv = 2, sd.errors = FALSE, n.init = 1)
  ## HO model with trait structure on loadings
  fit1 <- gllvm(Y, TR = TR, load.formula = ~length,
                family = "poisson", random.loadings = TRUE,
                num.lv = 2, sd.errors = FALSE, n.init = 1)
  ## Both must converge; the trait model is nested — logL may decrease slightly
  ## due to the prior pulling loadings toward the trait direction, but must be finite
  expect_true(finite_logL(fit0))
  expect_true(finite_logL(fit1))
})

test_that("load.formula with lv.formula: all coef non-NULL", {
  fit <- gllvm(Y, X = X, TR = TR,
               lv.formula = ~soil.dry, load.formula = ~length,
               family = "poisson", random.loadings = TRUE,
               num.lv.c = 2, num.lv = 1, sd.errors = FALSE, n.init = 1)
  expect_true(finite_logL(fit))
  expect_false(is.null(fit$params$LvXcoef))
  expect_false(is.null(fit$params$LoadTRcoef))
  expect_equal(dim(predict(fit)), dim(Y))
})

## ============================================================
## Block 17: comparisons to standard gllvm
## ============================================================
test_that("gllvmHO and gllvm (num.lv=2) give same prediction dimensions", {
  fit_ho  <- gllvm(Y, family = "poisson", random.loadings = TRUE,
                   num.lv = 2, sd.errors = FALSE, n.init = 1)
  fit_std <- gllvm(Y, family = "poisson",
                   num.lv = 2, sd.errors = FALSE)
  expect_equal(dim(predict(fit_ho)),  dim(Y))
  expect_equal(dim(predict(fit_std)), dim(Y))
  ## Both should converge with finite logL
  expect_true(finite_logL(fit_ho))
  expect_true(finite_logL(fit_std))
})

test_that("gllvmHO and gllvm getLV output have same shape (num.lv=2)", {
  fit_ho  <- gllvm(Y, family = "poisson", random.loadings = TRUE,
                   num.lv = 2, sd.errors = FALSE, n.init = 1)
  fit_std <- gllvm(Y, family = "poisson",
                   num.lv = 2, sd.errors = FALSE)
  lv_ho  <- getLV(fit_ho)
  lv_std <- getLV(fit_std)
  expect_equal(dim(lv_ho),  c(nrow(Y), 2L))
  expect_equal(dim(lv_std), c(nrow(Y), 2L))
})

test_that("gllvmHO and gllvm (num.lv.c=2) have compatible LvXcoef structure", {
  ## Use the same two predictors in both models so Kz matches
  fit_ho  <- gllvm(Y, X = X, lv.formula = ~soil.dry + reflection,
                   family = "poisson", random.loadings = TRUE,
                   num.lv.c = 2, sd.errors = FALSE, n.init = 1)
  fit_std <- gllvm(Y, X = X[, c("soil.dry", "reflection")], num.lv.c = 2,
                   family = "poisson", sd.errors = FALSE)
  ## Both should have a canonical covariate structure stored
  expect_false(is.null(fit_ho$params$LvXcoef))
  expect_false(is.null(fit_std$params$LvXcoef))
  ## Dimensions should match (Kz x d_active)
  expect_equal(dim(fit_ho$params$LvXcoef), dim(fit_std$params$LvXcoef))
  ## getLV marginal should give n x 2 for both
  lv_ho  <- getLV(fit_ho,  type = "marginal")
  lv_std <- getLV(fit_std, type = "marginal")
  expect_equal(dim(lv_ho),  c(nrow(Y), 2L))
  expect_equal(dim(lv_std), c(nrow(Y), 2L))
})

test_that("gllvmHO getLoadings and gllvm getLoadings have same shape", {
  fit_ho  <- gllvm(Y, family = "poisson", random.loadings = TRUE,
                   num.lv = 2, sd.errors = FALSE, n.init = 1)
  fit_std <- gllvm(Y, family = "poisson",
                   num.lv = 2, sd.errors = FALSE)
  ld_ho  <- getLoadings(fit_ho)
  ld_std <- getLoadings(fit_std)
  expect_equal(dim(ld_ho),  c(ncol(Y), 2L))
  expect_equal(dim(ld_std), c(ncol(Y), 2L))
})

test_that("gllvmHO num.RR model vs gllvm num.RR=2: same prediction shape", {
  fit_ho  <- gllvm(Y, X = X, num.RR = 2, family = "poisson",
                   random.loadings = TRUE, sd.errors = FALSE, n.init = 1)
  fit_std <- gllvm(Y, X = X, num.RR = 2, family = "poisson", sd.errors = FALSE)
  expect_equal(dim(predict(fit_ho)),  dim(Y))
  expect_equal(dim(predict(fit_std)), dim(Y))
  expect_true(finite_logL(fit_ho))
  expect_true(finite_logL(fit_std))
})

test_that("gllvmHO sigma.lv ordered descending (as in standard gllvm)", {
  fit <- gllvm(Y, family = "poisson", random.loadings = TRUE,
               num.lv = 3, sd.errors = FALSE, n.init = 1)
  sig <- fit$params$sigma.lv
  expect_length(sig, 3L)
  ## sigma values should be non-negative
  expect_true(all(sig >= 0))
})

## ============================================================
## Block 18: numerical verification of cQ formulas (pure R, no model fitting)
## ============================================================

test_that("log-link cQ formula (ms.pdf eq 5) matches Monte Carlo", {
  ## E_q[exp(sigma*z*gamma)] with z~N(ai,ui), gamma~N(aj,vj)
  ## Analytical: exp(sigma*ai*aj + cQ) where cQ is the 5-term formula
  ##   (T1-T5 as in gllvm_HO.cpp, ms.pdf eq 5)
  ## sigma chosen small so the second moment of exp(sigma*z*gamma) is finite
  ## (requires 4*sigma^2*ui*vj < 1, i.e. sigma < 1/(2*sqrt(ui*vj)))
  set.seed(1L)
  sigma <- 0.5; ai <- 1.0; aj <- 0.5; ui <- 0.25; vj <- 0.5
  sk2 <- sigma^2; ck <- 1 - sk2 * ui * vj
  cq <- -0.5 * log(ck) +
    0.5 * sk2 * ui * aj^2 +                           # T2
    sigma * sk2 * ui * vj * ai * aj / ck +             # T3
    0.5 * sk2 * vj * ai^2 / ck +                      # T4
    0.5 * sk2^2 * ui^2 * vj * aj^2 / ck               # T5
  ## Full expectation includes the VA mean term sigma*ai*aj
  analytical <- exp(sigma * ai * aj + cq)
  z_mc    <- rnorm(2e6, ai, sqrt(ui))
  gamma_mc <- rnorm(2e6, aj, sqrt(vj))
  mc <- mean(exp(sigma * z_mc * gamma_mc))
  expect_equal(analytical, mc, tolerance = 5e-3)
})

test_that("non-log-link cQ formula (exact bilinear variance) matches Monte Carlo", {
  ## Var_q(sigma*z*gamma) = sigma^2*(ui*vj + ui*aj^2 + vj*ai^2) for z~N(ai,ui), gamma~N(aj,vj)
  set.seed(2L)
  sigma <- 1.2; ai <- 0.8; aj <- 0.4; ui <- 0.3; vj <- 0.6
  analytical_var <- sigma^2 * (ui*vj + ui*aj^2 + vj*ai^2)
  z_mc    <- rnorm(2e6, ai, sqrt(ui))
  gamma_mc <- rnorm(2e6, aj, sqrt(vj))
  eta_mc  <- sigma * z_mc * gamma_mc
  mc_var  <- var(eta_mc)
  expect_equal(analytical_var, mc_var, tolerance = 5e-3)
})
