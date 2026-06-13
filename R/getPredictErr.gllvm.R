#' @title Extract prediction errors for latent variables from gllvm object
#' @description  Calculates the prediction errors for latent variables and random effects for gllvm model.
#'
#' @param object   an object of class 'gllvm'.
#' @param CMSEP logical, if \code{TRUE} conditional mean squared errors for predictions are calculated. If \code{FALSE}, prediction errors are based on covariances of the variational distributions for \code{method ="VA"} and \code{method ="EVA"}.
#' @param cov if \code{TRUE}, return as covariances/variances of predictions. Otherwise \code{FALSE} (default) return as standard errors of predictions.
#' @param ...	 not used
#'
#' @details 
#' Calculates conditional mean squared errors for predictions.
#' If variational approximation is used, prediction errors can be based on covariances 
#' of the variational distributions, and therefore they do not take into account 
#' the uncertainty in the estimation of (fixed) parameters. 
#'
#' @return Function returns following components:
#'  \item{lvs }{prediction errors for latent variables}
#'  \item{row.effects }{prediction errors for random row effects if included}
#'
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton
#'
#' @examples
#'\dontrun{
#'# Load a dataset from the mvabund package
#'data(antTraits, package = "mvabund")
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# prediction errors for latent variables:
#'getPredictErr(fit)
#'}
#'
#'@aliases getPredictErr getPredictErr.gllvm
#'@method getPredictErr gllvm
#'@export
#'@export getPredictErr.gllvm
getPredictErr.gllvm = function(object, CMSEP = TRUE, cov = FALSE, ...)
{
  if(!is.list(object$sd)){
    stop("Cannot calculate prediction errors without standard errors in the model.")
  }
  # backward compatibility
  
  if(is.null(object$params$row.params.random) && !inherits(object$row.eff, "formula") && object$row.eff == "random")object$params$row.params.random <- object$params$row.params
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X
  if(is.null(object$col.eff$col.eff))object$col.eff$col.eff <- FALSE
  
  # end backward compatibility
  
  n <- nrow(object$y)
  p <- ncol(object$y)
  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c
  num.RR <- object$num.RR
  
  if((num.lv.c+num.lv)==0&object$randomB==FALSE&is.null(object$params$row.params.random)&is.null(object$randomX)&object$col.eff$col.eff!="random"){
    stop("Cannot calculate prediction errors without random-effects in the model.")
  }
    
  out <- list()
  
  if(object$method == "LA"){
    if(cov){
      if((num.lv+num.RR+num.lv.c)>0) out$lvs <- object$prediction.errors$lvs
      if(!is.null(object$params$row.params.random)) out$row.effects <- lapply(object$prediction.errors$row.params,diag)
      if(object$col.eff$col.eff == "random") {
        out$Br <- object$prediction.errors$col.eff
        row.names(out$Br) <- row.names(object$params$Br)
        colnames(out$Br) <- colnames(object$y)
      }
      if(object$randomB!=FALSE) out$b.lv <- object$prediction.errors$Ab.lv
      if(!is.null(object$randomX)){
        out$Br  <- object$prediction.errors$Br
        row.names(out$Br) <- row.names(object$params$Br)
        colnames(out$Br) <- colnames(object$y)
      }
    } else {
      if((num.lv+num.RR+num.lv.c)>0) out$lvs <- sqrt(apply(object$prediction.errors$lvs,1,diag))
      if(!is.null(object$params$row.params.random)){
        out$row.effects <- object$prediction.errors$row.params
        out$row.effects <- lapply(out$row.effects,function(x)sqrt(diag(x)))
      }
      if(object$col.eff$col.eff=="random"){
      out$Br <- object$prediction.errors$Br
      out$Br <- sqrt(abs(out$Br))
      row.names(out$Br) <- row.names(object$params$Br)
      colnames(out$Br) <- colnames(object$y)
      }
      if(object$randomB!=FALSE) out$b.lv <- sqrt(abs(object$prediction.errors$Ab.lv))
      if(!is.null(object$randomX)){
        out$Br  <- sqrt(object$prediction.errors$Br)
        # out$Br  <- sqrt(apply(object$prediction.errors$Br,1,diag))
      }
    }
  }

  if((object$method %in% c("VA", "EVA"))){
    if(CMSEP) {
      sdb <- CMSEPf(object)

      # sdb<-sdA(object)
      if(object$num.lvcor >0){
        if((object$num.lvcor > 1) && (object$Lambda.struc %in% c("diagU","UNN","UU"))) {
          A<-array(diag(object$A[,,1]), dim = c(nrow(object$A[,,1]), object$num.lvcor,object$num.lvcor))
          for (i in 1:dim(A)[1]) {
            A[i,,]<-A[i,,]*object$AQ
          }
        } else if((object$num.lvcor > 0) & (object$corP$cstruclv !="diag")) {
          A<-array(0, dim = c(nrow(object$A[,,1]), object$num.lvcor,object$num.lvcor))
          if(all(dim(A) == dim(object$A))){
            A<- object$A
          } else {
            for (i in 1:object$num.lvcor) {
              A[,i,i]<- diag(object$A[,,i])
            }
          }
          if(object$num.lvcor==1) A <- matrix(A[,1,1])
        } else if((num.lv.c+num.lv)>0 & num.RR==0){
          A<-object$A
          if((num.lv.c+num.lv)==1) A <- A[,1,1]
        }
        
        if(object$num.lv.c > 0 |object$num.RR > 0){
          # if(NROW(A) != n) {
          if(inherits(object$lvCor,"formula")){
            if(length(dim(A)) <3) {
              object$A <- A <- as.matrix(object$TMBfn$env$data$dLV%*%A)
            } else {
              object$A <- array(0,dim=c(n,dim(A)[2:3]))
              for (k in 1:dim(A)[3]) {
                object$A[,,k] = as.matrix(object$TMBfn$env$data$dLV%*%A[,,k]) # !!!
              }
              A <- object$A
            }
          }
        }
        if(num.RR>0){
          #variational covariances but add 0s for RRR
          A <- array(0,dim=c(n,num.lv.c+num.RR+num.lv,num.lv.c+num.RR+num.lv))
          A[,-c((num.lv.c+1):(num.lv.c+num.RR)),-c((num.lv.c+1):(num.lv.c+num.RR))] <- object$A
        }
      } else if(num.RR>0){
        #variational covariances but add 0s for RRR
        A <- array(0,dim=c(n,num.lv.c+num.RR+num.lv,num.lv.c+num.RR+num.lv))
        A[,-c((num.lv.c+1):(num.lv.c+num.RR)),-c((num.lv.c+1):(num.lv.c+num.RR))] <- object$A
      } else if((num.lv.c+num.lv)>0 & num.RR==0){
        A<-object$A
        if((num.lv.c+num.lv)==1) A <- A[,1,1]
      }
      
      if(!is.null(object$params$row.params.random)){
        for(re in 1:ncol(object$TMBfn$env$data$trmsize))
        object$Ar[[re]]<-diag(sdb$Ar[[re]]+object$Ar[[re]])
      }
      if(object$col.eff$col.eff == "random" | !is.null(object$randomX)){
        if(object$col.eff$Ab.struct %in% c("diagonal", "blockdiagonal")){
          object$Ab <- matrix(diag(sdb$Ab+Matrix::bdiag(object$Ab)), ncol = p)
        }else if(object$col.eff$Ab.struct == "diagonalCL2"){
          # ordering is m independent blocks of p
          object$Ab <- matrix(diag(sdb$Ab+Matrix::bdiag(object$Ab)[order(rep(1:p,times=nrow(object$params$Br))),order(rep(1:p,times=nrow(object$params$Br)))]), ncol = p)
        }else if(object$col.eff$Ab.struct %in% c("unstructured")){
          object$Ab <- matrix(diag(sdb$Ab+object$Ab[[1]]), ncol = p)
        }else if(object$col.eff$Ab.struct %in% c("MNdiagonal", "MNunstructured")){
          # ordering is p blocks of m
          object$Ab <- matrix(diag(sdb$Ab + kronecker(cov2cor(object$Ab[[2]]), object$Ab[[1]])), ncol = p)
        }else if(object$col.eff$Ab.struct %in% c("diagonalCL1", "CL1", "CL2")){
          # ordering is m blocks of size p
          object$Ab <- matrix(diag(sdb$Ab),ncol=p)+matrix(diag(object$Ab),byrow=TRUE,ncol=p)
        }
      }

      if((num.lv+num.lv.c)>0){ object$A<-sdb$A+A} else{object$A <- sdb$A}
      if(num.RR>0&object$randomB!=FALSE){
       covsB <- as.matrix(Matrix::bdiag(lapply(seq(dim(object$Ab.lv)[1]), function(k) object$Ab.lv[k , ,])))
        
        for(i in 1:n){
          Q <- as.matrix(Matrix::bdiag(replicate(num.RR+num.lv.c,object$lv.X.design[i,,drop=F],simplify=F)))
          temp <- Q%*%covsB%*%t(Q) #variances and single dose of covariances
          temp[col(temp)!=row(temp)] <- 2*temp[col(temp)!=row(temp)] ##should be double the covariance
          A[i,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)] <- A[i,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)] + temp
        }
        object$A <- A
        
      }
    }else if(!CMSEP&(num.RR+num.lv.c)>0){
      sdb <- list(Ab_lv = 0)
      
      if(num.RR>0&object$randomB!=FALSE){
        covsB <- as.matrix(Matrix::bdiag(lapply(seq(dim(object$Ab.lv)[1]), function(q) object$Ab.lv[q , ,])))

        for(i in 1:n){
          Q <- as.matrix(Matrix::bdiag(replicate(num.RR+num.lv.c,object$lv.X.design[i,,drop=F],simplify=F)))
          temp <- Q%*%covsB%*%t(Q) #variances and single dose of covariances
          # temp[col(temp)!=row(temp)] <- 2*temp[col(temp)!=row(temp)] ##should be double the covariance
          A[i,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)] <- A[i,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)] + temp
        }
        object$A <- A
        
      }
    }
  
  
  r=0
  if(cov){
    if(!is.null(object$params$row.params.random)){
      # r=1
      out$row.effects <- (object$Ar)
    }
    if(object$col.eff$col.eff=="random"){
      out$Br <- object$Ab
      row.names(out$Br) <- row.names(object$params$Br)
      colnames(out$Br) <- colnames(object$y)
    }

    if(length(dim(object$A))==2){
      out$lvs <- (object$A[,1:(num.lv+num.lv.c+num.RR)+r])
    } else  if((num.lv+num.lv.c+num.RR)>0){
      if((num.lv+num.lv.c+num.RR) ==1) {
        out$lvs <- (as.matrix(object$A[,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.lv.c+num.RR)+r]))
      } else {
        out$lvs <- (object$A[,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.lv.c+num.RR)+r])
      }
    }
    
    if(!is.null(object$randomX)){
      out$Br <- object$Ab
      colnames(out$Br) <- colnames(object$y)
      row.names(out$Br) <- row.names(object$params$Br)
    }
    
    if(object$randomB!=FALSE){
      out$b.lv <- sdb$Ab_lv
      if(object$randomB%in%c("P","iid","single"))out$b.lv <- (abs(out$b.lv + simplify2array(lapply(seq(dim(object$Ab.lv)[1]), function(q) diag(object$Ab.lv[q , ,])))))
      if(object$randomB%in%c("LV"))out$b.lv <- (abs(out$b.lv + t(simplify2array(lapply(seq(dim(object$Ab.lv)[1]), function(q) diag(object$Ab.lv[q , ,]))))))
    }

  } else {
    if(!is.null(object$params$row.params.random)){
      # r=1
      out$row.effects <- object$Ar
      out$row.effects <- lapply(out$row.effects, sqrt)
    }
    if(object$col.eff$col.eff == "random"){
      out$Br <- sqrt(abs(object$Ab))
      row.names(out$Br) <- row.names(object$params$Br)
      colnames(out$Br) <- colnames(object$y)
    }

    if(length(dim(object$A))==2&(num.lv+num.lv.c+num.RR)>0){
      out$lvs <- sqrt(object$A[,1:(num.lv+num.lv.c+num.RR)+r])
    } else if((num.lv+num.lv.c+num.RR)>0){
      if((num.lv+num.lv.c+num.RR) ==1) {
        out$lvs <- sqrt(abs(as.matrix(object$A[,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.lv.c+num.RR)+r])))
      } else {
        out$lvs <- sqrt(abs(apply((object$A[,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.lv.c+num.RR)+r]),1,diag)))
      }
    }
    
    if(!is.null(object$randomX)){
      out$Br <- sqrt(abs(object$Ab))
      colnames(out$Br) <- colnames(object$y)
      row.names(out$Br) <- row.names(object$params$Br)
    }
    
    if(object$randomB!=FALSE){
      out$b.lv <- sdb$Ab_lv
      if(object$randomB%in%c("P","iid","single"))out$b.lv <- sqrt(abs(out$b.lv + simplify2array(lapply(seq(dim(object$Ab.lv)[1]), function(q) diag(object$Ab.lv[q , ,])))))
      if(object$randomB%in%c("LV"))out$b.lv <- sqrt(abs(out$b.lv + t(simplify2array(lapply(seq(dim(object$Ab.lv)[1]), function(q) diag(object$Ab.lv[q , ,]))))))    }
    
  }
  }

  if((num.lv+num.lv.c+num.RR) > 1 & is.matrix(out$lvs)) out$lvs <- t(out$lvs)
  
  return(out)
}

#' Hessian-based CMSEP correction for gllvmHO models
#'
#' Returns the Hessian-correction matrices to be ADDED to the VA variances to
#' give the full conditional MSEP.  Called by \code{getPredictErr.gllvmHO}.
#'
#' @param fit  A fitted \code{gllvmHO} object with \code{$Hess} populated.
#' @return A list with \code{$A} (n x d) and \code{$A_lv} (p x d) correction
#'   matrices (variances in the natural parameterisation, before sigma/alpha
#'   scaling).
#' @keywords internal
CMSEPf_HO <- function(fit) {
  if (is.null(fit$Hess))
    stop("No Hessian stored; refit model or ensure Hessian computation succeeded.")

  n <- nrow(fit$y); p <- ncol(fit$y)
  H    <- fit$Hess$Hess.full
  pnms <- rownames(H)

  incla_u <- pnms == "u"         # n*d_va_z entries
  incla_a <- pnms == "a_lv_sp"   # p*d_va_a entries

  ## VA dims from the actual Hessian block sizes, not from fit$num.lv
  ## (fit$num.lv is unconstrained-only; d_va_z includes constrained dims too)
  d_va_z <- if (any(incla_u)) sum(incla_u) %/% n else 0L
  d_va_a <- if (any(incla_a)) sum(incla_a) %/% p else 0L

  ## In HO both z_i and a_j are random VA effects; the meaningful CMSEPf correction
  ## propagates uncertainty in ONE set of VA means into the OTHER via
  ##   Var_correction(z_i) = H_{uu}^{-1} H_{ua} diag(A_lv) H_{au} H_{uu}^{-1}

  ## ---- CMSEP correction for site scores (n x d_va_z) ----
  if (d_va_z > 0L && d_va_a > 0L) {
    D_u  <- tryCatch(solve(H[incla_u, incla_u, drop = FALSE]),
                     error = function(e) MASS::ginv(H[incla_u, incla_u, drop = FALSE]))
    H_ua <- H[incla_u, incla_a, drop = FALSE]
    A_lv_vec <- as.vector(.ho_diag(fit$B))   # p * d_va_a
    Bw_u  <- H_ua * rep(sqrt(pmax(A_lv_vec, 0)), each = n * d_va_z)
    diag_u <- base::diag(D_u %*% tcrossprod(Bw_u) %*% t(D_u))
    A_sites <- matrix(0, n, d_va_z)
    for (k in seq_len(d_va_z)) A_sites[, k] <- diag_u[(k - 1L)*n + seq_len(n)]
  } else {
    A_sites <- matrix(0, n, max(d_va_z, 1L))
  }

  ## ---- CMSEP correction for species loadings (p x d_va_a) ----
  if (d_va_a > 0L && d_va_z > 0L) {
    D_a  <- tryCatch(solve(H[incla_a, incla_a, drop = FALSE]),
                     error = function(e) MASS::ginv(H[incla_a, incla_a, drop = FALSE]))
    H_au <- H[incla_a, incla_u, drop = FALSE]
    A_u_vec <- as.vector(.ho_diag(fit$A))    # n * d_va_z
    Bw_a  <- H_au * rep(sqrt(pmax(A_u_vec, 0)), each = p * d_va_a)
    diag_a <- base::diag(D_a %*% tcrossprod(Bw_a) %*% t(D_a))
    A_species <- matrix(0, p, d_va_a)
    for (k in seq_len(d_va_a)) A_species[, k] <- diag_a[(k - 1L)*p + seq_len(p)]
  } else {
    A_species <- matrix(0, p, max(d_va_a, 1L))
  }

  list(A = A_sites, A_lv = A_species, d_va_z = d_va_z, d_va_a = d_va_a)
}

#' CMSEP-corrected covariance for b_z (canonical covariate coefficients) in gllvmHO
#'
#' Uses the block-matrix inversion formula
#'   cov(b_z) = D^{-1} + D^{-1} C A_fixed B D^{-1}
#' where D = H_{bz,bz}, B = H_{fixed,bz}, C = H_{bz,fixed} and
#' A_fixed = cov.mat.mod (Schur-complement covariance of fixed effects).
#' D^{-1} is approximated by the block-diagonal Ab.lv matrix
#' (the VA posterior covariance, which is the optimum of the ELBO w.r.t. b_z).
#'
#' @param fit gllvmHO object with $Hess and $Ab.lv populated.
#' @return Full (Kz*d_c) x (Kz*d_c) CMSEP covariance matrix, or NULL if unavailable.
#' @keywords internal
CMSEPf_HO_bz <- function(fit) {
  if (is.null(fit$Hess) || is.null(fit$Ab.lv)) return(NULL)
  H    <- fit$Hess$Hess.full
  pnms <- rownames(H)
  incl <- fit$Hess$incl

  bz_mask <- pnms == "b_z"
  if (!any(bz_mask)) return(NULL)

  d_c <- dim(fit$Ab.lv)[1L]
  Kz  <- dim(fit$Ab.lv)[2L]

  ## VA posterior covariance for b_z as block-diagonal: bdiag(Ab.lv[1,,], ..., Ab.lv[d_c,,])
  ## TMB stores b_z column-major: [b_z[:,1], ..., b_z[:,d_c]], each column has Kz entries.
  ## matrix() wraps the slice to avoid R's drop=TRUE collapsing Kz=1 to a scalar.
  Ab_block <- as.matrix(Matrix::bdiag(lapply(seq_len(d_c), function(l) matrix(fit$Ab.lv[l, , ], Kz, Kz))))

  ## Cross-Hessian blocks between b_z and fixed-effect parameters
  C <- H[bz_mask, incl, drop = FALSE]   # (Kz*d_c) × n_fixed
  B <- H[incl, bz_mask, drop = FALSE]   # n_fixed  × (Kz*d_c)
  A <- fit$Hess$cov.mat.mod             # n_fixed  × n_fixed (cov of fixed effects)

  ## CMSEP correction: D^{-1} C A B D^{-1}  with D^{-1} ≈ Ab_block
  correction <- Ab_block %*% C %*% A %*% B %*% Ab_block
  Ab_block + correction
}

#' CMSEP-corrected covariance for b_gamma (trait coefficients) in gllvmHO
#'
#' Analogous to \code{CMSEPf_HO_bz} but for the b_gamma parameters.
#' @keywords internal
CMSEPf_HO_bgamma <- function(fit) {
  if (is.null(fit$Hess) || is.null(fit$Ab.load)) return(NULL)
  H    <- fit$Hess$Hess.full
  pnms <- rownames(H)
  incl <- fit$Hess$incl

  bg_mask <- pnms == "b_gamma"
  if (!any(bg_mask)) return(NULL)

  d_t <- dim(fit$Ab.load)[1L]
  Kt  <- dim(fit$Ab.load)[2L]

  Ab_block <- as.matrix(Matrix::bdiag(lapply(seq_len(d_t), function(l) matrix(fit$Ab.load[l, , ], Kt, Kt))))

  C <- H[bg_mask, incl, drop = FALSE]
  B <- H[incl, bg_mask, drop = FALSE]
  A <- fit$Hess$cov.mat.mod

  correction <- Ab_block %*% C %*% A %*% B %*% Ab_block
  Ab_block + correction
}

#' Joint CMSEP for (b_z, a_lv_sp) in gllvmHO
#'
#' Returns the Schur-complement CMSEP for the stacked VA parameter vector
#' (b_z; a_lv_sp), including the off-diagonal cross-covariance block
#' Cov_CMSEP(b_z, a_lv_sp).  This is needed for the full bilinear variance
#' formula in randomCoefplot.gllvmHO; the two parameters cannot be corrected
#' separately because the cross-block only appears when both are included
#' jointly in C_va = H[{b_z,a_lv_sp}, fixed].
#'
#' Ordering: b_z rows first (Kz*d_c), then a_lv_sp rows (p*d_va_a).
#' b_z column-major: position of b_z[k,l] = (l-1)*Kz + k.
#' a_lv_sp column-major: position of a_lv_sp[j,ia] = n_bz + (ia-1)*p + j.
#'
#' @return list: $joint (full matrix), $n_bz, $n_asp, $Kz, $d_c, $d_va_a
#' @keywords internal
CMSEPf_HO_bilinear <- function(fit) {
  if (is.null(fit$Hess) || is.null(fit$Ab.lv)) return(NULL)

  n <- nrow(fit$y); p <- ncol(fit$y)
  H    <- fit$Hess$Hess.full
  pnms <- rownames(H)
  incl <- fit$Hess$incl

  bz_mask  <- pnms == "b_z"
  asp_mask <- pnms == "a_lv_sp"
  if (!any(bz_mask)) return(NULL)

  d_c    <- dim(fit$Ab.lv)[1L]
  Kz     <- dim(fit$Ab.lv)[2L]
  d_va_a <- if (any(asp_mask)) sum(asp_mask) %/% p else 0L
  n_bz   <- Kz * d_c
  n_asp  <- p * d_va_a

  ## Block-diagonal VA posterior covariance: b_z first, then a_lv_sp
  Ab_bz <- as.matrix(Matrix::bdiag(lapply(seq_len(d_c), function(l)
    matrix(fit$Ab.lv[l, , ], Kz, Kz))))
  B_asp <- if (d_va_a > 0L && any(asp_mask))
    diag(as.vector(.ho_diag(fit$B)), n_asp)
  else
    matrix(numeric(0), 0L, 0L)
  D_va <- as.matrix(Matrix::bdiag(list(Ab_bz, B_asp)))

  ## Select VA rows from H in the same order as D_va (b_z then a_lv_sp)
  va_rows <- c(which(bz_mask), if (any(asp_mask)) which(asp_mask) else integer(0))
  C_va <- H[va_rows, incl, drop = FALSE]
  A    <- fit$Hess$cov.mat.mod

  joint <- D_va + D_va %*% C_va %*% A %*% t(C_va) %*% D_va
  list(joint = joint, n_bz = n_bz, n_asp = n_asp,
       Kz = Kz, d_c = d_c, d_va_a = d_va_a)
}

#'@export getPredictErr.gllvmHO
#'@method getPredictErr gllvmHO
getPredictErr.gllvmHO <- function(object, CMSEP = TRUE, cov = FALSE, ...) {
  num.RR   <- object$num.RR
  num.lvc  <- object$num.lv.c
  num.lv   <- object$num.lv   # unconstrained
  d_total  <- num.RR + num.lvc + num.lv
  Kz       <- if (!is.null(object$lv.X)) ncol(as.matrix(object$lv.X)) else 0L
  Kt       <- if (!is.null(object$TR))   ncol(as.matrix(object$TR))   else 0L
  n        <- nrow(object$y)
  p        <- ncol(object$y)

  ## VA dims per block.
  ## When Kz > 0, RR z is deterministic (no VA residual); VA covers only lvc+lv dims,
  ## which start at HO col num.RR+1.  When Kz == 0, all dims are VA starting at col 1.
  ## The unified formula ho_idx = seq(d_total - d_va + 1, d_total) handles both cases.
  d_va_z   <- ncol(.ho_diag(object$A))
  d_va_a   <- ncol(.ho_diag(object$B))

  ## Pad VA variances into full d_total HO-order matrices (0 for deterministic dims)
  A_sites_full   <- matrix(0.0, n, d_total)
  A_species_full <- matrix(0.0, p, d_total)
  if (d_va_z > 0L) {
    ho_idx_z <- seq(d_total - d_va_z + 1L, d_total)
    A_sites_full[, ho_idx_z] <- .ho_diag(object$A)
  }
  if (d_va_a > 0L) {
    ho_idx_a <- seq(d_total - d_va_a + 1L, d_total)
    A_species_full[, ho_idx_a] <- .ho_diag(object$B)
  }

  if (CMSEP) {
    if (!is.null(object$Hess)) {
      sdb <- CMSEPf_HO(object)
      if (d_va_z > 0L && sdb$d_va_z > 0L)
        A_sites_full[, ho_idx_z]   <- A_sites_full[, ho_idx_z]   + sdb$A[, seq_len(d_va_z), drop = FALSE]
      if (d_va_a > 0L && sdb$d_va_a > 0L)
        A_species_full[, ho_idx_a] <- A_species_full[, ho_idx_a] + sdb$A_lv[, seq_len(d_va_a), drop = FALSE]
    } else {
      warning("No Hessian in gllvmHO fit; using VA variances only (CMSEP skipped).")
    }
  }

  ## Add b_z VA covariance contribution: Var(z_il) += x_i^T Ab.lv[l,,] x_i
  ## b_z dims follow HO order (RR first, then lvc), so ho_col = l for all l = 1..d_c.
  ## RR dims have zero VA residual but non-zero b_z variance; lvc dims accumulate both.
  d_active <- num.RR + num.lvc
  if (!is.null(object$Ab.lv) && Kz > 0L && d_active > 0L) {
    lv_X_mat <- as.matrix(object$lv.X)
    d_c <- dim(object$Ab.lv)[1L]
    for (l in seq_len(d_c)) {
      cov_l   <- object$Ab.lv[l, , , drop = FALSE]
      dim(cov_l) <- c(Kz, Kz)
      for (i in seq_len(n)) {
        xi <- lv_X_mat[i, , drop = FALSE]
        A_sites_full[i, l] <- A_sites_full[i, l] + as.numeric(xi %*% cov_l %*% t(xi))
      }
    }
  }
  ## Add b_gamma VA covariance contribution: Var(gamma_jl) += t_j^T Ab.load[l,,] t_j
  if (!is.null(object$Ab.load) && Kt > 0L && d_active > 0L) {
    TR_mat <- as.matrix(object$TR)
    d_t <- dim(object$Ab.load)[1L]
    for (l in seq_len(d_t)) {
      cov_l   <- object$Ab.load[l, , , drop = FALSE]
      dim(cov_l) <- c(Kt, Kt)
      for (j in seq_len(p)) {
        tj <- TR_mat[j, , drop = FALSE]
        A_species_full[j, l] <- A_species_full[j, l] + as.numeric(tj %*% cov_l %*% t(tj))
      }
    }
  }

  out <- list()
  if (cov) {
    out$lvs      <- lapply(seq_len(n), function(i) diag(A_sites_full[i, ],   d_total))
    out$loadings <- lapply(seq_len(p), function(j) diag(A_species_full[j, ], d_total))
  } else {
    out$lvs      <- sqrt(A_sites_full)    # n × d_total SDs
    out$loadings <- sqrt(A_species_full)  # p × d_total SDs
  }
  out
}

#'@export getPredictErr
getPredictErr <- function(object, ...)
{
  UseMethod(generic = "getPredictErr")
}