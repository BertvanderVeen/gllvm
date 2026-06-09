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

  n <- nrow(fit$y); p <- ncol(fit$y); d <- fit$num.lv
  H    <- fit$Hess$Hess.full
  pnms <- rownames(H)

  ## In HO both z_i and a_j are random VA effects; the truly-fixed params {b, sigmaLV}
  ## have near-zero cross-Hessian with the VA means (verified numerically — the ELBO
  ## decouples them at convergence).  The meaningful CMSEPf correction propagates
  ## uncertainty in ONE set of VA means into the OTHER:
  ##   - site CMSEP : how uncertain a_j estimates make z_i estimates uncertain
  ##   - species CMSEP : how uncertain z_i estimates make a_j estimates uncertain
  ##
  ## Using the implicit-function-theorem derivative
  ##   dz_i / da_j = -H_{uu}^{-1} H_{ua}
  ## and treating a_j posterior variance as diag(A_lv):
  ##   Var_correction(z_i) = D_u * H_{u,a} * diag(A_lv_vec) * H_{a,u} * D_u

  incla_u <- pnms == "u"         # n*d params
  incla_a <- pnms == "a_lv_sp"   # p*d params

  ## ---- CMSEP correction for site scores: uncertainty from a_j estimates ----
  if (any(incla_u) && any(incla_a)) {
    D_u  <- tryCatch(solve(H[incla_u, incla_u, drop = FALSE]),
                     error = function(e) MASS::ginv(H[incla_u, incla_u, drop = FALSE]))
    H_ua <- H[incla_u, incla_a, drop = FALSE]   # (n*d) x (p*d)
    ## Scale columns by sqrt(A_lv) to form D_u * H_ua * diag(A_lv) * H_ua' * D_u
    A_lv_vec <- as.vector(if (!is.null(fit$A_lv_diag)) fit$A_lv_diag else fit$A_lv)  # p*d
    Bw_u <- H_ua * rep(sqrt(pmax(A_lv_vec, 0)), each = sum(incla_u))  # (n*d) x (p*d)
    diag_u <- base::diag(D_u %*% tcrossprod(Bw_u) %*% t(D_u))
    A_sites <- matrix(0, n, d)
    for (k in seq_len(d)) A_sites[, k] <- diag_u[(k - 1L)*n + seq_len(n)]
  } else {
    A_sites <- matrix(0, n, d)
  }

  ## ---- CMSEP correction for species loadings: uncertainty from z_i estimates
  if (any(incla_a) && any(incla_u)) {
    D_a  <- tryCatch(solve(H[incla_a, incla_a, drop = FALSE]),
                     error = function(e) MASS::ginv(H[incla_a, incla_a, drop = FALSE]))
    H_au <- H[incla_a, incla_u, drop = FALSE]   # (p*d) x (n*d)
    A_u_vec <- as.vector(if (!is.null(fit$A_diag)) fit$A_diag else fit$A)  # n*d
    Bw_a <- H_au * rep(sqrt(pmax(A_u_vec, 0)), each = sum(incla_a))  # (p*d) x (n*d)
    diag_a <- base::diag(D_a %*% tcrossprod(Bw_a) %*% t(D_a))
    A_species <- matrix(0, p, d)
    for (k in seq_len(d)) A_species[, k] <- diag_a[(k - 1L)*p + seq_len(p)]
  } else {
    A_species <- matrix(0, p, d)
  }

  list(A = A_sites, A_lv = A_species)
}

#'@export getPredictErr.gllvmHO
#'@method getPredictErr gllvmHO
getPredictErr.gllvmHO <- function(object, CMSEP = TRUE, cov = FALSE, ...) {
  d         <- object$num.lv
  ## Use diagonal summaries regardless of VA.struct (n x d, p x d)
  A_sites   <- if (!is.null(object$A_diag))    object$A_diag    else object$A
  A_species <- if (!is.null(object$A_lv_diag)) object$A_lv_diag else object$A_lv

  if (CMSEP) {
    if (!is.null(object$Hess)) {
      sdb       <- CMSEPf_HO(object)
      A_sites   <- A_sites   + sdb$A
      A_species <- A_species + sdb$A_lv
    } else {
      warning("No Hessian in gllvmHO fit; using VA variances only (CMSEP skipped).")
    }
  }

  out <- list()
  if (cov) {
    out$lvs      <- lapply(seq_len(nrow(A_sites)),
                           function(i) diag(A_sites[i, ],   d))
    out$loadings <- lapply(seq_len(nrow(A_species)),
                           function(j) diag(A_species[j, ], d))
  } else {
    out$lvs      <- sqrt(A_sites)    # n x d  SDs for site scores
    out$loadings <- sqrt(A_species)  # p x d  SDs for species loadings
  }
  out
}

#'@export getPredictErr
getPredictErr <- function(object, ...)
{
  UseMethod(generic = "getPredictErr")
}