## methods for extending emmeans to handle glmmTMB objects

## NOTE: methods are dynamically exported by emmeans utility -- see code in zzz.R

##' Downstream methods
##' 
##' @name downstream_methods
##' @aliases emmeans.gllvm
##' 
##' @description
##' Methods have been written that allow \code{glmmTMB} objects to be used with
##' several downstream packages that enable different forms of inference.
##' For some methods (\code{Anova} and \code{emmeans}, but \emph{not} \code{effects} at present),
##' set the \code{component} argument
##' to "cond" (conditional, the default), "zi" (zero-inflation) or "disp" (dispersion) in order to produce results
##' for the corresponding part of a \code{glmmTMB} model.
##' 
##' In particular,
##' \itemize{
##' \item \code{car::Anova} constructs type-II and type-III Anova tables
##' for the fixed effect parameters of any component
##' \item the \code{emmeans} package computes estimated marginal means (previously known as least-squares means)
##' for the fixed effects of any component
##' \item the \code{effects} package computes graphical tabular effect displays
##' (only for the fixed effects of the conditional component)
##' }
##' @param mod a glmmTMB model
##' @param component which component of the model to test/analyze ("cond", "zi", or "disp")
##' @param \dots Additional parameters that may be supported by the method.
##' @details While the examples below are disabled for earlier versions of
##' R, they may still work; it may be necessary to refer to private
##' versions of methods, e.g. \code{glmmTMB:::Anova.glmmTMB(model, ...)}.
##' @importFrom stats delete.response
##' @examples
##' warp.lm <- glmmTMB(breaks ~ wool * tension, data = warpbreaks)
##' salamander1 <- up2date(readRDS(system.file("example_files","salamander1.rds",package="glmmTMB")))
##' if (require(emmeans)) {
##'     emmeans(warp.lm, poly ~ tension | wool)
##'     emmeans(salamander1, ~ mined, type="response")
##'     emmeans(salamander1, ~ mined, component="zi", type="response")
##' }
##' if (getRversion() >= "3.6.0") {
##'    if (require(car)) {
##'        Anova(warp.lm,type="III")
##'        Anova(salamander1)
##'        Anova(salamander1, component="zi")
##'    }
##'    if (require(effects)) {
##'        plot(allEffects(warp.lm))
##'        plot(allEffects(salamander1))
##'    }
##' }
NULL  ## don't document the files here!


## recover_data method -- DO NOT export -- see zzz.R
## do not document either



recover_data.gllvm<- function(object, ...) {
  fcall <- getCall(object)
  if (!requireNamespace("emmeans"))
    stop("please install (if necessary) and load the emmeans package")
    dat <- data.frame(cbind(object$X,object$lv.X))
    if(!is.null(object$X)){
      formula <- object$formula
      if(inherits(formula,"formula")){
        formula <- deparse(formula)
      }
      if((object$num.lv.c+object$num.RR)>0){
        lv.formula <- object$lv.formula
        formula <- paste(formula,paste(all.vars(lv.formula),collapse="+"),sep="+")
      }
    }else if((object$num.lv.c+object$num.RR)>0){
        formula <- object$lv.formula
    }
    
    trms <- terms(as.formula(formula))
    attr(dat, "predictors") = all.vars(trms)
    attr(dat,"terms") <- trms
    return(dat)
}


##  emm_basis method -- Dynamically exported, see zzz.R
## don't document, causes confusion

## @rdname downstream_methods
## @aliases downstream_methods
## @param component which component of the model to compute emmeans for (conditional ("cond"), zero-inflation ("zi"), or dispersion ("disp"))
## vcov. user-specified covariance matrix
emm_basis.gllvm <- function (object, trms, xlev, grid, ...) {
  ## browser()
  L <- list(...)
  if (length(L)>0) {
    ## don't warn: $misc and $options are always passed through ...
    ## warning("ignored extra arguments to emm_basis.glmmTMB: ",
    ## paste(names(L),collapse=", "))
  }
    V <- vcov(object)

  misc <- list()
    dfargs = list(df = attr(logLik(object),"df"))
    dffun = function(k, dfargs) attr(logLik(object),"df")
  fam <- list()
  fam$family <- object$family
  if(object$family %in% c("poisson", "negative.binomial", "tweedie", "gamma", "exponential"))
    fam$link <- "log"
  if (object$family == "binomial" || object$family == 
      "beta") 
    fam$link <- binomial(link = object$link)$link
  if (object$family == "ordinal") 
    fam$link <- "probit"
  if (object$family == "ZIP") 
    fam4link <-  "other"
  if (object$family == "gaussian") 
    fam$link <- "identity"
  
  misc <- emmeans::.std.link.labels(fam, misc)
  ## (used to populate the reminder of response scale)
  contrasts <- attr(cbind(object$X.design,object$X.lv), "contrasts")
  ## keep only variables found in conditional fixed effects
  contrasts <- contrasts[names(contrasts) %in% c(all.vars(object$terms),all.vars(object$terms.lv))]
  m <- model.frame(trms, grid, na.action=na.pass, xlev=xlev)
  X <- model.matrix(trms, m, contrasts.arg=contrasts)
  bhat <- object$params$Xcoef
  
  if((object$num.RR+object$num.lv.c)>0){
    bhat <- cbind(bhat, object$params$theta %*% t(object$params$LvXcoef))
  }
 
    nbasis <- estimability::all.estble
  

  dfargs <- list(df=attr(logLik(object),"df"))
  dffun <- function(k, dfargs) dfargs$df
  namedList(X, bhat, nbasis, V, dffun, dfargs, misc)
}
namedList <- function (...)  {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(nm <- names(L))) 
    nm <- snm
  if (any(nonames <- nm == "")) 
    nm[nonames] <- snm[nonames]
  setNames(L, nm)
}

stop("doesnt work because vcov for reduced rank + mglm cannot be found")