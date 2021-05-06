#' Obtain complete set of LT or LD estimates
#'
#' @description When supplied with a model object that has fitted
#' dose-response lines for each of several levels of a factor,
#' `extractLR` calls the function `fieller` to calculate lethal time
#` or lethal dose or other such estimates.
#'
#' @details Fixed coefficients from `obj` must be for intercepts and
#' for slopes.  Starting the model formula with `0+` will commonly
#' do what is required. The coefficients `fixef(obj)[a]` are assumed
#' to specify line intercepts, while `fixef(obj)[b]` specify the
#' corresponding slopes.  These replace the arguments `nEsts`
#' (subscripts for intercepts were `1:nEsts`) and `slopeAdd`
#' (subscripts for slopes were `(nEsts+1):(nEsts+slopeAdd)`).
#'
#' @param p Target mortality proportion.
#' @param obj \code{merMod} object, created using \code{lmer()} or
#' \code{glmerMod} object, created using \code{glmer()}.
#' @param link Link function, specified as a character string.
#' For a folded power function, specify \code{"fpower"}.
#' @param logscale Logical.  Specify \code{TRUE}, if LT values are
#' to be back-transformed from a logarithmic scale.
#' @param a Subscripts for intercepts.
#' @param b Subscripts for corresponding slopes.
#' @param lambda (\code{extractLTpwr} only) Power for power function.
#' @param eps Replace \code{prob} by \code{prob+eps} before transformation.
#' @param offset Use to undo scaling of time or dose variable. This is
#' passed to the \code{\link{fieller}} function that \code{extractLT}
#' calls.
#' @param df.t Degrees of freedom for variance-covariance matrix.
#' If not supplied, it will be calculated internally, using a
#' conservative (low) value.  One is using the \code{t} distribution as
#' an approximation.  Getting a good approximation is not straightforward.
#' See `vignette('introduction-mixed-models', package='afex')`, page 19.
#'
#' @return Matrix holding LD or LD estimates.
#' @export
#'
#' @importFrom stats coef vcov
#' @importFrom lme4 fixef
#
#' @rdname extractLT
#' @export

extractLT <-
  function(p=0.99,
           obj,
           link="logit",
           logscale=FALSE,
           a=1:3,
           b=4:6,
           eps=0,
           offset=0,
           df.t=NULL){
  if(length(a)!=length(b))stop("`a` and `b` must be the same length")
  if(inherits(obj,'lm')){
    if(is.null(df.t))df.t <- summary(obj)$df.residual
    bfun <- coef
    varfun <- vcov
  } else if(class(obj)%in%c("lmerMod","glmerMod")) {
    ngrps <- summary(obj)$ngrps
    bfun <- fixef
    varfun <- vcov
  } else if(class(obj)=="glmmTMB") {
    ngrps <- summary(obj)$ngrps[['cond']]
    bfun <- function(x)fixef(x)[["cond"]]
    varfun <- function(x)vcov(x)[["cond"]]
  } else {
    cl <- class(obj)
    print(paste("No provision for object of class",cl))
    return()
  }
  LT99 <- matrix(0, nrow=length(a), ncol=4)
  pAdj <- (p+eps)/(1+2*eps)
  blmm <- bfun(obj)
  vlmm <- varfun(obj)
  rownames(LT99) <- names(blmm)[1:length(a)]
  if(is.null(df.t))df.t <- ngrps[length(ngrps)]-1
  for(i in 1:length(a)){
    ab <- c(a[i], b[i])
    bi <- blmm[ab]
    vii <- vlmm[ab,ab]
    LT99[i,] <- qra::fieller(pAdj, bi, vii, df.t=df.t, offset=offset,
                             logscale=logscale, link=link, eps=eps)[1:4]
  }
  colnames(LT99) <- c("est", "var", "lwr", "upr")
  LT99
}

#' @rdname extractLT
#' @export

extractLTpwr <-
  function(p=0.99,
                      obj,
                      link="fpower",
                      logscale=FALSE,
                      a=1:3,
                      b=1:3,
                      lambda=0,
                      eps=0.015,
                      offset=0,
                      df.t=NULL){
  if(length(a)!=length(b))stop("`a` and `b` must be the same length")
  LT99 <- matrix(0, nrow=length(a), ncol=4)
  rownames(LT99) <- names(fixef(obj))[1:length(a)]
  ngrps <- summary(obj)$ngrps
  if(is.null(df.t))df.t <- ngrps[length(ngrps)]-1
  for(i in 1:length(a)){
    ab <- c(a[i], b[i])
    blmm <- coef(summary(obj))[ab]
    vlmm <- vcov(obj)[ab,ab]
    LT99[i,] <- qra::fieller2(p, blmm,vlmm, link="fpower",
                              lambda=lambda,eps=eps, df.t=df.t,
                              logscale=logscale, offset=offset)[1:4]
  }
  colnames(LT99) <- c("est", "var", "lwr", "upr")
  LT99
}
