#' Obtain complete set of LT or LD estimates
#'
#' @param p Desired mortality proportion.
#' @param obj \code{merMod} object, created using \code{lmer()} or
#' \code{glmer()}.
#' @param link Link function, specified as a character string.
#' For a folded power function, specify \code{"fpower"}.
#' @param logscale Logical.  Specify \code{TRUE}, if LT values are
#' to be back-transformed from a logarithmic scale.
#' @param nEsts Number of groups. The first \code{nEsts} estimates
#' are assumed to be intercepts.
#' @param slopeAdd Add to the intercept number to get corresponding slope.
#' @param lambda (\code{extractLTpwr} only) Power for power function.
#' @param eps Replace \code{prob} by \code{prob+eps} before transformation.
#' @param scaling Use to undo scaling of time or dose variable.
#' @param df.t Degrees of freedom for variance-covariance matrix.
#' If not supplied, it will be calculated internally.
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
           nEsts=8,
           slopeAdd=8,
           eps=0,
           scaling=1,
           df.t=NULL){
  LT99 <- matrix(0, nrow=nEsts, ncol=4)
  pAdj <- (p+eps)/(1+2*eps)
  rownames(LT99) <- names(fixef(obj))[1:nEsts]
  if(is.null(df.t))df.t <- summary(obj)$ngrps-nEsts
  for(i in 1:nEsts){
    ab <- c(i, slopeAdd+i)
    blmm <- coef(summary(obj))[ab]
    vlmm <- vcov(obj)[ab,ab]
    LT99[i,] <- qra::fieller(pAdj, blmm,vlmm, df.t=16, offset=scaling,
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
                      nEsts=8,
                      slopeAdd=8,
                      lambda=0,
                      eps=0.015,
                      scaling=1,
                      df.t=NULL){
  LT99 <- matrix(0, nrow=nEsts, ncol=4)

  rownames(LT99) <- names(fixef(obj))[1:nEsts]
  if(is.null(df.t))df.t <- summary(obj)$ngrps-nEsts-1
  for(i in 1:nEsts){
    ab <- c(i, slopeAdd+i)
    blmm <- coef(summary(obj))[ab]
    vlmm <- vcov(obj)[ab,ab]
    LT99[i,] <- qra::fieller2(p, blmm,vlmm, link="fpower",
                              lambda=lambda,eps=eps, df.t=df.t,
                              logscale=logscale, offset=scaling)[1:4]
  }
  colnames(LT99) <- c("est", "var", "lwr", "upr")
  LT99
}
