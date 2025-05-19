#' Obtain complete set of LT or LD estimates
#'
#' @description When supplied with a model object that has fitted
#' dose-response lines for each of several levels of a factor,
#' \code{extractLT} calls the function \code{fieller} to calculate lethal time
#` or lethal dose or other such estimates.
#'
#' @details Fixed coefficients from \code{obj} must be for intercepts and
#' for slopes.  Starting the model formula with \code{0+} will commonly
#' do what is required. The coefficients \code{fixef(obj)[a]} are assumed
#' to specify line intercepts, while \code{fixef(obj)[b]} specify the
#' corresponding slopes.  These replace the arguments \code{nEsts}
#' (subscripts for intercepts were \code{1:nEsts)} and \code{slopeAdd}
#' (subscripts for slopes were \code{(nEsts+1):(nEsts+slopeAdd)}).
#'
#' @param obj \code{merMod} object, created using \code{lmer()} or
#' \code{glmerMod} object, created using \code{glmer()}.
#' @param a Subscripts for intercepts.
#' @param b Subscripts for corresponding slopes.
#' @param link Link function, for use with objects where no
#' link was specified in the function call, but it is required
#' to back-transform a transformation that was performed prior
#' to the function call.  Otherwise leave as \code{link=NULL}, and
#' the link function will be extracted as \code{family(obj)[['link']]}.
#' For a folded power function, with \code{extractLTpwr()}, the only
#' available link is \code{fpower}, and the exponent \code{lambda} must be
#' specified.
#' @param logscale Logical.  Specify \code{TRUE}, if LT values are
#' to be back-transformed from a logarithmic scale.
#' @param p Target response proportion.
#' @param lambda (\code{extractLTpwr} only) Power for power function.
#' @param eps Replace \code{prob} by \code{prob+eps} before transformation.
#' @param offset Use to undo scaling of time or dose variable. This is
#' passed to the \code{\link{fieller}} function that \code{extractLT}
#' calls.
#' @param df.t Degrees of freedom for a t-distribution approximation
#' for `t` or `z` statistics.  If NULL, a conservative (low) value will
#' be used.  For linear (but not generalized linear) models and mixed
#' models, approximations are implemented in the \pkg{afex} package.
#'  See \code{vignette('introduction-mixed-models', package="afex")}, page 19.
#'
#' @return Matrix holding LD or LD estimates.
#' @examples
#' pcheck <- suppressWarnings(requireNamespace("glmmTMB", quietly = TRUE) &&
#'                            packageVersion("glmmTMB") >= "1.1.2")
#' if(pcheck){
#' form <- cbind(Dead,Live)~0+trtGp/TrtTime+(1|trtGpRep)
#' HawMed <- droplevels(subset(HawCon, CN=="MedFly"&LifestageTrt!="Egg"))
#' HawMed <- within(HawMed,
#'                  {trtGp <- factor(paste0(CN,LifestageTrt, sep=":"))
#'                  trtGpRep <- paste0(CN,LifestageTrt,":",RepNumber)
#'                  scTime <- scale(TrtTime) })
#' HawMedbb.cll <- glmmTMB::glmmTMB(form, dispformula=~trtGp+splines::ns(scTime,2),
#'                                  family=glmmTMB::betabinomial(link="cloglog"),
#'                                  data=HawMed)
#' round(qra::extractLT(p=0.99, obj=HawMedbb.cll, link="cloglog",
#'                a=1:3, b=4:6, eps=0, df.t=NULL)[,-2], 2)} else
#' message("Example requires `glmmTMB` version >= 1.1.2: not available")
#'
#' @export
#'
#' @importFrom stats coef vcov family
#' @importFrom lme4 fixef
#' @importFrom splines ns
#
#' @rdname extractLT
#' @export

extractLT <-
  function(obj,
           a=1:3,
           b=4:6,
           link=NULL,
           logscale=FALSE,
           p=0.99,
           eps=0,
           offset=0,
           df.t=NULL){
  if(is.null(link))link <- family(obj)[['link']]
  if(inherits(obj,'lm')){
    if(is.null(df.t))df.t <- summary(obj)$df.residual
    bfun <- coef
    varfun <- vcov
  } else if(inherits(obj, "merMod")) {
    ngrps <- summary(obj)$ngrps
    bfun <- fixef
    varfun <- vcov
  } else if (inherits(obj, "glmmTMB")) {
    ngrps <- summary(obj)$ngrps[['cond']]
    bfun <- function(x)fixef(x)[["cond"]]
    varfun <- function(x)vcov(x)[["cond"]]
  } else {
    cl <- class(obj)
    stop(paste("No provision for object of class",cl))
  }
  pAdj <- (p+eps)/(1+2*eps)
  blmm <- bfun(obj)
  vlmm <- varfun(obj)
  if(any(is.na(vlmm)))stop("Variance covariance matrix contains NAs")
  check <- ifelse(length(a)!=length(b), 1,
                  ifelse(length(unique(b-a))!=1, 2,
                         ifelse(2*unique(b-a)!=length(blmm),2, 0)))
  if(check > 0){
    mess <- c("Lengths of a and b differ", "Intercepts do not match slopes")
    stop(mess[check])
  }
  if(is.null(df.t))df.t <- as.numeric(ngrps[length(ngrps)]-1)
  xCentile <- matrix(0, nrow=length(a), ncol=4,
                   dimnames=list(names(blmm)[1:length(a)],
                                 c("est", "var", "lwr", "upr")))
  attr(xCentile, 'Centile')  <- 100*p
  attr(xCentile, 'df.t')  <- df.t
  for(i in 1:length(a)){
    ab <- c(a[i], b[i])
    bi <- blmm[ab]
    vii <- vlmm[ab,ab]
    xCentile[i,] <- qra::fieller(pAdj, bi, vii, df.t=df.t, offset=offset,
                             logscale=logscale, link=link, eps=eps)[1:4]
  }
  xCentile
}

#' @rdname extractLT
#' @export

extractLTpwr <-
  function(obj,
           a=1:3,
           b=1:3,
           link="fpower",
           logscale=FALSE,
           p=0.99,
                      lambda=0,
                      eps=0.015,
                      offset=0,
                      df.t=NULL){
  check <- ifelse(length(a)!=length(b), 1,
                    ifelse(length(unique(b-a))!=1, 2,
                           ifelse(2*length(unique(b-a))!=length(blmm),2, 0)))
  if(check > 0){
      mess <- c("Lengths of a and b differ", "Intercepts do not match slopes")
      stop(mess[check])
  }
  blmm <- coef(summary(obj))
  vlmm <- vcov(obj)
  ngrps <- summary(obj)$ngrps
  if(is.null(df.t))df.t <- as.numeric(ngrps[length(ngrps)]-length(blmm))
  xCentile <- matrix(0, nrow=length(a), ncol=4,
                     dimnames=list(names(blmm)[1:length(a)],
                                   c("est", "var", "lwr", "upr")))
  attr(xCentile, 'Centile')  <- 100*p
  attr(xCentile, 'df.t')  <- df.t
  for(i in 1:length(a)){
    ab <- c(a[i], b[i])
    bi <- blmm[ab]
    vii <- vlmm[ab,ab]
    xCentile[i,] <- qra::fieller2(p, bi, vii, link="fpower",
                              lambda=lambda,eps=eps, df.t=df.t,
                              logscale=logscale, offset=offset)[1:4]
  }
  xCentile
  }

