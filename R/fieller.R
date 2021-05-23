#' Confidence Limits for Lethal Dose Estimate From Dose-response Line
#'
#' @description This uses Fieller's formula to calculate a confidence
#' interval for a specified mortality proportion, commonly
#' 0.50, or 0.90, or 0.99.  Here "dose" is a generic term for
#' any measure of intensity of a treatment that is designed
#' to induce insect death.
#'
#' @param phat Mortality proportion
#' @param b Length 2 vector of intercept and slope
#' @param vv Variance-covariance matrix for intercept and slope
#' @param df.t Degrees of freedom for variance-covariance
#'    matrix
#' @param offset Offset to be added to intercept. This can be of
#' length 2, in order to return values on the original scale,
#' in the case where \code{b} and \code{vv} are for values that
#' have been scaled by subtracting \code{offset[1]} and dividing by
#' \code{offset[2]}.
#' @param logscale Should confidence limits be back transformed
#' from log scale?
#' @param link Link function that transforms expected mortalities
#'    to the scale of the linear predictor
#' @param lambda The power \eqn{\lambda}, when using the
#'    \code{link="fpower"}. (This applies to \code{fieller2}
#'    only.)
#' @param eps If \code{eps>0} \code{phat} is replaced by
#'    \eqn{\frac{p+\epsilon}{1+2*\epsilon}} before applying
#'    the transformation.
#' @param type The default is to use Fieller's formula.  The
#'    Delta (\code{type="Delta"}) method, which relies on a first
#'    order Taylor series approximation to the variance, is
#'    provided so that it can be used for comparative purposes.
#'    It can be reliably used only in those cases where the
#'    interval is essentially the same as given by
#'    \code{type="Fieller"}.
#' @param maxg Maximum value of \code{g} for which a
#' confidence interval will be calculated. Must be \code{< 1}.
#'
#' @return A vector, with elements
#' \item{est}{Estimate}
#' \item{var}{Variance, calculated using the Delta method}
#' \item{lwr}{Lower bound of confidence interval}
#' \item{upr}{upper bound of confidence interval}
#' \item{g}{If \code{g} is close to 0 (perhaps \code{g < 0.05}),
#' confidence intervals will be similar to those calculated
#' using the Delta method, and the variance can reasonably
#' be used for normal theory inference.}
#'
#' @references Joe Hirschberg & Jenny Lye (2010) A Geometric
#' Comparison of the Delta and Fieller Confidence Intervals,
#' The American Statistician, 64:3, 234-241, DOI: 10.1198/ tast.2010.08130
#'
#' E C Fieller (1944). A Fundamental Formula in the Statistics
#' of Biological Assay, and Some Applications. Quarterly
#' Journal of Pharmacy and Pharmacology, 17, 117-123.
#'
#' David J Finney (1978). Statistical Method in Biological Assay (3rd ed.),
#' London, Charles Griffin and Company.
#'
#' @seealso \code{\link{varRatio}}
#'
#' @examples
#' redDel <- subset(DAAG::codling,
#'                  Cultivar=="Red Delicious"&year==1988)
#' redDel.glm <- glm(pobs~ct, data=redDel,
#'                   family=quasibinomial(link='cloglog'))
#' vv <- summary(redDel.glm)$cov.scaled
#' fieller(0.99, b=coef(redDel.glm), vv=vv, link='cloglog')
#'
#' @rdname fieller
#' @export
fieller <-
  function (phat, b, vv, df.t = Inf, offset = 0, logscale = FALSE,
            link = "logit", eps=0, type=c("Fieller","Delta"), maxg=0.99)
  {
    if(!type[1]%in%c("Fieller","Delta")){
      warning(paste0("Illegal interval type '",type[1],"': assuming Fieller"))
    }
    if(maxg>=1)maxg <- 0.99
    if(length(offset)==1) offset<- c(offset[1],1)
    offset <- as.vector(offset)
    unscale <- function(x, offset, logscale=FALSE){
      x <- x*offset[2]+offset[1]
      if(logscale)exp(x) else x
    }
    v11 <- vv[1, 1]
    v12 <- -vv[1, 2]
    v22 <- vv[2, 2]
    b <- as.vector(b)
    a <- stats::make.link(link)[["linkfun"]]((phat+eps)/(1+2*eps)) - b[1]
    m <- a/b[2]
    tau2 <- v11 - 2 * m * v12 + m^2 * v22
    v <- tau2*offset[2]^2/b[2]^2
    if (df.t == Inf)
      tt <- 1.96
    else tt <- stats::qt(0.975, df.t)
    if(type[1] == "Fieller") g <- (tt/b[2])^2 * v22 else
      if(type[1] == "Delta") g <- 0 else g <- (tt/b[2])^2 * v22
    if (g > maxg || g < 0) {
      m <- unscale(m, offset, logscale=logscale)
      return(c(estval = m, var = v, lower = NA, upper = NA,
                  g = g))
    }
    xhat0 <- (m - g*v12/v22)/(1 - g)
    Ix <- tt/b[2]/(1 - g) * sqrt(tau2 - g * v11 + g/v22 * v12^2)
    Ix <- abs(Ix)
    m <- unscale(m, offset=offset, logscale=logscale)
    lwr <- unscale(xhat0 - Ix, offset=offset, logscale=logscale)
    upr <- unscale(xhat0 + Ix, offset=offset, logscale=logscale)
    c(est = m, var = v, lwr = lwr, upr = upr, g = g)
  }

#' @rdname fieller
#' @export
fieller2 <-
  function (phat, b, vv, df.t = Inf, offset = 0, logscale = FALSE,
            link = "fpower", lambda=0, eps=0, type=c("Fieller","Delta"),
            maxg=0.99)
  {
    if(!type[1]%in%c("Fieller","Delta")){
      warning(paste0("Illegal interval type '",type[1],"': assuming Fieller"))
    }
      if(maxg>=1)maxg <- 0.99
      if(length(offset)==1) offset<- c(offset[1],1)
      offset <- as.vector(offset)
      unscale <- function(x, offset, logscale=FALSE){
        x <- x*offset[2]+offset[1]
        if(logscale)exp(x) else x
      }
      v11 <- vv[1, 1]
      v12 <- -vv[1, 2]
      v22 <- vv[2, 2]
      b <- as.vector(b)
      if(link!="fpower")warning(paste0("Illegal link ", link,
                                ". Using 'fpower'."))
      a <- qra::fpower(phat, lambda, eps) - b[1]
      m <- a/b[2]
      tau2 <- v11 - 2 * m * v12 + m^2 * v22
      v <- tau2*offset[2]^2/b[2]^2
      if (df.t == Inf)
        tt <- 1.96
      else tt <- stats::qt(0.975, df.t)
      if(type[1] == "Fieller") g <- (tt/b[2])^2 * v22 else
        if(type[1] == "Delta") g <- 0 else g <- (tt/b[2])^2 * v22
      if (g > maxg || g < 0) {
        m <- unscale(m, offset, logscale=logscale)
        return(c(estval = m, var = v, lower = NA, upper = NA,
                 g = g))
      }
      xhat0 <- (m - g*v12/v22)/(1 - g)
      Ix <- tt/b[2]/(1 - g) * sqrt(tau2 - g * v11 + g/v22 * v12^2)
      Ix <- abs(Ix)
      m <- unscale(m, offset=offset, logscale=logscale)
      lwr <- unscale(xhat0 - Ix, offset=offset, logscale=logscale)
      upr <- unscale(xhat0 + Ix, offset=offset, logscale=logscale)
      c(est = m, var = v, lwr = lwr, upr = upr, g = g)
  }
