#' Reproduce data for the linear model scale-location diagnostic plot
#'
#' The values returned are those used for \code{plot(x.lm, which=3)},
#' where \code{x.lm} is a linear model or a generalized linear model.
#' Plot the object returned to assess how successful the weights,
#' determined using the function \code{\link{scaleLocAdjust}}, have been
#' in adjusting for heterogenous variances.
#'
#' @param x Model fitted using \code{lm()} or \code{glm()}
#' @param span span parameter for use in smoothing the square
#' root of standardized deviance residuals.
#'
#' @return A data frame, with:
#' \item{linpred}{Predicted values, on the scale of the linear predictor}
#' \item{absrSmooth}{Smoothed values of the square roots of absolute
#' values of standardised deviance residuals.}
#'
#' @export
#'
#' @examples
#' royal <- subset(qra::codling1988, Cultivar=="ROYAL")
#' royal.glm <- glm(cbind(dead,total-dead)~ct, data=royal,
#'                  family=quasibinomial(link='cloglog'))
#' royalFix <- qra::scaleLocAdjust(royal.glm, lambda=2)
#' ## Check range of indicated prior weights
#' range(royalFix[[2]])
#' ## Range of updated dispersion estimates
#' range(summary(royalFix[[1]])[['dispersion']]/royalFix[[2]])
#' xy <- qra::checkDisp(royalFix[[1]])
#' plot(xy)
checkDisp <- function(x, span=0.75){
  dropInf <- function(x, h) {
    if (any(isInf <- h >= 1)) {
      x[isInf] <- NaN
    }
    x
  }
  r <- residuals(x)
  w <- weights(x)
s <- if (inherits(x, "rlm"))
    x$s
  else if (inherits(x,"glm"))
    sqrt(summary(x)$dispersion)
  else sqrt(deviance(x)/df.residual(x))
  hii <- influence(x, do.coef = FALSE)$hat
  r.w <- if (is.null(w))
    r
  else sqrt(w) * r
  rs <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
  sqrtabsr <- sqrt(abs(rs))
  hat <- predict(x)
  sqrtabsfit <- loess(sqrtabsr~hat, weights=w, family="symmetric", span=span)[["fitted"]]
  data.frame(linpred=hat, absrSmooth=sqrtabsfit)
}

