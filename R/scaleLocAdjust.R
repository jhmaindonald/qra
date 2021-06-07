#' Estimate dispersion as a function of predicted values
#'
#' A loess smooth is applied to the square roots of the standardized
#' deviance residuals. The inverses of values from the smooth, raised
#' to the power of \code{lambda}, are then used as prior weights to
#' update the model.  A value of \code{lambda} that is a little more
#' than 2.0 has often worked well.
#'
#' This function is primarily for experimental use, in investigating
#' possible ways to model a dispersion factor that varies with the
#' fitted value.
#'
#' @param x Model fitted using \code{glm} or, possibly \code{lm}
#' @param lambda Power of smooth of square roots of absolute
#' values of residuals, to try for values whose inverses will
#' be used as weights
#' @param span span parameter for use in smoothing the square
#' root of standardized deviance residuals.
#'
#' @return A list, with elements
#' \item{model}{Model updated to use the newly calculated weights}
#' \item{estDisp}{Estimated dispersions}
#'
#' @note The dispersion estimates that correspond to the updated
#' model are obtained by dividing the dispersion value given
#' by \code{summary()} for the updated model by the (prior) weights
#' supplied when the model was updated. The approach for obtaining
#' varying dispersion estimates is used because, empirically, it
#' has been found to work well for at least some sets of data.  In
#' particular, there seems no obvious theoretical basis for the
#' choice of \code{lambda}.  In the example given, used because the
#' data is publicly available, the method has limited success.
#'
#' @seealso \code{\link{checkDisp}}
#'
#' @importFrom stats influence loess predict residuals
#' update weights deviance df.residual
#'
#' @export
#'
#' @examples
#' ROYAL <- subset(qra::codling1988, Cultivar=="ROYAL")
#' ROYAL.glm <- glm(cbind(dead,total-dead)~ct, data=ROYAL,
#'                   family=quasibinomial(link='cloglog'))
#' ROYALFix <- qra::scaleLocAdjust(ROYAL.glm)
#' ## Check range of indicated prior weights
#' range(ROYALFix[[2]])
#' ## Range of updated dispersion estimates
#' range(summary(ROYALFix[[1]])[['dispersion']]/ROYALFix[[2]])
scaleLocAdjust <- function(x, lambda=2, span=0.75){
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
  sqrtabsfit <- loess(sqrtabsr~hat, weights=w,
                      family="symmetric", span=span)[["fitted"]]
  sqrtabsfit[sqrtabsfit<1] <- 1
  model <- update(x, weights=1/sqrtabsfit^lambda)
  snew <- if (inherits(model, "rlm"))
    model$s
  else if (inherits(model,"glm"))
    sqrt(summary(model)$dispersion)
  else sqrt(deviance(x)/df.residual(x))
  list(model, estDisp=snew^2*sqrtabsfit^lambda)
}
