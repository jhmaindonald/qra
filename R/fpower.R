#' Folded Power Transformation
#'
#' The name \dQuote{folded Power Transformation} is used because
#' this does for power transformations what Tukey's folded logarithm
#' does for the logarithmic tranformation.  The function calculates
#' \deqn{f(p, \lambda, \epsilon) = \frac{p+\epsilon}{1-p+\epsilon}^\lambda}
#' where \eqn{\lambda} is the power and \eqn{\epsilon} is a positive
#' offset that ensures that \eqn{\frac{p+\epsilon}{1-p+\epsilon}} is
#' greater than 0 and finite.
#'
#' @param p Mortality proportion
#' @param lambda Power \code{lambda} for the power transformation
#' @param eps If \code{eps>0} \code{phat} is replaced by
#'    \eqn{\frac{p+\epsilon}{1+\epsilon}} before applying
#'    the power transformation.
#'
#' @return The transformed values of \code{fpower(p)}.
#'
#' @export
#'
#' @examples
#' p <- (0:10)/10
#' ytrans <- fpower(p, lambda=0.25, eps=1/450)
fpower <-
  function (p, lambda, eps)
    if(lambda==0)log((p+eps)/(1-p+eps)) else
 sign(lambda)*((p+eps)/(1-p+eps))^lambda
