
#' Title Function to calculate ratio of \code{p+eps} to \code{1-p+eps}.
#'
#' @description
#' This is a convenience function that returns
#' \eqn{\frac{p+\epsilon}{1-p+\epsilon}}.  It calculates
#' the argument that is supplied to the \code{log}
#' function in Tukey's \sQuote{flog}.
#'
#' @param p Proportion
#' @param eps Offset. The choice \code{eps}=0.01 has the
#' same effect as replacing \eqn{\frac{r}{n-r}} by
#' \eqn{\frac{r+0.5}{n-r+0.5}} when \eqn{n=50}, or by
#' \eqn{\frac{r+1}{n-r+1}} when \eqn{n=100}
#'
#' @return \code{(p+eps)/(1-p+eps)}
#' @export
#'
#' @examples
#' foldp(c(0.2,0.75), 0)
foldp <-
function(p,eps){
  (p+eps)/(1-p+eps)
}
