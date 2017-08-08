#' Extract scaling coefficients from vector returned by \code{scale()}
#'
#' The function \code{scale()} replaces \code{x} by \code{(x-a)/b},
#' where \code{a} is \code{mean(x)} and \code{b} is \code{sd(x)}.
#' The quantities \code{a} and \code{b} are available as attributes
#' of the object that is returned.
#'
#' @param z Object returned by \code{scale()}
#'
#' @return A vector, whose elements are the scaling coefficients
#' \code{a} and \code{b}, or if \code{scale=FALSE} then \code{a}.
#'
#' @export
#'
#' @examples
#' z <- scale(1:10)
#' getScaleCoef(z)
getScaleCoef <- function(z){
  sc <- as.vector(unlist(attributes(z)[c("scaled:center","scaled:scale")]))
  if(is.na(sc[2]))sc <- sc[1]
  sc
}
