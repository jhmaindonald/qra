#' Use given vector to identify groups with specified categories
#'
#' Any one-dimensional object whose values distinguish groups
#' may be supplied.
#'
#' @param x One-dimensional object whose values distinguish
#' groups
#' @param f One-dimensional object or list of objects, the
#' combinations of whose values specify categories within
#' which groups are to be defined.
#'
#' @return Integer vector whose values, within each specified
#' category, run from  1 to the number of groups
#'
#' @export
#'
#' @examples
#' repnum <- with(qra::codling1988, gpsWithin(cultRep, Cultivar))
#' table(codling1988$Cultivar,repnum)
#'
gpsWithin <- function(x,f){
  z <- lapply(split(x,f), function(x)match(x,unique(x)))
  unsplit(z,f)}
