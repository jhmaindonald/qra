% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpsWithin.R
\name{gpsWithin}
\alias{gpsWithin}
\title{Use given vector to identify groups with specified categories}
\usage{
gpsWithin(x, f)
}
\arguments{
\item{x}{One-dimensional object whose values distinguish
groups}

\item{f}{One-dimensional object or list of objects, the
combinations of whose values specify categories within
which groups are to be defined.}
}
\value{
Integer vector whose values, within each specified
category, run from  1 to the number of groups
}
\description{
Any one-dimensional object whose values distinguish groups
may be supplied.
}
\examples{
repnum <- with(qra::codling1988, gpsWithin(cultRep, Cultivar))
table(codling1988$Cultivar,repnum)

}
