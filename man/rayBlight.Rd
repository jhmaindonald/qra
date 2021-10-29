\name{rayBlight}
\alias{rayBlight}
\docType{data}
\title{
Incidence of ray blight disease of pyrethrum
}
\description{
An assessment of the incidence of ray blight disease of pyrethrum in 62 sampling units, containing 6 plants each.
}
\usage{data("rayBlight")}
\format{
  The format is:
 int [1:62] 4 6 6 6 6 6 6 6 4 6 ...
}
\source{
\pkg{epiphy} package.
}
\references{
Pethybridge SJ, Esker P, Hay F, Wilson C, Nutter FW. 2005. Spatiotemporal description of epidemics caused by Phoma ligulicola in Tasmanian pyrethrum fields. Phytopathology 95, 648-658.
}
\examples{
data(rayBlight)
barplot(table(rayBlight))
}
\keyword{datasets}
