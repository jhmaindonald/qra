\name{malesINfirst12}
\alias{malesINfirst12}
\docType{data}
\title{
Number of males among first 12 in families of 13 children
}
\description{
The number of male children among the first 12 children of family size 13 in 6115 families taken from the hospital records in the nineteenth century Saxony (Lindsey (1995), p.59). The thirteenth child is ignored to assuage the effect of families non-randomly stopping when a desired gender is reached.
}
\usage{data("malesINfirst12")}
\format{
  A data frame with 13 observations on the following 2 variables.
  \describe{
    \item{\code{No_of_Males}}{a numeric vector}
    \item{\code{freq}}{a numeric vector}
  }
}
\details{
Data are available in the \pkg{fitODBOD} package.
}
\source{
\pkg{fitODBOD} package
}
\references{
Edwards, A. W. F. (1958). An analysis of Geissler's data on the human sex ratio. Annals of human genetics, 23(1), 6-15.

Geissler, A. (1889) Beiträge zur Frage des Geschlechtsverhältnisses der Geborenen. Z. Köngl. Sächs. Statist. Bur., 35,
1±24.

Lindsey, J. K., & Altham, P. M. E. (1998). Analysis of the human sex ratio by using overdispersion models. Journal of the Royal Statistical Society: Series C (Applied Statistics), 47(1), 149-157.
}
\examples{
data(malesINfirst12)
boxplot(freq ~ No_of_Males, data=malesINfirst12)
}
\keyword{datasets}
