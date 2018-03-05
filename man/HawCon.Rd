\name{HawCon}
\alias{HawCon}
\docType{data}
\title{
Hawaiian Contemporary Cold Treatment Dataset
}
\description{
The counts of live/dead were derived by injecting a known number of individuals of the target life stage into citrus fruits, subjecting them to treatment and then counting the number of individuals emerging.
}
\usage{data("HawCon")}
\format{
  A data frame with 106 observations on the following 10 variables.
  \describe{
    \item{\code{Species}}{Species of fruitfly}
    \item{\code{CN}}{Common name, in abbreviated form.
    MedFly is \sQuote{Mediterranean Fruit Fly}.  MelonFly is
    \sQuote{Melon Fly}}
    \item{\code{LifestageTrt}}{Lifestage treated}
    \item{\code{RepNumber}}{Replicate number}
    \item{\code{PropDead}}{Fraction dead}
    \item{\code{TrtTime}}{Treatment time (days)}
    \item{\code{Dead}}{a numeric vector}
    \item{\code{Live}}{a numeric vector}
    \item{\code{Total}}{a numeric vector}
  }
}
\details{
The help page for \code{HawCon} in the \pkg{ColdData} has further
details.
}
\source{
Dr Peter Follett
}
\references{
A paper is in the course of preparation.
}
\examples{
data(HawCon)
str(HawCon)
}
\keyword{datasets}
