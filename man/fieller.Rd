% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fieller.R
\name{fieller}
\alias{fieller}
\alias{fieller2}
\title{Confidence Limits for LT Estimate from Dose-response Line}
\usage{
fieller(phat, b, vv, df.t = Inf, offset = 0, logscale = FALSE,
  link = "logit", eps = 0, type = c("Fieller", "Delta"), maxg = 0.99)

fieller2(phat, b, vv, df.t = Inf, offset = 0, logscale = FALSE,
  link = "fpower", lambda = 0, eps = 0, type = c("Fieller", "Delta"),
  maxg = 0.99)
}
\arguments{
\item{phat}{Mortality proportion}

\item{b}{Length 2 vector of intercept and slope}

\item{vv}{Variance-covariance matrix for intercept and slope}

\item{df.t}{Degrees of freedom for variance-covariance
matrix}

\item{offset}{Offset to be added to intercept. This can be of
length 2, in order to return values on the original scale,
in the case where \code{b} and \code{vv} are for values that
have been scaled by subtracting \code{offset[1]} and dividing by
\code{offset[2]}.}

\item{logscale}{Should confidence limits be back transformed
from log scale?}

\item{link}{Link function that transforms expected mortalities
to the scale of the linear predictor}

\item{eps}{If \code{eps>0} \code{phat} is replaced by
\eqn{\frac{p+\epsilon}{1+2*\epsilon}} before applying
the transformation.}

\item{type}{The default is to use Fieller's formula.  A
generally less preferred alternative is to use the Delta
(\code{type="Delta"}) method, which relies on a first
order Taylor series approximation to the variance.}

\item{maxg}{Maximum value of \code{g} for which a
confidence interval will be calculated. Must be \code{< 1}.}

\item{lambda}{The power \eqn{\lambda}, when using the
\code{link="fpower"}. (This applies to \code{fieller2}
only.)}
}
\value{
A vector, with elements
\item{est}{Estimate}
\item{var}{Variance, calculated using the Delta method}
\item{lwr}{Lower bound of confidence interval}
\item{upr}{upper bound of confidence interval}
\item{g}{If \code{g} is close to 0 (perhaps \code{g < 0.05}),
confidence intervals will be similar to those calculated
using the Delta method, and the variance can reasonably
be used for normal theory inference.}
}
\description{
This uses Fieller's formula to calculate a confidence
interval for a specified mortality proportion, commonly
0.50, or 0.90, or 0.99
}
\examples{
redDel <- subset(DAAG::codling,
                 Cultivar=="Red Delicious"&year==1988)
redDel.glm <- glm(pobs~ct, data=redDel,
                  family=quasibinomial(link='cloglog'))
vv <- summary(redDel.glm)$cov.scaled
fieller(0.99, b=coef(redDel.glm), vv=vv, link='cloglog')

}
\references{
Joe Hirschberg & Jenny Lye (2010) A Geometric
Comparison of the Delta and Fieller Confidence Intervals,
The American Statistician, 64:3, 234-241, DOI: 10.1198/ tast.2010.08130

E C Fieller (1944). A Fundamental Formula in the Statistics
of Biological Assay, and Some Applications. Quarterly
Journal of Pharmacy and Pharmacology, 17, 117-123.

David J Finney (1978). Statistical Method in Biological Assay (3rd ed.),
London, Charles Griffin and Company.
}
\seealso{
\code{\link{varRatio}}
}
