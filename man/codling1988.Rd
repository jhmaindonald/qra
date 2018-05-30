\name{codling1988}
\alias{codling1988}
\docType{data}
\title{Dose-mortality data, for fumigation of codling moth with methyl bromide}
\description{
  Data are from trials that studied the mortality response of codling moth
  to fumigation with methyl bromide, for the year 1988 only
}
\usage{data(codling1988)}
\format{
  A data frame with 77 observations on the following 10 variables.
  \describe{
    \item{dose}{Injected dose of methyl bromide, in gm per cubic meter}
    \item{ct}{Concentration-time sum}    
    \item{total}{Number of insects in chamber}
    \item{dead}{Number of insects dying}
    \item{PropDead}{Proportion dying}
    \item{Cultivar}{a factor with levels \code{BRAEBURN} \code{FUJI} \code{GRANNY} \code{Gala} \code{ROYAL} \code{Red Delicious} \code{Splendour}}
    \item{rep}{replicate number, within \code{Cultivar}}
  }
}
\details{
    The research that generated these data was in part funded by New Zealand
  pipfruit growers.  The published analysis was funded by New Zealand
  pipfruit growers. See also \code{sorption}.
 }
\source{
  Maindonald, J.H.; Waddell, B.C.; Petry, R.J. 2001.
  Apple cultivar effects on codling moth (Lepidoptera: Tortricidae)
  egg mortality following fumigation with methyl bromide.
  Postharvest Biology and Technology 22: 99-110.
}
\keyword{datasets}
