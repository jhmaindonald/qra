#' qra: A package for calculations that relate to dose-mortality,
#' or time-mortality, or other such models.
#'
#' The qra package provides the functions:
#' checkDisp,	getRho, extractLT, getScaleCoef, scaleLocAdjust,
#' fieller,	gpsWithin, varRatio, foldp,	graphSum, fpower
#'
#' @importFrom stats make.link
#'
#' @import knitr lattice latticeExtra rmarkdown lme4
#'
#' @section qra functions:
#' \code{fieller}: Calculates lethal dose estimates, using Fieller's
#' formula to calculate 95% or other confidence intervals.  Several
#' of the functions noted below are useful ancillaries to \code{fieller}
#' and \code{fieller2}, notably \code{foldp}, \code{fpower}, \code{extractLT},
#' and \code{getScaleCoef}.
#'
#' \code{fieller2}: Use when a folded power link has been used.
#' See \code{\link{fpower}}.
#'
#' \code{extractLT}: Obtains complete set of LT or LD estimate, when it is
#' convenient to get results from several models at the same time.
#'
#' \code{foldp}: Calculates the ratio of \code{p+eps} to \code{1-code+eps}
#'
#' \code{getRho} Extracts estimates of the intra-class
#' correlation from a glmmTMB model object with betabinomial error.
#' See the vignette [timeMortality] for details of the parametization
#' used for the \code{betabinomial} error.
#'
#' \code{getScaleCoef}: Extracts the scale coefficients from a vector
#' that has been scaled using \code{scale}, as needed so that the scaling
#' can be undone.
#'
#' \code{gpsWithin} Renumbers group identifiers so that they run from
#' 1 to number of groups within for each level of the specified factor.
#'
#' \code{scaleLocAdjust}: Returns, for \code{glmmTMB} models with a
#' betabinomial error, dispersion factors (i.e., multipliers for the
#' binomial variance) as functions of predicted values.
#'
#' \code{varRatio}: Returns a first order approximation to the variance
#' of the $y$-ordinate to slope ratio.  This is used in the
#' \code{type="Delta"} approximation, for calculation of LT and LD
#' confidence intervals.  Primarily, this is provided for purposes
#' of comparison, to make it easy to show how poor the approximation
#' can be, and to warn against its general dewvuse!
#'
#' @details
#' Vignettes provide examples of the use of the functions.
#'
#' @keywords internal
"_PACKAGE"


