#' First order approximation to variance of y-ordinate to slope ratio
#'
#' In contexts where an LD99 estimate will be used as a data value
#' in a further analysis step, the inverse of the variance may be
#' used as a weight.  The y-ordinate is for the link function
#' transformed value of a specified mortality proportion, commonly
#' 0.50, or 0.90, or 0.99
#'
#' @param phat Mortality proportion
#' @param b Length 2 vector of intercept and slope
#' @param vv Variance-covariance matrix for intercept and slope
#' @param link Link function that transforms expected mortalities
#'    to the scale of the linear predictor
#'
#' @return A vector, with elements
#' \item{xhat}{Estimate}
#' \item{var}{Variance, calculated using the Delta method,  See
#' the help page for \code{\link{varRatio}} for further details
#' and references.}
#'
#' @export
#'
#' @examples
#' redDel <- subset(DAAG::codling,
#'                  Cultivar=="Red Delicious"&year==1988)
#' redDel.glm <- glm(pobs~ct, data=redDel,
#'                   family=quasibinomial(link='cloglog'))
#' vv <- summary(redDel.glm)$cov.scaled
#' varRatio(0.99, b=coef(redDel.glm), vv=vv, link="cloglog")
varRatio <-
  function (phat=0.99, b, vv, link = "cloglog")
  {
    v11 <- vv[1, 1]
    v12 <- -vv[1, 2]
    v22 <- vv[2, 2]
    g.link <- stats::make.link(link)[["linkfun"]]
    b <- as.vector(b)
    a <- g.link(phat) - b[1]
    if(b[2]>0){
      m <- a/b[2]
      v <- (v11 - 2 * m * v12 + m^2 * v22)/b[2]^2
    } else {m <- NA; v <- -1}
    c(xhat = m, var = v)
  }
