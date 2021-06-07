#' First order approximation to variance of y-ordinate to slope ratio
#'
#' In contexts where an LD99 estimate will be used as a data value
#' in a further analysis step, the inverse of the variance may be
#' used as a weight.  The y-ordinate is for the link function
#' transformed value of a specified mortality proportion, commonly
#' 0.50, or 0.90, or 0.99
#'
#' This function should only be used, in order to speed up
#' calculations that use the function \code{\link{fieller}}
#' (call \code{fieller} with (\code{type="Delta"})),
#' in a context where it is to be used many times,
#' and where a check has been made that its use leads to
#' confidence intervals that are a close approximation to those
#' given with the default argument (\code{type="Fieller"}).
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
#' the help page for \code{\link{fieller}} for further details
#' and references.}
#'
#' @export
#'
#' @examples
#' redDel <- subset(qra::codling1988, Cultivar=="Red Delicious")
#' redDel.glm <- glm(cbind(dead,total-dead)~ct, data=redDel,
#'                   family=quasibinomial(link='cloglog'))
#' vv <- summary(redDel.glm)$cov.scaled
#' qra::varRatio(0.99, b=coef(redDel.glm), vv=vv, link="cloglog")
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
