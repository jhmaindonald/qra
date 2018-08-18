#' Extract estimates of the intra-class correlation from a glmmTMB
#' model object with beta-binomial error.
#'
#' @description The intra-class correlation is calculated as
#'   \eqn{(1+exp(\theta))^{-1}}, where \eqn{\theta} is the
#'   estimate given by the formula specified in the argument
#'   \code{dispformula}.
#' @param obj \pkg{glmmTMB} model object with betabinomial error
#' @return Vector with as many elements as there are observations,
#' holding the estimates of \eqn{\rho}.
#' @examples{
#' codling1988 <- qra::codling1988
#' codling1988$gp <- with(codling1988, paste0(Cultivar,rep))
#' ge16xl.TMB <- glmmTMB::glmmTMB(formula=cbind(dead,total-dead)~0+Cultivar/dose+(1|gp),
#' dispformula=~Cultivar+poly(dose,2), family=glmmTMB::betabinomial(link='logit'),
#' data=subset(codling1988,dose>=16))
#' rho <- getRho(ge16xl.TMB)
#' }
#' @importFrom stats model.matrix
#' @export
getRho <- function(obj){
  mm <- model.matrix(obj$modelInfo$allForm$dispformula,
                     data=obj$frame)
  fixdisp <- fixef(obj)[["disp"]]
  1/(1+exp(mm%*%fixdisp))
}

