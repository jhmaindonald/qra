#' Extract estimates of the intra-class correlation from a glmmTMB
#' model object with beta-binomial error.
#'
#' @description The intra-class correlation is calculated as
#'   \eqn{(1+exp(\theta))^{-1}}, where \eqn{\theta} is the
#'   estimate given by the formula specified in the argument
#'   \code{dispformula}.
#' @details The variance for the betabinomial model is then
#' obtained by multiplying the binomial variance by
#' \eqn{1+(n-1)\rho}, where $n$ is the binomial `size`.
#' @param obj \pkg{glmmTMB} model object with betabinomial error,
#' and with a `dispformula` argument supplied.
#' @param varMult If \code{TRUE} return, in addition to \code{rho},
#' the factor \code{mult} by which the variance is inflated
#' relative to the binomial.
#' @return if \code{varMult==FALSE} return (as a vector) the estimates
#' \eqn{\rho}, else (\code{varMult==TRUE})  return
#' \code{list(rho, mult)}.
#' @examples
#' pcheck <- suppressWarnings(requireNamespace("glmmTMB", quietly = TRUE))
#' if(pcheck) pcheck & packageVersion("glmmTMB") >= "1.1.2"
#' if(pcheck){
#' form <- cbind(Dead,Live)~0+trtGp/TrtTime+(1|trtGpRep)
#' HawMed <- droplevels(subset(HawCon, CN=="MedFly"&LifestageTrt!="Egg"))
#' HawMed <- within(HawMed,
#'                  {trtGp <- factor(paste0(CN,LifestageTrt, sep=":"))
#'                  trtGpRep <- paste0(CN,LifestageTrt,":",RepNumber)
#'                  scTime <- scale(TrtTime) })
#' HawMedbb.TMB <- glmmTMB::glmmTMB(form, dispformula=~trtGp+splines::ns(scTime,2),
#'                                  family=glmmTMB::betabinomial(link="cloglog"),
#'                                  data=HawMed)
#' rho <- qra::getRho(HawMedbb.TMB)} else
#' message("Example requires `glmmTMB` version >= 1.1.2: not available")
#'
#' @importFrom stats model.matrix model.frame
#' @export
getRho <- function(obj, varMult=FALSE){
  if(class(obj)!="glmmTMB"){cl <- class(obj)
    stop(paste("No provision for object of class",cl))
  }
  mm <- model.matrix(obj$modelInfo$allForm$dispformula,
                     data=obj$frame)
  fixdisp <- fixef(obj)[["disp"]]
  rho <- 1/(1+exp(mm%*%fixdisp))
  if(varMult){
    total <- rowSums(model.frame(obj)[["cbind(Dead, Live)"]])
    list(rho=rho, mult=(1+(total-1)*rho))
  } else
    c("rho"=rho)
}

