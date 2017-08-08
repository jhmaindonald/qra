#' Draw plots f residuals from model fitted to insect mortality data
#'
#' Datasets that are in mind hold, for each replicate of
#' each combination of each of a several factors (e.g.,
#' species, lifestages, temperatures), model residuals for
#' each of a number of values of "dose".
#'
#' @param df Data frame from which data will be taken
#' @param subSet NULL, or an expression, such  as for example
#' \code{expression(LifeStage=='Eggs')}) that evaluates to a logical that
#' specifies the required data subset.  If not NULL then the subsetting
#' information is pasted on after the main title
#' @param logScale Logical, indicating whether the dose ($x$-variable)
#' is on a log scale.
#' @param dosevar Character; column holding "dose" values
#' @param Res Character; name of column holding residuals
#' @param Rep Character; NULL, or column holding replicate number, within panel
#' @param span Span for fitted smooth(s).  To omit smooth(s), set to zero
#' @param byFacet Graphics formula specifying factor combination that
#' determines panel
#' @param layout Graphics formula that can be supplied to \code{grid_facet}
#' @param maint Main title
#' @param pch Plot shape(s), by default 16
#' @param ptSize Pointsize, by default 2
#' @param xzeroOffsetFrac $x$-axis zero offset fraction, required when
#' scale is logarithmic, by default 0.08
#' @param xlab Expression specifying x-axis label
#'
#' @import ggplot2
#'
#' @export
graphResid <- function(df, subSet=NULL, logScale=FALSE,
                     dosevar="logCT", Res="Res", Rep="Rep",
                     span=0.85,
                     byFacet=~Species, layout=NULL,
                     maint="MeBr",
                     pch=16, ptSize=2,
                     xzeroOffsetFrac=0.08,
                     xlab=expression(bold("CT ")*"(gm.h."*m^{-3}*")")){
  ## Prepare data
  checkNames <- c(Res,dosevar) %in% names(df)
  if(any(!checkNames)){
    stop(paste("The following variables were not found:",
               paste(c(Res,dosevar)[!checkNames],collapse=" ")))
  }
  if(!is.null(subSet)){
    txt1 <- as.character(subSet[[1]])
    if(txt1[1]=="==")
      addtxt <- paste0(": ", txt1[2],"=",txt1[3]) else if(txt1[1]=="&")
      {
        txt2 <- as.character(subSet[[1]][[2]])
        txt3 <- as.character(subSet[[1]][[3]])
        addtxt <- paste0(": ", txt2[2],"=",txt2[3], "; ",txt3[2],"=",txt3[3])
      }
  tf <- eval(subSet,df)
  ds <- subset(df, tf)
  } else {
    ds <- df
    addtxt <- ""
  }
  ytxt <- "Residual"
  ran <- range(ds[[Res]])
  dose <- ds[[dosevar]]
  if(!is.null(Rep))ds[[Rep]] <- factor(ds[[Rep]])
  if(logScale){
    zeroDosePos <- min(dose[dose>-Inf])-
      xzeroOffsetFrac*diff(range(dose[dose>-Inf]))
    xtik2 <- round(pretty(log2(exp(dose[dose!=-Inf]))))
    xtik <- c(zeroDosePos,log(2^xtik2))
    xaxlab <- c("CTL", paste(2^xtik2))
    ds[[dosevar]][dose==-Inf] <- zeroDosePos
  } else {
    xtik <- pretty(dose)
    xaxlab <- paste(xtik)
    zeroDosePos <- min(dose[dose>0]) - xzeroOffsetFrac*diff(range(dose))
  }
  xUser <- c(zeroDosePos, max(dose))
  xUser <- xUser+diff(xUser)*c(-0.001,0.001)
  ## Create graph
  x <- rep(zeroDosePos,2)
  y <-ran
  ctlLine <- data.frame("x"=x,"y"=ran)
  if(!is.null(Rep)){
  gg0 <- ggplot(ds)+geom_point(aes_(as.name(dosevar), as.name(Res),
                               shape=~Rep, color=as.name(Rep)),
                               size=ptSize, na.rm=TRUE)
  if(length(pch)==1)pch <- rep(pch,length(unique(ds[[Rep]])))
  }else
  gg0 <- ggplot(ds)+geom_point(aes_(as.name(dosevar), as.name(Res)),
                size=ptSize, na.rm=TRUE)
  gg0 <- gg0 + scale_shape_manual(values=pch) +
    geom_line(aes(x,y), size=1.5*ptSize, alpha=0.4, color="gray",
              data=ctlLine)
if(!is.null(byFacet)) gg0 <- gg0+facet_grid(byFacet)
  if(span>0)
    gg0 <- gg0 + geom_smooth(aes_(as.name(dosevar), as.name(Res),
                                  color=as.name(Rep)),
                             na.rm=TRUE, alpha=0.4, size=0.5, lty=2,
                             span=span, se=FALSE)
gg0 <- gg0+
    xlab(xlab)+ ylab(ytxt)+
    ggtitle(paste0(maint, addtxt))+
    theme(axis.title=element_text(size=13,face="bold"),
          axis.text.x = element_text(color = c("brown",rep("black",length(xtik)-1))),
          plot.title=element_text(size=13)) +
  scale_x_continuous(breaks=xtik, labels=paste(xaxlab))
if(!is.null(layout))gg0 <- gg0 + facet_grid(layout)
gg0
}
##help
