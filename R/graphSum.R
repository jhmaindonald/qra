#' Draw graphs of insect mortality data
#'
#' Datasets that are in mind hold, for each replicate of
#' each combination of each of a several factors (e.g.,
#' species, lifestages, temperatures), mortalities for
#' each of a number of values of "dose".  See for example
#' the dataset help page \code{\link[DAAG]{codling}}.
#'
#' @param df Data frame from which data will be taken
#' @param subSet NULL, or an expression, such  as for example
#' \code{expression(LifeStage=='Eggs')}) that evaluates to a logical that
#' specifies the required data subset.  If not NULL then the subsetting
#' information is pasted on after the main title
#' @param link Link function. If character, obtain from \code{\link{make.link}}.
#'   Alternatively, a function may be supplied as argument.
#' @param logScale Logical, indicating whether the dose ($x$-variable)
#' is on a log scale.
#' @param dead Character; name of column holding number dead
#' @param tot Character; column holding total number
#' @param dosevar Character; column holding "dose" values
#' @param Rep Character; NULL, or column holding replicate number, within panel
#' @param fitRep Character; NULL, or column holding replicate fitted values
#' @param fitPanel Character; NULL, or column holding panel fitted values
#' @param byFacet Graphics formula specifying factor combination that
#' determines panel
#' @param layout Graphics formula that can be supplied to \code{grid_facet}
#' @param maint Main title
#' @param ptSize Pointsize, by default 2
#' @param xzeroOffsetFrac $x$-axis zero offset fraction, required when
#' scale is logarithmic
#' @param yzeroOneOffsets Length two vector, giving 0% mortality and
#' 100% mortality offsets, as fractions of the range for other
#' mortalities, on the scale of the link function.
#' @param yEps Fractional increase at bottom and top of $y$ user range
#' to accommodate points for mortalities of 0 and 1.
#' @param xlab Expression specifying x-axis label
#' @param ylabel If not \code{NULL}, $y$-axis label
#' @param ytiklab Place $y$ axis tiks and labels at these mortalities
#'
#' @import ggplot2
#'
#' @importFrom stats make.link
#'
#' @importFrom lattice xyplot
#'
#' @export
graphSum <- function(df, subSet=NULL,
                     link="cloglog", logScale=FALSE,
                     dead="Dead", tot="Tot", dosevar="logCT", Rep="Rep",
                     fitRep=NULL, fitPanel=NULL,
                     byFacet=~Species, layout=NULL,
                     maint="Codling Moth, MeBr", ptSize=2,
                     xzeroOffsetFrac=0.08,
                     yzeroOneOffsets = c(-0.08, 0.08), yEps=0.005,
                     xlab=expression(bold("CT ")*"(gm.h."*m^{-3}*")"),
                     ylabel=NULL,
                     ytiklab=c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.99)){
  ## Prepare data
  checkNames <- c(dead,tot,dosevar) %in% names(df)
  if(any(!checkNames)){
    stop(paste("The following variables were not found:",
               paste(c(dead,tot,dosevar)[!checkNames],collapse=" ")))
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
  if(is.character(link)){
    linkFun <- make.link(link)$linkfun
    if(is.null(ylabel))ylabel <- paste0("Mortality (", link, " scale)")
    } else if(is.function(link)){
      linkFun <- link
      if(is.null(ylabel))ylabel <- "Mortality"} else
      stop(paste("Invalid link",link))
  rantot <- range(ds[[tot]])
  pobs <- ds[[dead]]/ds[[tot]]
  numCheck <- c(any(pobs<0), any(pobs>1), rantot[1]<=1)
  if(any(numCheck)){
    issues <- c(paste("Some values of", dead, "are negative"),
                paste("Some values of", dead, "are greater than", tot),
                "Some totals are 0 or 1")
    stop(paste(issues[numCheck], collapse = "\n"))
  }
  ran <- range(0.16/rantot[2],linkFun(c(min(pobs[pobs>0]),
                                        (rantot[2]-0.16)/rantot[2])))
  ran <- ran + diff(ran)*yzeroOneOffsets
  yUser <- ran+diff(ran)*c(-yEps,yEps)
  pobsran <- range(pobs[ds[[dead]]>0&ds[[dead]]<ds[[tot]]])
  ytiklab <- ytiklab[ytiklab>=pobsran[1]&ytiklab<=pobsran[2]]
  ytik <- c(ran[1], linkFun(ytiklab),ran[2])
  ytiklab <- c(0,ytiklab,1)
  dose <- ds[[dosevar]]
  if(!is.null(Rep))ds[[Rep]] <- factor(ds[[Rep]])
  y <- linkFun(pobs)
  ds <- within(ds, pointType <- factor(ifelse(y==-Inf,"0",ifelse(y==Inf,"1","0<p<1")),
                       levels=c("1","0<p<1","0")))
  mortTab <- table(ds$pointType)
  plotChar <- c(24,20,25)[(1:3)[mortTab!=0]]
  ds$y <- ifelse(y==-Inf, ran[1], ifelse(y==Inf, ran[2], y))
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
  ctlLine <- data.frame(x=x,y=ran)
  if(!is.null(Rep))
  gg0 <- ggplot(ds)+geom_point(aes_(as.name(dosevar), ~y, shape=~pointType,
                                   color=as.name(Rep)), size=ptSize) else
  gg0 <- ggplot(ds)+geom_point(aes_(as.name(dosevar), ~y, shape=~pointType),
                               size=ptSize)
  gg0 <- gg0 + scale_shape_manual(values=plotChar) +
    geom_line(aes(x,y), size=1.5*ptSize, alpha=0.4, color="gray",
              data=ctlLine)
  if(!is.null(fitRep)){
    if(fitRep%in%names(ds)){
    gg0<-gg0+geom_line(aes_(as.name(dosevar),as.name(fitRep),
                            color=as.name(Rep)), alpha=0.4,na.rm=T)
    } else print(paste("Column",fitRep,"was not found"))
  }
  if(!is.null(fitPanel)){
    if(fitPanel%in%names(ds)){
    gg0 <- gg0+
      geom_line(aes_(as.name(dosevar),as.name(fitPanel)),
                color="gray",size=ptSize, alpha=0.4, na.rm=T)
    } else print(paste("Column",fitPanel,"was not found"))
}
if(!is.null(byFacet)) gg0 <- gg0+facet_grid(byFacet)
gg0 <- gg0+
    xlab(xlab)+ ylab(ylabel)+
    ggtitle(paste0(maint, addtxt))+
    theme(axis.title=element_text(size=13,face="bold"),
          axis.text.x = element_text(color = c("brown",rep("black",length(xtik)-1))),
          plot.title=element_text(size=13)) +
  scale_x_continuous(breaks=xtik, labels=paste(xaxlab))+
  scale_y_continuous(breaks=ytik, labels=paste(ytiklab)) +
  coord_cartesian(ylim=yUser)
if(!is.null(layout))gg0 <- gg0 + facet_grid(layout)
gg0
}

