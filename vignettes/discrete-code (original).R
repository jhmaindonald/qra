## ---- set-lattice ----
library(lattice)
sides <- list(tck = 0.6, pad1 = 0.75, pad2 = 0.75)
theme <- list(axis.line = list(alpha = 1, col = 'gray40',
                               fill = "transparent", lty = 1, lwd = 0.5),
              strip.border = list(alpha = 1, col = rep('gray40', 6),
                                  lty = rep(1, 6), lwd = rep(0.5, 6)),
              strip.shingle = list(alpha = 1, col = rep("gray80", 7)),
              box.3d = list(col = 'gray40'),
              axis.components = list(left = sides, top = sides,
                                     right = sides, bottom = sides),
              fontsize = list(text = 10, points = 6))
trellis.par.set(theme)

## ---- binom-dp ----
binom <- rbind(
paste(dbinom(x=0:3, size=3, prob=0.5)),
paste(pbinom(q=0:3, size=3, prob=0.5)))
colnames(binom) <- paste(0:3)
rownames(binom) <- c("dbinom(x=0:3, size=3, prob=0.5)",
                     "pbinom(q=0:3, size=3, prob=0.5)")
print(binom, quote=FALSE)

## ---- binom-q ----
binomq <- matrix(qbinom(p=c(0.125, 0.5, 0.875, 1.0),
                        size=3, prob=0.5), nrow=1)
rownames(binomq) <- "qbinom(p=c(0.125,0.5,0.875,1.0), size=3, prob=0.5)"
colnames(binomq) <- paste(c(0.125, 0.5, 0.875, 1.0))
print(binomq, quote=FALSE)

## ---- binom-gph1 ----
library(latticeExtra)
p <- c(pbinom(q=0:3, size=3, prob=0.5))  # Quantile change points
p1 <- c(0, p[1], NA, p[1]+0.001, p[2], NA, p[2]+0.001,
        p[3], NA, p[3]+0.001,p[4])
gph1 <- xyplot(p1 ~ qbinom(p1, size=3, prob=0.5), type='l',
               xlim=c(-0.01,3.01), ylim=c(0,1.003),
               xlab="Number of heads, in 3 tosses",
               scales=list(x=list(at=0:3), y=list(at=(0:5)*0.2)),
               ylab="Cumulative probability",
               main = list("A: (Cumulative) distribution function",
                           font=1, x=0, y=0.25, just="left", cex=1.25)) +
  latticeExtra::layer(panel.xyplot(x[c(2,5,8,11)], y[c(2,5,8,11)], pch=16))
gph1 <- update(gph1, axis=latticeExtra::axis.grid)

## ---- binom-gph2 ----
p <- c(pbinom(q=0:3, size=3, prob=0.5))  # Quantile change points
p1 <- c(0, p[1], NA, p[1]+0.001, p[2], NA, p[2]+0.001,
        p[3], NA, p[3]+0.001,p[4])
gph2 <- xyplot(qbinom(p1, size=3, prob=0.5) ~ p1, type='l',
               xlim=c(0,1.003), ylim=c(-0.01,3.01),
               scales=list(x=list(at=(0:5)*0.2),y=list(at=0:3)),
               xlab="(Cumulative) probability", ylab="Quantile",
               main = list("B: Quantile function", font=1, x=0, y=0.25,
                           just="left", cex=1.25))+
  latticeExtra::layer(panel.xyplot(x[c(2,5,8,11)], y[c(2,5,8,11)], pch=16))
#    latticeExtra::layer(panel.xyplot(x,y,type="S",lwd=2))
gph2 <- update(gph2, axis=latticeExtra::axis.grid)

## ---- quantile ----
probrange <- t(cbind("x = Number of heads"=0:3,
               "Interval"=
c("0, 0.125", ">0.125, 0.5", ">0.5, 0.875", ">0.875, 1.0")))
xtab <- xtable::xtable(probrange,align=c('l',rep('r',4)))
colnames(xtab) <- rep("",4)
print(xtab, type='html', hline.after=NULL)

## ---- u2q ----
round(qnorm(c(0.2,0.5,0.8)), 2)

## ---- gamlss ----
suppressPackageStartupMessages(library(gamlss))

## ---- cfDBI-BB ----
par(mar=c(2.6,4.1,2.6,1.6), mgp=c(2.5,.75,0), mfrow=c(1,2))
x <- 0:10
cols <- c('gray80',"deepskyblue","deepskyblue4","red")
dDBImax <- dDBI(5, mu=.5, sigma=0.7, bd=10)*1.025
plot(x, dBI(x, mu=.5, bd=10), type="h", lwd=16, col='gray80',
     xlab=paste("Number out of ",max(x)),
     ylab="Probability", lend=1, ylim=c(0, dDBImax*1.02), fg='gray')
lines(x-0.1, dDBI(x, mu=.5, sigma=0.702, bd=10), type="h",
      col=cols[2], lwd=2)
lines(x, dDBI(x, mu=.5, sigma=2, bd=10), type="h",
      col=cols[3],lwd=2)
# DBB
lines(x+0.1, dBB(x, mu=.5, sigma=.125, bd=10), type="h",
      col="red", lwd=2)
legend(x="topright", legend=c(expression("Binom","DBI: "*Phi==0.7,
                                         "DBI: "*Phi==2,"BB: "*Phi==2)),
       box.col='gray',
       col=cols, cex=0.85, lwd=c(6,2,2,2), border="gray40", bty="l")
mtext(side=3, line=0.5, adj=0, cex=1.25, 'A: Binomial vs alternatives')
plot(x, dBI(x, mu=.5, bd=10), type="h", lwd=16, ylim=c(0,dDBImax*1.02),
     xlab=paste("Number out of ",max(x)),
     ylab="Probability", col='gray80', lend=1, fg='gray')
lines(x, dDBI(x, mu=.5, sigma=8.27, bd=10), type="h",
      col="deepskyblue4",lwd=2)
# DBB
lines(x+0.1, dBB(x, mu=.5, sigma=5/7, bd=10), type="h",
      col=cols[4], lwd=2)
legend(x="topright",
       legend=c(expression("Binom", "DBI: "*Phi==4.75,
                           "BB: "*Phi==4.75)), col=cols[-2], cex=0.85,
                           lwd=c(8,2,2,2), border="gray40", bty="l")
mtext(side=3, line=0.5, adj=0, cex=1.25,
      'B: Binomial vs further alternatives')

## ---- wtd-var ----
wtd.var <- function(number,w){
  N <- sum(w)
  mu <- weighted.mean(number, w)
  var <- sum(w*(number-mu)^2)/sum(w)
  c(mean=mu,var=var, N=N)
}

## ---- fitVSobs ----
fitVSobs <- function(Number="number", Freq="freq", Fitprob="fitprob",
                     df=rugeFit, digits=3, size=10, rho=NULL,
                     distribution='poisson',
                     xlab="Count", ylab="Frequency",
                     atright=list(at=(1:5)*10, labels=(1:5)/20),
                     cex.title=1.1,
                     main=list("A: Radioactive decay counts",
                               font=1, x=0, y=0.25, just="left", cex=1.05)){
  number <- df[[Number]]
  freq <- df[[Freq]]
  fitted <- df[[Fitprob]]*sum(freq)
  mu <- weighted.mean(number,freq)
  var <- wtd.var(number,freq)['var']
  if(distribution == "binomial"){
    fitvar <- mu*(1-mu/size)
    titl <- paste("Data: Var =", round(var,digits),
                  "  Binomial fit: Var =", round(fitvar,digits),collapse="")
    main[[1]] <- paste0(main[[1]]," (Prob = ", round(mu/size,3),")", collapse="")
  } else
    if(distribution == "betabinomial"){
      fitvar <- mu*(1-mu/size)*(1+(size-1)*rho)
      titl <- paste("Data: Var =", round(var,digits),
                    "  Beta binomial fit: Var =", round(fitvar,digits))
      main[[1]] <- paste0(main[[1]]," (Prob = ", round(mu/size,3),")",
                          "rho =", round(rho,3),")", collapse="")
    } else if (distribution == "poisson"){ fitvar <- mu
    titl <- paste("Data: Var =", round(var,digits),
                  "  Poisson fit: Var =", round(fitvar,digits))
    }
  xlab <- paste(xlab, "(Mean = ", round(mu,2), ")")
  legend <- list(c("Data:", "Fitted:"))
  bcol <- trellis.par.get()$plot.polygon$col
  # if(is.null(atright))ylab.rt <- "" else
  #   ylab.rt <- "Relative frequency"
  gph <- barchart(freq~as.factor(number), data = df, box.ratio=3,
                  main=main,
                  par.settings=list(layout.widths=list(ylab.right=-1.5,
                                                       right.padding=1.5)),
                  ylim=c(0, 1.04*max(freq,fitted)),
                  xlab=xlab, ylab=ylab, ylab.right="",
                  scales=list(x=list(alternating=1),
                              y=list(tck=c(0.6,0))),
                  border='gray60', horiz=FALSE,
                  rectangles=TRUE,
                  key=list(space='top', columns=2, text=legend,
                           title=titl,cex.title=cex.title,
                           lines=list(lwd=c(10,2), col=c(bcol,1))))+
    latticeExtra::as.layer(xyplot(fitted~as.factor(number),
                                  type=c('h','g'), data = df, col='black',lwd=4))
  if(!is.null(atright))gph <- gph +
    latticeExtra::layer(panel.axis(side='right',at=at,labels=labels,
                                   half=FALSE), data=atright)
  gph
}

## ---- binData ----
htab <- table(factor(
  sapply(split(testDriveR::kerrich[,2],
               rep(1:200,rep(10,200))),sum), levels=0:10),
  dnn=list("A: Frequency of each of 0 to 10, in 10 tosses"))
heads <- data.frame(number=0:10, freq=as.numeric(htab))
muH <- with(heads, weighted.mean(number, w=freq))
heads$fitprob <- with(heads, dbinom(number, size=10, prob = muH/10))
tastab <- table(factor(epiphy::pyrethrum_ray_blight$i, levels=0:6),
                dnn=list("B: Frequency of each of 0 to 6, in 6 plants"))
diseased <- data.frame(
  number=0:6,
  freq=as.numeric(tastab)
)
mu <- with(diseased, weighted.mean(number, freq))
diseased$Binprob <- dbinom(0:6, size=6, prob=mu/6)

## ---- binAB ----
library(lattice)
gph1 <- fitVSobs(Number="number", Freq="freq", Fitprob="fitprob", df=heads,
                 digits=2, distribution="binomial", atright=list(at=(1:5)*10,
                                                                 labels=(1:5)/20),
                 main=list("A: No. of heads (10 tosses)",
                           cex=1.05, font=1, x=0, just="left"))
gph2 <- fitVSobs(Number="number", Freq="freq", Fitprob="Binprob", df=diseased,
                 digits=2, size=6, distribution="binomial",
                 atright=NULL, ylab="Frequency\nRelative Frequency",
                 main=list("B: No. of diseased plants (from 6)",
                           cex=1.05, font=1, x=0, just="left"))
gph2 <- gph2+latticeExtra::layer(panel.axis(side='left', at=(1:5)*10*62/50,
                                            labels=(1:5)*10/50, outside=FALSE, half=FALSE))
## Binomial variance:  5.451613*(1-.914) = 0.4688387
## BB fitted var = 0.4688387*2.613
## Sample variance = 1.169751 = 0.4688387*2.495

## ---- maleFit ----
muM <- with(fitODBOD::Male_Children,
            weighted.mean(No_of_Males, w=freq))
maleFit <- with(fitODBOD::Male_Children,
                data.frame(number=No_of_Males, freq=freq,
                           fitprob = dbinom(No_of_Males, size=12, prob = muM/12)))

## ---- xtra ----
## Not executed
plot(x, dBI(x, mu=.5, bd=10), type="h", lwd=16,
     xlab=paste("Number out of ",max(x)),
     ylab="Probability", col='gray80', lend=1)
lines(x, dDBI(x, mu=.5, sigma=8.27, bd=10), type="h",
      col="deepskyblue4",lwd=2)
# DBB
lines(x+0.1, dBB(x, mu=.5, sigma=0.714, bd=10), type="h",
      col="red", lwd=2)
legend(x="topright",
       legend=c(expression("Binom","DBI: "*phi==4.25,"BB: "*phi==4.25)),
       col=cols[-2], cex=0.85, lwd=c(8,2,2,2), border="gray40", bty="l")

## ---- cfFits ----
library(gamlss)
doBI <- gamlss(cbind(number, 6-number)~1, weights=freq,
               family=BI, data=diseased, trace=FALSE)
doBB <- gamlss(cbind(number, 6-number)~1, weights=freq,
               family=BB, data=diseased, trace=FALSE)
doDBI <- gamlss(cbind(number, 6-number)~1, weights=freq,
                family=DBI, data=diseased, n.cyc=100, trace=FALSE)

## ---- cfsim ----
par(mar=c(2.6,4.1,2.6,1.6), mgp=c(2.5,.75,0), mfrow=c(1,3))
## Binomial fits
rqres.plot(doBI, plot.type='all', type="QQ")
mtext(side=3, line=0.5, "A: Quantile-quantile", adj=0, cex=1.0)
res <- rqres.plot(doBI, plot.type='all', type="wp")
mtext(side=3, line=0.5, "B: Worm plot", adj=0, cex=1.0)
plot(density(res), xlab="Quantile residuals", main="")
mtext(side=3, line=0.5, "C: Density plot", adj=0, cex=1.0)

## ---- cfq2 ----
par(mar=c(2.6,4.1,2.6,1.6), mgp=c(2.5,.75,0), mfrow=c(1,3))
## Betabinomial fits
rqres.plot(doBB, plot.type='all', type="QQ")
mtext(side=3, line=0.5, "A: Quantile-quantile", adj=0, cex=0.85)
res <- rqres.plot(doBB, plot.type='all', type="wp")
mtext(side=3, line=0.5, "B: Worm plot", adj=0, cex=0.85)
plot(density(res), xlab="Quantile residuals", main="")
mtext(side=3, line=0.5, "C: Density plot", adj=0, cex=0.85)

## ---- cf-AIC ----
aicStat <- bbmle::AICtab(doBI, doBB, doDBI, mnames=
                           c("Binomial","Betabinomial","Double binomial"), base=TRUE)
print(aicStat, type='html', hline.after=NULL, right=TRUE)

## ---- CBin ----
## Not used
doCBin <- bbmle::mle2(fitODBOD::EstMLECorrBin, start=list(p=.8,cov=.05),
                      data = list(x=diseased$number,freq=diseased$freq))

## ---- Saxonymales ----
maleDF <- with(fitODBOD::Male_Children,
               data.frame(number=No_of_Males, freq=freq))
## Numbers of male children, out of 12
as.table(setNames(maleDF$freq, nm=0:12))

## ---- maleFit ----
doBI <- gamlss(cbind(number, 12-number)~1, weights=freq,
               family=BI, data=maleDF, trace=FALSE)
maleDF$binprob <- with(maleDF, dBI(number, bd=12, mu = fv(doBI)[1]))

## ---- alt ----
gph <- fitVSobs(Number="number", Freq="freq", Fitprob="binprob", df=maleDF,
                distribution="binomial", size=12,
                atright=list(at=(1:4)*305.75, labels=(1:4)*.05),
                main="B: No. of males out of 12 -- binomial fit")

## ---- cfMales ----
## Male children
doBI <- gamlss(cbind(number, 12-number)~1, weights=freq,
               family=BI, data=maleDF, trace=FALSE)
doBB <- gamlss(cbind(number, 12-number)~1, weights=freq,
                       family=BB, data=maleDF, trace=FALSE)
doDBI <- gamlss(cbind(number, 12-number)~1, weights=freq,
                       family=DBI, data=maleDF, trace=FALSE, n.cyc=100)

## ---- binFitProbs ----
maleDF$BBfit <- with(maleDF, dBB(number, bd=12,
                                 mu = fv(doBB, 'mu')[1],
                                 sigma=fv(doBB, 'sigma')[1]))
maleDF$DBIfit <- with(maleDF, dBB(number, bd=12,
                                  mu = fv(doDBI, 'mu')[1],
                                  sigma=fv(doDBI, 'sigma')[1]))

## ---- rqmales ----
par(mar=c(2.6,4.1,2.6,1.6), mgp=c(2.5,.75,0), mfrow=c(1,2))
rqres.plot(doBI, plot.type='all', type="wp", main="")
mtext(side=3, line=0.75, "A: Binomial model", adj=0, cex=1.25)
rqres.plot(doBB, plot.type='all', type="wp", main="", ylab='')
mtext(side=3, line=0.75, "B: Betabinomial model", adj=0, cex=1.25)

## ---- cf-AIC-males ----
aicStat <- AIC(doBI, doBB, doDBI)
rownames(aicStat) <-
  c(doBI="Binomial",doBB="Betabinomial",doDBI="Double binomial")[rownames(aicStat)]
aicStat$dAIC <- with(aicStat, round(AIC-AIC[1],1))
aicStat

## ---- countdata ----
## Radioactive count data
rugetab <- with(VGAM::ruge, setNames(as.table(counts), number),
                dnn=list("A: Numbers of scintillations in 2608 1/8 minute intervals"))
## Machinist accidents
machtab <- with(VGAM::machinists, setNames(as.table(ofreq), accidents),
                dnn=list("B: Numbers of accidents in three months"))

## ---- ruge ----
rugetab

## ---- accs ----
machtab

## ---- countFits ----
## Radioactive count data
muruge <- with(VGAM::ruge,
               weighted.mean(number, w=counts))
rugeFit <- with(VGAM::ruge, data.frame(number=number, freq=counts,
                                       fitprob = dpois(number, lambda = muruge)))
## Machinist accidents
mumach <- with(VGAM::machinists, weighted.mean(accidents, w=ofreq))
machFit <- with(VGAM::machinists,
                data.frame(number=accidents, freq=ofreq,
                           fitprob = dpois(accidents, lambda = mumach)))

## ---- poisAB ----
gph1 <- fitVSobs(Number="number", Freq="freq", Fitprob="fitprob", df=rugeFit,
                 distribution='poisson', digits=2,
                 xlab="Scintillations: 1/8 minute intervals",
                 main=list("A: Radioactive decay counts -- poisson fit",
                           x=0, cex=1.05, font=1, just="left"),
                 atright=list(at=(0:4)*0.05*2608, labels=(0:4)*0.05))
## Machinists
gph2 <- fitVSobs(Number="number", Freq="freq", Fitprob="fitprob", df=machFit,
                 distribution='poisson', digits=2,
                 xlab="Accident counts: 3-month intervals",
                 ylab="Frequency",
                 main=list("B: Machinist accidents data", x=0,
                           cex=1.05, font=1, just="left"),
                 atright=list(at=(0:4)*0.2*414, labels=(0:4)*0.2))

## ---- cfFits-poiss ----
dopoiss <- gamlss(number~1, weights=freq,
                  family=PO, data=machFit, trace=FALSE)
doNBI <- gamlss(number~1, weights=freq,
                family=NBI, data=machFit, trace=FALSE)
doPIG <- gamlss(number~1, weights=freq,
                family=PIG, data=machFit, trace=FALSE)
doZIP <- gamlss(number~1, weights=freq,
                  family=ZIP, data=machFit, trace=FALSE)
doZINBI <- gamlss(number~1, weights=freq,
                family=ZINBI, data=machFit, trace=FALSE)
doZIPIG <- gamlss(number~1, weights=freq,
                family=ZIPIG, data=machFit, trace=FALSE, n.cyc=25)

## ---- cf-AIC-counts ----
aicStat <- AIC(dopoiss, doNBI, doPIG, doZIP, doZINBI, doZIPIG)
aicStat$dAIC <- with(aicStat, round(AIC-AIC[1],1))
rownames(aicStat) <-
  c(dopoiss="Poisson",doNBI="Negative binomial",
                       doPIG="Poisson inverse gamma",
  doZIP="ZI Poisson",doZINBI="ZI Negative binomial",
  doZIPIG="ZI Poisson inverse gamma")[rownames(aicStat)]
aicStat

## ---- wormpois ----
par(mar=c(2.6,4.1,2.6,1.6), mgp=c(2.5,.75,0), mfrow=c(1,3))
rqres.plot(dopoiss, plot.type='all', type="wp")
mtext(side=3, line=1, "A: Worm plot, poisson fit", cex=0.85, adj=0)
rqres.plot(doZIP, plot.type='all', type="wp")
mtext(side=3, line=1, "B: Zero-inflated poisson", cex=0.85, adj=0)
rqres.plot(doNBI, plot.type='all', type="wp", main="Worm plot")
mtext(side=3, line=1, "C: Negative binomial", cex=0.85, adj=0)

## ---- EstFreqs ----
mach <- setNames(numeric(9), nm=0:8)
mach[names(machtab)] <- machtab
N <- sum(machtab)
NBIprob <- dNBI(0:8, mu=fv(doNBI)[1], sigma=fv(doNBI, 'sigma')[1])
PIGprob <- dPIG(0:8, mu=fv(doPIG)[1], sigma=fv(doPIG, 'sigma')[1])
ZIPprob <- dZIP(0:8, mu=fv(doZIP)[1], sigma=fv(doZIP, 'sigma')[1])
estFreqs <- rbind(NegBin=NBIprob, PIG=PIGprob, ZeroInfPoiss=ZIPprob)*N
rbind(Frequencies=mach, round(estFreqs,1))
