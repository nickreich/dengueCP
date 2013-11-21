## running Dengue analysis on Puerto Rico data
## Nicholas Reich and Michael Johannson
## March 2013

library(ggplot2)
library(reshape2)
library(doMC)
source("TSIR_Utils.R")

dat <- read.csv("../data/model1.csv")

#dat[, 4:7] <- dat[,4:7]+1
# dat <- dat[1:600,]

ddat <- melt(dat[,c("date",
                                     "den.cases.str1",
                                     "den.cases.str2",
                                     "den.cases.str3",
                                     "den.cases.str4")],
             id="date")

qplot(date, value, data=ddat, geom="line", group=variable,
      main="Bi-weekly dengue incidence over time by strain (in Puerto Rico)",
      ylab="case counts", xlab="time")+
        facet_grid(variable~.)



model.a.formula <- formula(log(It) ~ log(It.minus.1) + Zt.minus.1-1)
model.b.formula <- formula(log(It) ~ factor(str) + log(It.minus.1) + Zt.minus.1-1)
model.c.formula <- formula(log(It) ~ factor(biweeks) + log(It.minus.1) + Zt.minus.1-1)
model.d.formula <- formula(log(It) ~ factor(biweeks):factor(str)+log(It.minus.1)+Zt.minus.1-1)

DATE.STRING <- format(Sys.time(), "%m%d%Y")
DF <- 3
birth.lag <- 8

lambda.max <- 4 # in years
grid.size <- 10 #lambda.max*26+1
lambdas <- seq(0, lambda.max*26, length.out=grid.size)
k <- round(qexp(.75, 1/max(lambdas)))
starting.idx <- birth.lag + k + 1
analysis.start.date <- dat[starting.idx,"date"]
max.iter <- 100
deltas <- 1

## register the parallel backend
#registerDoMC()


## parallel loop for exp model fitting
delta.grid.size <- 11
delta.vec <- rep(seq(-1,1, length.out=delta.grid.size), each=length(lambdas))
lambdas.vec <- rep(lambdas, times=delta.grid.size)
params <- cbind(delta.vec, lambdas.vec)


exp.logliks <- foreach(i = 1:nrow(params), .combine=rbind) %dopar% {
        tmp <- get.cross.protect.data(dat,
                                      case.cols=paste("den.cases.str", 1:4, sep=""),
                                      delta=params[i,"delta.vec"], cr.protect.lag=k,
                                      lambda=params[i,"lambdas.vec"],
                                      analysis.start.time=analysis.start.date,
                                      birth.lag=birth.lag,
                                      n.iter=max.iter, tol=1e-5,
                                      smooth=TRUE, n.smooth=1,
                                      df.low=DF, df.high=DF)
        analysis.data <- tmp$analysis.data[tmp$nonzero.subset,]
        ll.Ea <- logLik(lm(model.a.formula, data=data.frame(analysis.data)))
        ll.Eb <- logLik(lm(model.b.formula, data=data.frame(analysis.data)))
        ll.Ec <- logLik(lm(model.c.formula, data=data.frame(analysis.data)))
        ll.Ed <- logLik(lm(model.d.formula, data=data.frame(analysis.data)))
        message(paste("round", i, "complete ::", Sys.time()))
        c(params[i,], ll.Ea, ll.Eb, ll.Ec, ll.Ed, length(tmp$nonzero.subset))
}

colnames(exp.logliks) <- c("delta", "lambda", "loglik.Ea", "loglik.Eb", "loglik.Ec", "loglik.Ed", "N.fit")

Ec.sum <- TSIR.post.estimation.expmodels(exp.logliks, ll.col="loglik.Ec",
                                         DF=3, k=k, model.name="model E_C",
                                         model.formula=model.c.formula,
                                         dat=dat,
                                         st.date=analysis.start.date)

Eb.sum <- TSIR.post.estimation.expmodels(exp.logliks, ll.col="loglik.Eb",
                                         DF=3, k=k, model.name="model E_B",
                                         model.formula=model.b.formula,
                                         dat=dat,
                                         st.date=analysis.start.date)
Ed.sum <- TSIR.post.estimation.expmodels(exp.logliks, ll.col="loglik.Ed",
                                         DF=3, k=k, model.name="model E_D",
                                         model.formula=model.d.formula,
                                         dat=dat,
                                         st.date=analysis.start.date)


Ec.dat <- get.dat.final(Ec.sum)
par(mfrow=c(4,1), mar=c(2, 2, 2, 1))
rge <- range(Ec.dat$rho.str1, Ec.dat$rho.str2, Ec.dat$rho.str3, Ec.dat$rho.str4)
with(Ec.dat, plot(date, rho.str1, ylim=rge, type="l"))
with(Ec.dat, points(date, rho.str2, ylim=rge, type="l"))
with(Ec.dat, points(date, rho.str3, ylim=rge, type="l"))
with(Ec.dat, points(date, rho.str4, ylim=rge, type="l"))


## plot raw logliklihood
ggplot(data.frame(exp.logliks)) + 
        geom_tile(aes(x=lambda, y=delta, fill=loglik.Ec))

## transmission parameter plots
trans.param.ests <- trans.param.plot(Eb.sum, Ec.sum, Ed.sum)


#tmp <- data.frame(exp.logliks)
#exp.logliks.sub <- subset(tmp, delta<.8)
#ggplot(exp.logliks.sub) + 
#        geom_tile(aes(x=lambda, y=delta, fill=loglik.Ec))
