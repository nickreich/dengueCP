library(lattice)
library(ggplot2)
library(doSMP)
library(xtable)
source("../TSIR_Utils.R")


n.files <- 12000
filenames <- paste("../../data/simulatedData_fromMI/simdata_5_local/simdata",1:n.files,".csv", sep="")

## calculate births for each dataset
BIRTHS.PER.YEAR <- 10000000*.02 ## = yearly number of births
births <- rep(BIRTHS.PER.YEAR/26, 26*40+1)

## calculate dates for each dataset
years <- 2000:2040
new.dates <- rep(as.Date("2000/01/01"), length=26*40+1)
idx <- 1
for(year in years) {
	if(year<2040) {
		new.dates[idx:(idx+25)] <- seq(as.Date(paste(year,"/1/1", sep="")),
					       by="2 weeks", length.out=26)
	} else {
		new.dates[idx] <- as.Date("2040/1/1")
	}
	idx <- idx+26
}

## model formulas
model.a.formula <- formula(log(It) ~ log(It.minus.1) + Zt.minus.1-1)
model.b.formula <- formula(log(It) ~ factor(str) + log(It.minus.1) + Zt.minus.1-1)
model.c.formula <- formula(log(It) ~ factor(biweeks) + log(It.minus.1) + Zt.minus.1-1)
model.d.formula <- formula(log(It) ~ factor(biweeks):factor(str)+log(It.minus.1)+Zt.minus.1-1)

##############
## EXP FITS ##
##############

## fitting parameters
DF <- 2 ## KEEP SAME AS IN ACTUAL DATA FITTING?
DATE.STRING <- format(Sys.time(), "%m%d%Y")
GLOBAL.INDEX <- TRUE
birth.lag <- 8
max.iter <- 100
grid.size <- 201
delta <- 1 #seq(-1,1, length.out=grid.size)
k <- 260
lambdas <- round(seq(0, 200, length.out=grid.size))
analysis.start.date <- as.Date("2011-01-01")

ll.store.list <- vector("list", n.files)
final.logliks <- matrix(NA, nrow=n.files, ncol=9)
N.fits <- matrix(NA, nrow=n.files, ncol=2)
colnames(N.fits) <- c("N.fit.N", "N.fit.Ns")

## parallel settings
w1 <- startWorkers(workerCount = 11)
registerDoSMP(w1)

for(j in 1:length(filenames)){
	filename <- filenames[j]
	dat <- read.csv(filename, row.names=1)

	## plot biweek data
	ddat <- melt(dat[,c("time", "y1", "y2", "y3", "y4")], id="time")
	##pdf(output.names[j])
	##gg1 <- qplot(time, value, data=ddat, geom="line", group=variable,
	##	     ylab="case counts", xlab="time") + #,main=fignames[j]
	##		     facet_grid(variable~.)
	##print(gg1)
	##dev.off()

	## tie all data together
	new.dat <- cbind(biweek=c(rep(1:26, times=40), 1),
			 date=new.dates,
			 dat,
			 births=births)

	## fit.model.Na
	tmp <- get.cross.protect.data(new.dat,
				      case.cols=paste("y", 1:4, sep=""),
				      delta=0, cr.protect.lag=0,
				      analysis.start.time=analysis.start.date,
				      birth.lag=birth.lag,
				      n.iter=max.iter, tol=1e-5,
				      smooth=TRUE, n.smooth=1,
				      df.low=DF, df.high=DF)

	global.idx <- tmp$nonzero.subset
	analysis.data <- tmp$analysis.data[global.idx,]
	bestfit.Nafd <- lm(model.a.formula, data=data.frame(analysis.data))
	fit.na <- logLik(bestfit.Nafd)
	bestfit.Nbfd <- lm(model.b.formula, data=data.frame(analysis.data))
	fit.nb <- logLik(bestfit.Nbfd)
	bestfit.Ncfd <- lm(model.c.formula, data=data.frame(analysis.data))
	fit.nc <- logLik(bestfit.Ncfd)
	bestfit.Ndfd <- lm(model.d.formula, data=data.frame(analysis.data))
	fit.nd <- logLik(bestfit.Ndfd)
	N.fits[j,2] <- length(tmp$nonzero.subset)


	## fit model N
	long.dat <- melt(new.dat,
			 id.vars=c("biweek", "date", "births", "time"),
			 variable_name="strain")
	outcome.idx <- which(long.dat$date>analysis.start.date)
	lag1.idx <- outcome.idx-1
	lag1.cases <- long.dat[lag1.idx,"value"]

	long.analysis.dat <- cbind(long.dat[outcome.idx,],
				   lag1=lag1.cases)

	if(GLOBAL.INDEX) {
		subs <- global.idx
	} else {
		subs <- (long.analysis.dat[,"lag1"]>0 & long.analysis.dat[,"value"]>0)
	}

	lag1.norm.fit <- glm(log(value)~log(lag1),
			     family=gaussian,
			     data=long.analysis.dat,
			     subset=subs)
	y.hat <- predict(lag1.norm.fit, type="response")
	cor(y.hat, long.analysis.dat[subs,"value"])
	fit.n <- logLik(lag1.norm.fit)
	N.fits[j,1] <- length(y.hat)


	## run optimization
	ll.store <- foreach(i = 1:length(lambdas), .combine=rbind) %dopar% {
		tmp <- get.cross.protect.data(new.dat,
					      case.cols=paste("y", 1:4, sep=""),
					      delta=delta, cr.protect.lag=k,
					      lambda=lambdas[i],
					      analysis.start.time=analysis.start.date,
					      birth.lag=birth.lag,
					      n.iter=max.iter, tol=1e-5,
					      smooth=FALSE, n.smooth=1,
					      df.low=DF, df.high=DF)
		if(GLOBAL.INDEX) {
			subs <- global.idx
		} else {
			subs <- tmp$nonzero.subset
		}
		analysis.data <- tmp$analysis.data[subs,]
		N.fit.Es <- nrow(analysis.data)

		## Exp models
		fit.Ea <- logLik(lm(model.a.formula,
				   data=data.frame(analysis.data)))
		fit.Eb <- logLik(lm(model.b.formula,
				   data=data.frame(analysis.data)))
		fit.Ec <- logLik(lm(model.c.formula,
				   data=data.frame(analysis.data)))
		fit.Ed <- logLik(lm(model.d.formula,
				   data=data.frame(analysis.data)))
		##if(i%%10==0)
		##	message(paste("    round", i, "complete ::", Sys.time()))
		c(fit.Ea, fit.Eb, fit.Ec, fit.Ed, N.fit.Es, lambdas[i])
	}
	colnames(ll.store) <- c("mod.a", "mod.b", "mod.c", "mod.d", "N.fit", "lambda")
	ll.store.list[[j]] <- ll.store

	final.logliks[j,] <- c(fit.n, fit.na, fit.nb, fit.nc, fit.nd,
			       apply(ll.store[,1:4], FUN=max, MAR=2))

	message(paste("dataset", j, "complete ::", Sys.time()))
}

stopWorkers(w1)

## check sample sizes were the same for all lambdas
for(i in sample(1:n.files, 10)) print(table(ll.store.list[[i]][,"N.fit"]))

## save modeled sample sizes for exp fits
N.fits.exp <- rep(NA, n.files)
for(i in 1:n.files) N.fits.exp[i] <- ll.store.list[[i]][1,"N.fit"]
N.fits.all <- cbind(N.fits, N.fits.exp)
N.fits.reps <- cbind(N.fits.all[,"N.fit.N"],
		     N.fits.all[,"N.fit.Ns"], N.fits.all[,"N.fit.Ns"], N.fits.all[,"N.fit.Ns"], N.fits.all[,"N.fit.Ns"],
		     N.fits.all[,"N.fits.exp"], N.fits.all[,"N.fits.exp"], N.fits.all[,"N.fits.exp"], N.fits.all[,"N.fits.exp"])
stable.data <- N.fits.exp == N.fits.all[,"N.fit.Ns"]

##################################
## summarize model GOF results  ##
##################################

## number of coefficients in models a->d + 1 for estimated variance
n.mod.params <- c(3, 7, 29, 106)

## multiplying DF (for smooth spline) by 4 b/c splines fit for each of 4 serotypes
## adding 1 in third c()'d value for lambda parameter
DFs <- matrix(c(3, n.mod.params+DF*4, n.mod.params+DF*4+1), byrow=TRUE, nrow=n.files, ncol=9)
BICs <- -2*final.logliks + DFs*log(N.fits.reps)
AICs <- -2*final.logliks + 2*DFs
final.loglik.diffs <- (final.logliks-final.logliks[,1])/(2*DFs)
colnames(AICs) <- c("N", "Na", "Nb", "Nc", "Nd", "Ea", "Eb", "Ec", "Ed")
rownames(AICs) <- paste("dataset", 1:n.files)

best.mods <- best.idx <- rep(NA, n.files)
for(i in 1:n.files) {
	n.mod <- ncol(final.logliks)
	## choose best model as one with lowest AIC.
	best.idx[i] <- which( AICs[i,] == min(AICs[i,]) )
	best.mods[i] <- colnames(AICs)[best.idx[i]]
}
final.loglik.mods <- data.frame(round(AICs, 1), best.mods)

good.idx <- which(final.logliks[,1]<0)

save.image(file=paste("../../data/simulationResults_", DATE.STRING, ".rda", sep=""))


## summarize point estimates from the model fits
chosen.ests <- matrix(NA, ncol=8, nrow=n.files) ## estimates chosen by our algorithm
colnames(chosen.ests) <- c("lambda", "lam.ci.l", "lam.ci.u", "loglik",
			   "rho1", "rho2", "rho3", "rho4")
correctModel.ests <- chosen.ests ## estimates from the "correct" model
lik.cut <- qchisq(.95, df=1)/2

for(j in 1:n.files){
	## get estimates from "chosen" model by chi-square tests
	max.lik.chosen <- final.logliks[j, best.idx[j]]
	chosen.ests[j,"loglik"] <- max.lik.chosen
	if (best.idx[j]<6){
		chosen.ests[j,"lambda"] <- 0
	} else {
		ll.col.num <- best.idx[j]-5
		tmp.idx <- which(ll.store.list[[j]][,ll.col.num] == max.lik.chosen)
		chosen.ests[j,"lambda"] <- ll.store.list[[j]][tmp.idx,"lambda"]
		lambda.ci.idx <- which(ll.store.list[[j]][,ll.col.num] > max.lik.chosen-lik.cut)
		chosen.ests[j,"lam.ci.l"] <- min(ll.store.list[[j]][lambda.ci.idx,"lambda"])
		chosen.ests[j,"lam.ci.u"] <- max(ll.store.list[[j]][lambda.ci.idx,"lambda"])
	}
	## refit best model
	tmp <- get.cross.protect.data(new.dat,
				      case.cols=paste("y", 1:4, sep=""),
				      delta=delta, cr.protect.lag=k,
				      lambda=chosen.ests[j,"lambda"],
				      analysis.start.time=analysis.start.date,
				      birth.lag=birth.lag,
				      n.iter=max.iter, tol=1e-5,
				      smooth=FALSE, n.smooth=1,
				      df.low=DF, df.high=DF)
	chosen.ests[j,c("rho1", "rho2", "rho3", "rho4")] <- tmp$curr.rho[1,]
	## correct model fits
	max.lik.correct <- final.logliks[j, 6]
	correctModel.ests[j,"loglik"] <- max.lik.correct
	tmp.idx <- which(ll.store.list[[j]][,1] == max.lik.correct)
	correctModel.ests[j,"lambda"] <- ll.store.list[[j]][tmp.idx,"lambda"]
	lambda.ci.idx <- which(ll.store.list[[j]][,1] > max.lik.correct-lik.cut)
	correctModel.ests[j,"lam.ci.l"] <- min(ll.store.list[[j]][lambda.ci.idx,"lambda"])
	correctModel.ests[j,"lam.ci.u"] <- max(ll.store.list[[j]][lambda.ci.idx,"lambda"])
	## refit best model
	tmp <- get.cross.protect.data(new.dat,
				      case.cols=paste("y", 1:4, sep=""),
				      delta=delta, cr.protect.lag=k,
				      lambda=correctModel.ests[j,"lambda"],
				      analysis.start.time=analysis.start.date,
				      birth.lag=birth.lag,
				      n.iter=max.iter, tol=1e-5,
				      smooth=FALSE, n.smooth=1,
				      df.low=DF, df.high=DF)
	correctModel.ests[j,c("rho1", "rho2", "rho3", "rho4")] <- tmp$curr.rho[1,]
}
chosen.ests.final <- data.frame(model=best.mods, chosen.ests)

write.csv(chosen.ests.final, file="../../data/chosenEsts_MIsims5_June2012.csv", row.names=F)
write.csv(correctModel.ests, file="../../data/correctEsts_MIsims5_June2012.csv", row.names=F)


## ENDED HERE
#save.image(file=paste("../../data/simulationResults_", DATE.STRING, ".rda", sep=""))


