library(xtable)
library(lattice)
library(ggplot2)



library("foreach")
if (.Platform$OS.type != "windows" && require("multicore")) {
  registerDoMC()

} else if (require("doSNOW")) {
  cl <- snow::makeCluster(2, type = "SOCK")
  on.exit(snow::stopCluster(cl), add = TRUE)
  registerDoSNOW(cl)
 
}
library(reshape2)

load("data/bkk.dengue.all.cases.rda")
source("code/TSIR_Utils_Hannah_crossprotection_serotype_specific.R")
n.strains<-4


############
## aset-up ##
############

Sys.time()

## model formulas
model.a.formula <- formula(log(It) ~ log(It.minus.1) + Zt.minus.1-1)
model.b.formula <- formula(log(It) ~ factor(str) + log(It.minus.1) + Zt.minus.1-1)
model.c.formula <- formula(log(It) ~ factor(biweeks) + log(It.minus.1) + Zt.minus.1-1)
model.d.formula <- formula(log(It) ~ factor(biweeks):factor(str)+log(It.minus.1)+Zt.minus.1-1)


DATE.STRING <- format(Sys.time(), "%m%d%Y")
DF <- 2
birth.lag <- 8

lambda.max <- 6.5
lambda.estimate.max<-5
# in years
grid.size<- 4
  #lambda.max*26+1
lambdas <- seq(0, lambda.estimate.max*26, length.out=grid.size)
lambdas1 <- seq(0, 4*26, length.out=grid.size)
lambdas2 <- seq(0, 4*26, length.out=grid.size)
lambdas3 <- seq(0*26, 4*26, length.out=grid.size)
lambdas4 <- seq(0*26, 4*26, length.out=grid.size)

kcpss<-array(dim=c(n.strains,n.strains))
kcpss[] <- max(round(qexp(.75, 1/(lambda.max*26))))


starting.idx <- birth.lag + kcpss[1,1] + 1
analysis.start.date <- bkk.dengue.all.cases[starting.idx,"date"]
max.iter <- 10
deltas <- 1


## register the parallel backend
#HC this is now above


################
## EXP models ##
################

Sys.time()

## model N_Exp: no cross-protection with and without strain-specific transmission

zeros<-array(0, dim=c(n.strains, n.strains))

tmp <- get.cross.protect.data(bkk.dengue.all.cases,
			      case.cols=paste("den.cases.str", 1:4, sep=""),
			      delta.cpss=zeros, k.cpss=zeros, lambda=NULL,
			      analysis.start.time=analysis.start.date,
			      birth.lag=birth.lag,
			      n.iter=10, tol=1e-5,
			      smooth=TRUE, n.smooth=1,
			      df.low=DF, df.high=DF)
analysis.data <- tmp$analysis.data[tmp$nonzero.subset,]
bestfit.NaExp <- lm(model.a.formula, data=data.frame(analysis.data))
bestfit.NbExp <- lm(model.b.formula, data=data.frame(analysis.data))
ll.store.NaExp <- logLik(bestfit.NaExp)
ll.store.NbExp <- logLik(bestfit.NbExp)


## parallel loop for exp model fitting
#HC so for now I am estimating two parameters for  lambda 
#HC fixing delta at 1 (seems to be what was done in the paper for the exponential model)
#HCdelta.grid.size <- 50
#HCnumberp<-delta.grid.size*(length(lambdas))^2
#HC#HCdelta1.vec <- rep(seq(0,1, length.out=delta.grid.size), each=numberp)
#HCdelta2.vec <- rep(seq(0,1, length.out=delta.grid.size), each=numberp)

# HC making the lambda grid to take from- depending on how many lambdas there are- the size of this grid will change and the repeating will change too
#HC length(lambdas)= gridsize

lambdas1.vec <- rep(lambdas1, times=(length(lambdas)^3))#HC so this changes every 1
lambdas2.vec <- rep(lambdas2, each=(length(lambdas)^3))#HC this changes every grid.size^(noofparameters-1)
lambdas3.vec <-rep(rep(lambdas3, each=(length(lambdas)^2),times=(length(lambdas))))#HC this changes every grid.size^(2)
lambdas4.vec <- rep(rep(lambdas4,times=(length(lambdas)^2) , each=(length(lambdas)))) #HC this changes every grid.size


paramsexp <- cbind( lambdas1.vec, lambdas2.vec,lambdas3.vec,lambdas4.vec)
exp.logliks<-array(dim=c(length(paramsexp[,1]), (length(paramsexp[1,])+4+1)))
                   

                   
                   
# array to put the lambdas in 
lambdacpss<-array(dim=c(n.strains,n.strains))

#foreach(p = 1:nrow(params), .combine=rbind) %dopar% {
  
for ( p in 1:nrow(paramsexp)){
  #foreach(p = 1:nrow(params), .combine=rbind) %dopar% {
  
 
  #HC This array is to be altered depending on which combination of lambdas is being considered
  
  lambdacpss[1,1:4]<-paramsexp[p,"lambdas1.vec"] 
  lambdacpss[2,1:4]<-paramsexp[p,"lambdas2.vec"] 
  lambdacpss[3,1:4]<-paramsexp[p,"lambdas3.vec"]  
  lambdacpss[4,1:4]<-paramsexp[p,"lambdas4.vec"] 


  
  lambdacpss[1,1]<-3000
  lambdacpss[2,2]<-3000
  lambdacpss[3,3]<-3000
  lambdacpss[4,4]<-3000
  
	tmp <- get.cross.protect.data(bkk.dengue.all.cases,
				      case.cols=paste("den.cases.str", 1:4, sep=""),
                                #HC k is now just the max for the cross protective period
				      delta.cpss=1, k.cpss=kcpss,
				      lambda.cpss=lambdacpss,
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
 	message(paste("round", p, "complete ::", Sys.time()))
exp.logliks[p,] <- c(signif(paramsexp[p,],3)
                     , signif(ll.Ea,3)
                     , signif(ll.Eb,3)
                     , signif(ll.Ec,3)
                     , signif(ll.Ed,3)
                     ,signif(length(tmp$nonzero.subset),3))
  
}
Sys.time()

colnames(exp.logliks) <- c("lambda1", "lambda2", "lambda3","lambda4", "loglik.Ea", "loglik.Eb", "loglik.Ec", "loglik.Ed", "N.fit")

exp.logliks[which(exp.logliks[,5]==max(exp.logliks[,5])),1:4]
exp.logliks[which(exp.logliks[,6]==max(exp.logliks[,6])),1:4]
exp.logliks[which(exp.logliks[,7]==max(exp.logliks[,7])),1:4]
exp.logliks[which(exp.logliks[,8]==max(exp.logliks[,8])),1:4]

save.image(file=paste("../../data/TSIR_maxlam75yr_smooth", DF, "df_", DATE.STRING, ".rda", sep=""))


###########################
## Fixed duration models ##
###########################

## parallel loop for fixed duration models
k.min <- 1
k.max <- lambda.max*26 #100
delta.grid.size <- 4
k.grid.size <- 4
  #k.max # 101


numberp<-delta.grid.size*(k.grid.size)^2
numberpk<-k.grid.size*(delta.grid.size)
delta1.vec <- rep(seq(0,1, length.out=delta.grid.size), each=numberp)
delta2.vec <- rep(seq(0,1, length.out=delta.grid.size), each=numberpk)
delta2.vec2<-rep(delta2.vec,times=k.grid.size)
k1.vec <- rep(round(seq(k.min, k.max, length.out=k.grid.size)), each=delta.grid.size)
k1.vec2<-rep(k1.vec,times=k.grid.size^2)
k2.vec <- rep(round(seq(k.min, k.max, length.out=k.grid.size)), times=numberp)


params <- cbind(delta1.vec, delta2.vec2, k1.vec2, k2.vec)

Sys.time()

fixeddur.lls<-array(dim=c(nrow(params),11))
colnames(fixeddur.lls)<-c("delta1","delta2", "k1","k2", "loglik.Fa", "loglik.Fb", "loglik.Fc", "loglik.Fd", "eps", "iters", "N.fit")
for(p in 1:nrow(params)){
  #                   , .combine=rbind) %dopar% {
 
  
  deltacpss[1:3,1:3]<-params[p,"delta1.vec"]
  deltacpss[4,1:3]<-params[p,"delta2.vec2"]
  deltacpss[1:3,4]<-params[p,"delta2.vec2"]
  deltacpss[1,1]<-0
  deltacpss[2,2]<-0
  deltacpss[3,3]<-0
  deltacpss[4,4]<-0
  
  kcpss[1:3,1:3]<-params[p,"k1.vec2"]
  kcpss[4,1:3]<-params[p,"k2.vec"]
  kcpss[1:3,4]<-params[p,"k2.vec"]
  kcpss[1,1]<-0
  kcpss[2,2]<-0
  kcpss[3,3]<-0
  kcpss[4,4]<-0
  
	tmp <- get.cross.protect.data(bkk.dengue.all.cases,
				      case.cols=paste("den.cases.str", 1:4, sep=""),
				      delta.cpss= deltacpss,
				      cr.protect.lag.cpss=kcpss,lambda.cpss=NULL
				      analysis.start.time=analysis.start.date,
				      birth.lag=birth.lag,
				      n.iter=max.iter, tol=1e-5,
				      smooth=TRUE, n.smooth=1,
				      df.low=DF, df.high=DF)
	analysis.data <- tmp$analysis.data[tmp$nonzero.subset,]
	ll.Fa <- logLik(lm(model.a.formula, data=data.frame(analysis.data)))
	ll.Fb <- logLik(lm(model.b.formula, data=data.frame(analysis.data)))
	ll.Fc <- logLik(lm(model.c.formula, data=data.frame(analysis.data)))
	ll.Fd <- logLik(lm(model.d.formula, data=data.frame(analysis.data)))
	message(paste("round", p, "complete ::", Sys.time()))
  fixeddur.lls[p,]<-	c(params[p,], ll.Fa, ll.Fb, ll.Fc, ll.Fd, tmp$eps, tmp$iters, length(tmp$nonzero.subset))
}



Sys.time()
fixeddur.logliks<-fixeddur.lls
colnames(fixeddur.logliks) <- c("delta1.vec","delta2.vec", "k1.vec","k2.vec", "loglik.Fa", "loglik.Fb", "loglik.Fc", "loglik.Fd", "eps", "iters", "N.fit")


save.image(file=paste("../../data/TSIR_maxlam75yr_smooth", DF, "df_", DATE.STRING, ".rda", sep=""))



## model N_fd: no cross-protection with and without strain-specific transmission
N.cp.data <- get.cross.protect.data(bkk.dengue.all.cases,
                                    case.cols=paste("den.cases.str", 1:4, sep=""),
                                    delta=0, cr.protect.lag=0,
                                    analysis.start.time=analysis.start.date,
                                    birth.lag=birth.lag,
                                    n.iter=max.iter, tol=1e-5,
                                    smooth=TRUE, n.smooth=1,
                                    df.low=DF, df.high=DF)
analysis.data <- N.cp.data$analysis.data[N.cp.data$nonzero.subset,]
N.fit.Ns <- nrow(analysis.data)
bestfit.Nafd <- lm(model.a.formula, data=data.frame(analysis.data))
ll.Nafd <- logLik(bestfit.Nafd)
bestfit.Nbfd <- lm(model.b.formula, data=data.frame(analysis.data))
ll.Nbfd <- logLik(bestfit.Nbfd)
bestfit.Ncfd <- lm(model.c.formula, data=data.frame(analysis.data))
ll.Ncfd <- logLik(bestfit.Ncfd)
bestfit.Ndfd <- lm(model.d.formula, data=data.frame(analysis.data))
ll.Ndfd <- logLik(bestfit.Ndfd)

## Model N_a (fixed duration version) metrics
alpha1.nafd <- coef(bestfit.Nafd)[1]
alpha1.nafd.ci <- confint(bestfit.Nafd)[1,]

yhat <- exp(predict(bestfit.Nafd))
rho.cols <- paste("rho.str", 1:4, sep="")
dat.cols <- paste("den.cases.str", 1:4, sep="")
rep.factors <- as.matrix(N.cp.data$dat.final[-1, rho.cols])
rep.factors <- as.numeric(rep.factors)[N.cp.data$nonzero.subset]
yhat.scaled <- yhat/rep.factors
obs.y <- as.matrix(N.cp.data$dat.final[-1,dat.cols])
obs.y <- as.numeric(obs.y)[N.cp.data$nonzero.subset]
spCorr.Na <- cor(obs.y, yhat.scaled, method="pearson")

## Model N_b (fixed duration version) metrics
alpha1.nbfd <- coef(bestfit.Nbfd)["log(It.minus.1)"]
alpha1.nbfd.ci <- confint(bestfit.Nbfd)["log(It.minus.1)",]

yhat <- exp(predict(bestfit.Nbfd))
rho.cols <- paste("rho.str", 1:4, sep="")
dat.cols <- paste("den.cases.str", 1:4, sep="")
rep.factors <- as.matrix(N.cp.data$dat.final[-1, rho.cols])
rep.factors <- as.numeric(rep.factors)[N.cp.data$nonzero.subset]
yhat.scaled <- yhat/rep.factors
obs.y <- as.matrix(N.cp.data$dat.final[-1,dat.cols])
obs.y <- as.numeric(obs.y)[N.cp.data$nonzero.subset]
spCorr.Nb <- cor(obs.y, yhat.scaled, method="pearson")

## Model N_c (fixed duration version) metrics
alpha1.ncfd <- coef(bestfit.Ncfd)["log(It.minus.1)"]
alpha1.ncfd.ci <- confint(bestfit.Ncfd)["log(It.minus.1)",]

yhat <- exp(predict(bestfit.Ncfd))
rho.cols <- paste("rho.str", 1:4, sep="")
dat.cols <- paste("den.cases.str", 1:4, sep="")
rep.factors <- as.matrix(N.cp.data$dat.final[-1, rho.cols])
rep.factors <- as.numeric(rep.factors)[N.cp.data$nonzero.subset]
yhat.scaled <- yhat/rep.factors
obs.y <- as.matrix(N.cp.data$dat.final[-1,dat.cols])
obs.y <- as.numeric(obs.y)[N.cp.data$nonzero.subset]
spCorr.Nc <- cor(obs.y, yhat.scaled, method="pearson")

## Model N_d (fixed duration version) metrics
alpha1.ndfd <- coef(bestfit.Ndfd)["log(It.minus.1)"]
alpha1.ndfd.ci <- confint(bestfit.Ndfd)["log(It.minus.1)",]

yhat <- exp(predict(bestfit.Ndfd))
rho.cols <- paste("rho.str", 1:4, sep="")
dat.cols <- paste("den.cases.str", 1:4, sep="")
rep.factors <- as.matrix(N.cp.data$dat.final[-1, rho.cols])
rep.factors <- as.numeric(rep.factors)[N.cp.data$nonzero.subset]
yhat.scaled <- yhat/rep.factors
obs.y <- as.matrix(N.cp.data$dat.final[-1,dat.cols])
obs.y <- as.numeric(obs.y)[N.cp.data$nonzero.subset]
spCorr.Nd <- cor(obs.y, yhat.scaled, method="pearson")


########################
## fit linear model N ##
########################

## set up lag 1 data
start.time <- analysis.start.date

tmp <- melt(bkk.dengue.all.cases,
	    id.vars=c("biweek", "date", "births"),
	    variable_name="strain")
outcome.idx <- which(tmp$date>start.time)
lag1.idx <- outcome.idx-1
lag1.cases <- tmp[lag1.idx,"value"]

dat <- cbind(tmp[outcome.idx,], lag1=lag1.cases)

subs <- dat[,"lag1"]>0 & dat[,"value"]>0
lag1.norm.fit <- glm(log(value)~log(lag1), family=gaussian, data=dat[subs,])
y.hat <- predict(lag1.norm.fit, type="response")
spCorr.N <- cor(y.hat, dat[subs,"value"])
ll.N <- logLik(lag1.norm.fit)

alpha1.n <- coef(lag1.norm.fit)[2]
alpha1.n.ci <- confint(lag1.norm.fit)[2,]
N.fit.N <- length(y.hat)

save.image(file=paste("../../data/TSIR_maxlam75yr_smooth", DF, "df_", DATE.STRING, ".rda", sep=""))


##################################
## get summaries from FD models ##
##################################

#pdf("../../figures/publicationReady/fixeddur_A.pdf", height=6.5)
Fa.sum <- TSIR.post.estimation(fixeddur.logliks,
			       delta1.range=c(0, .97),delta2.range=c(0, .97), k1.range=c(0, 195),k2.range=c(0, 195), DF=3,
			       model.frmla=model.a.formula, ll.col="loglik.Fa",
			       dat=bkk.dengue.all.cases,
			       st.date=analysis.start.date,
			       main="Fixed duration model A log-likelihood surface")
#dev.off()

#pdf("../../figures/publicationReady/fixeddur_B.pdf", height=6.5)
Fb.sum <- TSIR.post.estimation(fixeddur.logliks,
			       delta1.range=c(0, .97), delta2.range=c(0, .97), k1.range=c(0, 195),k2.range=c(0, 195), DF=3,
			       model.frmla=model.b.formula, ll.col="loglik.Fb",
			       dat=bkk.dengue.all.cases,
			       st.date=analysis.start.date,
			       main="Fixed duration model B log-likelihood surface")
#dev.off()

#pdf("../../figures/publicationReady/fixeddur_C.pdf", height=3.5)
Fc.sum <- TSIR.post.estimation(fixeddur.logliks,
			       delta1.range=c(0, 1), delta2.range=c(0, 1),k1.range=c(0, 200),k2.range=c(0, 200), DF=3,
			       model.frmla=model.c.formula, ll.col="loglik.Fc",
			       dat=bkk.dengue.all.cases,
			       st.date=analysis.start.date, smooth=FALSE,
			       main="Fixed duration model C log-likelihood surface",
			       delta.plot.range=c(0,1))
#dev.off()

#pdf("../../figures/publicationReady/fixeddur_D.pdf", height=6.5)
Fd.sum <- TSIR.post.estimation(fixeddur.logliks,
			       delta1.range=c(-.5, 1), delta2.range=c(-.5, 1), k1.range=c(0, 200),k2.range=c(0, 200), DF=3,
			       model.frmla=model.d.formula, ll.col="loglik.Fd",
			       dat=bkk.dengue.all.cases,
			       st.date=analysis.start.date,
			       main="Fixed duration model D log-likelihood surface")
#dev.off()


## FD models for simulating datasets
## fix delta==1
Fc.sum.d1 <- TSIR.post.estimation.delta1(fixeddur.logliks, delta.point=1,
					 k.range=c(0, 90), DF=3,
					 model.frmla=model.c.formula, ll.col="loglik.Fc",
					 dat=bkk.dengue.all.cases,
					 st.date=analysis.start.date,
					 main="Fixed duration model C log-likelihood surface",
					 delta.plot.range=c(0,1))

## fix delta==0
Fc.sum.d0 <- TSIR.post.estimation.delta1(fixeddur.logliks, delta.point=0,
					 k.range=c(0, 1), DF=3,
					 model.frmla=model.c.formula, ll.col="loglik.Fc",
					 dat=bkk.dengue.all.cases,
					 st.date=analysis.start.date,
					 main="Fixed duration model C log-likelihood surface",
					 delta.plot.range=c(0,1))

save.image(file=paste("../../data/TSIR_maxlam75yr_smooth", DF, "df_", DATE.STRING, ".rda", sep=""))

###################################
## get summaries from EXP models ##
###################################

#pdf("../../figures/publicationReady/Exp_A.pdf")
Ea.sum <- TSIR.post.estimation.expmodels(exp.logliks, ll.col="loglik.Ea",
					 DF=3, k=k, model.name="model E_A",
					 model.formula=model.a.formula,
					 dat=bkk.dengue.all.cases,
					 st.date=analysis.start.date)
#dev.off()

#pdf("../../figures/publicationReady/Exp_B.pdf")
Eb.sum <- TSIR.post.estimation.expmodels(exp.logliks, ll.col="loglik.Eb",
					 DF=3, k=k, model.name="model E_B",
					 model.formula=model.b.formula,
					 dat=bkk.dengue.all.cases,
					 st.date=analysis.start.date)
#dev.off()

#pdf("../../figures/publicationReady/Exp_C.pdf", height=4.5)
Ec.sum <- TSIR.post.estimation.expmodels(exp.logliks, ll.col="loglik.Ec",
					 DF=3, k=k, model.name="model E_C",
					 model.formula=model.c.formula,
					 dat=bkk.dengue.all.cases,
					 st.date=analysis.start.date)
#dev.off()

#pdf("../../figures/publicationReady/Exp_D.pdf")
Ed.sum <- TSIR.post.estimation.expmodels(exp.logliks, ll.col="loglik.Ed",
					 DF=3, k=k, model.name="model E_D",
					 model.formula=model.d.formula,
					 dat=bkk.dengue.all.cases,
					 st.date=analysis.start.date)
#dev.off()

save.image(file=paste("../../data/TSIR_maxlam75yr_smooth", DF, "df_", DATE.STRING, ".rda", sep=""))


####################################
## make summary tables for xtable ##
####################################

## table for main paper
Cp.type <- c("none", "", "", "", "", "fixed duration", "", "", "", "exponential", "", "", "")
model.name <- c("$N$", "$N_a$", "$N_b$", "$N_c$", "$N_d$", "$F_a$", "$F_b$", "$F_c$", "$F_d$", "$E_a$", "$E_b$", "$E_c$", "$E_d$")
rep.frac <- c("-", rep("$\\bullet$", 12))
serotype.trans <- c("-", rep(c("-", "$\\bullet$"), times=6))
time.trans <- c("-", rep(c("-", "-", "$\\bullet$", "$\\bullet$"), times=3))
## format for strings
ff <- "%#.2f"
ffa <- "%#.3f"
## lambdas
all.lambda.mles <- sprintf(c(Fa.sum$lambda.mle, Fb.sum$lambda.mle, Fc.sum$lambda.mle, Fd.sum$lambda.mle,
                           Ea.sum$lambda.mle, Eb.sum$lambda.mle, Ec.sum$lambda.mle, Ed.sum$lambda.mle), fmt=ff)
all.lambda.ci.lows <- sprintf(c(Fa.sum$lambda.ci[1], Fb.sum$lambda.ci[1], Fc.sum$lambda.ci[1], Fd.sum$lambda.ci[1],
                              Ea.sum$lambda.ci[1], Eb.sum$lambda.ci[1], Ec.sum$lambda.ci[1], Ed.sum$lambda.ci[1]),
			    fmt=ff)
all.lambda.ci.highs <- sprintf(c(Fa.sum$lambda.ci[2], Fb.sum$lambda.ci[2], Fc.sum$lambda.ci[2], Fd.sum$lambda.ci[2],
                               Ea.sum$lambda.ci[2], Eb.sum$lambda.ci[2], Ec.sum$lambda.ci[2], Ed.sum$lambda.ci[2]),
			     fmt=ff)
lambda.w.cis <- c("-", "-", "-", "-", "-",
                  paste("$\\cithree{", all.lambda.ci.lows,"}{", all.lambda.mles,"}{", all.lambda.ci.highs,"}$"))
## alphas
all.alpha1.mles <- sprintf(c(alpha1.n, alpha1.nafd, alpha1.nbfd, alpha1.ncfd, alpha1.ndfd,
                           Fa.sum$alpha1.mle, Fb.sum$alpha1.mle, Fc.sum$alpha1.mle, Fd.sum$alpha1.mle,
                           Ea.sum$alpha1.mle, Eb.sum$alpha1.mle, Ec.sum$alpha1.mle, Ed.sum$alpha1.mle), fmt=ffa)
all.alpha1.ci.lows <- sprintf(c(alpha1.n.ci[1], alpha1.nafd.ci[1], alpha1.nbfd.ci[1],
				alpha1.ncfd.ci[1], alpha1.ndfd.ci[1],
				Fa.sum$alpha1.ci[1], Fb.sum$alpha1.ci[1], Fc.sum$alpha1.ci[1], Fd.sum$alpha1.ci[1],
				Ea.sum$alpha1.ci[1], Eb.sum$alpha1.ci[1], Ec.sum$alpha1.ci[1], Ed.sum$alpha1.ci[1]),
			      fmt=ffa)
all.alpha1.ci.highs <- sprintf(c(alpha1.n.ci[2], alpha1.nafd.ci[2], alpha1.nbfd.ci[2],
				 alpha1.ncfd.ci[2], alpha1.ndfd.ci[2],
                               Fa.sum$alpha1.ci[2], Fb.sum$alpha1.ci[2], Fc.sum$alpha1.ci[2], Fd.sum$alpha1.ci[2],
                               Ea.sum$alpha1.ci[2], Eb.sum$alpha1.ci[2], Ec.sum$alpha1.ci[2], Ed.sum$alpha1.ci[2]),
			       fmt=ffa)
alpha1.w.cis <- paste("$\\cithree{", all.alpha1.ci.lows,"}{", all.alpha1.mles,"}{", all.alpha1.ci.highs,"}$")
## other metrics
all.logliks <- round(c(ll.N, ll.Nafd, ll.Nbfd, ll.Ncfd, ll.Ndfd,
                       Fa.sum$logLik, Fb.sum$logLik, Fc.sum$logLik, Fd.sum$logLik,
                       Ea.sum$logLik, Eb.sum$logLik, Ec.sum$logLik, Ed.sum$logLik), 2)
N.dfs <- c(attr(ll.Nafd, "df"), attr(ll.Nbfd, "df"), attr(ll.Ncfd, "df"), attr(ll.Ndfd, "df"))+12
all.DF <- c(3, N.dfs,
            Fa.sum$df, Fb.sum$df, Fc.sum$df, Fd.sum$df,
            Ea.sum$df, Eb.sum$df, Ec.sum$df, Ed.sum$df)
change.in.ll <- ll.N-all.logliks
test.stat <- round(-2*change.in.ll, 1)
corr <- round(c(spCorr.N, spCorr.Na, spCorr.Nb, spCorr.Nc, spCorr.Nd,
                Fa.sum$corr, Fb.sum$corr, Fc.sum$corr, Fd.sum$corr,
                Ea.sum$corr, Eb.sum$corr, Ec.sum$corr, Ed.sum$corr), 3)
result.table <- data.frame(CP=CP.type, M=model.name, R=rep.frac, S=serotype.trans, T=time.trans, lambda=lambda.w.cis, alpha1=alpha1.w.cis, loglik=all.logliks, DF=all.DF, cor=corr)
result.table$BIC <- -2*result.table$loglik + result.table$DF*log(2145)
result.table$AIC <- -2*result.table$loglik + 2*result.table$DF
result.table$AICc <- result.table$AIC + (2*result.table$DF*(result.table$DF+1))/(2145-result.table$DF-1)
xtab <- xtable(result.table, caption="main results table", label="tab:mainResults",
	       align=c("l", "l", rep("c", 12)),
	       digits=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 3, 1, 1, 1))
print(xtab, sanitize.text.function=identity, include.rownames=FALSE, caption.placement="top",
      file="../../manuscripts/mainTable.tex")


## alternate smaller table
CP.type <- c("naive lag-1", "no CP", "fixed duration CP", "exponential CP")
## lambdas
all.lambda.mles <- sprintf(c(Fc.sum$lambda.mle, Ec.sum$lambda.mle), fmt=ff)
all.lambda.ci.lows <- sprintf(c(Fc.sum$lambda.ci[1], Ec.sum$lambda.ci[1]), fmt=ff)
all.lambda.ci.highs <- sprintf(c(Fc.sum$lambda.ci[2], Ec.sum$lambda.ci[2]), fmt=ff)
lambda.w.cis <- c("-", "-", paste("$\\cithree{", all.lambda.ci.lows,"}{", all.lambda.mles,"}{", all.lambda.ci.highs,"}$"))
## alphas
all.alpha1.mles <- sprintf(c(alpha1.n, alpha1.ncfd, Fc.sum$alpha1.mle, Ec.sum$alpha1.mle), fmt=ffa)
all.alpha1.ci.lows <- sprintf(c(alpha1.n.ci[1], alpha1.ncfd.ci[1], Fc.sum$alpha1.ci[1], Ec.sum$alpha1.ci[1]), fmt=ffa)
all.alpha1.ci.highs <- sprintf(c(alpha1.n.ci[2], alpha1.ncfd.ci[2], Fc.sum$alpha1.ci[2], Ec.sum$alpha1.ci[2]), fmt=ffa)
alpha1.w.cis <- paste("$\\cithree{", all.alpha1.ci.lows,"}{", all.alpha1.mles,"}{", all.alpha1.ci.highs,"}$")
## other metrics
all.logliks <- round(c(ll.N, ll.Ncfd, Fc.sum$logLik, Ec.sum$logLik), 2)
N.dfs <- c(3, attr(ll.Ncfd, "df")+12)
all.DF <- c(N.dfs, Fc.sum$df, Ec.sum$df)
change.in.ll <- ll.N-all.logliks
corr <- round(c(spCorr.N, spCorr.Nc, Fc.sum$corr, Ec.sum$corr), 3)
result.table <- data.frame(CP=CP.type, lambda=lambda.w.cis, alpha1=alpha1.w.cis, loglik=all.logliks, DF=all.DF, cor=corr)
xtab <- xtable(result.table, caption="main results table", label="tab:mainResults",
	       align=c("l", "l", rep("c", 5)),
	       digits=c(1, 1, 1, 1, 1, 0, 3))
print(xtab, sanitize.text.function=identity, include.rownames=FALSE, caption.placement="top",
      file="../../manuscripts/mainTableAlt.tex")


##############################################
## create table with reporting rates        ##
##############################################

## exponential model
rr.Ec <- get.reporting.rates(Ec.sum)

## fixed duration model
rr.Fc <- get.reporting.rates(Fc.sum)

xtable(rbind(rr.Fc, rr.Fc))

##############################################
## simulate data from best-fitting fd model ##
##############################################

Nc.sum <- list(bestfit = bestfit.Ncfd,
	       bestfit.data = get.cross.protect.data(bkk.dengue.all.cases,
					       case.cols=paste("den.cases.str", 1:4, sep=""),
					       delta=0,
					       cr.protect.lag=k,
					       analysis.start.time=analysis.start.date,
					       birth.lag=birth.lag,
					       n.iter=100, tol=1e-5,
					       smooth=TRUE, n.smooth=1,
					       df.low=DF, df.high=DF),
	       delta.mle = 0,
	       k.mle = 0,
	       lambda.mle = 0,
	       alpha1.mle = alpha1.ncfd)


tt.n <- simulate.many.TSIR.datasets(nsim=1000, Nc.sum, all.data=bkk.dengue.all.cases)

tt <- simulate.many.TSIR.datasets(nsim=1000, Fc.sum, all.data=bkk.dengue.all.cases)

## plot datasets
par(mfrow=c(4,1), mar=c(2,1,1,1))
tmp.col <- rgb(t(col2rgb("blue"))/255, alpha=.01)
plot.idx <- which(bkk.dengue.all.cases$date == analysis.start.date):nrow(bkk.dengue.all.cases)
for(j in 1:4){
        str.col <- paste0("den.cases.str", j)
        str.col2 <- paste0("I", j)
        plot(bkk.dengue.all.cases[plot.idx,"date"], bkk.dengue.all.cases[plot.idx,str.col], type="l", ylim=c(0,80), col="red", lwd=2)
        for(i in 1:1000){
                points(bkk.dengue.all.cases[plot.idx,"date"], tt[plot.idx,str.col2,i], type="l", col=tmp.col)
        }
        points(bkk.dengue.all.cases[plot.idx,"date"], bkk.dengue.all.cases[plot.idx,str.col], type="l", col="red", lwd=2)
}

## plot spectra
layout(matrix(1:10, ncol=2))
ylim <- c(0,120)
xlim <- c(0, 14)
only.new.dat <- TRUE
for(str in 1:4)
	plot.multiple.spectra(tt.n, tsir.sum=Nc.sum,
			      all.data=bkk.dengue.all.cases, str=str,
			      col="red", ylim=ylim, xlim=xlim,
			      only.new.data=only.new.dat)
plot.multiple.spectra(tt.n, tsir.sum=Nc.sum,
		      all.data=bkk.dengue.all.cases, str="all",
		      only.new.data=only.new.dat, plot.histogram=TRUE,
		      plot.spectra=FALSE)
##plot.multiple.spectra(tt.n, tsir.sum=Nc.sum,
##		      all.data=bkk.dengue.all.cases, str="all",
##		      col="red", ylim=ylim, xlim=xlim, only.new.data=only.new.dat)
for(str in 1:4)
	plot.multiple.spectra(tt, tsir.sum=Fc.sum,
			      all.data=bkk.dengue.all.cases, str=str,
			      col="blue", ylim=ylim, xlim=xlim,
			      only.new.data=only.new.dat)
##plot.multiple.spectra(tt, tsir.sum=Fc.sum.d1,
##		      all.data=bkk.dengue.all.cases, str="all",
##		      col="blue", ylim=ylim, xlim=xlim, only.new.data=only.new.dat)
plot.multiple.spectra(tt, tsir.sum=Fc.sum,
		      all.data=bkk.dengue.all.cases, str="all",
		      only.new.data=only.new.dat, plot.histogram=TRUE,
		      plot.spectra=FALSE)


##############################################
## simulate forward for best-fitting models ##
##############################################

ny <- seq(40, 150, by=10)
for(ny in ny){
tmp <- TSIR.simulate.forward(Ec.sum, model.c.formula, type="exp",
			     seasonality=TRUE, sim.years=ny,
			     data.years=25, nskip=26*20,
			     bootstrap=FALSE)
title(ny)
}

## examples for paper
tmp <- TSIR.simulate.forward(Ec.sum, model.c.formula, type="exp",
			     seasonality=TRUE, sim.years=100,
			     data.years=25, nskip=26*60,
			     bootstrap=FALSE)

tmp <- TSIR.simulate.forward(Ec.sum, model.c.formula, type="exp",
			     seasonality=TRUE, sim.years=70,
			     data.years=25, nskip=26*45,
			     bootstrap=FALSE)





## "good" example
tmp <- TSIR.simulate.forward(Ec.sum, model.c.formula, type="exp",
			     seasonality=TRUE, sim.years=100,
			     data.years=25, nskip=26*60,
			     bootstrap=TRUE, n.boots=10)

# why does the shape change when we include more/less years
tmp <- TSIR.simulate.forward(Ec.sum, model.c.formula, type="exp",
			     seasonality=TRUE, sim.years=100,
			     data.years=25, nskip=26*75,
			     bootstrap=TRUE, n.boots=10)


Nc.sum <- list(bestfit.data=N.cp.data,
               lambda.mle=0,
               alpha1.mle=alpha1.ncfd)

tmp <- TSIR.simulate.forward(Nc.sum, model.c.formula, type="none",
			     seasonality=TRUE, sim.years=100, data.years=25)


#######################################
## make transmission parameter plot  ##
#######################################
trans.param.plot(Fb.sum, Fc.sum, Fd.sum)

mat <- matrix(c(1, 0, 2:5), byrow=TRUE, ncol=2)
layout(mat, widths=c(1.5, 1))

par(mar=c(4, 4, 3, 1))
trans.param.ests <- trans.param.plot(Eb.sum, Ec.sum, Ed.sum)

###################
## precipitation ##
###################

## organize precipitation data
precip.dat <- read.csv("../../data/bangkok.precip.data.in.mm.csv", na.strings=-9999)
precip.dat$month <- round((precip.dat[,1]-floor(precip.dat[,1]))*12)+1
precip.dat$year <- floor(precip.dat[,1])

## aggregate the precipitation data
precip.agg <- aggregate(precip.dat$KrungThep, by=list(precip.dat$month),
			FUN=function(x) quantile(x, p=c(.10, .5, .90), na.rm=TRUE))
precip.agg.formatted <- data.frame(month.num=1:12,
				   month.name=month.name,
				   precip.agg[,2])
colnames(precip.agg.formatted) <- c("month.num", "month.names", "p10", "p50", "p90")

short.month.names <- substr(month.name, start=0, stop=3)

## plot precip in ggplot
##col1="black"
##i1 <- ggplot(precip.agg.formatted, aes(x=month.num)) + theme_bw()
##i2 <- i1 + geom_ribbon(aes(ymin=p10, ymax=p90), alpha=.1, fill=col1) + geom_line(aes(y=p50), col=col1, lwd=1.5) + xlim(1, 12) + ylim(0, max(precip.agg)) + xlab("
##month") + ylab("average monthly rainfall (in mm)") + scale_x_discrete(breaks = 1:12, labels=month.name)
##print(i2)

## plot precipitation
plot(1:12, precip.agg.formatted[,"p50"], ylim=range(precip.agg.formatted[,3:5]),
     type="n", lwd=3, xlab="", bty="n", ylab="rainfall (in mm)",
     las=1, xaxt="n", xlim=c(1, 13), main="median monthly precipitation in Bangkok")
polygon(c(1:12, 12:1)+.5, c(precip.agg.formatted[,"p10"],
			    rev(precip.agg.formatted[,"p90"])),
	col="lightgray", density=-1, border=FALSE)
lines(1:12+.5, precip.agg.formatted[,"p50"],lwd=3)
axis(1, at=seq(1, 13, by=2), labels=short.month.names[c(1, 3, 5, 7, 9, 11, 1)])

## get correlation
biweek.dates <- seq(as.Date("2010-01-01"), length.out=26, by="14 days")+7
month.dates <- seq(as.Date("2010-01-15"), length.out=12, by="1 month")

interp.month.dates <- seq(as.Date("2009-12-15"), length.out=14, by="1 month")
interp.precip <- precip.agg.formatted[c(12, 1:12, 1), "p50"]
precip.biweeks <- spline(interp.month.dates, cumsum(interp.precip),
			 xout=biweek.dates, method="fmm")
precip.biweek.counts <- precip.biweeks$y - c(0, precip.biweeks$y[1:25])

precip.fit<- lm(trans.param.ests[,"mle"]~precip.biweek.counts)

plot(precip.biweek.counts, trans.param.ests[,"mle"], xlab="biweekly precipitation (interpolated, in mm)", ylab="estimated transmission parameters", bty="n", type="n", las=1)
abline(precip.fit, lty=2, lwd=3, col="gray")
text(precip.biweek.counts, trans.param.ests[,"mle"], 1:26, cex=.8)
text(x=100, y=-.1, paste("correlation =",round(sqrt(summary(precip.fit)$r.squared), 2)),
     font=2)
title("precipitation vs. estimated transmission")

#################
## temperature ##
#################

## organize temperature data
temp.dat <- read.csv("../../data/bangkok.temp.monthly.in.c.csv", na.strings=-9999)
temp.dat$month <- round((temp.dat[,1]-floor(temp.dat[,1]))*12)+1
temp.dat$year <- floor(temp.dat[,1])

## aggregate the temperature data
temp.agg <- aggregate(temp.dat$KrungThep, by=list(temp.dat$month),
			FUN=function(x) quantile(x, p=c(.10, .5, .90), na.rm=TRUE))
temp.agg.formatted <- data.frame(month.num=1:12,
				   month.name=month.name,
				   temp.agg[,2])
colnames(temp.agg.formatted) <- c("month.num", "month.names", "p10", "p50", "p90")

## plot temperature
plot(1:12, temp.agg.formatted[,"p50"], ylim=range(temp.agg.formatted[,3:5]),
     type="n", lwd=3, xlab="", bty="n", ylab="temperature (in C)",
     las=1, xaxt="n", xlim=c(1, 13))
polygon(c(1:12, 12:1)+.5, c(temp.agg.formatted[,"p10"],
			    rev(temp.agg.formatted[,"p90"])),
	col="lightgray", density=-1, border=FALSE)
lines(1:12+.5, temp.agg.formatted[,"p50"],lwd=3)
axis(1, at=seq(1, 13, by=2), labels=short.month.names[c(1, 3, 5, 7, 9, 11, 1)])
title("median monthly temperature in Bangkok")

## get correlation
interp.temp <- temp.agg.formatted[c(12, 1:12, 1), "p50"]
temp.biweeks <- spline(interp.month.dates, interp.temp,
		       xout=biweek.dates, method="fmm")
temp.biweek.counts <- temp.biweeks$y

temp.fit<- lm(trans.param.ests[,"mle"]~temp.biweek.counts)

plot(temp.biweek.counts, trans.param.ests[,"mle"], xlab="biweekly temperature (interpolated, in C)", ylab="estimated transmission parameters", bty="n", type="n", las=1)
abline(temp.fit, lty=2, lwd=3, col="gray")
text(temp.biweek.counts, trans.param.ests[,"mle"], 1:26, cex=.8)
text(x=27, y=.1, paste("correlation =",round(sqrt(summary(temp.fit)$r.squared), 2)),
     font=2)
title("temperature vs. estimated transmission")
