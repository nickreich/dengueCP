## Functions for fitting TSIR model with cross-protection term
## based on code from /research/chickenPoxNRL/code/TSIRUtilities.R

## Nicholas Reich
## October 2010

require(ggplot2)
require(zoo)

## function to fit a complete cross-protection model for a single strain
##
## dat = data.frame with following columns, (plus others,
##       specified below, for case counts)
##          "biweek" = an integer in 1-26, indicating the biweek
##          "date"   = a date indicating the time period, t
##          "births" = an integer indicating the birth count for t
## case.cols = character object matching exactly the column names for
##             the case counts for each strain
## delta = the cross-protection parameter
## cr.protect.lag = an integer indicating the time an infected should
##                  be excluded from susceptible counts, i.e. "k"
## birth.lag = an integer indicating the time delay between birth
##   and introduction to the suceptible class (assumed from
##   mother's immunity)
## analysis.start.time = the first time at which the case counts should be used in fitting the transmission model
## n.iter = the number of iterations for the EM-type algorithm to run
## tol = tolerance for the EM-type algorithm to declare convergence
## smooth = indicator of whether to use local linear smooths
## df.low, df.high = min and max values for the smoother degrees if freedom
## nsmooth = number of smoothness levels to investigate

get.cross.protect.data <- function(dat, case.cols,
				   delta, cr.protect.lag, lambda=NULL,
				   birth.lag,
				   analysis.start.time,
				   n.iter=100, tol=1e-5,
				   smooth=TRUE, n.smooth=10,
				   df.low=1000, df.high=5000) {
	## check to make sure we have early enough data
	idx.start.date <- which(dat[,"date"]==analysis.start.time)
	k <- cr.protect.lag
	if(idx.start.date-k-birth.lag<1)
		stop("not enough early data points.")

	## 1. Reconstruct susceptible dynamics assuming strain-
	##    specific rho_t but no offset.
	##    Use smooth spline for fits.
	##    This is done to obtain starting values for rho.t.
	n.strains <- length(case.cols)

	## get indices for cumulative sums
	case.idx <- (idx.start.date-k):nrow(dat)
	analysis.idx <- idx.start.date:nrow(dat)
	birth.idx <- analysis.idx - birth.lag

	## get cumulative birth count
	cumbirths.lagged <- cumsum(dat[birth.idx,"births"])

	## get cumulative case counts for each strain and fit first model
	cum.cases <- matrix(0, ncol=n.strains, nrow=length(analysis.idx))
	for(i in 1:n.strains){
		cum.cases[,i] <- cumsum(dat[analysis.idx,case.cols[i]])
	}

	## start loop for EM-like algorithm
	eps <- 1
	curr.rho <- matrix(0, ncol=n.strains, nrow=length(analysis.idx))
	Q.terms <- matrix(0, ncol=n.strains, nrow=length(analysis.idx))
	Q.sums <- matrix(0, ncol=n.strains, nrow=length(analysis.idx))
	fits <- vector("list", n.strains)
	iter <- 1
	while(eps>tol){

		## A. fit model for each strain
		for(i in 1:n.strains){
			fits[[i]] <- reconstruct.susc.class(cum.cases=
							    cum.cases[,i],
							    cum.births=
							    cumbirths.lagged,
							    delta=delta, k=k,
							    offset=Q.sums[,i],
							    smooth=smooth,
							    n.smooth=n.smooth,
							    df.low=df.low,
							    df.high=df.high,
							    plot=FALSE)
			curr.rho[,i] <- fits[[i]]$rho.t
		}
		#matplot(curr.rho, type="l", lty=1)

		## B. Calculate Q terms using fitted rho_ti

		## calculate the Q terms for each strain
		if(is.null(lambda)){
			for(i in 1:n.strains) {
				Q.terms[,i] <- calc.offset(case.counts=
							   dat[case.idx,case.cols[i]],
							   rho.t=curr.rho[,i],
							   k=k)
			}
		} else {
			for(i in 1:n.strains) {
				Q.terms[,i] <- calc.exp.offset(case.counts=
							   dat[case.idx,case.cols[i]],
							   rho.t=curr.rho[,i],
							   k=k,
							   lambda=lambda)
			}
		}


		## calculate the offsets as the sum of all Q terms for j!=i
		for(i in 1:n.strains){
			Q.sums[,i] <- apply(Q.terms[,-i],
					    MAR=1,
					    FUN=sum)
		}
		Q.sums <- Q.sums*delta

		## check stopping criteria
		if(iter==n.iter) {
			print(paste("delta = ", round(delta,3),
				    "k = ", k,
				    "final tol = ", round(eps, 5),
				    sep=" :: "))
			break
		}

		if(iter==1) {
			eps <- 1
		} else {
			eps <- max((curr.rho-old.rho)^2)
		}

		old.rho <- curr.rho
		iter <- iter+1
	}


	## 5. Use Z_ti and rho_ti from final converged model to fit
	##    transmission equation.

	## set up final matrix for analysis
	n.obs <- nrow(curr.rho)
	rhos.final <- curr.rho
	colnames(rhos.final) <- paste("rho.str", 1:n.strains, sep="")
	Zt.final <- matrix(NA, ncol=n.strains, nrow=n.obs)
	for(i in 1:n.strains) Zt.final[,i] <- fits[[i]]$Zt
	colnames(Zt.final) <- paste("zt.str", 1:n.strains, sep="")

	dat.final <- data.frame(tail(dat, n.obs), rhos.final, Zt.final)

	final.data <- get.TSIR.data(dat.final,
				    case.cols=case.cols,
				    rho.cols=colnames(rhos.final),
				    zt.cols=colnames(Zt.final),
				    analysis.start.time=analysis.start.time)

	return(list(curr.rho=curr.rho,
		    cum.cases=cum.cases,
		    fits=fits,
		    iters=iter,
		    eps=eps,
		    dat.final=dat.final,
		    analysis.data=final.data$analysis.data,
		    nonzero.subset=final.data$nonzero.subset))
}

## function to reconstruct the susceptible class based on births
## and cumulative observed cases
## cum.cases = cumulative cases across the time-periods
## cum.births = cumulative births across the time-periods
## delta, k = parameters for cross-protection, described above
## offset = the Q terms specified by the cross-protection terms,
##          serve as offset for linear fits
## smooth = indicator of whether to use local linear smooths
## df.low, df.high = min and max values for the df for smoother
## nsmooth = number of smoothness levels to investigate
## plot = logical, whether to show plots or not
## tol = the tolerance for smooth.spline. defaults to smooth.spline default

reconstruct.susc.class <- function(cum.cases, cum.births, delta, k,
				   offset,
				   smooth=TRUE, n.smooth=10,
				   df.low=5, df.high=50,
				   plot=TRUE){
	## check lengths
	n.obs <- length(cum.cases)
	if( !(n.obs==length(cum.births) & n.obs==length(offset)) )
		stop("case & birth vectors should have length = length(offset)")

	## relabel, and select only
	X.cumcases <- cum.cases
	Y.cumbirths <- cum.births
	Q.offset <- offset

	## fix error if IQR = 0
	tol <- ifelse(IQR(X.cumcases)==0,
		      1e-6*median(cum.cases),
		      1e-6*IQR(cum.cases))

	## GLOBAL LINEAR FIT
	## equation 8 from Finkenstadt & Grenfell (2000)
	global.fit <- lm(Y.cumbirths ~ X.cumcases + offset(Q.offset))
	Zt.global <- resid(global.fit)
	yhat.global <- predict(global.fit)

	if(smooth){
		## LOCAL LINEAR FIT
		## pick a set of df's
		df <- seq(df.low, df.high, length.out=n.smooth)

		## run model for each bandwidth, calculate SSE1, SSE2
		sse1 <- sse2 <- rep(0, length(df))
		for(i in 1:length(df)){
			## note the offset "cheat" in model formula
			fit <- smooth.spline(x=X.cumcases,
					     y=Y.cumbirths-Q.offset,
					     df=df[i], tol=tol)
			## need to adjust fitted values for offset "cheat"
			yhat.local <- predict(fit, X.cumcases)$y + Q.offset
			sse1[i] <- sum((Y.cumbirths-yhat.local)^2)
			sse2[i] <- sum((yhat.global-yhat.local)^2)
		}

		## final model is one with smallest |SSE1-SSE2|
		which.fit <- which(min(abs(sse1-sse2)) == abs(sse1-sse2))
		fit <- smooth.spline(x=X.cumcases,
				     y=Y.cumbirths-Q.offset,
				     df=df[which.fit], tol=tol)
		Zt.local <- residuals(fit)
		if(plot){
			plot(df, sse1, type="b", ylim=range(sse1, sse2), ylab="SSEs")
			points(df, sse2, type="b")
			print(paste("chosen df =", df[which.fit]))
		}
	}

	## return appropriate counts, for all cumulative cases, even ones not fit
	if(smooth){
		Zt <- Zt.local
		rho.t <- predict(fit, cum.cases, deriv=1)$y
		fit <- fit
	} else {
		Zt <- Zt.global
		rho.t <- rep(coef(global.fit)[2], n.obs)
		fit <- global.fit
	}

	## note that length(offset)==length(Zt) but
	##           length(cum.cases)==length(rho.t)==length(offset)+k
	out <- list(cum.cases=cum.cases,
		    cum.births=cum.births,
		    offset=offset,
		    Zt=Zt,
		    rho.t=rho.t,
		    k=k,
		    fit=fit)
	return(out)
}


## this function calculates the Q/offset terms given case counts, rho.t and k
calc.offset <- function(case.counts, rho.t, k){

	if(length(case.counts)!=length(rho.t)+k)
		stop("length(case.counts) should = length(rho.t) + k")

	## impute rho.t[1] for all t<1
	impute.length <- length(case.counts) - length(rho.t)
	new.rho <- rep(rho.t[1], impute.length)
	new.rho <- c(new.rho, rho.t)


	scaled.cases <- case.counts*new.rho
	nobs <- length(rho.t)

	## the two sequences that make up the kernel of the offset term
	seq1.scaled.cases <- scaled.cases[(k+1):(nobs+k)]
	seq2.scaled.cases <- scaled.cases[1:nobs]

	Q.tmp <- cumsum(seq1.scaled.cases-seq2.scaled.cases)


	## return a list of Q terms corresponding to the last
	## length(cum.cases)-k case counts.
	return(Q.sums=Q.tmp)
}

## this function calculates the Q/offset terms given case counts, rho.t, lambda and k
calc.exp.offset <- function(case.counts, rho.t, k, lambda){

	## in terms of the parameters in this function, we are calculating
	## \sum_{t'=0}^{k+1} \rho_{t-t',j}C_{t-t',j} Pr(t'<T|T<k+1)

	if(length(case.counts)!=length(rho.t)+k)
		stop("length(case.counts) should = length(rho.t) + k")

	## impute rho.t[1] for all t<1
	impute.length <- length(case.counts) - length(rho.t)
	new.rho <- rep(rho.t[1], impute.length)
	new.rho <- c(new.rho, rho.t)

	scaled.cases <- case.counts*new.rho
	nobs <- length(rho.t)

	## the first sequence that contributes to the Q term
	seq1.scaled.cases <- scaled.cases[(k+1):(nobs+k)]

	## calculations for the second sequence that contributes to the Q term
	## define vector of exponential probabilities (should be length k, from 0 to k-1)
	## these are pr(v<X<v+1|X<k) = [pr(X<v+1) - pr(X<v)]/pr(X<k)
	max.prob <- pexp(k, 1/lambda)
	exp.probs <- (pexp(1:k, 1/lambda) - pexp(0:(k-1), 1/lambda)) / max.prob

	## for each t, sum last k case counts, including current t
	seq2.scaled.cases <- rep(NA, nobs)
	for(i in 1:nobs){
		idx <- (i+k):(i+1) ## index on case.idx scale i=1 => time k+1
		seq2.scaled.cases[i] <- sum(exp.probs*scaled.cases[idx])
	}

	Q.tmp <- cumsum(seq1.scaled.cases-seq2.scaled.cases)

	## return a list of Q terms corresponding to the last
	## length(cum.cases)-k case counts.
	return(Q.sums=Q.tmp)
}


## function to fit TSIR model
## dat = a matrix of the sort returned by reconstruct.susc.class
## case.cols = the column names for the case counts at time t
## rho.cols = the column names for the rho.t
## zt.cols = the column names for the z.t
## NOTE: the previous three must line up, i.e. the first item in case.cols and rho.cols and zt.cols must always refer to the first strain, the second to the second strain, etc...
## analysis.start.time = the first time at which the case counts should be used in fitting the transmission model

get.TSIR.data <- function(dat,
			  case.cols,
			  rho.cols,
			  zt.cols,
			  analysis.start.time){

	fit.subset <- 2:nrow(dat)

	## calculate case counts
	It.all <- as.matrix(dat[, case.cols]*dat[,rho.cols])
	It <- It.all[fit.subset,]

	## get autoregressive case counts
	It.minus.1 <- It.all[fit.subset-1,]

	## get autoregressive residuals
	Zt.minus.1 <- as.matrix(dat[fit.subset-1, zt.cols])

	## get biweeks
	biweeks <- dat[fit.subset,"biweek"]

	## make data into long format
	analysis.data <- NULL
	for(i in 1:length(case.cols)){
		str.dat <- cbind(It[,i],
				 It.minus.1[,i],
				 Zt.minus.1[,i],
				 biweeks,
				 str=i)
		analysis.data <- rbind(analysis.data, str.dat)
	}
	colnames(analysis.data) <- c("It", "It.minus.1",
				     "Zt.minus.1", "biweeks", "str")
	rownames(analysis.data) <- 1:nrow(analysis.data)

	## make subset removing <=0 case counts
	nonzero.subset <- which(analysis.data[,"It"]>0 &
				analysis.data[,"It.minus.1"]>0)

	## fit model
	analysis.subset <- analysis.data[nonzero.subset,]

	return(list(analysis.data=analysis.data,
		    nonzero.subset=nonzero.subset))
}



TSIR.post.estimation <- function(ll.store, delta.range, k.range,
				 DF, model.frmla, ll.col, cor.method="pearson",
				 dat, st.date, plot.cor=FALSE,
				 delta.plot.range=c(-1,1), ...){
	require(lattice)

	## calculate MLEs
	liks <- ll.store[,ll.col]
	ll.tmp <- cbind(ll.store, liks=liks)

	## calculate MLEs
	idx.to.keep <- (ll.tmp[,"delta.vec"]>=delta.range[1] &
			ll.tmp[,"delta.vec"]<=delta.range[2] &
			ll.tmp[,"k.vec"]>=k.range[1] &
			ll.tmp[,"k.vec"]<=k.range[2] )
	max.lik <- max(liks[idx.to.keep], na.rm=TRUE)
	mle.idx <- which(liks==max.lik)
	delta.mle <- ll.tmp[mle.idx,"delta.vec"]
	k.mle <- ll.tmp[mle.idx,"k.vec"]
	N.fit <- ll.tmp[mle.idx, "N.fit"]

	## calculate univariate CIs
	ci.limit <- max.lik - qchisq(.95, 2)/2
	ci.idx <- which(liks>ci.limit & idx.to.keep)
	delta.limits <- range(ll.tmp[ci.idx, "delta.vec"])
	k.limits <- range(ll.tmp[ci.idx, "k.vec"])

	## fit model at MLE
	bestfit.data <- get.cross.protect.data(dat,
					       case.cols=paste("den.cases.str", 1:4, sep=""),
					       delta=delta.mle,
					       cr.protect.lag=k.mle,
					       analysis.start.time=st.date,
					       birth.lag=birth.lag,
					       n.iter=100, tol=1e-5,
					       smooth=TRUE, n.smooth=1,
					       df.low=DF, df.high=DF)
	analysis.data <- bestfit.data$analysis.data[bestfit.data$nonzero.subset,]
	bestfit <- lm(model.frmla, data=data.frame(analysis.data))

	alpha1.mle <- coef(bestfit)["log(It.minus.1)"]
	alpha1.ci <- confint(bestfit)["log(It.minus.1)",]

        ## get 95% CI for average duration, i.e. "lambda"
        lik.limit <- max.lik-qchisq(.95, 2)/2
        ## identify grid points within the 95% conf region
        in.ci <- which(liks >= lik.limit & idx.to.keep)
        lambdas <- ll.tmp[,"delta.vec"] * ll.tmp[,"k.vec"]
        lambda.ci <- range(lambdas[in.ci])
        lambda.mle <- delta.mle*k.mle


	## plots
	if(plot.cor) par(mfrow=c(1,2))

	#plot.lik.surface(ll.store, ll.col=ll.col, delta.range=delta.range, k.range=k.range, ...)
	## get delta plot range
	plotting.subset <- which(ll.store[,"delta.vec"] >= delta.plot.range[1] &
				 ll.store[,"delta.vec"] <= delta.plot.range[2])
	ggplot.lik.surface(ll.store[plotting.subset,], ll.col=ll.col, N.fit=N.fit)

	## get correlation between observed and predicted
	yhat <- exp(predict(bestfit))
	rho.cols <- paste("rho.str", 1:4, sep="")
	dat.cols <- paste("den.cases.str", 1:4, sep="")
	rep.factors <- as.matrix(bestfit.data$dat.final[-1, rho.cols])
	rep.factors <- as.numeric(rep.factors)[bestfit.data$nonzero.subset]
	yhat.scaled <- yhat/rep.factors
	obs.y <- as.matrix(bestfit.data$dat.final[-1,dat.cols])
	obs.y <- as.numeric(obs.y)[bestfit.data$nonzero.subset]

	spCorr <- cor(obs.y, yhat.scaled, method=cor.method)
	if(plot.cor) {
		plot(jitter(obs.y,2), yhat.scaled, pch=19, cex=.5,
		     ylim=range(obs.y, yhat.scaled), xlim=range(obs.y, yhat.scaled),
		     xlab="observed case count (jittered)", ylab="predicted case count")
		abline(a=0, b=1, col="gray")
		par(mfrow=c(1,1))
	}



	## return list
	return(list(bestfit=bestfit,
		    bestfit.data=bestfit.data,
		    delta.mle=delta.mle,
		    delta.ci.95=delta.limits,
		    k.mle=k.mle/26,
		    k.ci.95=k.limits/26,
                    lambda.mle=lambda.mle/26,
                    lambda.ci.95=lambda.ci/26,
		    alpha1.mle=alpha1.mle,
		    alpha1.ci=alpha1.ci,
		    corr=spCorr,
		    logLik=logLik(bestfit),
		    df=attr(logLik(bestfit), "df")+DF*4+2))

}


TSIR.post.estimation.delta1 <- function(ll.store, k.range, delta.point=1,
					DF, model.frmla, ll.col, cor.method="pearson",
					dat, st.date, plot.cor=FALSE,
					delta.plot.range=c(-1,1), ...){
	require(ggplot2)

	## fix delta range to be the point 1.
	delta.range <- rep(delta.point,2)

	## calculate MLEs
	liks <- ll.store[,ll.col]
	ll.tmp <- cbind(ll.store, liks=liks)

	## calculate MLEs
	idx.to.keep <- (ll.tmp[,"delta.vec"]>=delta.range[1] &
			ll.tmp[,"delta.vec"]<=delta.range[2] &
			ll.tmp[,"k.vec"]>=k.range[1] &
			ll.tmp[,"k.vec"]<=k.range[2] )
	max.lik <- max(liks[idx.to.keep], na.rm=TRUE)
	mle.idx <- which(liks==max.lik)
	delta.mle <- ll.tmp[mle.idx,"delta.vec"]
	k.mle <- ll.tmp[mle.idx,"k.vec"]
	N.fit <- ll.tmp[mle.idx, "N.fit"]

	## calculate univariate CIs
	ci.limit <- max.lik - qchisq(.95, 2)/2
	ci.idx <- which(liks>ci.limit & idx.to.keep)
	delta.limits <- range(ll.tmp[ci.idx, "delta.vec"])
	k.limits <- range(ll.tmp[ci.idx, "k.vec"])

	## fit model at MLE
	bestfit.data <- get.cross.protect.data(dat,
					       case.cols=paste("den.cases.str", 1:4, sep=""),
					       delta=delta.mle,
					       cr.protect.lag=k.mle,
					       analysis.start.time=st.date,
					       birth.lag=birth.lag,
					       n.iter=100, tol=1e-5,
					       smooth=TRUE, n.smooth=1,
					       df.low=DF, df.high=DF)
	analysis.data <- bestfit.data$analysis.data[bestfit.data$nonzero.subset,]
	bestfit <- lm(model.frmla, data=data.frame(analysis.data))

	alpha1.mle <- coef(bestfit)["log(It.minus.1)"]
	alpha1.ci <- confint(bestfit)["log(It.minus.1)",]

        ## get 95% CI for average duration, i.e. "lambda"
        lik.limit <- max.lik-qchisq(.95, 2)/2
        ## identify grid points within the 95% conf region
        in.ci <- which(liks >= lik.limit & idx.to.keep)
        lambdas <- ll.tmp[,"delta.vec"] * ll.tmp[,"k.vec"]
        lambda.ci <- range(lambdas[in.ci])
        lambda.mle <- delta.mle*k.mle


	## plots
	if(plot.cor) par(mfrow=c(1,2))

	liks.plot <- data.frame(ll.tmp[idx.to.keep, c("k.vec", "liks")])
	liks.plot$years <- liks.plot$k.vec/26
	with(liks.plot,
	     plot(years, liks, type = "l",
	     xlab="average length of cross-protection (in years)",
	     ylab="log-likelihood", xlim=c(0,7)))
	lines(x=rep(k.mle,2)/26, y=range(liks.plot[,"liks"]), lty=2)
	lines(x=k.limits/26, y=rep(ci.limit, 2), lty=2, col="gray")
	lines(x=rep(k.limits[1]/26,2), y=c(min(liks.plot$liks), ci.limit), lty=2, col="gray")
	lines(x=rep(k.limits[2]/26,2), y=c(min(liks.plot$liks), ci.limit), lty=2, col="gray")


	## get correlation between observed and predicted
	yhat <- exp(predict(bestfit))
	rho.cols <- paste("rho.str", 1:4, sep="")
	dat.cols <- paste("den.cases.str", 1:4, sep="")
	rep.factors <- as.matrix(bestfit.data$dat.final[-1, rho.cols])
	rep.factors <- as.numeric(rep.factors)[bestfit.data$nonzero.subset]
	yhat.scaled <- yhat/rep.factors
	obs.y <- as.matrix(bestfit.data$dat.final[-1,dat.cols])
	obs.y <- as.numeric(obs.y)[bestfit.data$nonzero.subset]

	spCorr <- cor(obs.y, yhat.scaled, method=cor.method)
	if(plot.cor) {
		plot(jitter(obs.y,2), yhat.scaled, pch=19, cex=.5,
		     ylim=range(obs.y, yhat.scaled), xlim=range(obs.y, yhat.scaled),
		     xlab="observed case count (jittered)", ylab="predicted case count")
		abline(a=0, b=1, col="gray")
		par(mfrow=c(1,1))
	}



	## return list
	return(list(bestfit=bestfit,
		    bestfit.data=bestfit.data,
		    delta.mle=delta.mle,
		    delta.ci.95=delta.limits,
		    k.mle=k.mle/26,
		    k.ci.95=k.limits/26,
                    lambda.mle=lambda.mle/26,
                    lambda.ci.95=lambda.ci/26,
		    alpha1.mle=alpha1.mle,
		    alpha1.ci=alpha1.ci,
		    corr=spCorr,
		    logLik=logLik(bestfit),
		    df=attr(logLik(bestfit), "df")+DF*4+2))

}


ggplot.lik.surface <- function(ll.store, ll.col, N.fit, ...) {
	require(ggplot2)
	require(RColorBrewer)

	ll.tmp <- data.frame(deltas = ll.store[,"delta.vec"],
			     ks = ll.store[,"k.vec"],
			     logliks = ll.store[,ll.col],
			     N.fit=ll.store[,"N.fit"])

	## calculate MLEs
	idx.to.keep <- which(ll.tmp[,"N.fit"]==N.fit)
		#which(ll.tmp[,"deltas"]>=delta.range[1] &
		#	     ll.tmp[,"deltas"]<=delta.range[2] &
		#	     ll.tmp[,"ks"]>=k.range[1] &
		#	     ll.tmp[,"ks"]<=k.range[2])
	max.lik <- max(ll.tmp[idx.to.keep, "logliks"], na.rm=TRUE)
	mle.idx <- which(ll.tmp[,"logliks"]==max.lik)
	delta.mle <- ll.tmp[mle.idx,"deltas"]
	k.mle <- ll.tmp[mle.idx,"ks"]

	## equation from BolkerBook p217
	conf.contours <- max.lik-qchisq(c(.90, .95), 2)/2
	at.levels <- c(Inf, ceiling(max.lik), conf.contours,
		       ceiling(max.lik) - c(5, 7, 9, 11, 13),
					#c(5, 10, 15, 20, 30),
		       -Inf)
	lab.ats <- rev(round(at.levels, 1))[-1]
	cats <- cut(ll.tmp[,"logliks"], breaks=at.levels, labels=lab.ats)
	cats <- factor(cats, levels=rev(levels(cats)))
	cats1 <- cats
	cats1[-idx.to.keep] <- NA

	dd <- data.frame(ll.tmp, cats1=factor(cats1))
	pal <- brewer.pal(9, "Blues")[2:9]
	p <- ggplot(dd, aes(ks, deltas, fill=cats1)) +
		theme_bw() + geom_tile() + scale_fill_manual(values=pal, name="log-lik")+
			geom_point(x=k.mle, y=delta.mle, show_guide=FALSE, shape=3) +
				scale_x_continuous("k, the duration of cross-protection (in years)", breaks=c(0, 26, 52, 78, 104), labels=0:4) + scale_y_continuous(expression(delta))
	print(p)
}




plot.lik.surface <- function(ll.store, ll.col,
			       delta.range, k.range, ...) {

	ll.tmp <- data.frame(deltas = ll.store[,"delta.vec"],
			     ks = ll.store[,"k.vec"]/26,
			     logliks = ll.store[,ll.col])

	## calculate MLEs
	k.range <- k.range/26
	idx.to.keep <- which(ll.tmp[,"deltas"]>=delta.range[1] &
			     ll.tmp[,"deltas"]<=delta.range[2] &
			     ll.tmp[,"ks"]>=k.range[1] &
			     ll.tmp[,"ks"]<=k.range[2])
	max.lik <- max(ll.tmp[idx.to.keep, "logliks"], na.rm=TRUE)
	mle.idx <- which(ll.tmp[,"logliks"]==max.lik)
	delta.mle <- ll.tmp[mle.idx,"deltas"]
	k.mle <- ll.tmp[mle.idx,"ks"]

	## equation from BolkerBook p217
	conf.contours <- max.lik-qchisq(c(.90, .95), 2)/2
	at.levels <- c(max.lik+1, conf.contours, conf.contours-5*1:10, -Inf)
	col.palette <- gray(25:100/100)

	## plot likelihood surface
	print(
	      levelplot(ll.tmp$logliks~ll.tmp$ks*ll.tmp$deltas,
			aspect="fill",
			xlim=c(-.1, 4.1),
			ylim=c(-1.1,1.1),
			col.regions=col.palette, at=at.levels,
			## columns of levelplot correspond to rows of matrix
			xlab=expression(paste("duration of cross protection, ",
				italic(k), " (in years)")),
			## rows of correspond to columns of matrix
			ylab=expression(paste("fraction of population with cross protection, ", delta)),
			panel=function(..., at, contour, region) {
				panel.levelplot(..., at=at.levels,
						contour = FALSE,
						region = TRUE)
				panel.contourplot(..., at=conf.contours, lwd=2,
						  contour = TRUE,
						  region = FALSE,
						  labels=list(labels=c("90%", "95%"),
						              cex=.6),
						  label.style="mixed")}, ...)
	      )
	trellis.focus("panel", 1, 1, highlight=FALSE)
	  lpoints(k.mle, delta.mle, col="black", pch=4)
	trellis.unfocus()

}




plot.lik.surface.expmodels <- function(ll.store) {

	## calculate MLEs
	max.lik <- max(ll.store, na.rm=TRUE)
	mle.idx <- which(ll.store==max.lik)

	mle.lambda <- lambdas[mle.idx]
	ci.95.limit <- max.lik - qchisq(.95, 1)/2
	ci.99.limit <- max.lik - qchisq(.99, 1)/2
	ci.95.idx <- which(ll.store>ci.95.limit)
	ci.99.idx <- which(ll.store>ci.99.limit)
	lambda.ci.95 <- lambdas[range(ci.95.idx)]
	lambda.ci.99 <- lambdas[range(ci.99.idx)]

	## plot likelihood graph
	plot(lambdas, ll.store, type="l", lwd=2,
	     xlab="average length of cross protection (in years)",
	     ylab="log-likelihood", las=1, bty="n", xaxt="n")
	max.year <- floor(max(lambdas)/26)
	axis(1, at=seq(0, max.year*26, by=26), labels=0:max.year)
	lines(x=rep(mle.lambda,2), y=range(ll.store), lty=2)
	lines(x=lambda.ci.95, y=rep(ci.95.limit, 2), lty=2, col="gray")
	lines(x=rep(lambda.ci.95[1],2), y=c(min(ll.store), ci.95.limit), lty=2, col="gray")
	lines(x=rep(lambda.ci.95[2],2), y=c(min(ll.store), ci.95.limit), lty=2, col="gray")
	# lines(x=lambda.ci.99, y=rep(ci.99.limit, 2))

	return(round(c(mle.lambda/26, lambda.ci.95/26, max.lik),3))
}

TSIR.post.estimation.expmodels <- function(logliks, lambdas, DF, model.formula,
					   cor.method="pearson",
					   dat=bkk.dengue.cases,
					   st.date, plot.cor=FALSE,
					   k, model.name,
					   dat.col.names="den.cases.str"){
	require(lattice)

	## calculate MLEs
	max.lik <- max(logliks, na.rm=TRUE)
	mle.idx <- which(logliks==max.lik)

	mle.lambda <- lambdas[mle.idx]
	ci.95.limit <- max.lik - qchisq(.95, 1)/2
	ci.99.limit <- max.lik - qchisq(.99, 1)/2
	ci.95.idx <- which(logliks>ci.95.limit)
	ci.99.idx <- which(logliks>ci.99.limit)
	lambda.ci.95 <- lambdas[range(ci.95.idx)]
	lambda.ci.99 <- lambdas[range(ci.99.idx)]

	## plot likelihood graph
	if(plot.cor) par(mfrow=c(1,2))
	plot(lambdas, logliks, type="l", lwd=2,
	     xlab="average length of cross protection (in years)",
	     ylab="log-likelihood", las=1, bty="n", xaxt="n",
	     main=model.name)
	max.year <- floor(max(lambdas)/26)
	axis(1, at=seq(0, max.year*26, by=26), labels=0:max.year)
	lines(x=rep(mle.lambda,2), y=range(logliks), lty=2)
	lines(x=lambda.ci.95, y=rep(ci.95.limit, 2), lty=2, col="gray")
	lines(x=rep(lambda.ci.95[1],2), y=c(min(logliks), ci.95.limit), lty=2, col="gray")
	lines(x=rep(lambda.ci.95[2],2), y=c(min(logliks), ci.95.limit), lty=2, col="gray")

	## fit model at MLE
	bestfit.data <- get.cross.protect.data(dat,
					       case.cols=paste(dat.col.names, 1:4, sep=""),
					       delta=1,
					       cr.protect.lag=k,
					       lambda=mle.lambda,
					       analysis.start.time=st.date,
					       birth.lag=birth.lag,
					       n.iter=50, tol=1e-5,
					       smooth=TRUE, n.smooth=1,
					       df.low=DF, df.high=DF)
	analysis.data <- bestfit.data$analysis.data[bestfit.data$nonzero.subset,]
	bestfit <- lm(model.formula, data=data.frame(analysis.data))

	alpha1.mle <- coef(bestfit)["log(It.minus.1)"]
	alpha1.ci <- confint(bestfit)["log(It.minus.1)",]

	## get correlation between observed and predicted
	yhat <- exp(predict(bestfit))
	rho.cols <- paste("rho.str", 1:4, sep="")
	dat.cols <- paste(dat.col.names, 1:4, sep="")
	rep.factors <- as.matrix(bestfit.data$dat.final[-1, rho.cols])
	rep.factors <- as.numeric(rep.factors)[bestfit.data$nonzero.subset]
	yhat.scaled <- yhat/rep.factors
	obs.y <- as.matrix(bestfit.data$dat.final[-1,dat.cols])
	obs.y <- as.numeric(obs.y)[bestfit.data$nonzero.subset]

	spCorr <- cor(obs.y, yhat.scaled, method=cor.method)

	if(plot.cor) {
		plot(jitter(obs.y,2), yhat.scaled, pch=19, cex=.5,
		     ylim=range(obs.y, yhat.scaled), xlim=range(obs.y, yhat.scaled),
		     xlab="observed case count (jittered)", ylab="predicted case count")
		abline(a=0, b=1, col="gray")
		par(mfrow=c(1,1))
	}

	## return list
	return(list(bestfit=bestfit,
		    bestfit.data=bestfit.data,
		    lambda.mle=mle.lambda/26,
		    lambda.ci.95=lambda.ci.95/26,
		    alpha1.mle=coef(bestfit)["log(It.minus.1)"],
		    alpha1.ci.95=alpha1.ci,
		    corr=spCorr,
		    logLik=logLik(bestfit),
		    df=attr(logLik(bestfit), "df")+DF*4+1,
		    used.dat=dat,
		    used.case.cols=paste(dat.col.names, 1:4, sep=""),
		    used.k=k,
		    used.analysis.start.time=st.date,
		    used.birth.lag=birth.lag,
		    used.DF=DF))

}

trans.param.plot <- function(Bsum, Csum, Dsum){
	require(dichromat)
	r.ests <- matrix(nrow=134, ncol=5)
	colnames(r.ests) <- c("strain", "biweek", "mle", "cilow", "cihigh")

	model.factor <- rep(LETTERS[5:7], times=c(4, 26, 104))
	B.idx <- 1:4
	C.idx <- 5:30
	D.idx <- 31:134

	r.ests[B.idx,"strain"] <- 1:4
	r.ests[C.idx,"biweek"] <- 1:26
	r.ests[D.idx,"strain"] <- rep(1:4, each=26)
	r.ests[D.idx,"biweek"] <- rep(1:26, times=4)

	B.ests <- coef(Bsum$bestfit)[1:4]
	r.ests[B.idx,"mle"] <- B.ests - mean(B.ests)
	r.ests[B.idx,c("cilow", "cihigh")] <- confint(Bsum$bestfit)[1:4,] - mean(B.ests)

	C.ests <- coef(Csum$bestfit)[1:26]
	r.ests[C.idx,"mle"] <- C.ests - mean(C.ests)
	r.ests[C.idx,c("cilow", "cihigh")] <- confint(Csum$bestfit)[1:26,] - mean(C.ests)

	D.ests <- coef(Dsum$bestfit)[3:106]
	r.ests[D.idx,"mle"] <- D.ests - mean(D.ests)
	r.ests[D.idx,c("cilow", "cihigh")] <- confint(Dsum$bestfit)[3:106,] - mean(D.ests)

	r.ests <- data.frame(model=model.factor, r.ests)

	mycols <- colorschemes$Categorical.12[c(2,8,10,12)]


	## MAKE PLOT
	rng <- range(r.ests[C.idx,c("cilow", "cihigh")], r.ests[,"mle"])
	##layout(matrix(1:2, nrow=1), widths=c(1,3))
	## model B coefs
	##plot(r.ests[B.idx,"mle"], pch=4, col=mycols,
	##    bty="n", ylim=rng, xlim=c(.5, 4), las=1,
	##     xlab="", ylab="MLEs of transmission parameters", xaxt="n")
	##title("Model B")
	##axis(1, at=1:4, labels=FALSE)
	##for(i in 1:4){
	##  lines(c(i,i), r.ests[i,c("cilow", "cihigh")], col=mycols[i])
	##  axis(1, at=i, label=paste("DENV", i, sep=""), las=2, col.axis=mycols[i])
	##}
	##abline(h=1, col="grey", lty=2)
	## model C & D
	pchs <- c(3, 4, 8, 19)
	dates <- seq(as.Date("2010-01-01"), length.out=26, by="14 days")+7
	plot(dates, r.ests[C.idx,"mle"], ylim=rng, type="l", lwd=3,
	     xlab="", bty="n", ylab="", las=1)
	polygon(c(dates, rev(dates)), c(r.ests[C.idx,"cilow"], rev(r.ests[C.idx,"cihigh"])), col="lightgray", density=-1, border=FALSE)
	points(dates, r.ests[C.idx,"mle"], type="l", lwd=3)
	title("Estimated transmission parameters from model Ec")
	for(str in 1:4){
		tmp <- subset(r.ests, model=="G"&strain==str)
		##	points(tmp[,"mle"], col=mycols[str], type="l", cols=mycols[str])
		points(dates, tmp[,"mle"], col=mycols[str], pch=pchs[str], cex=.7)
		points(smooth.spline(dates, tmp[,"mle"], df=10),
		       type="l", lwd=1, col=mycols[str], lty=str)
	}
	legend(x=as.Date("2010-11-01"), y=rng[2], legend=paste("DENV", 1:4, sep=""), lty=1:4, col=mycols, bty="n", pch=pchs)

	abline(h=0, col="grey", lty=2)
	## return model C parameter estimates
	r.ests[C.idx,]
}


##################################
## FORWARD SIMULATION FUNCTIONS ##
##################################

## function to simulate a single dataset of case counts over the period of observed data.
## based on conversation with Aaron, Pej and Sourya on 6/22/2012
simulate.many.TSIR.datasets <- function(nsim, tsir.sum, all.data, plot=FALSE) {
	## plot "actual" data
	if(plot) {
		used.dat <- tsir.sum$bestfit.data$dat.final
		rho.cols <- c("rho.str1", "rho.str2", "rho.str3", "rho.str4")
		rhos <- apply(used.dat[,rho.cols], MAR=2, FUN=mean)
		plot(all.data$date, (all.data$den.cases.str1+.1)*rhos[1], type="l", lwd=2)
	}

	## simulate data and plot it
	simulated.data <- array(NA, dim=c(nrow(all.data), 12, nsim))
	dimnames(simulated.data) <- list(NULL,
					 c("biweek", "date", "births", "I1", "I2", "I3", "I4", "Z1", "Z2", "Z3", "Z4", "all.I"),
					 NULL)
	for(i in 1:nsim){
		simulated.data[,,i] <- simulate.TSIR.dataset(tsir.sum=tsir.sum, all.data=all.data)
		if(plot)
			lines(all.data$date, simulated.data[,"I1",i], col=rgb(.5, .5, .5, alpha=.01))
	}

	## add line on top again,for good measure
	if(plot)
		lines(all.data$date, (all.data$den.cases.str1+.1)*rhos[1], type="l", lwd=2)

	return(simulated.data)
}

## calculate spectra from simulated data, assuming output is from simulate.many.TSIR.datasets
## str can be a number 1-4 or "all" indicating the sum of all strains
plot.multiple.spectra <- function(simulated.data, all.data, tsir.sum, str=1, only.new.data=FALSE, col="black", taper=.25, plot.spectra=TRUE, plot.histogram=FALSE, ...) {
	require(splines)
	nsim <- dim(simulated.data)[3]
	new.col <- rgb(t(col2rgb(col))/255, alpha=.01)

	if(only.new.data){
		## get data used in analysis
		used.dat <- tsir.sum$bestfit.data$dat.final
		## find index of data predicted
		initial.data.idx <- which(all.data[,"date"] <= used.dat[1,"date"])
		idx <- (max(initial.data.idx)+1):nrow(all.data)
	} else {
		idx <- 1:nrow(all.data)
	}


	if(str=="all"){
		## create a sum of all cases column in real data matrix
		real.data.cols <- paste("den.cases.str", 1:4, sep="")
		all.data$all.cases <- rowSums(all.data[,real.data.cols])

		## define the columns for the spectra
		real.data.col <- "all.cases"
		sim.data.col <- "all.I"
	} else {
		## for one serotype
		real.data.col <- paste("den.cases.str", str, sep="")
		sim.data.col <- paste("I", str, sep="")
	}

	## plot detrended & normalized data spectrum
	fm1 <- lm(all.data[idx,real.data.col]~ns(date, 3), data=all.data[idx,])
	detr.dat <- resid(fm1)
	detr.dat <- (detr.dat-mean(detr.dat))/sd(detr.dat)
	data.spec <- spectrum(detr.dat, plot=FALSE, taper=taper)
	if(plot.spectra) {
		par(mar=c(3,3,1,1))
		plot((1/data.spec$freq)/26, data.spec$spec,#/max(data.spec$spec),
		     type="l", lty=1, lwd=1, ylab="", las=1.5,
		     col=rgb(t(col2rgb("black"))/255, alpha=.5),
		     bty="n", xlab="years", ...)

		## plot simulated spectrum
		for(i in 1:nsim){
			fm1 <- lm(simulated.data[idx, sim.data.col, i]~ns(date, 3),
				  data=data.frame(simulated.data[idx,,i]))
			detr.dat <- resid(fm1)
			detr.dat <- (detr.dat-mean(detr.dat))/sd(detr.dat)
			sim.spec <- spectrum(detr.dat, plot=FALSE, taper=taper)
					#sim.spec <- spectrum(simulated.data[idx, sim.data.col, i], plot=FALSE)
			lines((1/sim.spec$freq)/26, sim.spec$spec,#/max(sim.spec$spec),
			      type="l", col=new.col)
		}

	## reassert original data
		lines((1/data.spec$freq)/26, data.spec$spec,#/max(data.spec$spec),
		      lwd=1.5, col=rgb(t(col2rgb("black"))/255, alpha=.5))
	}

	if(plot.histogram){
		## plot histogram of spectral maxima
		maxima <- matrix(NA, nrow=nsim, ncol=4)
		## data.max.idx <- which(data.spec$spec==max(data.spec$spec))
		## data.max <- 1/data.spec$freq[data.max.idx]/26
		## set up all strains
		for(j in 1:4){
			sim.data.col <- paste("I", j, sep="")
			for(i in 1:nsim){
				fm1 <- lm(simulated.data[idx, sim.data.col, i]~ns(date, 3),
					  data=data.frame(simulated.data[idx,,i]))
				detr.dat <- resid(fm1)
				detr.dat <- (detr.dat-mean(detr.dat))/sd(detr.dat)
				sim.spec <- spectrum(detr.dat, plot=FALSE, taper=taper)
				max.idx <- which(sim.spec$spec==max(sim.spec$spec))
				maxima[i,j] <- 1/sim.spec$freq[max.idx]/26
			}
		}
		hist(as.numeric(maxima), las=1, breaks=500, main="", xlim=c(0, 14), ylim=c(0,10), probability=TRUE)
	}
}

plot.spectral.max.histogram <- function() {
	## plot simulated spectrum
	for(i in 1:nsim){
		fm1 <- lm(simulated.data[idx, sim.data.col, i]~ns(date, 3),
			  data=data.frame(simulated.data[idx,,i]))
		detr.dat <- resid(fm1)
		detr.dat <- (detr.dat-mean(detr.dat))/sd(detr.dat)
		sim.spec <- spectrum(detr.dat, plot=FALSE, taper=taper)
		#sim.spec <- spectrum(simulated.data[idx, sim.data.col, i], plot=FALSE)
		lines((1/sim.spec$freq)/26, sim.spec$spec,#/max(sim.spec$spec),
		      type="l", col=new.col)
	}

}


## tsir.sum is the object that comes out of the TSIR.post.estimation.expmodels function
simulate.TSIR.dataset <- function(tsir.sum, all.data=bkk.dengue.all.cases) {
	n.times <- nrow(all.data)

	used.dat <- tsir.sum$bestfit.data$dat.final

	## create a matrix that can hold all the data we need, but fill in just what we know:
	##   births for all weeks, case counts (rho*C) for the first L biweeks
	##   colnames <- biweek, date, I1, I2, I3, I4, Z1, Z2, Z3, Z4, births
	new.dat <- matrix(NA, ncol=12, nrow=n.times)

	colnames(new.dat) <- c("biweek", "date", "births", "I1", "I2", "I3", "I4", "Z1", "Z2", "Z3", "Z4", "all.I")

	## get average reporting rates
	rho.cols <- c("rho.str1", "rho.str2", "rho.str3", "rho.str4")
	rhos <- apply(used.dat[,rho.cols], MAR=2, FUN=mean)

	## retrieve initial conditions data from the complete data set
	initial.data.idx <- which(all.data[,"date"] <= used.dat[1,"date"])
	new.dat[,"biweek"] <- all.data[,"biweek"]
	new.dat[,"date"] <- all.data[,"date"]
	new.dat[,"births"] <- all.data[,"births"]
	new.dat[initial.data.idx, "I1"] <- (all.data[initial.data.idx,"den.cases.str1"]+.1)*rhos[1]
	new.dat[initial.data.idx, "I2"] <- (all.data[initial.data.idx,"den.cases.str2"]+.1)*rhos[2]
	new.dat[initial.data.idx, "I3"] <- (all.data[initial.data.idx,"den.cases.str3"]+.1)*rhos[3]
	new.dat[initial.data.idx, "I4"] <- (all.data[initial.data.idx,"den.cases.str4"]+.1)*rhos[4]
	new.dat[max(initial.data.idx), "Z1"] <- used.dat[1, "zt.str1"]
	new.dat[max(initial.data.idx), "Z2"] <- used.dat[1, "zt.str2"]
	new.dat[max(initial.data.idx), "Z3"] <- used.dat[1, "zt.str3"]
	new.dat[max(initial.data.idx), "Z4"] <- used.dat[1, "zt.str4"]


	## parameters for transmission equation: alpha_1, zeta, r*_t, sigma^2_epsilon
	fm <- tsir.sum$bestfit
	coefs <- coef(fm)
	s2.eps <- summary(fm)$sigma^2
	zeta <- coefs["Zt.minus.1"]
	alpha1 <- coefs["log(It.minus.1)"]
	rt.names <- paste("factor(biweeks)", 1:26, sep="")
	r.t <- exp(coefs[rt.names])

	## parameters for susceptible accounting
	delta <- tsir.sum$delta.mle
	k <- tsir.sum$k.mle*26

	## iterate for rest of the values
	I.cols <- c("I1", "I2", "I3", "I4")
	Z.cols <- c("Z1", "Z2", "Z3", "Z4")
	start.idx <- length(initial.data.idx)+1
	for(i in start.idx:n.times){
		curr.r <- r.t[new.dat[i, "biweek"]]
		new.dat[i, I.cols] <- get.new.It(I.tminus1=new.dat[i-1, I.cols],
						 Z.tminus1=new.dat[i-1, Z.cols],
						 rt=curr.r,
						 alpha1=alpha1,
						 zeta=zeta,
						 s2.eps=s2.eps)
		new.dat[i, Z.cols] <- get.new.Zt(i=i,
						 new.dat=new.dat,
						 birth.lag=birth.lag,
						 delta=delta,
						 k=k)
	}

	## sample cases from all observed cases
	sampled.dat <- new.dat
	for(i in 1:4){
		n.rows <- nrow(new.dat)
		## choosing rbeta(, 4, 45) b/c it has mean ~1/500
		## or, get sampling fractions for each timepoint based on rhos
		short.rhos <- used.dat[,rho.cols[i]] ## rhos for all simulated data
		n.times.needed <- n.rows-length(short.rhos)
		long.rhos <- c(rep(short.rhos[1], n.times.needed),
			       short.rhos)
		case.obs.probs <- 1/long.rhos #rbeta(n.rows, 10, 4900)
		sampled.dat[,I.cols[i]] <- rbinom(n.rows,
						  round(new.dat[,I.cols[i]]),
						  case.obs.probs)
	}
	sampled.dat[,"all.I"] <- rowSums(sampled.dat[,c("I1", "I2", "I3", "I4")])
	return(sampled.dat)
}


## functions for iterating the system equations
get.new.It <- function(I.tminus1, Z.tminus1, rt, alpha1, zeta, s2.eps) {
	rt*(I.tminus1^alpha1)*exp(zeta*Z.tminus1)*exp(rnorm(1, sd=sqrt(s2.eps)))
}
get.new.Zt <- function(i, new.dat, birth.lag=8, delta=1, k){
	I.cols <- c("I1", "I2", "I3", "I4")
	Z.cols <- c("Z1", "Z2", "Z3", "Z4")
	Z.tminus1 <- new.dat[i-1, Z.cols]
	B.tminusd <- new.dat[i-birth.lag, "births"]
	I.t <- new.dat[i, I.cols]
	I.net.loss <- rep(NA, 4)
	for(j in 1:4){ ## the following two lines are taken from the summand in eq 8 supplementary materials
		I.tminus1.j <- sum(new.dat[i-1, I.cols[-j]])
		I.tminuskplus1.j <- sum(new.dat[i-k-1, I.cols[-j]])
		I.net.loss[j] <- I.tminus1.j-I.tminuskplus1.j
	}
	return(Z.tminus1 + B.tminusd + - I.t - delta*I.net.loss)
}

## dat = matrix with k+1 rows, 9 columns (biweek, I[1:4], Z[1:4])
## b.t is an integer representing the births in week t
## r.t is a vector (length 26) with the fitted r*s
## xi is the fitted xi parameter
## alpha is alpha_1, the fitted exponent of I_{t-1}
## delta is the fitted delta parameter
## lambda is the fitted lambda parameter
single.step.TSIR <- function(dat, b.t, r.t, xi, alpha, delta, lambda, s2.eps=0){
	## get new infected count (adapted equation 16 in Finkenstadt 2000)
	epsilon <- rnorm(4, 0, sd=sqrt(s2.eps))
	idx <- nrow(dat) ## equivalent to "t-1" index
	biweek <- ifelse(dat[idx,1]==26, 1, dat[idx,1]+1)
	r <- r.t[biweek]
	new.I <- r * dat[idx,2:5]^alpha * exp(xi*dat[idx,6:9]) * exp(epsilon)

	## get new "susceptible" count (adapted equation 17 in Finkenstadt 2000)
	new.case.counts <- rbind(dat[,2:5], new.I)
	Q.sums <- calc.exp.offset.for.sim(new.case.counts, delta, lambda)
#	print(Q.sums)
	new.Z <- b.t + dat[idx,6:9] - new.I - Q.sums
	return(c(biweek, new.I, new.Z))
}

## case.counts is a matrix with 4 columns (the last k I[1:4] counts)
## this function calculates the Q/offset terms given case counts, delta and lambda
calc.exp.offset.for.sim <- function(case.counts, delta, lambda){
	nobs <- nrow(case.counts)

	## the first sequences that contribute to the Q term (sum of case counts, but serotype)
	seq1.cases <- case.counts[nobs,]

	## Calculations for the second sequence that contributes to the Q term
	k <- nobs
	## define vector of exponential probabilities (should be length k, from 0 to k-1)
	## these are pr(v<X<v+1|X<k) = [pr(X<v+1) - pr(X<v)]/pr(X<k)

	max.prob <- pexp(k, 1/lambda)

	exp.probs <- (pexp(1:k, 1/lambda) - pexp(0:(k-1), 1/lambda)) / max.prob

	## for each t, sum last k case counts, including current t
	scaled.cases <- exp.probs*case.counts
	seq2.cases <- apply(scaled.cases, FUN=sum, MAR=2)

	Q.sums <- delta*(seq1.cases-seq2.cases)
	Q.final <- rep(NA, 4)
	for(i in 1:4)
		Q.final[i] <- sum(Q.sums[-i])

	return(Q.final)
}

## initial.dat =  matrix with k+1 rows, 9 columns (biweek, I[1:4], Z[1:4])
full.TSIR.sim <- function(niter, initial.dat, b.t, r.t, xi, alpha, delta, lambda, s2.eps=0) {
	n.start <- nrow(initial.dat)
	data <- matrix(NA, ncol=9, nrow=niter+n.start)
	colnames(data) <- c("biweek",
			    paste("I", 1:4, sep=""),
			    paste("Z", 1:4, sep=""))

	data[1:n.start,] <- initial.dat

	for(t in (1+n.start):(niter+n.start)){
		## print(t)
		idx <- (t-n.start):(t-1)
		tmp.dat <- data[idx,]
		data[t,] <- single.step.TSIR(dat=tmp.dat,
					     b.t=b.t,
					     r.t=r.t,
					     xi=xi,
					     alpha=alpha,
					     delta=delta,
					     lambda=lambda,
					     s2.eps=s2.eps)
	}
	return(data)
}


## tsir.sum is an object returned from one of the TSIR post.estimation functions
## tsir.data is the object returned from get.cross.protect.data()
## formula is the model formula used to fit the data
## type indicates whether the model used fixed duration or exponential cross-protection
## seasonality is binary indicator of inclusion in model
## data.years = number of years of data to use as initial condition
## sim.years = number of years to simulate
## s2.eps = variance of epsilon (simulation noise term)
## plot.log = plot on log scale or not
## nskip = number of biweeks to skip when calculating spectral density to adjust for transition
## bootstrap = indicator of whether bootstrapped intervals for spectra are calculated
## n.boots = number of boostrap iterations

TSIR.simulate.forward <- function(tsir.sum,
				  formula, type=c("fd", "exp", "none"),
				  seasonality, data.years=20,
				  sim.years=50, s2.eps=0, plot.log=FALSE,
				  plot=TRUE, plot.ts=FALSE,
				  dat.col.names="den.cases.str",
				  nskip=1,
				  bootstrap=FALSE, n.boots=10){
        op <- par(no.readonly = TRUE)
	## check for interaction models -- don't know how to deal with these yet
	if(any(grepl(":", formula)))
		stop("can't deal with interaction models yet.")

	type <- match.arg(type)
	tsir.data <- tsir.sum$bestfit.data
	cols <- paste(dat.col.names, 1:4, sep="")

	## fit the accounted data to the linear model
	data.used <- data.frame(tsir.data$analysis.data[tsir.data$nonzero.subset,])
	tsir.fit <- lm(formula, data=data.used)

	## pull out the coefficients
	coefs <- coef(tsir.fit)
	hat.alpha <- tsir.sum$alpha1.mle
	hat.xi <- coefs["Zt.minus.1"]
	if(seasonality){
		rt.names <- paste("factor(biweeks)", 1:26, sep="")
		hat.r.t <- coefs[rt.names]
	} else {
		hat.r.t <- rep(1, 26)
	}
	if(type=="exp") {
		hat.lambda <- tsir.sum$lambda.mle
		hat.delta <- 1
	}
	if(type=="none") {
		hat.lambda <- 0
		hat.delta <- 0
	}
        if(type=="fd") {
		hat.delta <- tsir.sum$delta.mle
		hat.k <- tsir.sum$k.mle
	}

	## get data from last K biweeks of real data
	births <- 10000
	K <- 26*data.years ## 20 years worth of data
	nn <- nrow(tsir.data$dat.final)
	idx <- (nn-K+1):nn
	initial.cases <- tsir.data$dat.final[idx,]
	I1 <- (initial.cases[,cols[1]]+1)*initial.cases$rho.str1
	I2 <- (initial.cases[,cols[2]]+1)*initial.cases$rho.str2
	I3 <- (initial.cases[,cols[3]]+1)*initial.cases$rho.str3
	I4 <- (initial.cases[,cols[4]]+1)*initial.cases$rho.str4
	Z1 <- initial.cases$zt.str1
	Z2 <- initial.cases$zt.str2
	Z3 <- initial.cases$zt.str3
	Z4 <- initial.cases$zt.str4

	## simulate forward
	initial.dat <- cbind(initial.cases$biweek, I1, I2, I3, I4, Z1, Z2, Z3, Z4)
	niter <- sim.years*26
	sim.data <- full.TSIR.sim(niter, initial.dat, r.t=hat.r.t,
				  xi=hat.xi, alpha=hat.alpha, b.t=births,
				  delta=hat.delta, lambda=hat.lambda, s2.eps=s2.eps)
	test <- data.frame(time=1:nrow(sim.data), sim.data)
	sim.data.noCP <- full.TSIR.sim(niter, initial.dat, r.t=hat.r.t,
				       xi=hat.xi, alpha=hat.alpha, b.t=births,
				       delta=0, lambda=0, s2.eps=s2.eps)
	test.noCP <- data.frame(time=1:nrow(sim.data), sim.data.noCP)

	if(plot){
	## plot time series
		if(plot.ts)
			display.mat <- matrix(c(1:10), nrow=5, byrow=FALSE)
		else
			display.mat <- matrix(c(1:5), nrow=5)
		par(mar=c(2,4,1,2))
		layout(display.mat)
		if(plot.log){
			plot.fun <- function(x) log(x)
		} else {
			plot.fun <- function(x) (x)
		}
		rge <- range(plot.fun(test$I1), plot.fun(test$I2), plot.fun(test$I3), plot.fun(test$I4))/1000
		I.sums <- rowSums(test[,c("I1", "I2", "I3", "I4")])
		I.sums.noCP <- rowSums(test.noCP[,c("I1", "I2", "I3", "I4")])
		if(plot.ts) {
			plot(test$time/26, plot.fun(test$I1)/1000, type="h", ylab="serotype 1", bty="n", las=1, ylim=rge, col="red")
			abline(v=nrow(initial.dat)/26, lty=2, col="gray")
			plot(test$time/26, plot.fun(test$I2)/1000, type="h", ylab="serotype 2", bty="n", las=1, ylim=rge, col="red")
			abline(v=nrow(initial.dat)/26, lty=2, col="gray")
			plot(test$time/26, plot.fun(test$I3)/1000, type="h", ylab="serotype 3", bty="n", las=1, ylim=rge, col="red")
			abline(v=nrow(initial.dat)/26, lty=2, col="gray")
			plot(test$time/26, plot.fun(test$I4)/1000, type="h", ylab="serotype 4", bty="n", las=1, ylim=rge, col="red")
			abline(v=nrow(initial.dat)/26, lty=2, col="gray")
			plot(test$time/26, I.sums/1000, type="h", ylab="all serotypes", bty="n", las=1, col="red")
			I.sums.noCP.w0 <- as.numeric(matrix(c(rep(0, length(I.sums.noCP)), I.sums.noCP/1000), ncol=length(I.sums.noCP), byrow=TRUE))
			lines(rep(test$time/26, each=2), I.sums.noCP.w0, col=rgb(0,0,1, alpha=.5))
			abline(v=nrow(initial.dat)/26, lty=2, col="gray")
		}

		## move bootstrapping so it only has to happen once
		## downsample the resulting simulations before calculating the spectra to add additional uncertainty.


		## plot and calculate spectral density
		data.idx <- 1:nrow(initial.dat)
		sim.idx <- (nrow(initial.dat)+nskip):nrow(test)
		par(mar=c(2,3,1,2))
		for(i in 1:4){
			col <- paste("I", i, sep="")
			data.spec <- spectrum(test[data.idx, col], plot=FALSE)
			sim.spec <- spectrum(test[sim.idx, col], plot=FALSE)
			sim.spec.noCP <- spectrum(test.noCP[sim.idx, col], plot=FALSE)
			y.range <- c(0, 1)
			x.range <- c(0, max(c(1/sim.spec.noCP$freq,
					      1/sim.spec$freq,
					      1/data.spec$freq)/26))
			plot((1/data.spec$freq)/26, data.spec$spec/max(data.spec$spec),
			     type="l", lty=1, lwd=2, ylab="", las=1,
			     xlim=x.range, ylim=y.range, col="black", bty="n")
			## add in bootstrapped spectra
			if(bootstrap){
				booted.spectra <- run.bootstrap.spectra(tsir.sum, formula, n.boots,
									births, nfreq=length(sim.spec$freq),
									seasonality, data.years, sim.years, s2.eps, nskip, j.col=col)
				for(j in 1:n.boots){
					lines((1/sim.spec$freq)/26, booted.spectra[,j],
					      type="l", lwd=1, lty=1, col="gray70")
				}
			}
			## add simulated spectra on top of other data
			lines((1/sim.spec$freq)/26, sim.spec$spec/max(sim.spec$spec),
			      type="l", lwd=2, lty=2, col="red")
			lines((1/sim.spec.noCP$freq)/26, sim.spec.noCP$spec/max(sim.spec.noCP$spec),
			      type="l", lwd=2, lty=3, col="blue")


		}
		data.spec <- spectrum(I.sums[data.idx], plot=FALSE)
		sim.spec <- spectrum(I.sums[sim.idx], plot=FALSE)
		sim.spec.noCP <- spectrum(I.sums.noCP[sim.idx], plot=FALSE )
		y.range <- c(0, 1)
		x.range <- c(0, max(c(1/sim.spec.noCP$freq,
				      1/sim.spec$freq,
				      1/data.spec$freq)/26))
		plot((1/data.spec$freq)/26, data.spec$spec/max(data.spec$spec),
		     type="l", lty=1, lwd=2, bty="n", las=1,
		     xlim=x.range, ylim=y.range, col="black", ylab="")
		lines((1/sim.spec$freq)/26, sim.spec$spec/max(sim.spec$spec),
		      type="l", lwd=2, lty=2, col="red")
		lines((1/sim.spec.noCP$freq)/26, sim.spec.noCP$spec/max(sim.spec.noCP$spec),
		      type="l", lwd=2, lty=3, col="blue")
	}
	par(op)

	return(list(initial.dat=initial.dat,
		    sim.data=sim.data,
		    sim.data.noCP=sim.data.noCP))
}


run.bootstrap.spectra <- function(tsir.sum, tsir.formula, n.boots, births, nfreq, seasonality, data.years, sim.years, s2.eps, nskip, j.col){
	## draw the values of lambda at which bootstraps will be run
	#lambda.range <- tsir.sum$lambda.ci.95
	#lambda.draw <- runif(n.boots, min=lambda.range[1], max=lambda.range[2])
	## turn lambda sampling off
	#lambda.draw <- rep(tsir.sum$lambda.mle, n.boots)

	## make storage
	booted.spectra <- matrix(NA, nrow=nfreq, ncol=n.boots)
	colnames(booted.spectra) <- paste("rep", 1:n.boots, sep="")

	used.dat <- tsir.sum$used.dat

	## run the bootstraps
	for(j in 1:n.boots){
		## bootstrap data
		case.cols <- paste("den.cases.str", 1:4, sep="")

		## line to activate bootstrapping
		boot.data <- bootstrap.serotyped.cases(used.dat[,case.cols])

		## line to deactivate bootstrapping
		#boot.data <- used.dat

		new.dat <- data.frame(biweek=used.dat$biweek,
				      date=used.dat$date,
				      den.cases.str1=boot.data[, "den.cases.str1"],
				      den.cases.str2=boot.data[, "den.cases.str2"],
				      den.cases.str3=boot.data[, "den.cases.str3"],
				      den.cases.str4=boot.data[, "den.cases.str4"],
				      births=used.dat$births)
		## fit model at a sequence of 10 lambdas
		n.lambdas <- 10
		lambda.seq <- seq(1, 3, length.out=n.lambdas)*26
		bestfits <- vector("list", n.lambdas)
		logliks <- rep(NA, n.lambdas)
		for(i in 1:n.lambdas){
			bestfits[[i]] <- get.cross.protect.data(new.dat,
								case.cols=tsir.sum$used.case.cols,
								delta=1,
								cr.protect.lag=tsir.sum$used.k,
								lambda=lambda.seq[i],
								analysis.start.time=tsir.sum$used.analysis.start.time,
								birth.lag=tsir.sum$used.birth.lag,
								n.iter=50, tol=1e-5,
								smooth=TRUE, n.smooth=1,
								df.low=tsir.sum$used.DF, df.high=tsir.sum$used.DF)
			analysis.data <- bestfits[[i]]$analysis.data[bestfits[[i]]$nonzero.subset,]
			## FIT MODEL WITH MCMC, FOR POSTERIOR DIST'N?
			tsir.boot.fit <- lm(tsir.formula, data=data.frame(analysis.data))
			logliks[i] <- logLik(tsir.boot.fit)
		}
		## fit model at the best lambda
		best.idx <- which(logliks==max(logliks))
		best.lambda <- lambda.seq[best.idx]
		message(paste("best lambda =", round(best.lambda/26,1)))
		bestfit.data <- bestfits[[best.idx]]
		analysis.data <- bestfit.data$analysis.data[bestfit.data$nonzero.subset,]
		## FIT MODEL WITH MCMC, FOR POSTERIOR DIST'N?
		tsir.boot.fit <- lm(tsir.formula, data=data.frame(analysis.data))

		## pull out the coefficients
		coefs <- coef(tsir.boot.fit)
		hat.alpha <- coefs["log(It.minus.1)"]
		hat.xi <- coefs["Zt.minus.1"]
		if(seasonality){
			rt.names <- paste("factor(biweeks)", 1:26, sep="")
			hat.r.t <- coefs[rt.names]
		} else {
			hat.r.t <- rep(1, 26)
		}

		K <- 26*data.years ## 20 years worth of data
		nn <- nrow(bestfit.data$dat.final)
		idx <- (nn-K+1):nn
		cols <- tsir.sum$used.case.cols
		initial.cases <- bestfit.data$dat.final[idx,]
		## create bootstrap samples of initial case counts
		## currently not bootstrapping these.
		I1.boot <- initial.cases[,cols[1]]
		I2.boot <- initial.cases[,cols[2]]
		I3.boot <- initial.cases[,cols[3]]
		I4.boot <- initial.cases[,cols[4]]
		I1 <- (I1.boot)*initial.cases$rho.str1
		I2 <- (I2.boot)*initial.cases$rho.str2
		I3 <- (I3.boot)*initial.cases$rho.str3
		I4 <- (I4.boot)*initial.cases$rho.str4
		Z1 <- initial.cases$zt.str1
		Z2 <- initial.cases$zt.str2
		Z3 <- initial.cases$zt.str3
		Z4 <- initial.cases$zt.str4
		initial.dat <- cbind(initial.cases$biweek, I1, I2, I3, I4, Z1, Z2, Z3, Z4)
		niter <- sim.years*26

		sim.boot.data <- full.TSIR.sim(niter, initial.dat, r.t=hat.r.t,
					       xi=hat.xi, alpha=hat.alpha, b.t=births,
					       delta=1, lambda=best.lambda, s2.eps=s2.eps)
		## account for reporting rate
		n.times <- nrow(sim.boot.data)
		I1.sampling.prob <- 1/mean(initial.cases$rho.str1)
		I2.sampling.prob <- 1/mean(initial.cases$rho.str2)
		I3.sampling.prob <- 1/mean(initial.cases$rho.str3)
		I4.sampling.prob <- 1/mean(initial.cases$rho.str4)
		I1.downsampled <- rbinom(n.times, round(sim.boot.data[,"I1"]), I1.sampling.prob)
		I2.downsampled <- rbinom(n.times, round(sim.boot.data[,"I2"]), I2.sampling.prob)
		I3.downsampled <- rbinom(n.times, round(sim.boot.data[,"I3"]), I3.sampling.prob)
		I4.downsampled <- rbinom(n.times, round(sim.boot.data[,"I4"]), I4.sampling.prob)
		test <- data.frame(time=1:nrow(sim.boot.data),
				   I1=I1.downsampled,
				   I2=I2.downsampled,
				   I3=I3.downsampled,
				   I4=I4.downsampled)
		## end of reporting rate accounting

		## old code, before reporting rate accounting added
		## test <- data.frame(time=1:nrow(sim.boot.data), sim.boot.data)
		sim.idx <- (nrow(initial.dat)+nskip):nrow(test)
		sim.spec <- spectrum(test[sim.idx, j.col], plot=FALSE)
		booted.spectra[,j] <- sim.spec$spec/max(sim.spec$spec)
		print(paste(j, "done", Sys.time()))
	}
	return(booted.spectra)
}

## this function takes a bootstrap sample with replacement of cases from a time series
## cases is a vector of integer case counts
bootstrap.sample.time.series <- function(cases) {
	n.times <- length(cases)
	disagg.cases <- rep(1:n.times, times=cases)
	sampled.cases <- sample(disagg.cases, size=length(disagg.cases), replace=TRUE)
	agg.cases <- rep(0, n.times)

	for(i in 1:n.times) agg.cases[i] <- sum(sampled.cases==i)
	return(agg.cases)
}

## this function takes a bootstrap sample with replacement of cases from a time series
## cases is a matrix with four (ordered) columns of integer case counts
bootstrap.serotyped.cases <- function(cases) {
	n.times <- nrow(cases)
	## disaggregate serotype-specific case counts
	disagg.cases.sty1 <- rep(1:n.times, times=cases[,1])
	disagg.cases.sty2 <- rep(1:n.times, times=cases[,2])
	disagg.cases.sty3 <- rep(1:n.times, times=cases[,3])
	disagg.cases.sty4 <- rep(1:n.times, times=cases[,4])
	## set up matrix to sample from
	sty.case.counts <- apply(cases, FUN=sum, MAR=2)
	stys <- rep(1:4, times=sty.case.counts)
	disagg.cases.sty.all <- c(disagg.cases.sty1, disagg.cases.sty2,
				  disagg.cases.sty3, disagg.cases.sty4)
	disagg.cases.all <- cbind(stys, disagg.cases.sty.all)
	sampled.rows <- sample(1:nrow(disagg.cases.all), size=nrow(disagg.cases.all), replace=TRUE)
	sampled.cases <- disagg.cases.all[sampled.rows,]
	agg.cases <- table(sampled.cases[,2], sampled.cases[,1])
	## include times with no cases, using "zoo" package
	agg.cases.zoo <- zoo(agg.cases, as.integer(rownames(agg.cases)))
	agg.cases.full <- merge(zoo(,1:n.times), agg.cases.zoo)
	agg.cases.full[is.na(agg.cases.full)] <- 0
	colnames(agg.cases.full) <- colnames(cases)
	return(agg.cases.full)
}

#####################
### OLD FUNCTIONS  ##
#####################

TSIR.post.estimation.old <- function(ll.store, ks, deltas,
				 good.cols=1:ncol(ll.store),
				 good.rows=1:nrow(ll.store),
				 DF, model.frmla, cor.method="pearson",
				 dat=bkk.dengue.cases, st.date=as.Date("1979-01-01")){
	require(lattice)

	## calculate MLEs
	aa <- max(ll.store[good.rows,good.cols], na.rm=TRUE)
	mle.idx <- which(ll.store==aa, arr.ind=TRUE)
	delta.mle <- deltas[mle.idx[1,2]]
	k.mle <- ks[mle.idx[1,1]]

	## calculate univariate CIs
	ci.limit <- aa - qchisq(.95, 2)/2
	ci.idx <- which(ll.store>ci.limit, arr.ind=TRUE)
	ci.idx <- subset(ci.idx, ci.idx[,"col"] %in% good.cols)
	ci.idx <- subset(ci.idx, ci.idx[,"row"] %in% good.rows)
	delta.limits <- deltas[range(ci.idx[,2])]
	k.limits <- ks[range(ci.idx[,1])]

	## fit model at MLE
	bestfit.data <- get.cross.protect.data(dat,
					       case.cols=paste("den.cases.str", 1:4, sep=""),
					       delta=delta.mle,
					       cr.protect.lag=k.mle,
					       analysis.start.time=st.date,
					       birth.lag=birth.lag,
					       n.iter=50, tol=1e-5,
					       smooth=TRUE, n.smooth=1,
					       df.low=DF, df.high=DF)
	analysis.data <- bestfit.data$analysis.data[bestfit.data$nonzero.subset,]
	bestfit <- lm(model.frmla, data=data.frame(analysis.data))

	## get correlation between observed and predicted
	yhat <- exp(predict(bestfit))
	rho.cols <- paste("rho.str", 1:4, sep="")
	dat.cols <- paste("den.cases.str", 1:4, sep="")
	rep.factors <- as.matrix(bestfit.data$dat.final[-1, rho.cols])
	rep.factors <- as.numeric(rep.factors)[bestfit.data$nonzero.subset]
	yhat.scaled <- yhat/rep.factors
	obs.y <- as.matrix(bestfit.data$dat.final[-1,dat.cols])
	obs.y <- as.numeric(obs.y)[bestfit.data$nonzero.subset]

	spCorr <- cor(obs.y, yhat.scaled, method=cor.method)
	plot(jitter(obs.y,2), yhat.scaled, pch=19, cex=.5,
	     ylim=range(obs.y, yhat.scaled), xlim=range(obs.y, yhat.scaled),
	     xlab="observed case count (jittered)", ylab="predicted case count")
	abline(a=0, b=1, col="gray")


	## return list
	return(list(bestfit=bestfit,
		    delta.mle=delta.mle,
		    delta.ci.limits=delta.limits,
		    k.mle=k.mle,
		    k.ci.limits=k.limits,
		    alpha1.mle=coef(bestfit)["log(It.minus.1)"],
		    corr=spCorr,
		    logLik=logLik(bestfit),
		    df=attr(logLik(bestfit), "df")+DF*4+2))

}


ggplot.lik.surface.old <- function(ll.store,
			       good.cols=1:ncol(ll.store),
			       good.rows=1:nrow(ll.store),
			       tile.width=1, ...) {

	## melt data
	colnames(ll.store) <- deltas
	rownames(ll.store) <- ks
	ll.tmp <- melt(ll.store)
	ll.tmp <- cbind(ll.tmp, tile.width)
	names(ll.tmp) <- c("k", "delta", "ll", "tile.width")

	## calculate MLEs
	max.lik <- max(ll.store[good.rows,good.cols], na.rm=TRUE)
	mle.idx <- which(ll.store==max.lik, arr.ind=TRUE)
	delta.mle <- deltas[mle.idx[1,2]]
	k.mle <- ks[mle.idx[1,1]]

	## equation from BolkerBook p217
	conf.contours <- max.lik-qchisq(c(.90, .95), 2)/2
	at.levels <- round(max.lik) - 5*0:10
#	at.levels <- c(max.lik, conf.contours, conf.contours-5*1:10)
	lim <- c(min(10 * floor(at.levels/10)), max(10*ceiling(at.levels/10)))
	col.palette <- gray(seq(0, 100, length.out=5)/100)

	## plot likelihood surface
	p <- ggplot(ll.tmp, aes(k, delta, fill=ll, z=ll)) +
		theme_bw() + geom_tile(aes(width=tile.width)) +
			scale_fill_gradientn(col.palette, breaks=at.levels,
					     limits=lim, name="log-likelihood", labels=round(at.levels,2)) +
				stat_contour(breaks=conf.contours, colour="black") +
					geom_point(x=k.mle, y=delta.mle, legend=FALSE) +
						scale_x_continuous("duration of cross-protection (in years)", breaks=c(0, 26, 52, 78, 104), labels=0:4) +
							scale_y_continuous(expression(delta))

	print(p)
}



plot.lik.surface.old <- function(ll.store,
			     good.cols=1:ncol(ll.store),
			     good.rows=1:nrow(ll.store), ...) {

	## calculate MLEs
	max.lik <- max(ll.store[good.rows,good.cols], na.rm=TRUE)
	mle.idx <- which(ll.store==max.lik, arr.ind=TRUE)
	delta.mle <- deltas[mle.idx[1,2]]
	k.mle <- ks[mle.idx[1,1]]

	## equation from BolkerBook p217
	conf.contours <- max.lik-qchisq(c(.90, .95), 2)/2

	at.levels <- c(max.lik+1, conf.contours, conf.contours-5*1:10, -Inf)

	col.palette <- gray(25:100/100)

	## plot likelihood surface
	print(
	      levelplot(ll.store, row.values=ks, column.values=deltas, aspect="fill",
			ylim=c(-1.05,1.05),
			col.regions=col.palette, at=at.levels,
			xlab="k", ## columns of levelplot correspond to rows of matrix
			ylab=expression(delta), ## rows of correspond to columns of matrix
			panel=function(..., at, contour, region) {
				panel.levelplot(..., at=at.levels,
						contour = FALSE,
						region = TRUE)
				panel.contourplot(..., at=conf.contours, lwd=2,
						  contour = TRUE,
						  region = FALSE,
						  labels=list(labels=c("90%", "95%"),
						              cex=.6),
						  label.style="mixed")},
			...)
	      )
	trellis.focus("panel", 1, 1, highlight=FALSE)
	  lpoints(k.mle, delta.mle, col="black", pch=4)
	trellis.unfocus()

}
