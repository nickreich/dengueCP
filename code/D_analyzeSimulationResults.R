## analyze and process simulation results from MI data
## Nicholas Reich
## February 2012

library(ggplot2)

## starting from scratch
n.files <- 12000
chosen.ests.final <- read.csv("../data/chosenEsts_MIsims5_June2012.csv")
key <- read.csv("../data/simulatedData_fromMI/key_5.csv", row.names=1)
true.lambdas <- 1/key[,"delta"]

true.lambdas.verb <- rep("3 years", n.files)
param.cats <- rep("", n.files)
param.cats.level <- rep(0, n.files)
subparam.cats <- rep("", n.files)
subparam.cats.level <- param.cats.level
ci.coverage <- rep(NA, n.files)
for(i in 1:n.files){
	if(true.lambdas[i]==2) true.lambdas.verb[i] <- "2 years"
	if(true.lambdas[i]>1 & true.lambdas[i]<2) true.lambdas.verb[i] <- "1.5 years"
	if(true.lambdas[i]==1) true.lambdas.verb[i] <- "1 year"
	if(true.lambdas[i]==.5) true.lambdas.verb[i] <- "6 months"
	if(true.lambdas[i]<.5) true.lambdas.verb[i] <- "1 day"
	## get CI coverage probabilitites
	## assumed that if model N is chosen with 1 day CP = successful coverage
	## if chosen model is N, then ci coverage successful if true duration==
	if(substr(chosen.ests.final[i, "model"], 1, 1) =="N"){
		ci.coverage[i] <- true.lambdas[i]<.5
	} else {
		ci.coverage[i] <- true.lambdas[i]<=(chosen.ests.final[i,"lam.ci.u"]/26) & true.lambdas[i]>=(chosen.ests.final[i,"lam.ci.l"]/26)
	}
	enh <- ifelse(key[i,1]==1, "No E", "E")
	rr <- key[i,3]*100
	param.cats[i] <- paste(true.lambdas.verb[i], ", ", enh, ", ", rr, "%", sep="")
	## create re-ordering of factors
	param.cats.level[i] <- switch(param.cats[i],
				      "1 day, No E, 1%" = 1,
				      "1 day, No E, 10%" = 2,
				      "1 day, E, 1%" = 3,
				      "1 day, E, 10%" = 4,
				      "6 months, No E, 1%" = 5,
				      "6 months, No E, 10%" = 6,
				      "6 months, E, 1%" = 7,
				      "6 months, E, 10%" = 8,
				      "1 year, No E, 1%" = 9,
				      "1 year, No E, 10%" = 10,
				      "1 year, E, 1%" = 11,
				      "1 year, E, 10%" = 12,
				      "1.5 years, No E, 1%" = 13,
				      "1.5 years, No E, 10%" = 14,
				      "1.5 years, E, 1%" = 15,
				      "1.5 years, E, 10%" = 16,
				      "2 years, No E, 1%" = 17,
				      "2 years, No E, 10%" = 18,
				      "2 years, E, 1%" = 19,
				      "2 years, E, 10%" = 20,
				      "3 years, No E, 1%" = 21,
				      "3 years, No E, 10%" = 22,
				      "3 years, E, 1%" = 23,
				      "3 years, E, 10%" = 24)
	subparam.cats[i] <- paste(enh, ", ", rr, "%", sep="")
	subparam.cats.level[i] <- switch(subparam.cats[i],
					 "No E, 1%" = 3,
					 "No E, 10%" = 1,
					 "E, 1%" = 4,
					 "E, 10%" = 2)
}
param.cats <- factor(param.cats)
param.cats.ordered <- reorder(param.cats, param.cats.level)
subparam.cats <- factor(subparam.cats)
subparam.cats.ordered <- reorder(subparam.cats, subparam.cats.level)

## relevel and add to data
chosen.ests.final <- cbind(chosen.ests.final, true.lambdas.verb, param.cats=param.cats.ordered, subparam.cats=subparam.cats.ordered)
chosen.ests.final$true.lambdas.verb <- relevel(chosen.ests.final$true.lambdas.verb, "6 months")
chosen.ests.final$true.lambdas.verb <- relevel(chosen.ests.final$true.lambdas.verb, "1 day")

## calculate means by parameter set
means <- aggregate(chosen.ests.final$lambda, by=list(param.cats.ordered), FUN=function(x) mean(x)/26)
colnames(means)[1] <- "param.cats"
colnames(means)[2] <- "mean.ests"
subparam.cats <- rep("", nrow(means))
subparam.cats.level <- rep(0, nrow(means))
true.lambdas.verb <- rep("", nrow(means))
true.means <- rep(0, nrow(means))
for(i in 1:nrow(means)){
	tmp.params <- strsplit(as.character(means[i,"param.cats"]), ", ")[[1]]
	true.lambdas.verb[i] <- tmp.params[1]
	true.means[i] <- switch(true.lambdas.verb[i],
				"1 day" = 1/365,
				"6 months" = .5,
				"1 year" = 1,
				"1.5 years" = 1.5,
				"2 years" = 2,
				"3 years" = 3)
	subparam.cats[i] <- paste(tmp.params[2:3], collapse=", ")
	subparam.cats.level[i] <- switch(subparam.cats[i],
					 "No E, 1%" = 3,
					 "No E, 10%" = 1,
					 "E, 1%" = 4,
					 "E, 10%" = 2)
}
subparam.cats.ordered <- reorder(subparam.cats, subparam.cats.level)
true.lambdas.verb <- reorder(true.lambdas.verb, true.means)
means <- cbind(means, true.means, subparam.cats.ordered, true.lambdas.verb)
chosen.ests.final$subparam.cats.ordered <- subparam.cats.ordered


## plot histograms
qplot(lambda/26, data=chosen.ests.final, geom="histogram") + facet_grid(subparam.cats.ordered ~ true.lambdas.verb) + theme_bw() + geom_vline(aes(xintercept=mean.ests), means, col="blue") + geom_vline(aes(xintercept=true.means), means, col="red") + xlab("average duration of cross protection")

## summarize model choice
model.choice <- with(chosen.ests.final, table(param.cats.ordered, model))
model.choice.sums <- cbind(E=rowSums(model.choice[, 1:4]), ## sum the E model cols
			   N=rowSums(model.choice[, 5:8]))/500 ## model Nd never chosen

## summarize model choice
model.choice.sums <- data.frame(model.choice.sums, subparam.cats=subparam.cats.ordered, true.means, true.lambdas.verb)

## calculate metrics
pct.avg.bias <- (means$mean.ests-means$true.means)/means$true.means
ses <- aggregate(chosen.ests.final$lambda/26, by=list(param.cats.ordered), FUN=sd)
colnames(ses)[2] <- "ses"
mse <- (means$mean.ests-means$true.means)^2 + ses$ses^2
avg.bias <- means$mean.ests-means$true.means
ci.cov.sum <- aggregate(ci.coverage, by=list(param.cats.ordered), FUN=sum)
colnames(ci.cov.sum)[2] <- "ci.coverage"

## bind all results together
param.cat.split <- strsplit(as.character(means$subparam.cats), ", ")
enh <- rep("", nrow(means))
rr <- rep("", nrow(means))
for(i in 1:nrow(means)){
	enh[i] <- param.cat.split[[i]][1]
	rr[i] <- param.cat.split[[i]][2]
}
enh <- as.factor(enh)
rr <- as.factor(rr)
rr <- relevel(rr, "10%")
all.results <- cbind(means, mse, ses, pct.avg.bias, mod.N.chosen=model.choice.sums$N, avg.bias, ci.coverage=ci.cov.sum$ci.coverage/500, enh, rr)


## plot metrics
qplot(x=true.means, y=mse, data=all.results, geom="line", group=subparam.cats, colour=enh, lty=rr) + theme_bw() + xlab("average duration of cross-protection") + ylab("MSE")

qplot(x=true.means, y=mod.N.chosen, data=all.results, geom="line", group=subparam.cats, colour=enh, lty=rr) + theme_bw() + xlab("average duration of cross-protection") + ylab("fraction of datasets choosing a model with no CP")+ ylim(0, 1)

qplot(x=true.means, y=ses, data=all.results, geom="line", group=subparam.cats, colour=enh, lty=rr) + theme_bw() + xlab("average duration of cross-protection") + ylab("standard error of estimates")

qplot(x=true.means, y=pct.avg.bias, data=all.results, geom="line", group=subparam.cats, colour=enh, lty=rr) + theme_bw() + xlab("average duration of cross-protection") + ylab("% average bias")

qplot(x=true.means, y=avg.bias, data=all.results, geom="line", group=subparam.cats, colour=enh, lty=rr) + theme_bw() + xlab("average duration of cross-protection") + ylab("average bias")

qplot(x=true.means, y=ci.coverage, data=all.results, geom="line", group=subparam.cats, colour=enh, lty=rr) + theme_bw() + xlab("average duration of cross-protection") + ylab("95% CI coverage") + geom_hline(aes(yintercept=.95), col="gray60", lty=2) + ylim(.5, 1)
