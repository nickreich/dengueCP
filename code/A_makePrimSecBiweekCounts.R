## read in raw data and make files for analysis
## September 2010
## Nicholas Reich
library(ggplot2)

## data coming from one file
##  ../data/D73-10.csv

## constants
ANALYSIS.YEARS <- 1973:2010
FILENAME1 <- "../../data/D73-10.csv"
TODAY <- Sys.Date()


################################
## data read-in and cleaning  ##
################################


## read in data and add dates
dat <- read.csv(FILENAME1, header=TRUE)

## make two corrections in data years, based on visual inspection
bad.year.1 <- which(dat$S1YR==58)
dat[bad.year.1, "S1YR"] <- dat[bad.year.1-1, "S1YR"]
bad.year.2 <- which(dat$S1YR==9)
dat[bad.year.2, "S1YR"] <- dat[bad.year.2-1, "S1YR"]

## make all years into 4 digit years (19xx years have no 19)
pre.2k.idx <- which(dat[,"S1YR"]<100)
dat[pre.2k.idx, "S1YR"] <- dat[pre.2k.idx,"S1YR"]+1900

## change one NA in dat$ISOTYP values to a "9"
na.idx <- which(is.na(dat$ISOTYP))
dat[na.idx, "ISOTYP"] <- 9

## make indicators for serology and type of infection
prim.den.1 <- dat$ISOTYP==1 & dat$SEROLO==1
sec.den.1 <- dat$ISOTYP==1 & dat$SEROLO==2
ind.den.1 <- dat$ISOTYP==1 & dat$SEROLO==3
prim.den.2 <- dat$ISOTYP==2 & dat$SEROLO==1
sec.den.2 <- dat$ISOTYP==2 & dat$SEROLO==2
ind.den.2 <- dat$ISOTYP==2 & dat$SEROLO==3
prim.den.3 <- dat$ISOTYP==3 & dat$SEROLO==1
sec.den.3 <- dat$ISOTYP==3 & dat$SEROLO==2
ind.den.3 <- dat$ISOTYP==3 & dat$SEROLO==3
prim.den.4 <- dat$ISOTYP==4 & dat$SEROLO==1
sec.den.4 <- dat$ISOTYP==4 & dat$SEROLO==2
ind.den.4 <- dat$ISOTYP==4 & dat$SEROLO==3

sero.dat <- cbind(prim.den.1, sec.den.1, ind.den.1,
		  prim.den.2, sec.den.2, ind.den.2,
		  prim.den.3, sec.den.3, ind.den.3,
		  prim.den.4, sec.den.4, ind.den.4)

## aggregate counts by time
agg.dat <- aggregate(sero.dat, by=list(month=dat$S1MO, year=dat$S1YR), FUN=sum)
## NOTE: 129 records have dat$S1MO==0. (as of 10/16/2011)
## They disappear in date loop below, b/c we don't have a month to assign them to.

###
## clean up dates, add zero counts for months with no cases
###

## make date sequence for months with non-zero counts
date.seq <- paste(agg.dat[,"year"], "-",
		  agg.dat[,"month"], "-",
		  "01",
		  sep="")
dates <- as.Date(date.seq)

## make date sequence for all months in analysis years
full.date.seq <- as.Date(paste(rep(ANALYSIS.YEARS, each=12), "-",
			       rep(1:12, times=length(ANALYSIS.YEARS)), "-",
			       "01",
			       sep=""))

## fix empty complete dataset
complete.dat <- data.frame(full.date.seq,
			   prim.den.1=0, sec.den.1=0, ind.den.1=0,
			   prim.den.2=0, sec.den.2=0, ind.den.2=0,
			   prim.den.3=0, sec.den.3=0, ind.den.3=0,
			   prim.den.4=0, sec.den.4=0, ind.den.4=0)
colnames(complete.dat)[1] <- "dates"

## loop through and add data in non-zero years
for(i in 1:length(full.date.seq)){
	if( full.date.seq[i] %in% dates ) {
		idx <- match(full.date.seq[i], dates)
		complete.dat[i,2:13] <- agg.dat[idx, 3:14]
	}
}


#####################################################################
## make combined (all types of infections lumped together) dataset ##
#####################################################################

dengue.1.all <- complete.dat[,"prim.den.1"]+complete.dat[,"sec.den.1"]+complete.dat[,"ind.den.1"]
dengue.1.all.cum <- cumsum(dengue.1.all)
dengue.1.prim.cum <- cumsum(complete.dat[,"prim.den.1"])
dengue.1.sec.cum <- cumsum(complete.dat[,"sec.den.1"])

dengue.2.all <- complete.dat[,"prim.den.2"]+complete.dat[,"sec.den.2"]+complete.dat[,"ind.den.2"]
dengue.2.all.cum <- cumsum(dengue.2.all)
dengue.2.prim.cum <- cumsum(complete.dat[,"prim.den.2"])
dengue.2.sec.cum <- cumsum(complete.dat[,"sec.den.2"])

dengue.3.all <- complete.dat[,"prim.den.3"]+complete.dat[,"sec.den.3"]+complete.dat[,"ind.den.3"]
dengue.3.all.cum <- cumsum(dengue.3.all)
dengue.3.prim.cum <- cumsum(complete.dat[,"prim.den.3"])
dengue.3.sec.cum <- cumsum(complete.dat[,"sec.den.3"])

dengue.4.all <- complete.dat[,"prim.den.4"]+complete.dat[,"sec.den.4"]+complete.dat[,"ind.den.4"]
dengue.4.all.cum <- cumsum(dengue.4.all)
dengue.4.prim.cum <- cumsum(complete.dat[,"prim.den.4"])
dengue.4.sec.cum <- cumsum(complete.dat[,"sec.den.4"])

bkk.cases <- data.frame(complete.dat,
			dengue.1.all, dengue.1.all.cum, dengue.1.prim.cum, dengue.1.sec.cum,
			dengue.2.all, dengue.2.all.cum, dengue.2.prim.cum, dengue.2.sec.cum,
			dengue.3.all, dengue.3.all.cum, dengue.3.prim.cum, dengue.3.sec.cum,
			dengue.4.all, dengue.4.all.cum, dengue.4.prim.cum, dengue.4.sec.cum)


#################################
## make all infections dataset ##
#################################


## make bi-week data
## 1.  fit spline to cumulative monthly case data to provide a tight fit.
## 2.  interpolate cumulative values for casecounts every 14 days.
## 2a.  use one 15 day bi-week (randomly chosen from the 26) in non-leap years
## 2b.  use two 15 day bi-weeks (randomly chosen from the 26) in leap years
## 3.  back out incident counts from cumulative ones

## import births
births.raw <- read.csv("../../../ThaiPopData/ThaiBirthData_1973_2010.csv",
		       skip=1, header=FALSE, as.is=TRUE)

## make interp.dates, the dates on which to find bi-weeks
years <- ANALYSIS.YEARS
start.date <- as.Date(paste(years[1], "-01-01", sep=""))
interp.dates <- rep(start.date, 26*length(ANALYSIS.YEARS))
births <- rep(0, 26*length(ANALYSIS.YEARS))
curr.idx <- 1
for(year in years){
	## determine leap year and times for bi-week startings
	if((year%%4 == 0) & (year%%100 != 0) | (year%%400 == 0)){
		## IS A leap year
		interp.seq <- rep(14, 25)
		interp.seq[sample(1:25, 2)] <- 15
		tot.days <- 366
	} else{
		## IS NOT A leap year
		interp.seq <- rep(14, 25)
		interp.seq[sample(1:25, 1)] <- 15
		tot.days <- 365
	}
	cum.interp.seq <- cumsum(interp.seq)
	interp.dates[curr.idx] <- as.Date(paste(year,"-01-01", sep=""))
	interp.dates[(curr.idx+1):(curr.idx+25)] <-
		interp.dates[curr.idx]+cum.interp.seq

	## calculate births, assuming const. daily rate
	birth.col.idx <- which(births.raw[1,]==year)
	birth.row.idx <- which(births.raw[,1]==" Bangkok ")
	births.for.year <- as.numeric(births.raw[birth.row.idx,birth.col.idx])
	births[curr.idx:(curr.idx+24)] <- round(births.for.year*interp.seq/tot.days)
	births[curr.idx+25] <- births.for.year-sum(births[curr.idx:(curr.idx+24)])

	curr.idx <- curr.idx+26
}


## fit spline to cumulative monthly case data for tight fit, by strain
new.interp.dates <- c(start.date-15, interp.dates)
fit.str1 <- smooth.spline(bkk.cases[,"dates"], bkk.cases[,"dengue.1.all.cum"],
			  all.knots = TRUE)
interp.cases.str1 <- predict(fit.str1, as.numeric(new.interp.dates))
biweek.case.counts.str1 <- diff(interp.cases.str1$y)
biweek.case.counts.str1[which(biweek.case.counts.str1<0)] <- 0

fit.str2 <- smooth.spline(bkk.cases[,"dates"], bkk.cases[,"dengue.2.all.cum"],
			  all.knots = TRUE)
interp.cases.str2 <- predict(fit.str2, as.numeric(new.interp.dates))
biweek.case.counts.str2 <- diff(interp.cases.str2$y)
biweek.case.counts.str2[which(biweek.case.counts.str2<0)] <- 0

fit.str3 <- smooth.spline(bkk.cases[,"dates"], bkk.cases[,"dengue.3.all.cum"],
			  all.knots = TRUE)
interp.cases.str3 <- predict(fit.str3, as.numeric(new.interp.dates))
biweek.case.counts.str3 <- diff(interp.cases.str3$y)
biweek.case.counts.str3[which(biweek.case.counts.str3<0)] <- 0

fit.str4 <- smooth.spline(bkk.cases[,"dates"], bkk.cases[,"dengue.4.all.cum"],
			  all.knots = TRUE)
interp.cases.str4 <- predict(fit.str4, as.numeric(new.interp.dates))
biweek.case.counts.str4 <- diff(interp.cases.str4$y)
biweek.case.counts.str4[which(biweek.case.counts.str4<0)] <- 0


## finalize and save data
biweeks <- rep(1:26, length(years))
bkk.dengue.all.cases <- data.frame(biweek=biweeks,
				   date=interp.dates,
				   den.cases.str1=round(biweek.case.counts.str1),
				   den.cases.str2=round(biweek.case.counts.str2),
				   den.cases.str3=round(biweek.case.counts.str3),
				   den.cases.str4=round(biweek.case.counts.str4),
				   births=births)
save(bkk.dengue.all.cases, file="../../data/bkk.dengue.all.cases.rda")



## plot data
ddat <- melt(bkk.dengue.all.cases[,c("date",
				 "den.cases.str1",
				 "den.cases.str2",
				 "den.cases.str3",
				 "den.cases.str4")],
	     id="date")

qplot(date, value, data=ddat, geom="line", group=variable,
      main="Bi-weekly dengue incidence over time by strain (at one hospital in Bangkok)",
      ylab="case counts", xlab="time")+
	facet_grid(variable~.)

## for bar-plot
plot.den.ts <- function(dat, col) {
	v.seq <- seq(as.Date("1973-01-01"), as.Date("2011-01-01"), by="1 years")
	label.seq <- seq(1975, 2010, by=5)
	label.loc.seq <- as.numeric(as.Date(paste(label.seq, "-06-01", sep="")))+20
	plot(dat[,"dates"], dat[,col], type="h", ylim=rge, xaxt="n", ylab="")
	abline(v=v.seq, lty=2, col="gray")
	axis(1, at=v.seq, labels=FALSE)
	axis(1, at=label.loc.seq, tick=FALSE, labels=label.seq)
}

display.mat <- matrix(c(1:4), nrow=4, ncol=1)
par(mar=c(2,4,1,2))
d.tmp <- bkk.cases[,c("dates", "dengue.1.all", "dengue.2.all", "dengue.3.all", "dengue.4.all")]
rge <- range(d.tmp[,2:4])
layout(display.mat)
plot.den.ts(d.tmp, "dengue.1.all")
plot.den.ts(d.tmp, "dengue.2.all")
plot.den.ts(d.tmp, "dengue.3.all")
plot.den.ts(d.tmp, "dengue.4.all")


## bar-plot for monthly data
ddat.monthly <- melt(bkk.cases[,c("dates",
				  "dengue.1.all",
				  "dengue.2.all",
				  "dengue.3.all",
				  "dengue.4.all")],
		     id="dates")

ggplot(ddat.monthly, aes(dates, ymin=-0.5, ymax=value)) + geom_linerange() +
	facet_grid(variable~.) +
	opts(title="Bi-weekly dengue incidence over time by strain (at Queen Sirikit hospital in Bangkok)") +
	scale_y_continuous("case counts")


####################################
## make primary infection dataset ##
####################################

## fit spline to cumulative monthly case data for tight fit, by strain
new.interp.dates <- c(start.date-15, interp.dates)
fit.str1 <- smooth.spline(bkk.cases[,"dates"], bkk.cases[,"dengue.1.prim.cum"],
			  all.knots = TRUE)
interp.cases.str1 <- predict(fit.str1, as.numeric(new.interp.dates))
biweek.case.counts.str1 <- diff(interp.cases.str1$y)
biweek.case.counts.str1[which(biweek.case.counts.str1<0)] <- 0

fit.str2 <- smooth.spline(bkk.cases[,"dates"], bkk.cases[,"dengue.2.prim.cum"],
			  all.knots = TRUE)
interp.cases.str2 <- predict(fit.str2, as.numeric(new.interp.dates))
biweek.case.counts.str2 <- diff(interp.cases.str2$y)
biweek.case.counts.str2[which(biweek.case.counts.str2<0)] <- 0

fit.str3 <- smooth.spline(bkk.cases[,"dates"], bkk.cases[,"dengue.3.prim.cum"],
			  all.knots = TRUE)
interp.cases.str3 <- predict(fit.str3, as.numeric(new.interp.dates))
biweek.case.counts.str3 <- diff(interp.cases.str3$y)
biweek.case.counts.str3[which(biweek.case.counts.str3<0)] <- 0

fit.str4 <- smooth.spline(bkk.cases[,"dates"], bkk.cases[,"dengue.4.prim.cum"],
			  all.knots = TRUE)
interp.cases.str4 <- predict(fit.str4, as.numeric(new.interp.dates))
biweek.case.counts.str4 <- diff(interp.cases.str4$y)
biweek.case.counts.str4[which(biweek.case.counts.str4<0)] <- 0


## finalize and save data
biweeks <- rep(1:26, length(years))
bkk.dengue.prim.cases <- data.frame(biweek=biweeks,
				   date=interp.dates,
				   den.cases.str1=round(biweek.case.counts.str1),
				   den.cases.str2=round(biweek.case.counts.str2),
				   den.cases.str3=round(biweek.case.counts.str3),
				   den.cases.str4=round(biweek.case.counts.str4),
				   births=births)
save(bkk.dengue.prim.cases, file="../../data/bkk.dengue.prim.cases.rda")


## plot primary incidence data
ddat <- melt(bkk.dengue.prim.cases[,c("date",
				      "den.cases.str1",
				      "den.cases.str2",
				      "den.cases.str3",
				      "den.cases.str4")],
	     id="date")
qplot(date, value, data=ddat, geom="line", group=variable,
      main="Bi-weekly dengue primary infection incidence over time by strain (at one hospital in Bangkok)",
      ylab="case counts", xlab="time")+
	facet_grid(variable~.)
