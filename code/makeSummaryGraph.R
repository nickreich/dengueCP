require(RCurl)

## download CI data

myCsv <- getURL("https://docs.google.com/spreadsheet/pub?hl=en_US&hl=en_US&key=0AqbVcVVyI5k8dHVic0pxcjNIZnNkdWtkX3hHbEFuUlE&single=true&gid=1&output=csv")
d <- read.csv(textConnection(myCsv))
d <- d[c(3:1, 6:4),]

pdf("../../figures/publicationReady/cpSummary.pdf")
plot(d[,"cp.duration"], 1:6, pch="|", xlim=c(0, max(d[,"cp.duration.high"])),
     xlab="duration of cross-protection, in years",
     ylab="", yaxt="n", bty="n", ylim=c(1, 6.5))
text(d[,"cp.duration"], 1:6+.2, labels=d[,"cp.duration"], cex=.8)
for(i in 1:nrow(d)){
	lines(x = c(d[i,"cp.duration.low"], d[i,"cp.duration.high"]),
	      y = rep(i, 2),
	      lwd=3)
}
abline(h=3.5, lty=2, col="gray40")
axis(2, at=1:6, labels=d$model.name, las=2)
dev.off()
