## load bangkok case data
load("bkk.cases.long.rda")
library(ggplot2)
library(reshape2)

shinyServer(function(input, output) {
        
        time.min  <- as.Date("1979-01-01")
        time.max  <- as.Date("2011-01-01")
        stys <- 1:3
        sty.names <- paste("den.cases.str", stys, sep="")
        ddat.to.plot <- subset(ddat, variable %in% sty.names)
        
        output$plot_cases <- reactivePlot(function() {
                p <- ggplot(ddat.to.plot, aes(date, ymin=-0.5, ymax=value)) + 
                        geom_linerange() +
                        ggtitle("Bi-weekly dengue incidence over time by strain (at one hospital in Bangkok)") +
                        scale_y_continuous("case counts") + 
                        scale_x_date(limits=c(time.min, time.max)) +
                        facet_grid(variable~.)        
                print(p)
                })        
})