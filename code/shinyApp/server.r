## load bangkok case data
load("bkk.dengue.all.cases.rda")

shinyServer(function(input, output) {
        
        output$plot_sero1 <- reactivePlot(function() {
                with(bkk.dengue.all.cases,
                     plot(date, den.cases.str1,
                          type=input$chart.type))
        }) 
        output$plot_sero2 <- reactivePlot(function() {
                with(bkk.dengue.all.cases,
                     plot(date, den.cases.str2,
                          type=input$chart.type))
        }) 
        output$plot_sero3 <- reactivePlot(function() {
                with(bkk.dengue.all.cases,
                     plot(date, den.cases.str3,
                          type=input$chart.type))
        }) 
        output$plot_sero4 <- reactivePlot(function() {
                with(bkk.dengue.all.cases,
                     plot(date, den.cases.str4,
                          type=input$chart.type))
        })
        
})