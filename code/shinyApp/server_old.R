if (!require(quantmod)) {
        stop("This app requires the quantmod package. To install it, run 'install.packages(\"quantmod\")'.\n")
}

# Download data for a stock, if needed
require_symbol <- function(symbol) {
        if (!exists(symbol))
                getSymbols(symbol)
}


shinyServer(function(input, output) {
        
        # Make a chart for a symbol, with the settings from the inputs
        make_chart <- function(symbol) {
                require_symbol(symbol)
                
                chartSeries(get(symbol),
                            name = symbol,
                            type = input$chart_type,
                            subset = paste("last", input$time_num, input$time_unit),
                            log.scale = input$log_y,
                            theme = "white")
        }
        
        output$plot_aapl <- reactivePlot(function() { make_chart("AAPL") })
        output$plot_msft <- reactivePlot(function() { make_chart("MSFT") })
        output$plot_ibm <- reactivePlot(function() { make_chart("IBM") })
        output$plot_goog <- reactivePlot(function() { make_chart("GOOG") })
        output$plot_yhoo <- reactivePlot(function() { make_chart("YHOO") })
})