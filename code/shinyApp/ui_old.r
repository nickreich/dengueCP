shinyUI(pageWithSidebar(
        headerPanel("Stocks"),
        
        sidebarPanel(
                wellPanel(
                        p(strong("Stocks")),
                        checkboxInput(inputId = "stock_aapl", label = "Apple (AAPL)", value = TRUE),
                        checkboxInput(inputId = "stock_msft", label = "Microsoft (MSFT)", value = FALSE),
                        checkboxInput(inputId = "stock_ibm", label = "IBM (IBM)", value = FALSE),
                        checkboxInput(inputId = "stock_goog", label = "Google (GOOG)", value = TRUE),
                        checkboxInput(inputId = "stock_yhoo", label = "Yahoo (YHOO)", value = FALSE)
                ),
                
                selectInput(inputId = "chart_type",
                            label = "Chart type",
                            choices = c("Candlestick" = "candlesticks",
                                        "Matchstick" = "matchsticks",
                                        "Bar" = "bars",
                                        "Line" = "line")
                ),
                
                wellPanel(
                        p(strong("Date range (back from present)")),
                        sliderInput(inputId = "time_num",
                                    label = "Time number",
                                    min = 1, max = 24, step = 1, value = 6),
                        
                        selectInput(inputId = "time_unit",
                                    label = "Time unit",
                                    choices = c("Days" = "days",
                                                "Weeks" = "weeks",
                                                "Months" = "months",
                                                "Years" = "years"),
                                    selected = "Months")
                ),
                
                checkboxInput(inputId = "log_y", label = "log y axis", value = FALSE)
        ),
        
        mainPanel(
                conditionalPanel(condition = "input.stock_aapl",
                                 plotOutput(outputId = "plot_aapl")),
                
                conditionalPanel(condition = "input.stock_msft",
                                 plotOutput(outputId = "plot_msft")),
                
                conditionalPanel(condition = "input.stock_ibm",
                                 plotOutput(outputId = "plot_ibm")),
                
                conditionalPanel(condition = "input.stock_goog",
                                 plotOutput(outputId = "plot_goog")),
                
                conditionalPanel(condition = "input.stock_yhoo",
                                 plotOutput(outputId = "plot_yhoo"))
        )
))