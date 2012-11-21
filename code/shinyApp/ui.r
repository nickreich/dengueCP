shinyUI(pageWithSidebar(
        headerPanel("Dengue incidence in Bangkok"),
        
        sidebarPanel(
                ## select which data to plot
                wellPanel(
                        p(strong("Serotypes")),
                        checkboxInput(inputId = "sero1", label = "Serotype 1", value = TRUE),
                        checkboxInput(inputId = "sero2", label = "Serotype 2", value = FALSE),
                        checkboxInput(inputId = "sero3", label = "Serotype 3", value = FALSE),
                        checkboxInput(inputId = "sero4", label = "Serotype 4", value = TRUE)
                ),
                
                ## select which type of graph
                selectInput(inputId = "chart.type",
                            label = "Type of graph",
                            choices = c("Matchstick" = "h",
                                        "Line" = "l",
                                        "Points" = "p")
                )),
        
        mainPanel(                
                 conditionalPanel(condition = "input.sero1",
                                  plotOutput(outputId = "plot_sero1")),
                 conditionalPanel(condition = "input.sero2",
                                  plotOutput(outputId = "plot_sero2")),
                 conditionalPanel(condition = "input.sero3",
                                  plotOutput(outputId = "plot_sero3")),
                 conditionalPanel(condition = "input.sero4",
                                  plotOutput(outputId = "plot_sero4"))
        )
))