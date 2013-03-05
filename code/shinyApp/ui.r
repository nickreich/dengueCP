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
                
                # Specification of range within an interval
                sliderInput("range", "Range:",
                            min = 1, max = 1000, value = c(1,500)),
                
                # Provide a custom currency format for value display, with basic animation
                sliderInput("format", "Custom Format:", 
                            min = 0, max = 10000, value = 0, step = 2500,
                            format="$#,##0", locale="us", animate=TRUE)
        ),
                        
        mainPanel(
                  plotOutput(outputId = "plot_cases")
        )
))