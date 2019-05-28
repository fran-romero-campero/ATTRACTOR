library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Old Faithful Geyser Data"),

      
      # Show a plot of the generated distribution
      mainPanel(
        fluidRow(
            column(width = 2,
                   img(src='dribbbble.jpg', align = "center", width=200),
                   selectInput(inputId = "selected_gene",
                                                      label = "Select a gene", 
                                                      choices = c("gen1","gen2","gen3")),
                   radioButtons(inputId = "selected_analysis", label = "Choose your analysis",
                                choices = c("Peak Visualizer" = "peakvis", "Network Visualizer" = "networkvis"))),
            column(width = 8, 
                   tabsetPanel(type = "tabs",
                               tabPanel(title = "Network Visualizer", plotOutput(outputId = "plot")),
                               tabPanel(title = "Peak Visualizer", plotOutput(outputId = "plot3")))
                   ),
            column(width = 2, plotOutput(outputId = "plot2"))
        ),
        plotOutput("distPlot")
      )
   )


# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

