## Load libraries
library(shiny)
library(ggplot2)

## Load network
network.data <- read.table(file="data/attractor_network_topological_parameters.tsv",header = TRUE,as.is=TRUE)
head(network.data)

selected.colors <- c("blue4","blue","deepskyblue","gold","firebrick","gray47")
peak.times <- c("peak20","peak0","peak4","peak8","peak12","peak16")
names(selected.colors) <- peak.times
node.colors <- selected.colors[network.data$cluster.classification]
names(node.colors) <- NULL

## Extract gene ids
genes <- sort(network.data$name)


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Old Faithful Geyser Data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("bins",
                     "Number of bins:",
                     min = 1,
                     max = 50,
                     value = 30)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("networkPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ## Initial/default visualization of ATTRACTOR
  output$networkPlot <- renderPlot({
    ggplot(network.data, aes(x.pos,y.pos)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=node.colors,size=1)
  })
  
   
  
  
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

