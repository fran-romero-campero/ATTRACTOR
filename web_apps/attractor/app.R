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
  titlePanel(tags$b("ATTRACTOR, an Arabidopsis Thaliana TRanscriptionAl Circadian neTwORk")),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        
        tags$h3(tags$b("Gene Selection:")),
        
        selectizeInput(inputId = "selected.genes",
                       label = "Gene ID",
                       choices = genes,
                       selected = "AT1G22770",
                       multiple = TRUE),
        ## Button to trigger selections based on gene ID
        actionButton(inputId = "button_gene_id",label="Select Genes"),
        width = 3 
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("networkPlot"),
         
         width = 9
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
  
  selected_gene_id <- eventReactive(input$button_gene_id, {
    print("aquí llego 0")
    selected.genes.df <- subset(network.data, names %in% input$selected.genes)
    print(selected.genes.df)
#    return(selected.genes.df)
  })
  
  
  ## Visualization of selected genes by ID
  observeEvent(input$button_gene_id, {
    print("aquí llego 1")
    
    selected.genes.df <- subset(network.data, names %in% input$selected.genes)
    selected.nodes.colors <- selected.colors[selected.genes.df$cluster.classification]
    
    print(selected.genes.df)
    print(selected.nodes.colors)
    
    output$networkPlot <- renderPlot({
      ggplot(network.data, aes(x.pos,y.pos)) + 
        theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank()) + 
        geom_point(color=node.colors,size=1) +
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=3, fill=selected.nodes.colors,colour="black",pch=21)
      })
  })

}

# Run the application 
shinyApp(ui = ui, server = server)

