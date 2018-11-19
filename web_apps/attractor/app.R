## Load libraries
library(shiny)
library(ggplot2)
library(org.At.tair.db)


## Load network
network.data <- read.table(file="data/attractor_network_topological_parameters.tsv",header = TRUE,as.is=TRUE)

## Tranforming coordinates for a better visualization
x.coord <- as.numeric(network.data$y.pos)
y.coord <- as.numeric(network.data$x.pos)

network.data$x.pos <- x.coord
network.data$y.pos <- y.coord

pos.data <- t(matrix(data=c(x.coord,y.coord),ncol=2))

rotation.angle <- -pi/2
rotation.matrix <- matrix(data = c(cos(rotation.angle),sin(rotation.angle),-sin(rotation.angle),cos(rotation.angle)),nrow = 2,ncol = 2)
rotated.pos <- t(rotation.matrix %*% pos.data)

network.data$x.pos <- rotated.pos[,1]
network.data$y.pos <- rotated.pos[,2]
 
selected.colors <- c("blue4","blue","deepskyblue","gold","firebrick","gray47")
peak.times <- c("peak20","peak0","peak4","peak8","peak12","peak16")
names(selected.colors) <- peak.times
node.colors <- selected.colors[network.data$cluster.classification]
names(node.colors) <- NULL

## Extract gene ids
genes <- sort(network.data$name)

## Load all and circadian genes
my.key <- keys(org.At.tair.db, keytype="ENTREZID")
my.col <- c("SYMBOL", "TAIR")
columns(org.At.tair.db)
alias2symbol.table <- select(org.At.tair.db, keys=my.key, columns=my.col, keytype="ENTREZID")
alias2symbol.table <- subset(alias2symbol.table, genes %in% TAIR)
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 
genes.selectize <- paste(names(alias), alias, sep=" - ")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
  titlePanel(tags$b("ATTRACTOR, an Arabidopsis Thaliana TRanscriptionAl Circadian neTwORk")),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        
        tags$h3(tags$b("Gene Selection:")),
        ## Select a few genes
        selectizeInput(inputId = "selected.genes",
                       label = "Gene ID",
                       choices = genes.selectize,
                       selected = "AT1G22770",
                       multiple = TRUE),
        ## Button to trigger selections based on gene ID
        actionButton(inputId = "button_gene_id",label="Select Genes"),
        
        textAreaInput(inputId = "gene.list", label= "Set of genes", width="90%", 
                      height = "200px",placeholder = "Insert set of genes",
                      value= "AT2G23290
AT2G40900
AT2G40890
AT2G47450
AT2G45190
AT2G45400
AT2G33230
AT5G12440
AT4G17245
AT4G16780"),
        actionButton(inputId = "button_select_gene_list",label = "Select Genes"),
        
        
        sliderInput(inputId = "degree_range", label = h3("Degree Range"), min = 0, 
                    max = 11, value = c(2, 4)),
        actionButton(inputId = "button_degree",label = "Select Genes"),
        
        checkboxGroupInput(inputId = "selected.tfs",
                           label = "Select Transcription Factors:",
                           choices = list("CCA1","LHY", "TOC1", "PRR5", "PRR7", "PRR9", "PHYA","PHYB",
                                          "CRY2","FHY1","LUX","PIF3","PIF4","PIF5","ELF3","ELF4"),
                           inline = TRUE,width = "100%"),
        actionButton(inputId = "button_tfs",label = "Select Genes"),
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
  },height = 700)
  
  selected_gene_id <- eventReactive(input$button_gene_id, {
    print("aquí llego 0")
    selected.genes.agi <- as.vector(unlist(as.data.frame(strsplit(input$selected.genes," - "))[1,]))
    selected.genes.df <- subset(network.data, names %in% selected.genes.agi)
    print(selected.genes.df)
#    return(selected.genes.df)
  })
  
  
  ## Visualization of selected genes by ID
  observeEvent(input$button_gene_id, {
    print("aquí llego 1")
    selected.genes.agi <- as.vector(unlist(as.data.frame(strsplit(input$selected.genes," - "))[1,]))
    selected.genes.df <- subset(network.data, names %in% selected.genes.agi)
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
      },height = 700)
  })


  ## Visualization of selected gene list
  observeEvent(input$button_select_gene_list, {
    print("aquí llego 2")
    selected.gene.list <- as.vector(unlist(
      strsplit(input$gene.list, split="\n",
               fixed = TRUE)[1]))
    print(selected.gene.list)
    selected.genes.df <- subset(network.data, names %in% selected.gene.list)
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
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)
  })
  
  ## Visualization of selected genes according to their degree
  observeEvent(input$button_degree, {
    print("aquí llego 3")
    
    node.degree <- network.data$indegree
    degree.values <- input$degree_range
    selected.genes.df <- subset(network.data, indegree >= degree.values[1] & indegree <= degree.values[2])
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
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)
  })

  ## Visualization of selected genes according to their degree
  observeEvent(input$button_tfs, {
    print("aquí llego 4")

    if(length(input$selected.tfs) == 1)
    {
      gene.selection <- network.data[,input$selected.tfs] == 1
      sum(gene.selection)
    } else if(length(input$selected.tfs) > 1)
    {
      gene.selection <- rowSums(network.data[,input$selected.tfs]) == length(input$selected.tfs)
      
    }
    selected.genes.df <- network.data[gene.selection,]
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
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    },height = 700)
  })
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

