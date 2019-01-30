## Load libraries
library(shiny)
library(DT)
library(ggplot2)
library(org.At.tair.db)
library(SuperExactTest)

#Auxiliary functions
intersectSets <- function(tf1,tf2,set.of.genes, alias){
  intersection.data <- list()
  sets <- list(tf1, tf2, set.of.genes)
  #names(sets) <- c("cca1", "lhy", "peakZT0")
  results <- supertest(x = sets, n = 5778)
  results.table <- summary(results)
  p.value <- tail(results.table$P.value, n=1) #Get the last p-value
  enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]
  intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]
  intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]
  
  intersection.genes.agi <- intersection.genes
  intersection.genes.primary.symbol <- alias[intersection.genes]
  names(intersection.genes.primary.symbol) <- NULL
  gene.table <- matrix(nrow=length(intersection.genes), ncol=3)
  colnames(gene.table) <- c("AGI", "Primary Symbol", "Description")
  gene.table[,1] <- intersection.genes.agi
  gene.table[,2] <- intersection.genes.primary.symbol
  #  gene.table[,3] <- description
  
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment

  intersection.data[[3]] <- intersection.genes
  names(intersection.data) <- c("p-value", "enrichment", "genes")

  return(intersection.data)
}

## Function to generate output table
create.output.table <- function(input.gene.df,alias,tfs.names)
{
  output.selected.genes.df <- data.frame(matrix(nrow=nrow(input.gene.df), ncol=6))
  colnames(output.selected.genes.df) <- c("AGI ID", "Gene Name", "Gene Description", "Regulators","Expression Peak Time","Expression Trough Time")
  output.selected.genes.df$`Gene Description` <- input.gene.df$description
  
  for(i in 1:nrow(output.selected.genes.df))
  {
    tair.link <- paste0("https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",input.gene.df[i,1])
    output.selected.genes.df[i,1] <- paste(c("<a href=\"",
                                             tair.link,
                                             "\" target=\"_blank\">",
                                             input.gene.df[i,1], "</a>"),
                                           collapse="")
    output.selected.genes.df[i,2] <- alias[input.gene.df[i,1]]
    output.selected.genes.df[i,4] <- paste(tfs.names[which(input.gene.df[i,tfs.names] == 1)],collapse=", ")
    output.selected.genes.df[i,5] <-paste0("ZT",substr(input.gene.df[i,"peak.zt"],start=5,stop=nchar(input.gene.df[i,"peak.zt"])))
    output.selected.genes.df[i,6] <-paste0("ZT",substr(input.gene.df[i,"trough.zt"],start=7,stop=nchar(input.gene.df[i,"trough.zt"])))
  }
  
  return(output.selected.genes.df)
}

## Function to extract TF target from network representation
extract.targets <- function(tf.name, network.specification)
{
  return(network.specification$names[which(network.specification[,tf.name] == 1)])
}

## Function to extract set of circadian genes
extract.circadian.genes <- function(peak.time, trough.time, network.specification)
{
  if(peak.time == "Any" && trough.time == "Any")
  {
    res.circadian.genes <- network.specification$names
  } else if(peak.time == "Any")
  {
    res.circadian.genes <- network.specification$names[network.specification$trough.zt == paste0("trough",substr(trough.time,start=3,stop=nchar(trough.time)))]
  } else if(trough.time == "Any")
  {
    res.circadian.genes <- network.specification$names[network.specification$peak.zt == paste0("peak",substr(peak.time,start=3,stop=nchar(peak.time)))]
  } else
  {
    res.circadian.genes <- network.specification$names[network.specification$peak.zt == paste0("peak",substr(peak.time,start=3,stop=nchar(peak.time))) &
                                                       network.specification$trough.zt == paste0("trough",substr(trough.time,start=3,stop=nchar(trough.time)))   ]
  }
  
  return(res.circadian.genes)
}


## Load network
#network.data <- read.table(file="data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")
network.data <- read.table(file="data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")

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
node.colors <- selected.colors[network.data$peak.zt]
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

## Transcription factors AGI ids and names
tfs.names <- c("CCA1","LHY", "TOC1", "PRR5", "PRR7", "PRR9", "PHYA","PHYB",
         "CRY2","FHY1","LUX","PIF3","PIF4","PIF5","ELF4","ELF3")

tf.ids <- c("AT2G46830", "AT1G01060", "AT5G61380", "AT5G24470", "AT5G02810", "AT2G46790",
             "AT1G09570", "AT2G18790", "AT1G04400", "AT2G37678", "AT3G46640", "AT1G09530",
             "AT2G43010", "AT3G59060", "AT2G40080", "AT2G25930")

names(tf.ids) <- tfs.names

## Extract gene descriptions
gene.description <- network.data$description
names(gene.description) <- network.data$names

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
  titlePanel(tags$b("ATTRACTOR, an Arabidopsis Thaliana TRanscriptionAl Circadian neTwORk")),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        
        tags$h3(tags$b("Gene Selection:")),
        
        radioButtons(inputId = "gene_selection_mode",
                     label = "Gene Selection Mode", 
                     choices = c("Individual Genes", 
                                 "Gene List", 
                                 "Common TF target genes",
                                 "Topological parameter"),
                     selected = "Individual Genes"),
        
        ## Dynamic panel for selecting single genes
        conditionalPanel(condition = "input.gene_selection_mode == 'Individual Genes'",
          ## Select a few genes
          selectizeInput(inputId = "selected.genes",
                         label = "Gene ID",
                         choices = genes.selectize,
                         selected = "AT1G22770",
                         multiple = TRUE),
          ## Button to trigger selections based on gene ID
          actionButton(inputId = "button_gene_id",label="Select Genes")
        ),
        
        ## Dynamic panel for selecting gene list
        conditionalPanel(condition = "input.gene_selection_mode == 'Gene List'",
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
                         actionButton(inputId = "button_select_gene_list",label = "Select Genes")
        ),

        conditionalPanel(condition = "input.gene_selection_mode == 'Common TF target genes'",        
          checkboxGroupInput(inputId = "selected.tfs",
                             label = "Select Transcription Factors:",
                             choices = list("CCA1","LHY", "TOC1", "PRR5", "PRR7", "PRR9", "PHYA","PHYB",
                                            "CRY2","FHY1","LUX","PIF3","PIF4","PIF5","ELF3","ELF4"),
                             inline = TRUE,width = "100%"),
          checkboxInput(inputId =  "edges",label = "Visualize Edges",value = FALSE),
          actionButton(inputId = "button_tfs",label = "Select Genes")
        ),

        conditionalPanel(condition = "input.gene_selection_mode == 'Topological parameter'",
          sliderInput(inputId = "degree_range", label = h3("Degree Range"), min = 0, 
                      max = 11, value = c(2, 4)),
          actionButton(inputId = "button_degree",label = "Select Genes")
        ),
        
        
          tags$h3(tags$b("Multiset Intersections:")),
          
          selectInput(inputId = "tf1", label="Transcription Factor 1", 

                      choices = tfs.names, selected = NULL,
                      multiple = FALSE, selectize = TRUE),
          selectInput(inputId = "tf2", label="Transcription Factor 2", 
                      choices = tfs.names, selected = NULL,
                      multiple = FALSE, selectize = TRUE),
          tags$b("Cluster of Circadian Genes"),
          selectInput(inputId = "peak", label="Peak", 
                      choices = c("Any",paste("ZT",seq(from=0,to=20,by=4),sep="")), selected = NULL,

                      multiple = FALSE, selectize = TRUE),
          selectInput(inputId = "trough", label="Trough", 
                      choices = c("Any",paste("ZT",seq(from=0,to=20,by=4),sep="")), selected = NULL,
                      multiple = FALSE, selectize = TRUE),
        
        actionButton(inputId = "button_intersect", label = "Test"),
          
          

        width = 3 
      ),
      
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("networkPlot"),
         tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
         tags$br(),tags$br(),tags$br(),
         htmlOutput(outputId = "outputText"),
         DT::dataTableOutput("mytable"),
         
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
    selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
    
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
    selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
    
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
    
    
    ## Output table with gene info
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
    
  })
  
  ## Visualization of selected genes according to their degree
  observeEvent(input$button_degree, {
    print("aquí llego 3")
    
    node.degree <- network.data$indegree
    degree.values <- input$degree_range
    selected.genes.df <- subset(network.data, indegree >= degree.values[1] & indegree <= degree.values[2])
    selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
    
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
    
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)

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
    selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
    
    print(selected.genes.df)
    print(selected.nodes.colors)
    
    network.representation <- ggplot(network.data, aes(x.pos,y.pos)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=node.colors,size=1) +
     # geom_point(data = selected.tfs.df, size=8, fill=selected.tfs.df$color,colour="black",pch=21) +
      geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)

    if(input$edges)
    {
      for(i in 1:length(input$selected.tfs))
      {
        tf.xpos <- subset(network.data, names == tf.ids[input$selected.tfs[i]])[["x.pos"]]
        tf.ypos <- subset(network.data, names == tf.ids[input$selected.tfs[i]])[["y.pos"]]
        network.representation <- network.representation +
          annotate("segment",
                   x=rep(tf.xpos,nrow(selected.genes.df)),
                   y=rep(tf.ypos,nrow(selected.genes.df)),
                   xend=selected.genes.df$x.pos,
                   yend=selected.genes.df$y.pos, 
                   color="grey", arrow=arrow(type="closed",length=unit(0.1, "cm")))
      }
    }

    output$networkPlot <- renderPlot({
      network.representation
    },height = 700)

    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
  })
  
  ##Visualization of text and table with p.value, enrichment and genes.
  
  observeEvent(input$button_intersect, {
    print("Aquí llega Pedro")
    # tf1.filename <- paste0(input$tf1, "_targets_in_network.txt")
    # tf2.filename <- paste0(input$tf2, "_targets_in_network.txt")
    # set.of.genes.filename <- paste0(input$set, ".txt")
    
    ## Extract TF targets
    tf1.targets <- extract.targets(tf.name = input$tf1, network.specification = network.data)
    tf2.targets <- extract.targets(tf.name = input$tf2, network.specification = network.data)
    
    ## Extract circadian set of genes
    circadian.genes.set <- extract.circadian.genes(peak.time = input$peak,trough.time = input$trough,network.specification = network.data)
    
    #tf1 <- read.table(file=paste0("data/intersections/",tf1.filename), header = TRUE, as.is=TRUE)
    #tf2 <- read.table(file=paste0("data/intersections/",tf2.filename), header = TRUE, as.is=TRUE)
    #set.of.genes <- read.table(file=paste0("data/intersections/",set.of.genes.filename), header = TRUE, as.is=TRUE)
    
    #Apply the function intersectSets


    result <- intersectSets(tf1 = tf1.targets, tf2 = tf2.targets, set.of.genes = circadian.genes.set)

    p.value <- result[1][[1]]
    enrichment <- result[2][[1]]
    gene.table <- result[3][[1]]
    
    selected.genes.df <- subset(network.data, names %in% intersect.genes)
    selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
    
    print(selected.genes.df)
    print(selected.nodes.colors)
    
    network.representation <- ggplot(network.data, aes(x.pos,y.pos)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=node.colors,size=1) +
      #geom_point(data = selected.tfs.df, size=8, fill=selected.tfs.df$color,colour="black",pch=21) +
      geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)

    
    output$networkPlot <- renderPlot({
      network.representation
    },height = 700)
    
    ## Visualization of text with p value and enrichment

    


    

    if(p.value < 0.01)
    {
      text.intersection.result <- paste0("<b>The intersection between the targets of ", input$tf1, " and ", input$tf2,
                                         " is significant with a p-value of ", p.value,
                                         " and an enrichment of ", round(x = enrichment,digits = 2),
                                         "<b>") 
    } else
    {
      text.intersection.result <- paste0("<b>The intersection between the targets of ", input$tf1, " and ", input$tf2,
                                         " is NOT significant with a p-value of ", round(x=p.value, digits = 2),
                                         " and an enrichment of ", round(x = enrichment,digits = 2),
                                         "<b> <br> <br>") 
    }

    

    output$outputText <- renderText(expr = text.intersection.result, quoted = FALSE)
    ## Visualization of a table with genes in the intsersections
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

