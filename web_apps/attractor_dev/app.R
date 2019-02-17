## Load libraries
library(shiny)
library(shinycssloaders)
library(ggplot2)
library(org.At.tair.db)
library(SuperExactTest)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(ChIPseeker)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

#Auxiliary functions
intersectSets <- function(tf1,tf2,set.of.genes){
  intersection.data <- list()
  sets <- list(tf1, tf2, set.of.genes)
  #names(sets) <- c("cca1", "lhy", "peakZT0")
  results <- supertest(x = sets, n = 5778)
  results.table <- summary(results)
  p.value <- tail(results.table$P.value, n=1) #Get the last p-value
  enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]
  intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]
  intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment
  intersection.data[[3]] <- intersection.genes
  names(intersection.data) <- c("p-value", "enrichment", "genes")
  return(intersection.data)
}

intersect2sets <- function(set1, set2, alias, gene.descriptions){
  intersection.data <- list()
  sets <- list(set1, set2)
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
  
  
  
  
  intersection.genes.description <- gene.descriptions[intersection.genes]
  names(intersection.genes.description) <- NULL
  
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment
  intersection.data[[3]] <- data.frame(intersection.genes,intersection.genes.primary.symbol,intersection.genes.description,stringsAsFactors = F)
  
  names(intersection.data) <- c("p-value", "enrichment", "gene.table")
  return(intersection.data)
  
}

#intersectBed function 
intersectBed <- function(peaks.set1, peaks.set2)
{
  intersection <- matrix(ncol = 3, nrow=0 )
  current.intersection <- matrix(ncol = 3 )
  for (i in 1:nrow(peaks.set1))
  {
    #Set the current peak values of set1
    current.chr <- as.numeric(peaks.set1[i,1])
    current.start <- peaks.set1[i,2]
    current.end <- peaks.set1[i,3]
    #Checking if there is intersection between the current peak and any peak of set2
    option1 <- nrow(subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.start & peaks.set2[,3]>=current.start))
    option2 <- nrow(subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end))
    
    # print(i)
    
    if(option1+option2 > 0)
    {
      # print("HIT")
      if(option1>0)
      {
        hit.peak2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.start & peaks.set2[,3]>=current.start)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- current.start
        current.intersection[1,3] <- hit.peak2[1,3]
        
      }else
      {
        hit.peak2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- hit.peak2[1,2]
        current.intersection[1,3] <- current.end
      }
      
      intersection <- rbind(intersection, current.intersection)
    }
  }
  return(intersection)
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

#Bed files
bed.files <- c("../peak_visualizer/data/bed_files/PHYA_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/PHYB_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/PRR5_1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/TOC1_1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/CCA1_ZT02_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/CCA1_ZT14_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/LHY_1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/CRY2_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/FHY1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/LUX_ZT10_1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/LUX_ZT12_1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/PIF3_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/PIF4_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/PIF5_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/PRR7_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/PRR9_1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/ELF3_ZT0_1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/ELF3_ZT4_1_peaks.narrowPeak",
               "../peak_visualizer/data/bed_files/ELF4_1_peaks.narrowPeak")

chromosomes.length <- read.table(file="../peak_visualizer/data/bed_files/atha_chr_lengths.txt",as.is=T)[[1]]


bed.names  <-c("PHYA ZT00", "PHYB ZT00" ,"PRR5 ZT10", "TOC1 ZT15","CCA1 ZT02","CCA1 ZT14","LHY ZT02","CRY2 ZT08","FHY1 ZT04","LUX ZT10", "LUX ZT12", "PIF3 ZT08","PIF4 ZT04","PIF5 ZT04","PRR7 ZT12","PRR9 ZT??","ELF3 ZT00", "ELF3 ZT04", "ELF4 ZT10")
names(bed.files) <- bed.names
number.randomisation <- 5 #provisional

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
      
      
      
      
      tags$h3(tags$b("Topological Intersections")),
      selectInput(inputId = "peak_top", label="Peak",
                  choices = c("Any",paste("ZT",seq(from=0,to=20,by=4),sep="")), selected = NULL,
                  multiple = FALSE, selectize = TRUE), 
      selectInput(inputId = "trough_top", label = "Trough", 
                  choices = c("Any",paste("ZT",seq(from=0,to=20,by=4),sep="")), selected = NULL,
                  multiple = FALSE, selectize = TRUE), 
      selectInput(inputId = "topological_parameter", label = "Topological parameter", 
                  choices = c("Degree","Betweeness", "Closeness", "Eccentricity","Transitivity"), selected = NULL,
                  multiple = FALSE, selectize = TRUE),
      selectInput(inputId = "threshold", label = "Parameter threshold",
                  choices = c(0.75,0.90,0.95), selected = NULL, multiple = FALSE, selectize = TRUE),
                  
      
      actionButton(inputId = "top_intersect", label = "Test"),
      
      
      tags$h3(tags$b("Intersections of binding regions:")),
      
      selectInput(inputId = "tf1_bed", label="Transcription Factor 1", 
                  choices = bed.names, selected = NULL,
                  multiple = FALSE, selectize = TRUE),
      selectInput(inputId = "tf2_bed", label="Transcription Factor 2", 
                  choices = bed.names, selected = NULL,
                  multiple = FALSE, selectize = TRUE),
      numericInput(inputId = "number_random", label = "Number of randomisations", 
                   value = 1000, min = 5, max = 1000000, step = NA,
                   width = NULL),
      actionButton(inputId = "button_bed", label = "Test"),
      
      
      width = 3 
      ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("networkPlot"),
      tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
      tags$br(),tags$br(),tags$br(),
      htmlOutput(outputId = "outputText"),
      dataTableOutput(outputId = "outputTable"),
      
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
    intersect.genes <- result[3][[1]]
    
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
    
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
    
  })
  
  ##Visualization and intersection between topological parameters and cluster genes
  observeEvent(input$top_intersect, {
    print("Test top intersection")
    gene.names <- network.data$names
    
    if (input$topological_parameter == "Degree")
    {
      attractor.degree <- network.data$indegree + network.data$outdegree
      degree.threshold <- quantile(attractor.degree, prob=as.numeric(input$threshold))
      top.genes <- gene.names[attractor.degree > degree.threshold]
    } else if (input$topological_parameter == "Transitivity")
    {
      network.data$transitivity[is.na(network.data$transitivity)] <- 0
      attractor.trans <- network.data$transitivity
      trans.threshold <- quantile(attractor.trans, prob=as.numeric(input$threshold))
      top.genes <- gene.names[attractor.trans > trans.threshold]
    } else if (input$topological_parameter == "Closeness")
    {
      attractor.closeness <- network.data$closeness
      closeness.threshold <- quantile(attractor.closeness, prob=as.numeric(input$threshold))
      top.genes <- gene.names[attractor.closeness > closeness.threshold]
    } else if (input$topological_parameter == "Betweeness")
    {
      attractor.bet <- network.data$betweeness
      bet.threshold <- quantile(attractor.bet, prob=as.numeric(input$threshold))
      top.genes <- gene.names[attractor.bet > bet.threshold]
    } else if (input$topological_parameter == "Eccentricity")
    {
      attractor.eccen <- network.data$eccentricity
      eccen.threshold <- quantile(attractor.eccen, prob=as.numeric(input$threshold))
      top.genes <- gene.names[attractor.eccen > eccen.threshold]
    } 
    
    if (input$peak_top == "Any")
    {
      if (input$trough_top == "Any") 
      {
        zt.genes <- network.data$names
      } else
      {
        
        zt.genes <- subset(network.data, trough.zt == paste0("trough", substr(x = input$trough_top, start = 3, stop = nchar(input$trough_top))))$names
      }
      
    } else 
    {
      if (input$trough_top == "Any")
      {
        peak.selection <- paste0("peak", substr(x = input$peak_top, start = 3, stop = nchar(input$peak_top)))
        zt.genes <- subset(network.data, peak.zt == peak.selection)$names
      } else 
      {
        trough.selection <- paste0("trough", substr(x = input$trough_top, start = 3, stop = nchar(input$trough_top)))
        peak.selection <- paste0("peak", substr(x = input$peak_top, start = 3, stop = nchar(input$peak_top)))
        zt.genes <- subset(network.data, trough.zt == trough.selection & peak.zt == peak.selection)$names
      }
      
    }
    
    result <- intersect2sets(set1 = top.genes, set2 = zt.genes, alias = alias, gene.descriptions = gene.description)
    p.value <- result[1][[1]]
    enrichment <- result[2][[1]]
    intersect.genes <- result[3][[1]]$intersection.genes
    
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
    if (length(top.genes) == 0)
    {
      text.intersection.result <- "<b>There is no genes with this restriction<b>"
    } else 
    {
      if(p.value < 0.01)
      {
        text.intersection.result <- paste0("<b>The intersection between the genes with high ", input$topological_parameter,
                                           " and genes that show a ", input$peak_top, " peak and ", input$trough_top, " trough ",
                                           " is significant with a p-value of ", p.value,
                                           " and an enrichment of ", round(x = enrichment,digits = 2),
                                           "<b>") 
        
      } else
      {
        text.intersection.result <- paste0("<b>The intersection between the genes with high ", input$topological_parameter,
                                           " and genes that show a ", input$peak_top, " peak and ", input$trough_top, " trough ",
                                           " is NOT significant with a p-value of ", round(x=p.value, digits = 2),
                                           " and an enrichment of ", round(x = enrichment,digits = 2),
                                           "<b> <br> <br>") 
      }
      
    }
    
    
    output$outputText <- renderText(expr = text.intersection.result, quoted = FALSE)
    
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
    
  })
  
  
  ##Intersection between bed files 
  observeEvent(input$button_bed, {
    print("Test bed intersection")
    
    
    number.randomisation <- input$number_random
    bed1 <- bed.files[input$tf1_bed] #Set the bed to read
    bed2 <- bed.files[input$tf2_bed] #Set the bed to read
    peaks1 <- read.table(file=bed1,header = F, as.is = T) #Read the selected bed
    peaks2 <- read.table(file=bed2,header = F, as.is = T) #Read the selected bed
    
    real.intersection <- intersectBed(peaks.set1 = peaks1, peaks.set2 = peaks2)
    if (nrow(real.intersection) > 0)
    {
      print("Hay intersección")
      random.intersections <- vector(mode = "numeric",length=number.randomisation) #Creating vector
      for(j in 1:number.randomisation)
      {
        print(j)
        random.peaks2 <- matrix(nrow=nrow(peaks2),ncol=3) #Matriz con 3 columnas, una para el cromosoma, otra para el comienzo y otra para el final de la región aleatoria.
        for(k in 1:nrow(peaks2))
        {
          current.chr <- peaks2[k,1][[1]] #Chr de la iésima marca real
          current.start <- peaks2[k,2] #Start de la iésima marca real
          current.end <- peaks2[k,3] #End de la iésima marca real
          current.length <- current.end - current.start #Longitud de la iésima marca real
          
          chr.length <- chromosomes.length[current.chr] #Length del actual cromosoma
          #Ahora genero los mismos datos para regiones aleatorias
          random.start <- floor(runif(n = 1,min = 1,max = chr.length))
          random.end <- random.start + current.length
          
          random.peaks2[k,1] <- current.chr
          random.peaks2[k,2] <- random.start
          random.peaks2[k,3] <- random.end
        }
        
        
        random.intersections[j] <- nrow(intersectBed(peaks.set1 = peaks1, peaks.set2 = random.peaks2 )) 
        
        
      }
      
      p.value <- sum(random.intersections > nrow(real.intersection)) / number.randomisation
      if( p.value == 0)
      {
        p.value <- 1/number.randomisation
      }
      
      colnames(real.intersection) <- c("chromosome", "start", "end")
      granges.intersection <- makeGRangesFromDataFrame(real.intersection,
                                                       keep.extra.columns=FALSE,
                                                       ignore.strand=FALSE,
                                                       seqinfo=NULL,
                                                       seqnames.field="chromosome",
                                                       start.field="start",
                                                       end.field="end",
                                                       starts.in.df.are.0based=FALSE)
      
      
      
      peakAnno <- annotatePeak(granges.intersection, tssRegion=c(-2000, 2000),
                               TxDb=txdb, annoDb="org.At.tair.db")
      
      annot.peaks <- as.data.frame(peakAnno)
      target.genes <- subset(annot.peaks, distanceToTSS >= 2000 | distanceToTSS <= -2000)$geneId
      # target.genes <- paste(target.genes, collapse = ",")
      

      text.bed <- paste0("The estimated p-value is ", p.value)
    } else 
    {
      text.bed <- "No intersection"
    }
    
    
    selected.genes.df <- subset(network.data, names %in% target.genes)
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
      # req(input$button_bed)
      # Sys.sleep(time = 5)
      network.representation
    },height = 700)
    
    
    output$outputText <- renderText(expr = text.bed, quoted = FALSE)
    
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
    
  })
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

