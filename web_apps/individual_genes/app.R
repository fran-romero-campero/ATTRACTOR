# Authors: Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
#          Francisco J. Romero-Campero
# 
# Contact: Francisco J. Romero-Campero - fran@us.es 
# Date: June 2019

# Load neccesary libraries
library(shiny)
library(ChIPpeakAnno)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(Biostrings)
library(seqinr)
library(org.At.tair.db)
library(igraph)
library(ggplot2)

##Load the network data
network.data <- read.table(file="data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")
rownames(network.data) <- network.data$names

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

## Transcription factors AGI ids and names
tfs.names <- c("CCA1","LHY", "TOC1", "PRR5", "PRR7", "PRR9", "PHYA","PHYB",
               "CRY2","FHY1","LUX","PIF3","PIF4","PIF5","ELF4","ELF3")

tf.ids <- c("AT2G46830", "AT1G01060", "AT5G61380", "AT5G24470", "AT5G02810", "AT2G46790",
            "AT1G09570", "AT2G18790", "AT1G04400", "AT2G37678", "AT3G46640", "AT1G09530",
            "AT2G43010", "AT3G59060", "AT2G40080", "AT2G25930")

names(tf.ids) <- tfs.names

## Extract gene ids
genes <- sort(network.data$name)

## Load and extract Arabidopsis thaliana annotation regarding genes, exons and cds 
txdb <- TxDb.Athaliana.BioMart.plantsmart28
genes.data <- subset(genes(txdb), seqnames %in% c("1","2","3","4","5")) ## only nuclear genes are considered
genes.data <- as.data.frame(genes.data)
exons.data <- as.data.frame(exons(txdb))
cds.data <- as.data.frame(cds(txdb))

## Load all and circadian genes
alias2symbol.table <- AnnotationDbi::select(org.At.tair.db, 
                                            keys=keys(org.At.tair.db, keytype="ENTREZID"), 
                                            columns=c("SYMBOL", "TAIR"), keytype="ENTREZID")
# alias2symbol.table <- subset(alias2symbol.table, genes %in% TAIR)
agis <- alias2symbol.table$TAIR
agis <- intersect(agis, genes)
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias <- alias[agis]
alias[is.na(alias)] <- "" 
genes.selectize <- paste(names(alias), alias, sep=" - ")


## Setting conversion between alias and agis
agis <-alias2symbol.table$TAIR
names(agis) <- alias2symbol.table$SYMBOL
agis[is.na(agis)] <- ""

## Color vectors
line.colors <- c("blue","red", "darkgreen","black","#663300","#99003d","#b3b300","#4d0039","#4d2600","#006666","#000066","#003300","#333300","#660066")
area.colors <- c("skyblue","salmon", "lightgreen","lightgrey","#ffcc99","#ff99c2","#ffffb3","#ffe6f9","#ffe6cc","#80ffff","#b3b3ff","#99ff99","#e6e600","#ffb3ff")

## Load chromosome sequences
chr1 <- getSequence(read.fasta(file = "data/athaliana_genome/chr1.fa",seqtype = "AA"))[[1]]
chr2 <- getSequence(read.fasta(file = "data/athaliana_genome/chr2.fa",seqtype = "AA"))[[1]]
chr3 <- getSequence(read.fasta(file = "data/athaliana_genome/chr3.fa",seqtype = "AA"))[[1]]
chr4 <- getSequence(read.fasta(file = "data/athaliana_genome/chr4.fa",seqtype = "AA"))[[1]]
chr5 <- getSequence(read.fasta(file = "data/athaliana_genome/chr5.fa",seqtype = "AA"))[[1]]

## Function to compute the reverse complement
reverse.complement <- function(dna.sequence)
{
  return(c2s(comp(rev(s2c(dna.sequence)),forceToLower = FALSE)))
}

## Load Position Weight Matrices
## Open file connection
con <- file("data/jaspar_motifs/pfm_plants_20180911.txt",open = "r")

## Empty list for storing PWM
motifs.pwm <- vector(mode="list",length = 453)
motif.ids <- vector(mode="character",length=453)
motif.names <- vector(mode="character",length=453)

## Load 64 PWM
for(j in 1:453)
{
  ## First line contains motif id and name
  first.line <- readLines(con,1)
  
  motif.ids[j] <- strsplit(first.line,split=" ")[[1]][1]
  motif.names[j] <- strsplit(first.line,split=" ")[[1]][2]
  
  ## Next four line contians probabilites for each nucleotide
  a.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  c.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  g.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  t.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  
  ## Construct PWM
  motif.pwm <- matrix(nrow = 4,ncol=length(a.row))
  
  motif.pwm[1,] <- a.row
  motif.pwm[2,] <- c.row 
  motif.pwm[3,] <- g.row
  motif.pwm[4,] <- t.row
  
  rownames(motif.pwm) <- c("A","C","G","T")
  
  motifs.pwm[[j]] <- prop.table(motif.pwm,2)
}

## Close file connection
close(con)

## Naming list with PWM
names(motifs.pwm) <- motif.names
names(motif.ids) <- motif.names

## Load bigwig files for each transcription factor
bigwig.files <- c("data/bw_files/PHYA.bw",
                  "data/bw_files/PHYB_FLAG_27_1.bw",
                  "data/bw_files/PRR5_1.bw",
                  "data/bw_files/TOC1.bw",
                  "data/bw_files/CCA1_ZT02.bw",
                  "data/bw_files/CCA1_ZT14_1.bw",
                  "data/bw_files/LHY_1.bw",
                  "data/bw_files/CRY2.bw",
                  "data/bw_files/FHY1.bw",
                  "data/bw_files/LUX_ZT10.bw",
                  "data/bw_files/LUX_ZT12.bw",
                  "data/bw_files/PIF3.bw",
                  "data/bw_files/PIF4.bw",
                  "data/bw_files/PIF5.bw",
                  "data/bw_files/PRR7.bw",
                  "data/bw_files/PRR9_1.bw",
                  "data/bw_files/ELF3_ZT0.bw",
                  "data/bw_files/ELF3_ZT4.bw",
                  "data/bw_files/ELF4.bw")

names(bigwig.files) <- c("PHYA ZT00", "PHYB ZT00" ,"PRR5 ZT10", "TOC1 ZT15","CCA1 ZT02","CCA1 ZT14","LHY ZT02","CRY2 ZT08","FHY1 ZT04","LUX ZT10", "LUX ZT12","PIF3 ZT08","PIF4 ZT04","PIF5 ZT04","PRR7 ZT12","PRR9 ZT04","ELF3 ZT00", "ELF3 ZT04", "ELF4 ZT10")

## Load bed files for each transcription factor
bed.files <- c("data/bed_files/PHYA_peaks.narrowPeak",
               "data/bed_files/PHYB_peaks.narrowPeak",
               "data/bed_files/PRR5_1_peaks.narrowPeak",
               "data/bed_files/TOC1_1_peaks.narrowPeak",
               "data/bed_files/CCA1_ZT02_peaks.narrowPeak",
               "data/bed_files/CCA1_ZT14_peaks.narrowPeak",
               "data/bed_files/LHY_1_peaks.narrowPeak",
               "data/bed_files/CRY2_peaks.narrowPeak",
               "data/bed_files/FHY1_peaks.narrowPeak",
               "data/bed_files/LUX_ZT10_1_peaks.narrowPeak",
               "data/bed_files/LUX_ZT12_1_peaks.narrowPeak",
               "data/bed_files/PIF3_peaks.narrowPeak",
               "data/bed_files/PIF4_peaks.narrowPeak",
               "data/bed_files/PIF5_peaks.narrowPeak",
               "data/bed_files/PRR7_peaks.narrowPeak",
               "data/bed_files/PRR9_1_peaks.narrowPeak",
               "data/bed_files/ELF3_ZT0_1_peaks.narrowPeak",
               "data/bed_files/ELF3_ZT4_1_peaks.narrowPeak",
               "data/bed_files/ELF4_1_peaks.narrowPeak")

names(bed.files) <- c("PHYA ZT00", "PHYB ZT00" ,"PRR5 ZT10", "TOC1 ZT15","CCA1 ZT02",
                      "CCA1 ZT14","LHY ZT02","CRY2 ZT08","FHY1 ZT04","LUX ZT10", 
                      "LUX ZT12", "PIF3 ZT08","PIF4 ZT04","PIF5 ZT04","PRR7 ZT12",
                      "PRR9 ZT04","ELF3 ZT00", "ELF3 ZT04", "ELF4 ZT10")

## TF binding sites colors and symbol shapes
symbol.shapes <- c(17, 18, 19, 15)
symbol.color <- c("blue", "red", "darkgreen", "magenta")

## Colors used to represent repression/activation/neutrality in clock visualizer
repressor.color <- "firebrick1"
activator.color <- "seagreen3"
neutral.color <- "lightgrey"

## Colors to represent gene expression profiles in clock visualizer
selected.colors <- c("blue4","blue","deepskyblue","gold","firebrick","gray47")
peak.times <- c("peak20","peak0","peak4","peak8","peak12","peak16")
names(selected.colors) <- peak.times

## Node colors to represent in the global transcriptional network
node.colors <- selected.colors[network.data$peak.zt]
names(node.colors) <- NULL

## Auxiliary function to determine surrounding ZTs in clock visualizer
zts <- c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20")
zts.to.consider <- function(zt.point, zts=zts)
{
  zts.numeric <- seq(from=0,to=20,by=4)
  
  if(zt.point %in% zts.numeric)
  {
    return(c(zt.point,zt.point))
  } else
  {
    next.zt <- zts.numeric[which(zts.numeric >= zt.point)[1]]
    previous.zt <- zts.numeric[which(zts.numeric >= zt.point)[1] - 1]
    return(c(previous.zt, next.zt))
  }
}

# Circle and profile parameters for clock visualizer
radius.1 <- 100 #Outer circle radius
height <- 4 ## highest point in ylim for profile plot

#Function for radian conversion
radian.conversion <- function(alpha)
{
  rad <- (alpha*pi/180)
  return(rad)
}

## Generate coordinates for inner and outer circle in the clock representation for 
## clock visualizer
angle <- seq(from=0, to=2*pi, by=0.01)
x.circle.1 <- radius.1*sin(angle)
y.circle.1 <- radius.1*cos(angle)

radius.2 <- radius.1 - radius.1/12
x.circle.2 <- radius.2 * sin(angle)
y.circle.2 <- radius.2 * cos(angle)

## Define vectrors for location of transcription factors in the clock representation 
## in clock visualizer
agi.tfs <- c("AT2G46830", "AT1G01060", "AT5G61380", "AT5G24470", "AT5G02810", 
             "AT2G46790","AT1G09570", "AT2G18790", "AT1G04400", "AT2G37678", "AT3G46640", 
             "AT1G09530", "AT2G43010", "AT3G59060", "AT2G40080", "AT2G25930")
name.tfs <- c("CCA1", "LHY",  "TOC1", "PRR5", "PRR7", "PRR9", "PHYA", "PHYB", 
              "CRY2", "FHY1", "LUX", "PIF3", "PIF4", "PIF5", "ELF4", "ELF3")
agi.tfs.zts <- list(c("ZT02","ZT14"),
                    c("ZT02"),c("ZT15"),c("ZT10"),c("ZT12"),c("ZT04"),c("ZT00"),c("ZT00"),
                    c("ZT08"),c("ZT04"),c("ZT10","ZT12"),c("ZT08"),c("ZT04"),c("ZT04"),
                    c("ZT10"),c("ZT00","ZT04"))

## Pasting transcription factors with ZTs
tfs.with.zts <- c()
for (i in 1:length(name.tfs))
{
  tfs.with.zts <- c(tfs.with.zts, paste(name.tfs[i], agi.tfs.zts[[i]], sep= "_") )
}

## Determine the number of ZTs points for each transcription factor to represent
## this multiplicity in the clock for clock visualizer
agi.tfs.zts.multiplicity <- sapply(agi.tfs.zts,length)
names(agi.tfs.zts) <- agi.tfs
names(agi.tfs.zts.multiplicity) <- agi.tfs
names(name.tfs) <- agi.tfs

## Extract adjacency matrix
adj.global.matrix <- as.matrix(network.data[,tfs.with.zts])
rownames(adj.global.matrix) <- network.data$names

## Extract expression profiles
mean.expression <- as.matrix(network.data[,zts])
rownames(mean.expression) <- network.data$names

## Computing angles for each transcription factor according to its ZT and
## the muber of transcription factors in that ZT to represent tfs avoiding overlapping. 
splitted.tfs.names <- strsplit(tfs.with.zts,split="_")
tfs.angles <- vector(mode="numeric",length=length(tfs.with.zts))
tfs.zts <- vector(mode="numeric",length=length(tfs.with.zts))

for(i in 1:length(splitted.tfs.names))
{
  tfs.zts[i] <- substr(x=splitted.tfs.names[i][[1]][2],start = 3,stop=nchar(splitted.tfs.names[i][[1]][2]))
  tfs.angles[i] <- radian.conversion(15*as.numeric(tfs.zts[i]))
}

zt.multiplicity <- table(tfs.zts)
## Compute coordinates for each transcription factor setting a radius 
## and height to avoid the overlap
radius.to.multiply <- vector(mode="numeric",length=length(splitted.tfs.names))
height.to.multiply <- vector(mode="numeric",length=length(splitted.tfs.names))
node.labels <- vector(mode="numeric",length=length(splitted.tfs.names))
for(i in 1:length(splitted.tfs.names))
{
  node.labels[i] <- splitted.tfs.names[i][[1]][1]
  current.zt <- substr(x=splitted.tfs.names[i][[1]][2],start=3,stop=nchar(splitted.tfs.names[i][[1]][2]))
  current.multiplicity <- zt.multiplicity[current.zt]
  radius.to.multiply[i] <- (1 - (0.16*current.multiplicity))*radius.1
  height.to.multiply[i] <- (1 - (0.1*current.multiplicity))*height
  zt.multiplicity[current.zt] <- zt.multiplicity[current.zt] - 1
}

get.first <- function(my.vector)
{
  return(my.vector[[1]])
}

# names(height.to.multiply) <- tfs.names
names(height.to.multiply) <- tfs.with.zts


## Set the x.y coordinates for the positions 
tfs.x <- radius.to.multiply * sin(tfs.angles)
tfs.y <- radius.to.multiply * cos(tfs.angles)

## Generating a positions matrix 
matrix.pos <- matrix(data = c(tfs.x, tfs.y), nrow = length(tfs.x), ncol = 2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title, introductory text and sidebar navigation
  fluidRow(
    column(
      width = 2,
      img(src='attractor_logo_2.jpg', align = "center", width=200),
      tags$br(),
      radioButtons(inputId = "navigation_bar", width="100%",selected="home",
                   label="",
                   choices=c(
                      "Home" = "home",
                     "Individual gene analysis" = "individual_gene",
                     "Multiple transcription factor analysis" = "multiple_gene",
                     "Tutorials" = "tutorials",
                     "GitHub repository" = "github",
                     "Citation" = "citation"
                   ))),
    column(
      width = 8,
      tags$div(align = "center", 
               tags$h1(tags$b("ATTRACTOR,"),tags$i("Arabidopsis Thaliana"), "TRanscriptionAl Circadian neTwORk")),
      tags$br(),tags$br(),
      conditionalPanel(condition = "input.navigation_bar == 'home'",
        tags$div(align="justify", "The", tags$b("circadian clock"), "and", tags$b("light signalling"), "play central roles in", 
          tags$i("plant physiology"), "and", tags$i("development."), "As a consequence, massive amounts of",
          tags$b("omics data"), " have been generated to characterise their individual components. Nonetheless, 
          these data remains fragmented and researchers who want to explore the joint regulation exherted by the 
          circadian clock and light signalling need to consult different papers and resources making imperative 
          the used of", tags$b("molecular systems biology"), "techniques to integrate and make easily accesible 
          all the generated information."),
        tags$div(align="justify", tags$b("ATTRACTOR"),", is a web based tool for the analysis of the synergistic transcriptional control 
          exherted by the circadian clock and light signalling over genes exhibiting ryhtmic expression profiles in the model plant ", 
          tags$i(tags$b("Arabidopsis thaliana.")), tags$b("ATTRACTOR"), ", consists of a ", tags$b("transcriptional network"), 
          " that integrates transcriptomic data collected over diurnal cycles with 12 hours of light and 12 hours of darkness 
          with cistromic data generated using ChIP-seq for key transcriptional factors and regulators in the circadian clock 
          and light signalling. Specifically, our network is composed of 5778 nodes or genes with diurnal rythmic expression profiles and
          14529 edges or transcriptional regulations. The transcription factors and regulators included in our network comprise the
          components of the morning and central loops CCA1, LHY, the pseudo response regulator family members TOC1,
          PRR5, PRR7 and PRR9; as well as some components of the evening loop such as LUX, ELF3 and ELF4. In order to capture
          synergistic regulations with light signalling we added the light sensors and transcriptional regulators phytochromes
          PHYA and PHYB, the cryptochrome CRY2 as well as the light transcriptional factors from the phytochrome interacting factor
          family PIF5, PIF4 and PIF3. Finally, the phytochrome interacting transcriptional factor FHY1 (Far-red elongated Hypocotyl 1)
          is also included in our network."),
        tags$div(align="justify","Use the navigation bar on the left to check the different utilities in ATTRACTOR.")
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'individual_gene'",
        tags$div(align="justify", tags$b("ATTRACTOR"), "allows researchers to explore the coordinated regulation of several 
                 transcription factors or regulators over a selected individual gene as well as the effect observed in its
                 expression profile. Follow the following steps:"),
        tags$div(align="justify", "lalala" )
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'multiple_gene'",
                       tags$div(align="justify", tags$b("ATTRACTOR"), "allows researchers to explore the coordinated regulation of several 
                 transcription factors or regulators over their common targets."),
                       tags$div(align="justify", "lalala" )
      )
      
    ),
    column(
      width = 2,
      img(src='logo_ibvf.jpg', align = "center", width=100),
      img(src='logo_us.png', align = "center", width=100),
      tags$br(),tags$br(),tags$br(),
      img(src='logo_csic.jpg', align = "center", width=100)
    )
  ),
  
  ## Separation for different tools
  tags$br(),tags$br(),
  
  ## Conditional panel for individual gene analysis
  conditionalPanel(condition = "input.navigation_bar == 'individual_gene'",
                   fluidRow(
                     column(width = 3,
                            ## Select target gene to study
                            selectizeInput(inputId = "target.gene",
                                           label = "Gene",
                                           choices = genes.selectize,
                                           multiple = FALSE),
                            ## Check box for the TF peaks to represent
                            checkboxGroupInput(inputId = "selected.tfs",
                                               label = "Select Transcription Factors:",
                                               choices = list("PHYA ZT00", "PHYB ZT00", "ELF3 ZT00", "CCA1 ZT02", "LHY ZT02", "ELF3 ZT04","FHY1 ZT04", "PIF4 ZT04", "PIF5 ZT04", "PRR9 ZT04", "CRY2 ZT08", "PIF3 ZT08", "ELF4 ZT10", "PRR5 ZT10", "LUX ZT10", "PRR7 ZT12", "LUX ZT12","CCA1 ZT14", "TOC1 ZT15"),
                                               inline = TRUE,width = "100%")
                     ),
                     
                     column(width = 9,
                            tabsetPanel(type = "tabs",
                                        tabPanel(title = "Clock Visualizer", 
                                                 plotOutput(outputId = "clock",width = 600, height=600)),
                                        tabPanel(title = "Peak Visualizer",
                                                 column(wellPanel(
                                                   ## Numeric input for promoter length
                                                   numericInput(inputId = "promoter.length",
                                                                label = "Promoter Length",
                                                                value = 2000,
                                                                min = 500,
                                                                max = 2000,
                                                                step = 100),
                                                   ## Numeric input for 5' length
                                                   numericInput(inputId = "fiveprime.length",
                                                                label = "5' Length",
                                                                value = 500,
                                                                min = 100,
                                                                max = 500,
                                                                step = 100)),width=3),
                                                 column(wellPanel(
                                                   ## Selectize to choose target gene to represent
                                                   selectizeInput(inputId = "selected.motifs",
                                                                  label = "Select Motifs",
                                                                  choices = motif.names,
                                                                  multiple = TRUE),
                                                   
                                                   ## Checkbox to select all available motifs
                                                   checkboxInput(inputId = "all.motifs",
                                                                 label = "Select All Motifs:",
                                                                 value = FALSE),
                                                   
                                                   ## Numeric input for PWM score
                                                   numericInput(inputId = "min.score.pwm", 
                                                                label = "Motif Identification Score:",
                                                                value = 100, 
                                                                min = 80,
                                                                max = 100,
                                                                step = 5)),width=9),
                                                 
                                                 actionButton(inputId = "go",label = "GO"),
                                                 
                                                 fluidRow(
                                                   column(
                                                     plotOutput(outputId = "plot.to.chulo.to.wapo"),
                                                     width=12)
                                                 )),
                                        tabPanel(title = "Expression Visualizer", 
                                                 plotOutput(outputId = "expression",width = 600, height=600)))
                     )
                   )
  ),
  
  
  ## Conditional panel for multiple transcription factors and regulators analysis
  conditionalPanel(condition = "input.navigation_bar == 'multiple_gene'",
                   fluidRow(
                     column(width = 3,
                            ## Check box for the TFs to analyse
                            checkboxGroupInput(inputId = "selected.multiple.tfs",
                                               label = "Select Transcription Factors:",
                                               choices = list("PHYA ZT00", "PHYB ZT00", "ELF3 ZT00", "CCA1 ZT02", 
                                                              "LHY ZT02", "ELF3 ZT04","FHY1 ZT04", "PIF4 ZT04", 
                                                              "PIF5 ZT04", "PRR9 ZT04", "CRY2 ZT08", "PIF3 ZT08", 
                                                              "ELF4 ZT10", "PRR5 ZT10", "LUX ZT10", "PRR7 ZT12", 
                                                              "LUX ZT12","CCA1 ZT14", "TOC1 ZT15"),
                                               inline = TRUE,width = "100%"),
                            tags$b("Select a specific rythmic gene expression pattern with peak at:"),
                            selectInput(inputId = "peak", label="", 
                                        choices = c("Any ZT",paste("ZT",seq(from=0,to=20,by=4),sep="")), selected = NULL,
                                        multiple = FALSE, selectize = TRUE),
                            selectInput(inputId = "trough", label="and trough at:", 
                                        choices = c("Any ZT",paste("ZT",seq(from=0,to=20,by=4),sep="")), selected = NULL,
                                        multiple = FALSE, selectize = TRUE),
                            checkboxInput(inputId =  "edges",label = "Visualize Edges",value = FALSE),
                            actionButton(inputId = "go_multiple",label = "GO")
                     ),
                     
                     column(width = 9,
                            tabsetPanel(type = "tabs",
                                        tabPanel(title = "Network Visualization",
                                                 plotOutput("networkPlot")
                                        ),
                                        tabPanel(title = "Gene Table",
                                                 dataTableOutput(outputId = "outputTable"),
                                                 uiOutput(outputId = "download_ui_for_table")
                                         ),
                                         tabPanel(title = "GO Enrichment",
                                            tabsetPanel(type = "tabs",
                                                       tabPanel(title = "GO map"),
                                                       tabPanel(title = "GO barplot"),
                                                       tabPanel(title = "GO concept network"),
                                                       tabPanel(title = "GO table")
                                            )
                                         ),
                                        tabPanel(title = "Pathway Enrichment"),
                                        tabPanel(title = "TFBS Enrichment")
                            )
                     )
                                                 
                   )
  )
                     

)

## ATTRACTOR server
server <- function(input, output) {

  ## clock visualizer code
  output$clock <- renderPlot({
    
    ## Error message for the user
    validate(
      need(input$selected.tfs, "Please select some transcription factor")
    )
    
    ## Extracting agi ID for selected gene
    target.agi <- strsplit(x = input$target.gene, split = " - ")[[1]][1]
    
    #Plot circle
    par(mar=c(0,0,0,0))
    plot(x.circle.1,y.circle.1, type = "l", lwd=3, axes=FALSE, xlab = "", ylab="",xlim=c(-1.2 * radius.1, 1.2 * radius.1),ylim=c(-1.2 * radius.1, 1.2 * radius.1))
    lines(x.circle.2, y.circle.2, lwd=3)
    x.polygon <- c(sin(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
                   sin(seq(from=-pi, to=0, by=0.01))* radius.1)
    y.polygon <-c(cos(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
                  cos(seq(from=-pi, to=0, by=0.01))*radius.1)
    polygon(x = x.polygon, y = y.polygon, col = "black")
    for (i in 0:5)
    {
      angle.zt <- radian.conversion(alpha = 60*i)
      zt <- 4*i
      current.zt <- paste("ZT", zt,  sep = "")
      text(x = (radius.1 + radius.1/6)*sin(angle.zt), y = (radius.1 + radius.1/6)*cos(angle.zt), labels = current.zt,cex = 1.5,font=2)
      lines(x = c(radius.1 * sin(angle.zt), (radius.1 + radius.1/20)* sin(angle.zt)), 
            y = c(radius.1 * cos(angle.zt), (radius.1 + radius.1/20)* cos(angle.zt)), lwd=2)
    }
    
    ## Get agis and alias of selected tfs and extracting data from adj.global.matrix 
    ## for clock visualizer
    selected.tfs.alias <- sapply(X=strsplit(input$selected.tfs,split=" "),FUN = get.first)
    selected.tfs.agi <- agis[selected.tfs.alias]
    selected.tfs.zts.pasted <- vector(mode = "character", length = length(input$selected.tfs))
    for (i in 1:length(input$selected.tfs))
    {
      selected.tfs.zts.pasted[i] <- paste(strsplit(input$selected.tfs[i], split=" ")[[1]], collapse = "_")
    }
    
    
    ## Set transcription factors to keep in the matrix
    to.keep <- rep(FALSE,ncol(adj.global.matrix))
    for(i in 1:length(input$selected.tfs))
    {
      to.keep <- (to.keep | grepl(selected.tfs.zts.pasted[i],colnames(adj.global.matrix)))
    }
    
    ## Update adjacency matrix to create the transcriptional network 
    ## represented in the clock visualizer
    updated.adj.matrix.to.represent <- adj.global.matrix[target.agi,to.keep]
    for (i in 1:length(input$selected.tfs))
    {
      updated.adj.matrix.to.represent <- rbind(rep(0,length(input$selected.tfs)),
                                               updated.adj.matrix.to.represent)
    }
    
    updated.adj.matrix.to.represent <- cbind(updated.adj.matrix.to.represent, rep(0,length(input$selected.tfs)+1))
    colnames(updated.adj.matrix.to.represent) <- c(selected.tfs.zts.pasted, alias[target.agi])
    rownames(updated.adj.matrix.to.represent) <- c(selected.tfs.zts.pasted, alias[target.agi])
    
    ## Modify adj.matrix and matrix.pos to add the target.gene
    gene.peak.str <- subset(network.data, names == target.agi)$peak.zt 
    gene.peak <- as.numeric(substr(x=gene.peak.str,start=5,stop=nchar(gene.peak.str)))
    target.color <- selected.colors[paste0("peak",gene.peak)]
    
    ## Transpose the matrix to create graph adjacency with igraph
    new.matrix <- t(updated.adj.matrix.to.represent)
    
    ## Generating the complete network
    tfs.network <- graph.adjacency(adjmatrix = new.matrix, mode = "directed",weighted = TRUE)
    edge.weights <- E(tfs.network)$weight
    
    if (length(input$selected.tfs) > 1)
    {
      ## Edge colors
      for(k in 1:length(edge.weights))
      {
        if(edge.weights[k] == 1)
        {
          E(tfs.network)$color[k] <- activator.color #"darkgreen"
        } else if(edge.weights[k] == -1)
        {
          E(tfs.network)$color[k] <- repressor.color  #"darkred"
        } else if(edge.weights[k] == 2)
        {
          print(1)
          E(tfs.network)$color[k] <- neutral.color
        }
      }
      
    }
        
    ## Vertex colors
    node.colors <- vector(mode="character",length=nrow(new.matrix))
    for(k in 1:(length(node.colors)-1))
    {
      if(new.matrix[k,ncol(new.matrix)] == 1)
      {
        node.colors[k] <- activator.color
      } else if(new.matrix[k,ncol(new.matrix)] == -1)
      {
        node.colors[k] <- repressor.color
      } else #if (new.matrix[k,ncol(new.matrix)] == 0)
      {
        node.colors[k] <- neutral.color
      }
    }
    
    node.colors[length(node.colors)] <- target.color
    
    ## Modify the angles, radius, and labels positions to keep only the selected tfs
    tfs.angles <- tfs.angles[to.keep]
    radius.to.multiply <- radius.to.multiply[to.keep]
    node.labels <- node.labels[to.keep]
    
    ## Modify the angles, the radius and the positions to add the new node
    new.tfs.angles <- c(tfs.angles, radian.conversion(gene.peak*15))
    new.multiply <- c(radius.to.multiply, radius.1*0.3)
    tfs.x <- new.multiply * sin(new.tfs.angles)
    tfs.y <- new.multiply * cos(new.tfs.angles)
    
    new.matrix.pos <- matrix(data = c(tfs.x, tfs.y), nrow = nrow(new.matrix), ncol = 2)
    new.node.labels <- c(node.labels, alias[target.agi])
    
    ## Plot the network
    plot.igraph(tfs.network, layout=new.matrix.pos, add = TRUE, rescale=FALSE, vertex.size=radius.1*13,
                vertex.color = node.colors, vertex.label=new.node.labels, edge.arrow.size = 0.8, 
                edge.arrow.width=2, edge.curved= TRUE, edge.width = 4, vertex.label.dist = 0,
                vertex.label.cex=1, vertex.label.font=2,vertex.label.color="black",label.font=2)
    
  })
  
  ## Express visualizer code
  output$expression <- renderPlot({
    
    ## Error message for the user
    validate(
        need(input$selected.tfs, "Please select some transcription factor")
      )
    
    
    ## Get agis and alias of selected tfs and extracting data from adj.global.matrix 
    ## for clock visualizer
    selected.tfs.alias <- sapply(X=strsplit(input$selected.tfs,split=" "),FUN = get.first)
    selected.tfs.agi <- agis[selected.tfs.alias]
    selected.tfs.zts.pasted <- vector(mode = "character", length = length(input$selected.tfs))
    for (i in 1:length(input$selected.tfs))
    {
      selected.tfs.zts.pasted[i] <- paste(strsplit(input$selected.tfs[i], split=" ")[[1]], collapse = "_")
    }
    
    ## Get expression data for the selected gene
    target.agi <- strsplit(x = input$target.gene, split = " - ")[[1]][1]
    gene.expression <- as.vector(scale(mean.expression[target.agi,]))
    gene.expression <- c(gene.expression, gene.expression[1])
    extended.gene.expression <- approx(x = seq(from=0,to=24,by=4), y = gene.expression, xout=c(0,2,4,8,10,12,14,15,16,20,24))
    extended.gene.expression.values <- extended.gene.expression$y
    names(extended.gene.expression.values) <- c("ZT00", "ZT02", "ZT04", "ZT08", "ZT10", "ZT12", "ZT14", "ZT15", "ZT16", "ZT20", "ZT24")
    line.color <- selected.colors[network.data[target.agi, "peak.zt"]]
    
    ## Plot the initial visualization of the expression profile
    plot(x=seq(from=0,to=24,by=4),gene.expression,
         type="o",lwd=5,cex=1.5,
         ylim=c(-2.5,height),xlim=c(0,24),
         col=line.color,axes=FALSE,xlab="",ylab="", 
         main=paste(target.agi, alias[target.agi],sep=" - "))
    
    ## Add TFs to expression profile
    for(i in 1:length(selected.tfs.agi))
    {
      current.tf.name <- names(selected.tfs.agi[i])
      current.zt <- strsplit(x = selected.tfs.zts.pasted[i], split="_")[[1]][2]
      current.time.point <- as.numeric(substr(x = current.zt, start = 3, stop = nchar(current.tf.zt)))
      # current.time.point <- as.numeric(substr(x = current.tf.zts[j], start = 3, stop = nchar(current.tf.zts[j])))
      current.regulation <- adj.global.matrix[target.agi,selected.tfs.zts.pasted[i]]
        
      if(current.regulation == 1)
      {
        point.color <- activator.color #"seagreen3"#"darkgreen"
        arrow.angle <- 45
        draw.tf <- TRUE
      } else if (current.regulation == -1)
      {
        point.color <- repressor.color # "firebrick1"
        arrow.angle <- 90
        draw.tf <- TRUE
      } else if (current.regulation == 2)
      {
        point.color <- neutral.color
        arrow.angle <- 45
        draw.tf <- TRUE
      } else if (current.regulation == 0)
      {
        draw.tf <- FALSE
      }
      
      if(draw.tf)
      {
        arrows(x0 = current.time.point, y0 = height.to.multiply[selected.tfs.zts.pasted[i]],
               x1 = current.time.point ,y1= extended.gene.expression.values[current.zt] + 0.3,lwd=4,angle=arrow.angle,length=0.05,col=point.color)
        points(x = current.time.point,y=height.to.multiply[selected.tfs.zts.pasted[i]],lwd=4,cex=4, col=point.color, pch = 19)  
        text(x = current.time.point,y=height.to.multiply[selected.tfs.zts.pasted[i]],labels = current.tf.name, cex=1.2, font=2 )
      }
    }
    
    polygon(x=c(0,12,12,0),y=c(-2,-2,-2.3,-2.3),lwd=2)
    polygon(x=c(12,24,24,12),y=c(-2,-2,-2.3,-2.3),col = "black",lwd=2)
    
    axis(side = 2,at = -2:2,labels = FALSE,lwd=2)
    mtext("Normalized Gene Expression",side = 2,line = 1.3,cex = 1.5,at = 0)
    axis(side = 1,at=seq(from=0,to=24,by=2),line=-1,las=2,labels = paste("ZT",seq(from=0,to=24,by=2),sep=""),lwd=2)
    
  })
  
  ## Peak visualizer code
  observeEvent(input$go,{
    
    ## Extract target gene annotation 
    gene.name <-  strsplit(input$target.gene,split=" - ")[[1]][1]
    
    target.gene.body <- genes.data[gene.name,]
    target.gene.chr <- as.character(target.gene.body$seqnames)
    target.gene.start <- target.gene.body$start
    target.gene.end <- target.gene.body$end
    
    target.gene.strand <- as.character(target.gene.body$strand)
    
    ## Extract cds annotation
    cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
    
    ## Extract exons annotation
    exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
    
    ## Determine the genome range to plot including promoter, gene body and 5' UTR
    ## This depends on whether the gene is on the forward or reverse strand
    range.to.plot <- target.gene.body
    
    if(target.gene.strand == "+")
    {
      range.to.plot$start <- range.to.plot$start - input$promoter.length
      range.to.plot$end <- range.to.plot$end + input$fiveprime.length
    } else if (target.gene.strand == "-")
    {
      range.to.plot$end <- range.to.plot$end + input$promoter.length
      range.to.plot$start <- range.to.plot$start - input$fiveprime.length
    }
    
    ## Compute the length of the genome range to represent
    current.length <- range.to.plot$end - range.to.plot$start
    
    ## Determine upper limit of the graph
    number.tfs <- length(input$selected.tfs)
    upper.lim <- 25 * number.tfs
    
    ## Draw DNA strand
    gene.height <- -25
    cord.x <- 1:current.length
    output$plot.to.chulo.to.wapo <- renderPlot({
      
      ## Sanity checks
      validate(
        need(length(input$selected.tfs) > 0 , "Please select a set of transcription factors")
      )
      
      validate(
        need((length(input$selected.motifs) > 0 || input$all.motifs), "Please select a set of DNA motifs")
      )
      
      plot(cord.x, rep(gene.height,length(cord.x)),type="l",col="black",lwd=3,ylab="",
           cex.lab=2,axes=FALSE,xlab="",main="",cex.main=2,
           ylim=c(-30,upper.lim),xlim=c(-3000,max(cord.x)))
      
      ## Extract exons for target gene
      exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Transform exon coordinates to current range
      min.pos <- min(exons.data.target.gene$start)
      
      if(target.gene.strand == "+")
      {
        exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$promoter.length
        exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$promoter.length
      } else if(target.gene.strand == "-")
      {
        exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$fiveprime.length
        exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$fiveprime.length
      }
      
      ## Represent exons
      exon.width <- 2
      for(i in 1:nrow(exons.data.target.gene))
      {
        # Determine start/end for each exon
        current.exon.start <- exons.data.target.gene$start[i]
        current.exon.end <- exons.data.target.gene$end[i]
        
        ## Determine coordinates for each exon polygon and represent it
        exon.x <- c(current.exon.start,current.exon.end,current.exon.end,current.exon.start)
        exon.y <- c(gene.height + exon.width, gene.height + exon.width, gene.height - exon.width, gene.height - exon.width)
        
        polygon(x = exon.x, y = exon.y, col = "blue",border = "blue")
      }
      
      ## Extract cds for target gene
      cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
      
      ## Transform cds coordinates to current range
      if(target.gene.strand == "+")
      {
        cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$promoter.length
        cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$promoter.length
      } else if (target.gene.strand == "-")
      {
        cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$fiveprime.length
        cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$fiveprime.length
      }
      
      cds.width <- 3
      for(i in 1:nrow(cds.data.target.gene))
      {
        # Determine current cds start/end
        current.cds.start <- cds.data.target.gene$start[i]
        current.cds.end <- cds.data.target.gene$end[i]
        
        # Determine curret cds coordinates for the polygon and represent it
        cds.x <- c(current.cds.start,current.cds.end,current.cds.end,current.cds.start)
        cds.y <- c(gene.height + cds.width, gene.height + cds.width, gene.height - cds.width, gene.height - cds.width)
        
        polygon(x = cds.x, y = cds.y, col = "blue",border = "blue")
      }
      
      ## Draw arrow to represent transcription direction 
      if(target.gene.strand == "+")
      {
        lines(c(input$promoter.length,input$promoter.length,input$promoter.length+100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
        lines(c(input$promoter.length+50,input$promoter.length+100),y=c(gene.height+6,gene.height+5),lwd=3)
        lines(c(input$promoter.length+50,input$promoter.length+100),y=c(gene.height+4,gene.height+5),lwd=3)
      } else if (target.gene.strand == "-")
      {
        lines(c(current.length - input$promoter.length, current.length - input$promoter.length, current.length - input$promoter.length-100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
        lines(c(current.length - input$promoter.length-50, current.length - input$promoter.length - 100),y=c(gene.height + 6, gene.height + 5),lwd=3)
        lines(c(current.length - input$promoter.length-50, current.length - input$promoter.length - 100),y=c(gene.height + 4, gene.height + 5),lwd=3)
      }
      
      ## Draw promoter range
      if(target.gene.strand == "+")
      {
        axis(side = 1,labels = c(- input$promoter.length, - input$promoter.length / 2,"TSS"),at = c(1,input$promoter.length/2,input$promoter.length),lwd=2,cex=1.5,las=2,cex=2)
      } else if(target.gene.strand == "-")
      {
        axis(side = 1,labels = c("TSS",- input$promoter.length / 2,- input$promoter.length),at = c(current.length-input$promoter.length,current.length-input$promoter.length/2, current.length),lwd=2,cex=1.5,las=2,cex=2)
      }
      
      selected.bigwig.files <- bigwig.files[input$selected.tfs]
      selected.bed.files <- bed.files[input$selected.tfs]
      
      ## Since ChIPpeakAnno needs more than one region to plot our region
      ## is duplicated 
      regions.plot <- GRanges(rbind(range.to.plot,range.to.plot))
      
      ## Import signal from the bigwig files
      cvglists <- sapply(selected.bigwig.files, import, 
                         format="BigWig", 
                         which=regions.plot, 
                         as="RleList")
      
      names(cvglists) <- input$selected.tfs
      
      ## Compute signal in the region to plot
      chip.signal <- featureAlignedSignal(cvglists, regions.plot, 
                                          upstream=ceiling(current.length/2), 
                                          downstream=ceiling(current.length/2),
                                          n.tile=current.length) 
      
      ## Compute mean signal 
      chip.signal.means <- matrix(nrow=number.tfs, ncol=ncol(chip.signal[[1]]))
      
      for(i in 1:length(input$selected.tfs))
      {
        if(target.gene.strand == "+")
        {
          chip.signal.means[i, ] <- colMeans(chip.signal[[i]],na.rm = TRUE)
        } else if (target.gene.strand == "-")
        {
          chip.signal.means[i, ] <- rev(colMeans(chip.signal[[i]],na.rm = TRUE))
        }
      }
      
      ## Draw peak regions for each TF and determing TF binding sequences
      
      ## Determine TFBS motifs to search for
      if(input$all.motifs)
      {
        selected.motifs.pwm <- motifs.pwm
      } else
      {
        selected.motifs.pwm <- motifs.pwm[input$selected.motifs]
      }
      
      selected.motif.names <- names(selected.motifs.pwm)
      selected.motif.ids <- motif.ids[selected.motif.names]
      
      ## Initialize data frame containing TF binding sequences in the peak regions
      df.hits <- data.frame(0,0,"","","")
      colnames(df.hits) <- c("tf_number","position","id","name","seq")
      
      ## Width of the rectangule representing the peak reagion
      peak.width <- 1
      for(i in 1:length(input$selected.tfs))
      {
        ## Extract bed file name 1 and read it
        current.bed.file <- selected.bed.files[i]
        current.peaks <- read.table(file=current.bed.file,header = F, as.is = T)
        peak.coordinates <- subset(current.peaks, V1 == range.to.plot$seqnames & V2 >= range.to.plot$start & V3 <= range.to.plot$end) 
        current.peaks.to.plot <- peak.coordinates[,2:3]
        
        ## Transform coordinates 
        current.peaks.to.plot <- current.peaks.to.plot - range.to.plot$start
        
        ## Check if there are peaks for the target gene
        if(nrow(current.peaks.to.plot) > 0)
        {
          ## Normalization
          chip.signal.means[i, ] <- 10 * chip.signal.means[i, ] / max(chip.signal.means[i, ])
          
          #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
          for(j in 1:nrow(current.peaks.to.plot))
          {
            ## Extract start and end point of each peak region
            current.peak.start <- current.peaks.to.plot[j,1]
            current.peak.end <- current.peaks.to.plot[j,2]
            
            ## Computer coordinates for polygon and draw it
            peak.x <- c(current.peak.start,current.peak.end,
                        current.peak.end,current.peak.start)
            peak.y <- c(25*(i - 1) - 5 + peak.width, 25*(i - 1) - 5 + peak.width, 
                        25*(i - 1) - 5 - peak.width, 25*(i - 1) - 5 - peak.width)  
            
            polygon(x = peak.x, y = peak.y, col = area.colors[i], border = line.colors[i],lwd=2)
            
            ## Identify TF binding DNA motifs 
            peak.chr <- peak.coordinates[j, 1]
            peak.start <- peak.coordinates[j, 2]
            peak.end <- peak.coordinates[j, 3]
            
            ## Extract peak sequence
            if(peak.chr == "1")
            {
              peak.sequence <- c2s(chr1[peak.start:peak.end])
            } else if(peak.chr == "2")
            {
              peak.sequence <- c2s(chr2[peak.start:peak.end])
            } else if(peak.chr == "3")
            {
              peak.sequence <- c2s(chr3[peak.start:peak.end])
            } else if(peak.chr == "4")
            {
              peak.sequence <- c2s(chr4[peak.start:peak.end])
            } else if(peak.chr == "5")
            {
              peak.sequence <- c2s(chr5[peak.start:peak.end])
            }
            
            peak.rev.comp.sequence <- reverse.complement(peak.sequence)
            
            for(k in 1:length(selected.motifs.pwm))
            {
              motif.pwm <- selected.motifs.pwm[[k]]
              
              hits.fw <- matchPWM(motif.pwm, peak.sequence, 
                                  min.score = paste0(input$min.score.pwm,"%"))
              hits.fw.seqs <- as.data.frame(hits.fw)[[1]]
              hits.fw <- as(hits.fw, "IRanges")
              hits.fw.start <- start(hits.fw)
              hits.fw.end <- end(hits.fw)
              
              if(length(hits.fw.start) > 0)
              {
                df.hits.fw <- data.frame(rep(i,length(hits.fw.start)),
                                         ((hits.fw.start+hits.fw.end)/2) + current.peak.start,
                                         rep(selected.motif.ids[k],length(hits.fw.start)),
                                         rep(selected.motif.names[k],length(hits.fw.start)),
                                         hits.fw.seqs)
                colnames(df.hits.fw)  <- c("tf_number","position","id","name","seq")
                df.hits <- rbind(df.hits,df.hits.fw)
              }
              
              hits.rev <- matchPWM(motif.pwm, peak.rev.comp.sequence, 
                                   min.score = paste0(input$min.score.pwm,"%"))
              hits.rev.seqs <- as.data.frame(hits.rev)[[1]]
              hits.rev.seqs <- sapply(hits.rev.seqs,reverse.complement)
              names(hits.rev.seqs) <- NULL
              
              hits.rev <- as(hits.rev, "IRanges")
              hits.rev.start <- nchar(peak.sequence) - end(hits.rev) + 1
              hits.rev.end <- nchar(peak.sequence) - start(hits.rev) + 1
              
              if(length(hits.rev.start) > 0)
              {
                df.hits.rev <- data.frame(rep(i,length(hits.rev.start)),
                                          ((hits.rev.start+hits.rev.end)/2) + current.peak.start,
                                          rep(selected.motif.ids[k],length(hits.rev.start)),
                                          rep(selected.motif.names[k],length(hits.rev.start)),
                                          hits.rev.seqs)
                colnames(df.hits.rev)  <- c("tf_number","position","id","name","seq")
                df.hits <- rbind(df.hits,df.hits.rev)
              }
              
            }
            
          }
        }
      }
      
      ## Remove first line of the data frame added just for technical reason
      df.hits <- df.hits[-1,]
      
      ## Draw TF binding sites
      detected.tfbs <- unique(as.vector(df.hits$name))
      
      number.of.shapes <- ceiling(length(detected.tfbs) / length(symbol.color))
      
      necessary.shapes <- rep(symbol.shapes[1:number.of.shapes],each = length(detected.tfbs)/number.of.shapes)
      necessary.colors <- rep(symbol.color,number.of.shapes)
      
      if(length(detected.tfbs) > 0)
      {
        for(i in 1:length(detected.tfbs))
        {
          current.tfbs <- detected.tfbs[i]
          current.shape <- necessary.shapes[i]
          current.color <- necessary.colors[i]
          
          positions <- subset(df.hits, name == current.tfbs)
          
          for(j in 1:nrow(positions))
          {
            tf.to.draw <- positions$tf_number[j]
            pos.to.draw <- positions$position[j]
            
            points(x = pos.to.draw, y = 25*(tf.to.draw - 1) - 5 - 5*peak.width,
                   pch = current.shape, col = current.color, cex = 1)
          }
        }
        
        ## Add legend for TFBS
        legend.step <- 10
        for(i in 1:length(detected.tfbs))
        {
          points(x = -3000, y = upper.lim - (i-1)*legend.step, 
                 pch=necessary.shapes[i], col = necessary.colors[i],cex = 1)
          
          
          current.seq <- as.character(subset(df.hits,name == detected.tfbs[i])[["seq"]][[1]])
          current.label <- paste(c(detected.tfbs[i], "  -  ", current.seq ),collapse="")
          
          text(x = -2900, y = upper.lim - (i-1)*legend.step, labels = current.label,
               adj = 0,cex = 0.7)
        }
      }
      
      ## Draw profiles for TF binding
      for(i in 1:number.tfs)
      {
        ## Compute base line for current TF
        current.base.line <- 25 * (i - 1)
        
        ## Represent signal from the current TF
        lines(chip.signal.means[i,]+current.base.line,type="l",col=line.colors[i],lwd=3)
        
        ## Determine polygon coordinates and represent it
        cord.y <- c(current.base.line,chip.signal.means[i,]+current.base.line,current.base.line)
        cord.x <- 1:length(cord.y)
        
        polygon(cord.x,cord.y,col=area.colors[i])
        
        text(x = -50,y = 25*(i-1) + 12,labels = input$selected.tfs[i],adj = 1,col = line.colors[i],font = 2)
      }
      
    })
    
  })
  
  ## Multiple transcription factor code
  
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
  
  ## Determine common targets and perfomr analysis when button is clicked
  observeEvent(input$go_multiple, {

    ## Determine targets of selected TFs
    selected.only.tfs <- sapply(X=strsplit(input$selected.multiple.tfs,split=" "),FUN = get.first)
    #input$selected.multiple.tfs
    
    gene.selection <- rowSums(network.data[,selected.only.tfs]) == length(selected.only.tfs)
    selected.genes.df <- network.data[gene.selection,]    
    
    ## Node colors for representation
    selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
    
    ## Target gene representation on the network
    network.representation <- ggplot(network.data, aes(x.pos,y.pos)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=node.colors,size=1) +
      geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    
    ## Add edges to the network when selected
    if(input$edges)
    {
      for(i in 1:length(input$selected.multiple.tfs))
      {
        tf.xpos <- subset(network.data, names == tf.ids[selected.only.tfs[i]])[["x.pos"]]
        tf.ypos <- subset(network.data, names == tf.ids[selected.only.tfs[i]])[["y.pos"]]
        network.representation <- network.representation +
          annotate("segment",
                   x=rep(tf.xpos,nrow(selected.genes.df)),
                   y=rep(tf.ypos,nrow(selected.genes.df)),
                   xend=selected.genes.df$x.pos,
                   yend=selected.genes.df$y.pos, 
                   color="grey", arrow=arrow(type="closed",length=unit(0.1, "cm")))
      }
      
      network.representation <- network.representation + 
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
    }
    
    ## Update network representation on the app
    output$networkPlot <- renderPlot({
      network.representation
    },height = 700)
    
    
  })
  

}

# Run the application 
shinyApp(ui = ui, server = server)

