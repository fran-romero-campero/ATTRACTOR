library(shiny)
## Load libraries
library(shiny)
library(ChIPpeakAnno)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(Biostrings)
library(seqinr)
library(org.At.tair.db)

## Load and extract Arabidopsis thaliana annotation regarding genes, exons and cds 
txdb <- TxDb.Athaliana.BioMart.plantsmart28
genes.data <- subset(genes(txdb), seqnames %in% c("1","2","3","4","5")) ## only nuclear genes are considered
genes.data <- as.data.frame(genes.data)
exons.data <- as.data.frame(exons(txdb))
cds.data <- as.data.frame(cds(txdb))

## Load all and circadian genes
my.key <- keys(org.At.tair.db, keytype="ENTREZID")
my.col <- c("SYMBOL", "TAIR")
alias2symbol.table <- select(org.At.tair.db, keys=my.key, columns=my.col, keytype="ENTREZID")
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 
genes <- paste(names(alias), alias, sep=" - ")

circadian.genes <- read.table(file="data/genes_info/gene_info.txt",header=FALSE,col.names = c("gene.name","symbol"),fill = TRUE,as.is=TRUE)
circadian.genes <- paste(circadian.genes$gene.name, circadian.genes$symbol, sep = " - ")

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

names(bed.files) <- c("PHYA ZT00", "PHYB ZT00" ,"PRR5 ZT10", "TOC1 ZT15","CCA1 ZT02","CCA1 ZT14","LHY ZT02","CRY2 ZT08","FHY1 ZT04","LUX ZT10", "LUX ZT12", "PIF3 ZT08","PIF4 ZT04","PIF5 ZT04","PRR7 ZT12","PRR9 ZT04","ELF3 ZT00", "ELF3 ZT04", "ELF4 ZT10")

## TF binding sites colors and symbol shapes
symbol.shapes <- c(17, 18, 19, 15)
symbol.color <- c("blue", "red", "darkgreen", "magenta")


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
                   # Only show this panel if the peak visualizer is selected
                   conditionalPanel(condition = "input$selected_analysis = 'peakvis'",
                                    tags$head(
                                      tags$style(HTML("
                                                      .shiny-output-error-validation {
                                                      color: red;
                                                      font-size: 20px;
                                                      }
                                                      "))
                                      ),
                                    
                                    #shinythemes::themeSelector(),
                                    # Application title
                                    titlePanel("Peak Visualizer"),
                                    
                                    # Sidebar with a slider input for number of bins 
                                    sidebarLayout(
                                      fluidRow(
                                        
                                        column(wellPanel(
                                          ## Selectize to choose target gene to represent
                                          selectizeInput(inputId = "target.gene",
                                                         label = "Gene",
                                                         choices = genes,
                                                         multiple = FALSE),
                                          
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
                                                       step = 100)),width=4),
                                        
                                        column(wellPanel(
                                          ## Check box for the TF peaks to represent
                                          checkboxGroupInput(inputId = "names.tfs",
                                                             label = "Select Transcription Factors:",
                                                             choices = list("PHYA ZT00", "PHYB ZT00", "ELF3 ZT00", "CCA1 ZT02", "LHY ZT02", "ELF3 ZT04","FHY1 ZT04", "PIF4 ZT04", "PIF5 ZT04", "PRR9 ZT04", "CRY2 ZT08", "PIF3 ZT08", "ELF4 ZT10", "PRR5 ZT10", "LUX ZT10", "PRR7 ZT12", "LUX ZT12","CCA1 ZT14", "TOC1 ZT15"),
                                                             
                                                             
                                                             #list("CCA1 ZT02", "CCA1 ZT14","LHY", "TOC1", "PRR5", "PRR7", "PRR9", "PHYA","PHYB",
                                                             #               "CRY2","FHY1","LUX","PIF3","PIF4","PIF5","ELF3 ZT00","ELF3 ZT04","ELF4"),
                                                             inline = TRUE,width = "100%")),width=4),
                                        
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
                                                       label = "Minimum Score for Motif Identification:",
                                                       value = 100, 
                                                       min = 80,
                                                       max = 100,
                                                       step = 5)),width=4),
                                        actionButton(inputId = "go",label = "GO")
                                      ),
                                      
                                      # fluidRow(
                                      #   column(
                                      #     ## Action button to trigger plot drawing
                                      #     actionButton(inputId = "go",label = "GO"),width=12)),
                                      
                                      # Show a plot of the generated distribution
                                      fluidRow(
                                        column(
                                          plotOutput(outputId = "plot.to.chulo.to.wapo"),
                                          width=12)
                                      )
                                    )
                                      )
                                    

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

