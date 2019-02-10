# R script for pre-processing 
# Copyright (C) 2018  Francisco J. Romero-Campero, Pedro de los Reyes Rodríguez
# Ana Belén Romero Losada
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public
# License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Authors: Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
#          Francisco J. Romero-Campero
# 
# Contact: Francisco J. Romero-Campero - fran@us.es 
# Date: September 2018

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

names(bed.files) <- c("PHYA ZT00", "PHYB ZT00" ,"PRR5 ZT10", "TOC1 ZT15","CCA1 ZT02","CCA1 ZT14","LHY ZT02","CRY2 ZT08","FHY1 ZT04","LUX ZT10", "LUX ZT12", "PIF3 ZT08","PIF4 ZT04","PIF5 ZT04","PRR7 ZT12","PRR9 ZT??","ELF3 ZT00", "ELF3 ZT04", "ELF4 ZT10")
                    
## TF binding sites colors and symbol shapes
symbol.shapes <- c(17, 18, 19, 15)
symbol.color <- c("blue", "red", "darkgreen", "magenta")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
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

# Define server logic required to draw a histogram
server <- function(input, output) {

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
    number.tfs <- length(input$names.tfs)
    upper.lim <- 25 * number.tfs
    
    ## Draw DNA strand
    gene.height <- -25
    cord.x <- 1:current.length
    output$plot.to.chulo.to.wapo <- renderPlot({
      
      ## Sanity checks
      validate(
        need(length(input$names.tfs) > 0 , "Please select a set of transcription factors")
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
      
      selected.bigwig.files <- bigwig.files[input$names.tfs]
      selected.bed.files <- bed.files[input$names.tfs]

      ## Since ChIPpeakAnno needs more than one region to plot our region
      ## is duplicated 
      regions.plot <- GRanges(rbind(range.to.plot,range.to.plot))
      
      ## Import signal from the bigwig files
      cvglists <- sapply(selected.bigwig.files, import, 
                         format="BigWig", 
                         which=regions.plot, 
                         as="RleList")
      
      names(cvglists) <- input$names.tfs
      
      ## Compute signal in the region to plot
      chip.signal <- featureAlignedSignal(cvglists, regions.plot, 
                                          upstream=ceiling(current.length/2), 
                                          downstream=ceiling(current.length/2),
                                          n.tile=current.length) 
      
      ## Compute mean signal 
      chip.signal.means <- matrix(nrow=number.tfs, ncol=ncol(chip.signal[[1]]))
      
      for(i in 1:number.tfs)
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
      for(i in 1:number.tfs)
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
        
        text(x = -50,y = 25*(i-1) + 12,labels = input$names.tfs[i],adj = 1,col = line.colors[i],font = 2)
      }

    })

  })
      

  
}

# Run the application 
shinyApp(ui = ui, server = server)

