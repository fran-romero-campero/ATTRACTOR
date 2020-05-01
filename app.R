# Authors: Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
#          Francisco J. Romero-Campero
# 
# Contact: Francisco J. Romero-Campero - fran@us.es 

## Input to test 
## input <- list(selected.multiple.tfs = c("CCA1 ZT02", "PRR5 ZT10"), peak = "peak0", trough = "trough12")
## input <- list(selected.multiple.tfs = c("CCA1 ZT02", "PRR5 ZT10"), peak = "peak8", trough = "any")
## input <- list(selected.multiple.tfs = c("CCA1 ZT02", "CCA1 ZT14", "PRR5 ZT10"), peak = "peak8", trough = "any")
## input <- list(selected.multiple.tfs = c("CCA1 ZT02", "CCA1 ZT14", "PRR5 ZT10", "TOC1 ZT15"), peak = "peak8", trough = "any")

# Load neccesary libraries
library(shiny)
library(shinythemes)
#library(ChIPpeakAnno)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(Biostrings)
library(seqinr)
library(org.At.tair.db)
library(igraph)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(pathview)
library(shinycssloaders)
library(shinyWidgets)
library(shinyjs)
library(SuperExactTest)
library(VennDiagram)
library(Rsamtools)

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

## Load normalized gene expression data for animation
norm.data <- read.table(file = "data/normalized_gene_expression.txt",header = T,sep = "\t")

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
ath.universe <- alias2symbol.table$TAIR
alias2symbol.table <- subset(alias2symbol.table, TAIR %in% genes)
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
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
chr1 <- FaFile(file = "data/athaliana_genome/chr1.fa") #getSequence(read.fasta(file = "data/athaliana_genome/chr1.fa",seqtype = "AA"))[[1]]
chr2 <- FaFile(file = "data/athaliana_genome/chr2.fa") #getSequence(read.fasta(file = "data/athaliana_genome/chr2.fa",seqtype = "AA"))[[1]]
chr3 <- FaFile(file = "data/athaliana_genome/chr3.fa") #getSequence(read.fasta(file = "data/athaliana_genome/chr3.fa",seqtype = "AA"))[[1]]
chr4 <- FaFile(file = "data/athaliana_genome/chr4.fa") #getSequence(read.fasta(file = "data/athaliana_genome/chr4.fa",seqtype = "AA"))[[1]]
chr5 <- FaFile(file = "data/athaliana_genome/chr5.fa") #getSequence(read.fasta(file = "data/athaliana_genome/chr5.fa",seqtype = "AA"))[[1]]

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

names(bigwig.files) <- c("PHYA ZT00", "PHYB ZT00" ,"PRR5 ZT10", "TOC1 ZT15","CCA1 ZT02",
                         "CCA1 ZT14","LHY ZT02","CRY2 ZT08","FHY1 ZT04","LUX ZT10", 
                         "LUX ZT12","PIF3 ZT08","PIF4 ZT04","PIF5 ZT04","PRR7 ZT12",
                         "PRR9 ZT04","ELF3 ZT00", "ELF3 ZT04", "ELF4 ZT10")

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

## Function to generate output table to download
create.downloadable.output.table <- function(input.gene.df,alias,tfs.names)
{
  output.selected.genes.df <- data.frame(matrix(nrow=nrow(input.gene.df), ncol=6))
  colnames(output.selected.genes.df) <- c("AGI ID", "Gene Name", "Gene Description", "Regulators","Expression Peak Time","Expression Trough Time")
  output.selected.genes.df$`Gene Description` <- input.gene.df$description
  
  for(i in 1:nrow(output.selected.genes.df))
  {
    output.selected.genes.df[i,1] <- input.gene.df[i,1]
    output.selected.genes.df[i,2] <- alias[input.gene.df[i,1]]
    output.selected.genes.df[i,4] <- paste(tfs.names[which(input.gene.df[i,tfs.names] == 1)],collapse=", ")
    output.selected.genes.df[i,5] <-paste0("ZT",substr(input.gene.df[i,"peak.zt"],start=5,stop=nchar(input.gene.df[i,"peak.zt"])))
    output.selected.genes.df[i,6] <-paste0("ZT",substr(input.gene.df[i,"trough.zt"],start=7,stop=nchar(input.gene.df[i,"trough.zt"])))
  }
  
  return(output.selected.genes.df)
}

## Auxiliary function to compute enrichments for GO table
compute.enrichments <- function(gene.ratios, bg.ratios)
{
  gene.ratios.eval <- sapply(parse(text=gene.ratios),FUN = eval)
  bg.ratios.eval <- sapply(parse(text=bg.ratios),FUN = eval)
  enrichments <- round(x=gene.ratios.eval/bg.ratios.eval,digits = 2)
  enrichments.text <- paste(enrichments, " (", gene.ratios, "; ", bg.ratios, ")",sep="")
  
  return(enrichments.text)  
}

#GO links and tair link functions
go.link <- function(go.term)
{
  link <- paste0("http://amigo.geneontology.org/amigo/term/", go.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           go.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

gene.link.function <- function(gene.name)
{
  tair.link <- paste(c("https://www.arabidopsis.org/servlets/TairObject?name=",
                       gene.name,
                       "&type=locus"),collapse="")
  gene.link <- paste(c("<a href=\"",
                       tair.link,
                       "\" target=\"_blank\">",
                       gene.name, "</a>"),
                     collapse="")
  return(gene.link)
}

## KEGG pathway link
kegg.pathway.link <- function(kegg.pathway)
{
  link <- paste0("https://www.genome.jp/kegg-bin/show_pathway?",kegg.pathway)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           kegg.pathway, "</a>"),
                         collapse = "")
  return(complete.link)
}

## KEGG module link
kegg.module.link <- function(kegg.module)
{
  link <- paste0("https://www.genome.jp/kegg-bin/show_module?",kegg.module)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           kegg.module, "</a>"),
                         collapse = "")
  return(complete.link)
}

## TFBS jaspar link
tfbs.link <- function(motif.id)
{
  link <- paste0("http://jaspar.genereg.net/matrix/",motif.id)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           motif.id, "</a>"),
                         collapse = "")
  return(complete.link)
}


## Red gradient for animation
##red.gradient <- colorRampPalette(c("red", "white"))
##current.red.gradient <- c(red.gradient(5),rep("#FFFFFF",15))

# Define UI for application that draws a histogram
ui <- fluidPage(
  ##shinythemes::themeSelector(),
  theme = shinytheme("sandstone"),
  ##theme = shinytheme("simplex"),
  ##theme = shinytheme("journal"),
  
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
      
      conditionalPanel(condition = "input.navigation_bar == 'github'",
        tags$div(align = "justify", tags$b("ATTRACTOR,"), "is entirely developed using 
        the R package", tags$b( tags$a(href="https://shiny.rstudio.com/", "shiny.")), "The 
        source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
        tags$b("GitHub."), "If you experience any problem using ATTRACTOR please create an", tags$b(tags$a(href="https://github.com/fran-romero-campero/ATTRACTOR/issues","issue")), "in GitHub and we will address it."),
        tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/ATTRACTOR", "ATTRACTOR at GitHub"))))
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'citation'",
                       tags$div(align = "justify", "We are strongly committed to", tags$b("open access software"), 
                       "and", tags$b("open science."),"Following our philosophy we have deposited our GitHub code 
                       into", tags$a(href="https://zenodo.org/record/3780022#.XqwzNvlS9uQ", target="_blank",tags$b("Zenodo")), ", a
                       general-purpose open-access repository developed under the", 
                       tags$a(href="https://www.openaire.eu/", target="_blank", tags$b("European OpenAIRE program.")), "Meanwhile we publish 
                       our work in a journal if you find", tags$b("ATTRACTOR"), "useful in your research we would be most grateful if you cite 
                       our GitHub repository with a,", tags$b("DOI"),  "as follows:",
                       tags$br(),
                       tags$br(),
                       tags$div(tags$h4(tags$b("de los Reyes, P., Romero-Losada, A.B., Romero-Campero, F.J. (2020) ATTRACTOR, 
                       Arabidopsis Thaliana TRanscriptionAl Circadian neTwORk v1.0, Zenodo, doi:10.5381/zenodo.3780022")))),
                       
                       tags$br(),
                       tags$br(),
                       tags$div(align="center", img(src='smiley.png', align = "center", width=200,hight=200)),
                       tags$br()
                       
      ),
      
      
      conditionalPanel(condition = "input.navigation_bar == 'home'",
        tags$div(align="justify", "The", tags$b("circadian clock"), "and", tags$b("light signalling"), "play central roles in", 
          tags$i("plant physiology"), "and", tags$i("development."), "As a consequence, massive amounts of",
          tags$b("omics data"), " have been generated to characterize their individual components. Nonetheless, 
          these data remain fragmented and researchers who want to explore the synergistic regulation exerted by the 
          circadian clock and light signalling need to consult different papers and resources making imperative 
          the use of", tags$b("molecular systems biology"), "techniques to integrate and make easily accesible 
          all the generated information."),
        tags$div(align="justify", tags$b("ATTRACTOR,"),"is a web based tool for the analysis of the synergistic transcriptional control 
          exerted by the circadian clock and light signalling over genes exhibiting rhythmic expression profiles in the model plant ", 
          tags$i(tags$b("Arabidopsis thaliana.")), tags$b("ATTRACTOR,"), "consists of a ", tags$b("transcriptional network"), 
          " that integrates transcriptomic data collected over diurnal cycles with 12 hours of light and 12 hours of darkness 
          with cistromic data generated using ChIP-seq for key transcriptional factors and regulators in the circadian clock 
          and light signalling. Specifically, our network is composed of 5778 nodes or genes with diurnal rhythmic expression profiles and
          14529 edges or transcriptional regulations. The transcription factors and regulators included in our network comprise the
          components of the morning and central loops", tags$b("CCA1, LHY,"), "the pseudo response regulator family members", 
          tags$b("TOC1,PRR5, PRR7"), "and ", tags$b("PRR9;"), "as well as some components of the evening loop such as", tags$b("LUX, ELF3"), "and", tags$b("ELF4."),
          "In order to capture synergistic regulations with light signalling we added the light sensors and transcriptional regulators phytochromes",
          tags$b("PHYA"), "and", tags$b("PHYB,"),"the cryptochrome ", tags$b("CRY2"), "as well as the light transcriptional factors from the phytochrome interacting factor
          family", tags$b("PIF5, PIF4"), " and ", tags$b("PIF3."),"Finally, the phytochrome interacting transcriptional factor", tags$b("FHY1"), "(Far-red elongated Hypocotyl 1)
          is also included in our network."),
        tags$div(align="justify","Use the navigation bar on the left to explore the different utilities in ATTRACTOR or alternatively",
                 tags$a(href="https://www.youtube.com/watch?v=8o2otN-DY4c&t=1220s", target="_blank", tags$b("view our video tutorial."))),
        tags$br(), tags$br()#,
        # actionButton("run", "Run Animation"),
        # fluidRow(column(width = 9,
        #                 plotOutput("networkAnimation")),
        #          column(width = 3,
        #                 tags$br(), tags$br(),tags$br(), tags$br(),tags$br(), tags$br(),tags$br(), tags$br(),
        #                 plotOutput("clockAnimation")))
        # plotOutput("networkAnimation"),
        # tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(),
        # plotOutput("clockAnimation")
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'individual_gene'",
        tags$div(align="justify", tags$b("ATTRACTOR"), "allows researchers to explore the coordinated regulation of several 
                 transcription factors or regulators over an", tags$b("individually selected"), "gene as well as the effect observed in its
                 expression profile. Follow the steps below:",
                 tags$ol(
                   tags$li("Select a specific gene from our network using the", tags$b("Target Gene"),
                           "dropdown menu on the left below. You can enter either the AGI identifier or primary symbol for the",
                           "gene of interest."), 
                   tags$li("Select several transcription factors (TFs) to explore their regulation over the
                           previously selected target gene using the",tags$b("Select Transcription Factors"), "checkboxes on the left below."), 
                   tags$li("Results will be depicted on the tabs below. The", tags$b("Clock Visualizer tab"), "shows a circular 
                            representation of the diurnal cycle with the selected TFs and target gene located at
                            the time point of their expression peak. Binding of TFs to the target gene and 
                            the observed effect over its expression are represented using colored arrows. The", tags$b("Expression 
                            Visualizer tab"),  "shows the target gene expression profile in a linear representation with the selected
                            TFs located over the time point of their expression peak. The binding of 
                            the selected TFs to the target gene promoter and the observed effect is represented using colored arrows.
                            The", tags$b("Peak Visualizer tab"), "shows the selected target gene genomic location and the peaks detected
                            in our analysis of the correspondig TFs ChIP-seq data. Here you can specify the length of the gene promoter and 
                            3' region as well as DNA TFs binding motifs to search for in the detected peak regions with the specified score
                            for identity.")
                   # tags$li("On the tab", tags$b("Clock Visualizer"), "a circular representation of the diurnal cycle under study
                   #         will be depicted. The selected transcription factors will be located close to the edge of the circle at the
                   #         specific time point when the corresponding ChIP-seq data were generated. The selected target gene will be 
                   #         located towards the center of the circle pointing at the specific time point where its expression is highest. 
                   #         An arrow will be drawn from a transcription factor to the target gene when the corresponding transcription factor
                   #         binds to the promoter of the target gene. The node of the transcription factors and arrows will be colored red to represent a decrease in the target gene expression
                   #         after binding of the transcription factor, green to represent an increase or grey when no change is observed.")
                 ))
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'multiple_gene'",
                       tags$div(align="justify", tags$b("ATTRACTOR"), "allows researchers to explore the coordinated regulation of several 
                                transcription factors or regulators over their common targets. Follow the steps below:", 
                                
                                tags$ol(
                                  tags$li("Select your TFs of interest using the", tags$b("Select Transcription Factors"),
                                          "checkbox menu on the left below."),
                                  tags$li("You can also select a cluster of genes exhibiting a specific rhythmic pattern of expression by chosing
                                          the time point for their peak and trough using the dropdown menus", 
                                          tags$b("Select a specific rhythmic gene expression pattern with peak at: ... and trough at:")),
                                  tags$li("Check the box", tags$b("Visualize Edges"), "when you want to depict arrows from TFs to their target genes."),
                                  tags$li("Click on the", tags$b("SELECT GENES"), "to visualize your selected TFs common target genes exhibiting the 
                                           specified rhythmic expression pattern in our transcriptional network. Explore the different tabs to 
                                           download a table with the selected genes, perform a signficance analysis of the overlap between the selected 
                                           TFs targets and the specified expression pattern as weel as", tags$b("GO term, pathways 
                                           and binding motifs enrichment analysis"), ".")
                                )
                                )
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
                                           label = "Target Gene:",
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
                                                 tags$div(align="center",
                                                          plotOutput(outputId = "clock",width = 600, height=600))),
                                        tabPanel(title = "Expression Visualizer", 
                                                 tags$div(align="center",
                                                          plotOutput(outputId = "expression",width = 600, height=600))),
                                        tabPanel(title = "Peak Visualizer",
                                                 column(wellPanel(
                                                   ## Numeric input for promoter length
                                                   numericInput(inputId = "promoter.length",
                                                                label = "Promoter Length",
                                                                value = 2000,
                                                                min = 500,
                                                                max = 2000,
                                                                step = 100),
                                                   ## Numeric input for 3' length
                                                   numericInput(inputId = "threeprime.length",
                                                                label = "3' Length",
                                                                value = 500,
                                                                min = 100,
                                                                max = 500,
                                                                step = 100)),width=3),
                                                 column(wellPanel(
                                                   ## Selectize to choose target gene to represent
                                                   selectizeInput(inputId = "selected.motifs",
                                                                  label = "Select Motifs",
                                                                  selected ="EE",
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
                                                 
                                                 #actionButton(inputId = "go",label = "GO"),
                                                 
                                                 fluidRow(
                                                   column(
                                                     plotOutput(outputId = "peak_plot"),
                                                     width=12)
                                                 ))
                                        )
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
                            tags$b("Select a specific rhythmic gene expression pattern with peak at:"),
                            selectInput(inputId = "peak", label="", 
                                        choices = c("Any ZT" = "any", "ZT0" = "peak0", "ZT4" = "peak4", "ZT8" = "peak8",
                                                                      "ZT12" = "peak12", "ZT16" = "peak16", "ZT20" = "peak20"),
                                                  selected = NULL, multiple = FALSE, selectize = TRUE),
                            selectInput(inputId = "trough", label="and trough at:", 
                                        choices = c("Any ZT" = "any", "ZT0" = "trough0", "ZT4" = "trough4", "ZT8" = "trough8",
                                                    "ZT12" = "trough12", "ZT16" = "trough16", "ZT20" = "trough20"),
                                        selected = NULL,
                                        multiple = FALSE, selectize = TRUE),
                            checkboxInput(inputId =  "edges",label = "Visualize Edges",value = FALSE),
                            actionButton(inputId = "go_multiple",label = "Select Genes")
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
                                        tabPanel(title = "Overlap Significance",
                                                 tags$br(),
                                                 tags$div(align="justify", "In this section, we present the results of a significance analysis of the
                                                           overlap between the targets of the selected transcription factors and gene
                                                           with a specific expresion pattern."),
                                                 tags$br(),
                                                 textOutput("overlap.message"),
                                                 textOutput("overlap.significance.text"),
                                                 tags$br(),
                                                 tags$br(),
                                                 tags$div(align="center",
                                                          plotOutput("venn.diagram.plot"))
                                        ),
                                        tabPanel(title = "Functional Enrichment",
                                                 tabsetPanel(type = "tabs",
                                                             tabPanel(title = "GO Enrichment",
                                                                      tags$br(),
                                                                      tags$div(align="justify", "In this section you can perform a GO term
                                                                               enrichment analysis over the selected genes. First
                                                                               of all, you need to choose the background set of genes between
                                                                               the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in ATTRACTOR:"),
                                                                      tags$br(),
                                                                      radioButtons(inputId = "go.background", width="100%",selected="allgenome",
                                                                                   label="",
                                                                                   choices=c(
                                                                                     "Complete genome" = "allgenome",
                                                                                     "Genes in network" = "onlynet"
                                                                                   )), tags$br(),
                                                                      actionButton(inputId = "goterm",label = "GO terms analysis"),
                                                                      tags$br(),
                                                                      tags$br(),
                                                                      shinyjs::useShinyjs(),
                                                                      hidden(div(id='loading.div',h3('Please be patient, computing GO enrichment ...'))),
                                                                      tags$br(),
                                                                      tabsetPanel(type = "tabs",
                                                                                  tabPanel(title = "GO table",
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "textGOTable"),
                                                                                           tags$br(), tags$br(),
                                                                                           dataTableOutput(outputId = "output_go_table"),
                                                                                           htmlOutput(outputId = "revigo"),
                                                                                           uiOutput(outputId = "download_ui_for_go_table")
                                                                                  ),
                                                                                  tabPanel(title = "GO map",
                                                                                           # withSpinner(ui_element =
                                                                                           # plotOutput(outputId = "gomap"),type = 4)),
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "gomap_text"),
                                                                                           tags$br(),
                                                                                           div(style= "overflow:scroll; height:500px;",
                                                                                               addSpinner(plotOutput("gomap"), spin = "circle", color = "#E41A1C"))),
                                                                                  tabPanel(title = "GO barplot",
                                                                                           tags$br(),
                                                                                           htmlOutput(outputId = "barplot_text"),
                                                                                           tags$br(),
                                                                                           # withSpinner(ui_element = 
                                                                                           #   plotOutput(outputId = "bar.plot"),type = 4)),#,inline=TRUE))),
                                                                                           addSpinner(plotOutput("bar.plot"), spin = "circle", color = "#E41A1C")),
                                                                                  tabPanel(title = "GO Enrichment Map",
                                                                                           tags$br(), 
                                                                                           htmlOutput(outputId = "emapplot_text"),
                                                                                           tags$br(),
                                                                                           addSpinner(plotOutput(outputId = "emap.plot",inline=TRUE))),
                                                                                  tabPanel(title = "GO concept network",
                                                                                           tags$br(), 
                                                                                           htmlOutput(outputId = "cnetplot_text"),
                                                                                           tags$br(),
                                                                                           addSpinner(plotOutput(outputId = "cnet.plot",inline=TRUE)))
                                                                      )
                                                                      ),
                                                             tabPanel(title = "KEGG Pathway Enrichment",
                                                                      tags$br(),
                                                                      tags$div(align="justify", "In this section you can perform a KEGG pathways and modules
                                                                               enrichment analysis over the selected genes. First
                                                                               of all, you need to choose the background set of genes between
                                                                               the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in ATTRACTOR:"),
                                                                      tags$br(),
                                                                      radioButtons(inputId = "pathway_background", width="100%",selected="allgenome",
                                                                                   label="",
                                                                                   choices=c(
                                                                                     "Complete genome" = "allgenome",
                                                                                     "Genes in network" = "onlynet"
                                                                                   )),
                                                                      actionButton(inputId = "pathway_button",label = "KEGG pathway analysis"),
                                                                      tags$br(),
                                                                      tags$br(),
                                                                      shinyjs::useShinyjs(),
                                                                      hidden(div(id='loading.div.kegg',h3('Please be patient, computing KEGG pathway enrichment ...'))),
                                                                      tags$br(),
                                                                      tabsetPanel(type = "tabs",
                                                                                  tabPanel(title = "Enriched Pathway Table",
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "no_kegg_enrichment"),
                                                                                           dataTableOutput(outputId = "output_pathway_table"),
                                                                                           uiOutput(outputId = "download_ui_for_kegg_table")
                                                                                  ),
                                                                                  tabPanel(title = "Enriched Pathway Visualization",
                                                                                           uiOutput(outputId = "kegg_selectize"),
                                                                                           imageOutput("kegg_image"),
                                                                                           br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                                                                           br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                                                                           br(), br(), br(), br(), br()
                                                                                  ),
                                                                                  tabPanel(title = "Enriched Module Table",
                                                                                           htmlOutput(outputId = "text_module_kegg"),
                                                                                           br(), br(),
                                                                                           dataTableOutput(outputId = "output_module_table")
                                                                                  )
                                                                      )
                                                                      )
                                                 )
                                        ),
                                        tabPanel(title = "TFBS Enrichment",
                                                 tags$br(),
                                                 tags$div(align="justify", "In this section you can perform a 
                                                          Transcription Factor Binding Sites (TFBS) enrichment analysis 
                                                          over the promoters of the selected genes. First of all, you 
                                                          have to set the some required parameters: the background and 
                                                          the length of what will be considered the promoter of each gene"),
                                                 column(wellPanel(
                                                   ## Select the background
                                                   tags$div(align="justify", "Please, choose the background set of genes between
                                                          the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in ATTRACTOR:"),
                                                   # tags$br(),
                                                   radioButtons(inputId = "tfbs_background", width="100%",selected="onlynet",
                                                                label="",
                                                                choices=c(
                                                                  "Complete genome" = "allgenome",
                                                                  "Genes in network" = "onlynet")
                                                                )),width = 4),
                                                   column(wellPanel(
                                                     tags$div(align="justify", "Promoter length upstream of the start codon:"),
                                                     radioButtons(inputId = "up_promoter", width="100%",selected="2000",
                                                                  label="",
                                                                  choices=c(
                                                                    "500 bp" = "500",
                                                                    "1000 bp" = "1000",
                                                                    "1500 bp" = "1500",
                                                                    "2000 bp" = "2000")
                                                                  )),width = 3),
                                                 column(wellPanel(
                                                   tags$div(align="justify", "Promoter length downstream of the start codon:"),
                                                   radioButtons(inputId = "down_promoter", width="100%",selected="500",
                                                                label="",
                                                                choices=c(
                                                                  "0 bp" = "0",
                                                                  "200 bp" = "200",
                                                                  "500 bp" = "500")
                                                                )),width = 3),
                                                 column(wellPanel(
                                                   tags$div(align="justify", "Motif Identification Score:"),
                                                   radioButtons(inputId = "score", width="100%",selected="90",
                                                                label="",
                                                                choices=c(
                                                                  "80" = "80",
                                                                  "85" = "85",
                                                                  "90" = "90",
                                                                  "95" = "95")
                                                   )),width = 2),
                                                 tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
                                                 tags$br(),
                                                 tags$br(),
                                                 tags$br(),
                                                 actionButton(inputId = "tfbs_button",label = "TFBS enrichment analysis"),
                                                 tags$br(),
                                                 tags$br(),
                                                 shinyjs::useShinyjs(),
                                                 hidden(div(id='loading.div.tfbs',h3('Please be patient, computing TFBS enrichment ...'))),
                                                 tags$br(),
                                                 dataTableOutput(outputId = "output_tfbs_table"),
                                                 uiOutput(outputId = "download_ui_tfbs_table")
                                                 )
                            )
                     )
                                                 
                   )
  )

)

## ATTRACTOR server

#library(ChIPpeakAnno)

server <- function(input, output, session) {

  ## Animation in main page with node size
  # rv <- reactiveValues(i = 0)
  # 
  # increase.step <- 2
  # max.steps <- 1000
  # 
  # output$networkAnimation <- renderPlot( {
  #   ggplot(network.data, aes(x.pos,y.pos)) + 
  #     theme(panel.background = element_blank(), 
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank(),
  #           axis.title = element_blank(),
  #           axis.text = element_blank(),
  #           axis.ticks.y = element_blank()) + 
  #     geom_point(fill=node.colors,size=1.65^norm.data[[(rv$i %% 48)+1]],pch=21)
  # },height = 600)
  # 
  # observeEvent(input$run, {
  #   rv$i <- 0
  #   observe({
  #     isolate({
  #       rv$i <- rv$i + increase.step
  #       print(rv$i)
  #     })
  #     
  #     if(rv$i < max.steps) {
  #       invalidateLater(5, session)
  #     }
  #   })
  # })
  
  ## Animation in main page with red gradient
  # rv <- reactiveValues(i = 0)
  # 
  # increase.step <- 2
  # increase.step.sec <- 0.2
  # max.steps <- 1000
  # 
  # output$networkAnimation <- renderPlot( {
  #   ggplot(network.data, aes(x.pos,y.pos)) + 
  #     theme(panel.background = element_blank(), 
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank(),
  #           axis.title = element_blank(),
  #           axis.text = element_blank(),
  #           axis.ticks.y = element_blank()) + 
  #     geom_point(fill=current.red.gradient[ceiling(1.7^norm.data[[(rv$i %% 48)+1]])],size=5,pch=21)
  # },height = 600, width = 600)
  # 
  # output$clockAnimation <- renderPlot({
  #   #Plot circle
  #   par(mar=c(0,0,0,0))
  #   plot(x.circle.1,y.circle.1, type = "l", lwd=3, axes=FALSE, xlab = "", ylab="",xlim=c(-1.2 * radius.1, 1.2 * radius.1),ylim=c(-1.2 * radius.1, 1.2 * radius.1))
  #   lines(x.circle.2, y.circle.2, lwd=3)
  #   x.polygon <- c(sin(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
  #                  sin(seq(from=-pi, to=0, by=0.01))* radius.1)
  #   y.polygon <-c(cos(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
  #                 cos(seq(from=-pi, to=0, by=0.01))*radius.1)
  #   polygon(x = x.polygon, y = y.polygon, col = "black")
  #   for (j in 0:5)
  #   {
  #     angle.zt <- radian.conversion(alpha = 60*j)
  #     zt <- 4*j
  #     current.zt <- paste("ZT", zt,  sep = "")
  #     text(x = (radius.1 + radius.1/6)*sin(angle.zt), y = (radius.1 + radius.1/6)*cos(angle.zt), labels = current.zt,cex = 1.5,font=2)
  #     lines(x = c(radius.1 * sin(angle.zt), (radius.1 + radius.1/20)* sin(angle.zt)), 
  #           y = c(radius.1 * cos(angle.zt), (radius.1 + radius.1/20)* cos(angle.zt)), lwd=2)
  #   }
  #   
  #   radio.flecha <- 80
  #   angle.zt <- radian.conversion(alpha = 8*rv$i)
  #   arrows(x0 = 0, y0 = 0, x1 = sin(angle.zt)*radio.flecha, y1 = cos(angle.zt)*radio.flecha,lwd = 5)
  #   
  # }, height = 300, width = 300)
  # 
  # observeEvent(input$run, {
  #   rv$i <- 0
  #   observe({
  #     isolate({
  #       rv$i <- rv$i + increase.step
  #       print(rv$i)
  #     })
  #   
  #     
  #     if(rv$i < max.steps) {
  #       invalidateLater(5, session)
  #     }
  #   })
  #   
  # })

    
  ## clock visualizer code
  output$clock <- renderPlot({
    
    ## Error messages for the user
    validate(
      need(input$selected.tfs, "Please select some transcription factor"),
      need(input$target.gene, "Please select a target gene")
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
    
    ## Edge colors
    if (length(edge.weights) == 0)
    {
      # E(tfs.network)$color[k] <- NULL
      # print("no edge")
    }else 
    {
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
        need(input$selected.tfs, "Please select some transcription factor"),
        need(input$target.gene, "Please select a target gene")
      )
    
    
    ## Get agis and alias of selected tfs and extracting data from adj.global.matrix 
    ## for expression visualizer
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
      current.time.point <- as.numeric(substr(x = current.zt, start = 3, stop = nchar(selected.tfs.zts.pasted[i])))
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
  output$peak_plot <- renderPlot({
      
      ## Sanity checks
      validate(
        need(length(input$selected.tfs) > 0 , "Please select a set of transcription factors"),
        need(input$target.gene, "Please select a target gene")
      )
      
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
      
      ## Determine the genome range to plot including promoter, gene body and 3' UTR
      ## This depends on whether the gene is on the forward or reverse strand
      range.to.plot <- target.gene.body
      
      if(target.gene.strand == "+")
      {
        range.to.plot$start <- range.to.plot$start - input$promoter.length
        range.to.plot$end <- range.to.plot$end + input$threeprime.length
      } else if (target.gene.strand == "-")
      {
        range.to.plot$end <- range.to.plot$end + input$promoter.length
        range.to.plot$start <- range.to.plot$start - input$threeprime.length
      }
      
      ## Compute the length of the genome range to represent
      current.length <- range.to.plot$end - range.to.plot$start
      
      ## Determine upper limit of the graph
      number.tfs <- length(input$selected.tfs)
      upper.lim <- 25 * length(input$selected.tfs)
      
      ## Draw DNA strand
      gene.height <- -25
      cord.x <- 1:current.length
      
      plot(cord.x, rep(gene.height,length(cord.x)),type="l",col="black",lwd=3,ylab="",
           cex.lab=2,axes=FALSE,xlab="",main="",cex.main=2,
           ylim=c(-30,25 * length(input$selected.tfs)),
           xlim=c(-3000,max(cord.x)))
      
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
        exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$threeprime.length
        exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$threeprime.length
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
        cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$threeprime.length
        cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$threeprime.length
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
      chip.signal <- ChIPpeakAnno::featureAlignedSignal(cvglists, regions.plot, 
                                                        upstream=ceiling(current.length/2), 
                                                        downstream=ceiling(current.length/2),
                                                        n.tile=current.length) 
      
      ## Compute mean signal 
      chip.signal.means <- matrix(nrow=length(input$selected.tfs), ncol=ncol(chip.signal[[1]]))
      
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
      
      ## Determine TFBS motifs to search for and Selecting an
      ## example motif if the user does not select any of them
      if(input$all.motifs)
      {
        selected.motifs.pwm <- motifs.pwm
      } else if(length(input$selected.motifs)==0)
      {
        selected.motifs.pwm <- motifs.pwm["EE"]
      }else
      {
        selected.motifs.pwm <- motifs.pwm[input$selected.motifs]
      }
      
      selected.motif.names <- names(selected.motifs.pwm)
      selected.motif.ids <- motif.ids[selected.motif.names]
      
      ## Initialize data frame containing TF binding sequences in the peak regions
      df.hits <- data.frame(0,0,"","","")
      colnames(df.hits) <- c("tf_number","position","id","name","seq")
      
      ## Width of the rectangle representing the peak region
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
              peak.sequence <- toString(getSeq(x = chr1)[[1]][peak.start:peak.end]) #c2s(chr1[peak.start:peak.end])
            } else if(peak.chr == "2")
            {
              peak.sequence <- toString(getSeq(x = chr2)[[1]][peak.start:peak.end]) #c2s(chr2[peak.start:peak.end])
            } else if(peak.chr == "3")
            {
              peak.sequence <- toString(getSeq(x = chr3)[[1]][peak.start:peak.end]) #c2s(chr3[peak.start:peak.end])
            } else if(peak.chr == "4")
            {
              peak.sequence <- toString(getSeq(x = chr4)[[1]][peak.start:peak.end]) #c2s(chr4[peak.start:peak.end])
            } else if(peak.chr == "5")
            {
              peak.sequence <- toString(getSeq(x = chr5)[[1]][peak.start:peak.end]) #c2s(chr5[peak.start:peak.end])
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
      for(i in 1:length(input$selected.tfs))
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
  
  ## Determine common targets and perform analysis when button is clicked
  observeEvent(input$go_multiple, {

    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }

    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    if(input$peak != "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, (peak.zt == input$peak & trough.zt == input$trough))
    } else if(input$peak == "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, trough.zt == input$trough)
    } else if(input$peak != "any" && input$trough == "any")
    {
      selected.genes.df <- subset(selected.genes.df, peak.zt == input$peak)
    }
    
    ## Node colors for representation
    selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
    
    ## Determine significance of overlap
    ## Number of sets to overlap
    if(input$peak != "any" || input$trough != "any")
    {
      number.of.sets <- length(input$selected.multiple.tfs) + 1
    } else
    {
      number.of.sets <- length(input$selected.multiple.tfs)
    }
    
    list.of.sets <- vector(mode = "list",length=number.of.sets)
    
    ## TFs targets
    for(i in 1:length(selected.tfs.with.zts))
    {
      list.of.sets[[i]] <- rownames(network.data)[which(network.data[,selected.tfs.with.zts[[i]]] != 0)]
    }
    
    ## Gene cluster with specified expression pattern
    if(input$peak != "any" && input$trough != "any")
    {
      list.of.sets[[length(list.of.sets)]] <- rownames(subset(network.data, (peak.zt == input$peak & trough.zt == input$trough)))
      
      expression.gene.set.name <- paste0(paste0("peak at ZT",paste(substr(start = 5,stop=nchar(input$peak),x = input$peak))),
                                         paste0(" and trough at ZT",paste(substr(start = 5,stop=nchar(input$trough),x = input$trough))))
      name.of.sets <- c(input$selected.multiple.tfs, expression.gene.set.name)
    } else if(input$peak == "any" && input$trough != "any")
    {
      list.of.sets[[length(list.of.sets)]] <- rownames(subset(network.data, trough.zt == input$trough))
      
      expression.gene.set.name <- paste0(paste0("trough at ZT",paste(substr(start = 5,stop=nchar(input$trough),x = input$trough))))
      name.of.sets <- c(input$selected.multiple.tfs, expression.gene.set.name)
    } else if(input$peak != "any" && input$trough == "any")
    {
      list.of.sets[[length(list.of.sets)]] <- rownames(subset(network.data, peak.zt == input$peak))
      
      expression.gene.set.name <- paste0(paste0("peak at ZT",paste(substr(start = 5,stop=nchar(input$peak),x = input$peak))))
      name.of.sets <- c(input$selected.multiple.tfs, expression.gene.set.name)
    } else
    {
      name.of.sets <- input$selected.multiple.tfs
    }
    
    names(list.of.sets) <- name.of.sets
    
    ## Compute overlap p value and enrichment
    overlap.output <- supertest(x = list.of.sets, n = nrow(network.data))
    overlap.results <- summary(overlap.output)
    overlap.p.value <- tail(overlap.results$P.value, n=1)
    overlap.p.value <- round(overlap.p.value,digits = -log10(overlap.p.value)+2)
    overlap.enrichment <- round((overlap.results$Table)[["FE"]][nrow(overlap.results$Table)],digits=2)
    
    ## Ouput text for the overlap
    overlap.message <- "The overlap between the targets of the transcription factors"
    text.first.tfs <- paste(input$selected.multiple.tfs[1:(length(input$selected.multiple.tfs)-1)],collapse=",")
    text.all.tfs <- paste(c(text.first.tfs, "and", input$selected.multiple.tfs[length(input$selected.multiple.tfs)]),collapse=" ")
    overlap.message <- paste(overlap.message,text.all.tfs,sep=" ")
    
    if(input$peak != "any" && input$trough != "any")
    {
      overlap.message <- paste(overlap.message, paste(c("and the genes with expression peaks at ZT", substr(x=input$peak, start = 5,stop = nchar(input$peak)), 
      "and expression troughs at", substr(x=input$trough, start = 5,stop = nchar(input$trough))), collapse=" "),sep= " ")
    } else if(input$peak == "any" && input$trough != "any")
    {
      overlap.message <- paste(overlap.message, paste(c("and the genes with expression troughs at", input$trough), collapse=" "), sep = " ")
    } else if(input$peak != "any" && input$trough == "any")
    {
      overlap.message <- paste(overlap.message, paste(c("and the genes with expression peaks at ZT", substr(x=input$peak, start = 5,stop = nchar(input$peak))), collapse=" "), sep = " ")
    }
    
    if(overlap.p.value < 0.05)
    {
      overlap.message <- paste(overlap.message, paste0(paste(c("is significant with a p-value of",overlap.p.value,"and an enrichment of", overlap.enrichment),collapse=" "),"."), sep=" ")
    } else
    {
      overlap.message <- paste(overlap.message, paste0(paste(c("is NOT significant according to a p-value of",overlap.p.value),collapse=" "),"."), sep = " ")
    }
      
    
    output$overlap.significance.text <- renderText(expr = {
      overlap.message
    })
    
    ## Draw venn diagram
    if(number.of.sets == 1)
    {
      output$overlap.significance.text <- renderText(expr = {
        "A single transcription factor or set of genes with specific expression pattern was selected. The analysis of 
        overlap significance is not applicable."
      })
    }
    if(number.of.sets == 2)
    {
      output$venn.diagram.plot <- renderPlot(expr = {
        grid.newpage()
        draw.pairwise.venn(area1 = length(list.of.sets[[1]]),
                           area2 = length(list.of.sets[[2]]),
                           cross.area = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                           category = name.of.sets,
                           scaled = TRUE,lwd = 2,col = "black",
                           fill = c("blue","red"),alpha = 0.7,
                           cex = 2,cat.cex = 1.5,cat.pos = c(0,0))
      },width = 600,height = 600)
    } else if (number.of.sets == 3)
    {
      output$venn.diagram.plot <- renderPlot(expr = {
        grid.newpage()
        draw.triple.venn(area1 = length(list.of.sets[[1]]),
                         area2 = length(list.of.sets[[2]]),
                         area3 = length(list.of.sets[[3]]),
                         n12 = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                         n23 = length(intersect(list.of.sets[[2]],list.of.sets[[3]])),
                         n13 = length(intersect(list.of.sets[[1]],list.of.sets[[3]])), 
                         n123 =  length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[3]])),
                         category = name.of.sets,
                         scaled = TRUE,lwd = 2,col = "black",
                         fill = c("blue","red","darkgreen"),alpha = 0.7,
                         cex = 1.5,cat.cex = 1.5,cat.pos = c(-30,30,180))
      }, width = 600, height = 600)
    } else if (number.of.sets == 4)
    {
      output$venn.diagram.plot <- renderPlot(expr = {
        grid.newpage()
        draw.quad.venn(area1 = length(list.of.sets[[1]]),
                       area2 = length(list.of.sets[[2]]),
                       area3 = length(list.of.sets[[3]]),
                       area4 = length(list.of.sets[[4]]),
                       n12 = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                       n13 = length(intersect(list.of.sets[[1]],list.of.sets[[3]])), 
                       n14 = length(intersect(list.of.sets[[1]],list.of.sets[[4]])),
                       n23 = length(intersect(list.of.sets[[2]],list.of.sets[[3]])),
                       n24 = length(intersect(list.of.sets[[2]],list.of.sets[[4]])),
                       n34 = length(intersect(list.of.sets[[3]],list.of.sets[[4]])),
                       n123 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[3]])), 
                       n124 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[4]])),
                       n134 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[3]]),list.of.sets[[4]])),
                       n234 = length(intersect(intersect(list.of.sets[[2]],list.of.sets[[3]]),list.of.sets[[4]])),
                       n1234 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]), intersect(list.of.sets[[3]],list.of.sets[[4]]))),
                       category = name.of.sets,
                       scaled = TRUE,lwd = 2,col = "black",
                       fill = c("blue","red","darkgreen","orange"),alpha = 0.7,
                       cex = 1.5,cat.cex = 1.5,cat.pos = c(0,0,0,0))
      }, width = 600, height = 600)
    } else if (number.of.sets > 4)
    {
      output$venn.diagram.plot <- renderPlot(expr = {
            plot(overlap.output, Layout = "landscape")
      })
    }
    
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
    
    ## Output table with gene info
    output$outputTable <- renderDataTable({
      create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
    },escape=FALSE)
    
    ## Generate UI to download table and creating a downlodable table
    output$download_ui_for_table<- renderUI(
      tagList(downloadButton(outputId= "downloadData", "Download Selected Genes"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
    )
    
    output$downloadData<- downloadHandler(
      filename= function() {
        paste0(paste(input$selected.tfs,collapse = "_"), ".tsv")
      },
      content= function(file) {
        write.table(create.downloadable.output.table(input.gene.df=selected.genes.df,alias,tfs.names), 
                    file=file, 
                    sep = "\t", 
                    quote = FALSE,
                    row.names = FALSE)
      })
    
    
    
  })
  
  ##Perform GO terms enrichment analysis when button is clicked
  observeEvent(input$goterm,{
    
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )

    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    ## GO enrichment analysis
    
    ## Set the background to perform the GO terms enrichment analysis depending on the user selection
    if (input$go.background == "allgenome")
    {
      go.universe <- ath.universe
    } else
    {
      go.universe <- network.data$name
    }
    
    
    ## Show element when GO term enrichment analysis starts
    shinyjs::showElement(id = 'loading.div')

    enrich.go <- enrichGO(gene          = selected.genes.df$name,
                          universe      = go.universe,
                          OrgDb         = org.At.tair.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE,
                          keyType = "TAIR")

    ## Hide loading element when GO term enrichment is done
    shinyjs::hideElement(id = 'loading.div')
    
    ## Generate ouput table
    enrich.go.result <- as.data.frame(enrich.go)
    
    if(nrow(enrich.go.result) > 0)
    {
      ## GO term Description P-value Q-value Enrichment (SetRatio, BgRatio) Genes
      go.term.enrichments <- compute.enrichments(gene.ratios = enrich.go.result$GeneRatio,
                                                 bg.ratios = enrich.go.result$BgRatio)
      
      go.result.table <- data.frame(enrich.go.result$ID, enrich.go.result$Description,
                                    enrich.go.result$pvalue, enrich.go.result$qvalue,
                                    go.term.enrichments, 
                                    gsub(pattern = "/",replacement = " ",x = enrich.go.result$geneID),
                                    stringsAsFactors = FALSE)
      
      colnames(go.result.table) <- c("GO ID", "Description", "p-value", "q-value",
                                     "Enrichment (Target Ratio; BG Ration)","Genes")
      
      go.result.table.with.links <- go.result.table
      ## Add links to the genes
      genes.in.go.enrichment <- go.result.table$Genes
      
      ## Add link to genes
      for(i in 1:length(genes.in.go.enrichment))
      {
        go.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(genes.in.go.enrichment[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
      }
      
      ## Add links to GO ids
      go.result.table.with.links[["GO ID"]] <- sapply(X = go.result.table.with.links[["GO ID"]], FUN = go.link)
      
      ## Introductory text for GO enrichment table
      go.table.text <- "The table below summarizes the result of the GO term
      enrichment analysis. Each row represents a GO term significantly enriched in the selected
      gene set. The first column represents the GO term
      identifier. The second column contains a human readable description. For more details on the
      corresponding GO term, click on the identifier in the first column. The third and fourth
      column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
      of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
      n is the number of genes with annotation from the target set, N is the number of genes with
      annotation from the gene universe, m is the number of genes from the target set annotated with the
      corresponding GO term and M is the number of genes from the gene universe annotated with
      the GO term associated with the corresponding row. The enrichment is then computed as
      E = (m/n) / (M/N). Finally, the last column, contains the genes from the target set
      annotated with the GO term represented in the corresponding row."
      
      output$textGOTable <- renderText(expr = go.table.text)
      
      ## Output table with GO enrichment result
      output$output_go_table <- renderDataTable({
        ## Error message for the user
        validate(
          need(input$selected.multiple.tfs, "Please select some transcription factor")
        )
        go.result.table.with.links #go.result.table
      },escape=FALSE,options =list(pageLength = 5)) 
      
      output$download_ui_for_go_table<- renderUI(
        tagList(downloadButton(outputId= "downloadGOData", "Download GO Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
      )
      
      output$downloadGOData<- downloadHandler(
        filename= function() {
          paste0(paste(c("GO_enrichment",input$selected.tfs),collapse = "_"), ".tsv")
        },
        content= function(file) {
          write.table(go.result.table, 
                      file=file, 
                      sep = "\t", 
                      quote = FALSE,
                      row.names = FALSE)
        })
      
      ## Download result
      output$downloadData<- downloadHandler(
        filename= function() {
          paste("godata" , ".csv", sep="")
        },
        content= function(file) {
          write.csv(go.result.table,
                    file,
                    row.names=TRUE
          )
        })
      
      ## Link to REVIGO 
      revigo.data <- paste(revigo.data <- apply(go.result.table[,c("GO ID", "q-value")], 1, paste, collapse = " "), collapse="\n")
      
      url1 <- tags$a("here", href="#", onclick="document.revigoForm.submit();")
      url2 <- tags$form(
        name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank",
        tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText",
                      style="visibility: hidden", revigo.data)
      )
      
      output$revigo<- renderUI(
        tagList("The enriched GO terms above may be redundant. Visualize these results in REViGO in order to remove redundancy. Click", url1, url2)
      )
      
      
      output$barplot_text <- renderText("In the following barplot each bar represents a significantly enriched 
GO term. The length of the bar corresponds to the number of genes in the
                                        target set annotated with the given GO term. The bar color captures the level
                                        of significance from blue, less significant, to red, more significant.")
      
      
      ## GO map
      output$gomap_text <- renderText("The following figure corresponds to a subgraph
                                      induced by most significant GO terms from 
                                      Biological Process subcategory. Enriched terms 
                                      are colored and the color depends on the 
                                      adjusted p-value according to the Benjamini & Hochberg 
                                      method, increasing the p-value from purple to red")
      
      output$gomap <- renderPlot(
        width     = 1040,
        height    = 1000,
        res       = 120,
        expr = {
          ## Error message for the user
          validate(
            need(input$selected.multiple.tfs, "Please select some transcription factor")
          )
          goplot(enrich.go,showCategory = 10)
        })
      
      ## Barplot
       output$bar.plot <- renderPlot(
        width     = 870,
        height    = 600,
        res       = 120,
        expr = {
          ## Error message for the user
          validate(
            need(input$selected.multiple.tfs, "Please select some transcription factor")
          )
          barplot(enrich.go,drop=TRUE,showCategory = 10)
        })
       
       ##EMAP plot
       output$emapplot_text <- renderText("The following figure consists of an enrichment map where nodes represent enriched GO terms. The
        size of a node is proportional to the number of genes annotated with the corresponding GO term in the target set.
The node colors represent the level of significance from less signficant in blue to more significant in red. Edges are drawn
between two nodes when the corresponding GO terms are semantically related.")
       
       output$emap.plot <- renderPlot(
         width     = 870,
         height    = 600,
         res       = 120,
         expr = {
           emapplot(enrich.go)
         })

       ## CNET plot
       output$cnetplot_text <- renderText("The following figure corresponds to a gene-concept network. The beige
nodes represents GO terms and the grey nodes genes. An edge is drawn from a gene to a GO term when the gene is annotated
with the corresponding gene. The size of nodes representing GO terms is proportional to the number of genes annotated
with the corresponding GO term.")

       output$cnet.plot <- renderPlot(
         width     = 870,
         height    = 600,
         res       = 120,
         expr = {
           cnetplot(enrich.go)
         })
    }
  })
  
  ##Perform KEGG pathway enrichment analysis when button is clicked
  observeEvent(input$pathway_button,{
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )
    
    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    ## Set the background to perform pathway enrichment analysis depending on the user selection
    if (input$pathway_background == "allgenome")
    {
      pathway.universe <- ath.universe
    } else
    {
      pathway.universe <- network.data$name
    }
    
    ## Show element when kegg pathway enrichment starts
    shinyjs::showElement(id = 'loading.div.kegg')

    ## Compute KEGG pathway enrichment
    pathway.enrichment <- enrichKEGG(gene = selected.genes.df$name, 
                                     organism = "ath", 
                                     keyType = "kegg",
                                     universe = pathway.universe,
                                     qvalueCutoff = 0.05)
    pathway.enrichment.result <- as.data.frame(pathway.enrichment)
    
    ## Hide loading element when KEGG pathway enrichment is done
    shinyjs::hideElement(id = 'loading.div.kegg')
    
    ## Generate output table
    if(nrow(pathway.enrichment.result) > 0)
    {
      pathways.enrichment <- compute.enrichments(gene.ratios = pathway.enrichment.result$GeneRatio,
                                                 bg.ratios = pathway.enrichment.result$BgRatio)
      
      ## Separate genes with blank spaces
      kegg.enriched.genes <- pathway.enrichment.result$geneID
      for(i in 1:length(kegg.enriched.genes))
      {
        kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
      }
      
      ## Generate data frame with output table
      pathways.result.table <- data.frame(pathway.enrichment.result$ID, pathway.enrichment.result$Description,
                                          pathway.enrichment.result$pvalue, pathway.enrichment.result$qvalue,
                                          pathways.enrichment, 
                                          kegg.enriched.genes,
                                          stringsAsFactors = FALSE)
      
      colnames(pathways.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                           "Enrichment (Target Ratio; BG Ration)","Genes")
      
      kegg.result.table.with.links <- pathways.result.table
      
      ## Add links to genes
      for(i in 1:length(kegg.enriched.genes))
      {
        kegg.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(kegg.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
      }
      
      ## Add links to kegg pathways
      kegg.result.table.with.links[["KEGG ID"]] <- sapply(X=kegg.result.table.with.links[["KEGG ID"]],FUN = kegg.pathway.link)

      output$output_pathway_table <- renderDataTable({
        kegg.result.table.with.links
      },escape=FALSE,options =list(pageLength = 5)) 
      
      
      output$download_ui_for_kegg_table<- renderUI(
        tagList(downloadButton(outputId= "downloadKEGGData", "Download KEGG Pathway Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
      )
      
      output$downloadKEGGData<- downloadHandler(
        filename= function() {
          paste0(paste(c("KEGG_enrichment",input$selected.tfs),collapse = "_"), ".tsv")
        },
        content= function(file) {
          write.table(pathways.result.table, 
                      file=file, 
                      sep = "\t", 
                      quote = FALSE,
                      row.names = FALSE)
        })
      
      ## Download result
      output$downloadData<- downloadHandler(
        filename= function() {
          paste("keggdata" , ".csv", sep="")
        },
        content= function(file) {
          write.csv(pathways.result.table,
                    file,
                    row.names=TRUE
          )
        })
    } else
    {
      output$no_kegg_enrichment <- renderText(expr = "No enriched KEGG pathway was detected in the selected genes.")
    }

    ## Visualization of specific enriched pathways
    genes.pathway <- rep(0, length(pathway.universe))
    names(genes.pathway) <- pathway.universe
    
    genes.pathway[selected.genes.df$name] <- 1
    
    pathways.for.select <- paste(pathways.result.table[["KEGG ID"]], pathways.result.table[["Description"]], sep=" - ")
    
    output$kegg_selectize <- renderUI({
      selectInput(inputId = "kegg_pathway", 
                  label="Choose Pathway for Representation",
                  multiple = FALSE,
                  selected = pathways.for.select[1],
                  choices=pathways.for.select)
    })
    
    ## Enriched pathway image
    output$kegg_image <- renderImage({
      pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
               pathway.id = strsplit(input$kegg_pathway,split=" - ")[[1]][1],
               species = "ath",
               limit = list(gene=max(abs(genes.pathway)), cpd=1),
               gene.idtype ="kegg")
      
      list(src = paste(c(strsplit(input$kegg_pathway,split=" - ")[[1]][1],"pathview","png"), collapse="."),
           contentType="image/png",width=1200,height=900)
    },deleteFile = T)
    

    ## KEGG module enrichment analysis
    modules.enrichment <- enrichMKEGG(gene = selected.genes.df$name, 
                                      universe = pathway.universe, 
                                      organism = "ath", 
                                      keyType = "kegg",
                                      minGSSize = 4)
    
    modules.enrichment.result <- as.data.frame(modules.enrichment)
    if(nrow(modules.enrichment.result) > 0)
    {
      modules.enrichment <- compute.enrichments(gene.ratios = modules.enrichment.result$GeneRatio,
                                                bg.ratios = modules.enrichment.result$BgRatio)
      
      modules.enriched.genes <- modules.enrichment.result$geneID
      for(i in 1:length(modules.enriched.genes))
      {
        modules.enriched.genes[i] <- paste(strsplit(modules.enriched.genes[i],split="/")[[1]],collapse=" ")
      }
      
      modules.result.table <- data.frame(modules.enrichment.result$ID, modules.enrichment.result$Description,
                                         modules.enrichment.result$pvalue, modules.enrichment.result$qvalue,
                                         modules.enrichment, 
                                         modules.enriched.genes,
                                         stringsAsFactors = FALSE)
      
      colnames(modules.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                          "Enrichment (Target Ratio; BG Ration)","Genes")
      
      modules.result.table.with.links <- modules.result.table
      
      ## Add links to genes
      for(i in 1:length(modules.enriched.genes))
      {
        modules.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(modules.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
      }
      
      ## Add links to kegg pathways
      modules.result.table.with.links[["KEGG ID"]] <- sapply(X=modules.result.table.with.links[["KEGG ID"]],FUN = kegg.module.link)
      
      ## Generate output table
      output$output_module_table <- renderDataTable({
        modules.result.table.with.links
      },escape=FALSE,options =list(pageLength = 5)) 
    } else
    {
      output$text_module_kegg <- renderText(expr = "No enriched KEGG module was detected in the selected genes.")
    }
      
    
    

    
  })
  
  ##Perform TFBS enrichment analysis when button is clicked
  observeEvent(input$tfbs_button,{
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )
    
    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]  
    
    if(input$peak != "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, (peak.zt == input$peak & trough.zt == input$trough))
    } else if(input$peak == "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, trough.zt == input$trough)
    } else if(input$peak != "any" && input$trough == "any")
    {
      selected.genes.df <- subset(selected.genes.df, peak.zt == input$peak)
    }
    
    ## Set the background to perform TFBS enrichment analysis depending on the user selection
    if (input$tfbs_background == "allgenome")
    {
      tfbs.universe <- 33323
    } else
    {
      tfbs.universe <- 5778
    }
    
    # ## Show element when TFBS enrichment starts
    # shinyjs::showElement(id = 'loading.div.tfbs')
    
    file.precomputed <- paste0(c("precomputed_",input$up_promoter,"_",
                                 input$down_promoter, "_",
                                 input$score, "_",tfbs.universe,".tsv"),collapse="")
    
    ## Load file with precomputed results (background) and compute m and n
    precomputed.result <- read.table(file=paste0("data/precomputed_results_tfbs/",file.precomputed),header = T)
    m <- colSums(precomputed.result > 0) 
    n <- nrow(precomputed.result) - m
    
    ## Compute selection size (k) and number of ocurrences (x).
    # target.genes <- read.table(file = "peak_ZT0_trough_ZT12.txt",as.is = T)[[1]]
    target.genes <- intersect(rownames(precomputed.result),selected.genes.df$names)
    
    
    k <- length(target.genes)
    x <- colSums(precomputed.result[target.genes,] > 0)
    
    ## Compute p-values for enrichment according to a hypergeometric distribution
    p.values <- vector(mode="numeric", length=length(x))
    names(p.values) <- colnames(precomputed.result)
    
    for(i in 1:length(x))
    {
      p.values[i] <- phyper(q = x[i] - 1, m = m[i], n = n[i], k = k, lower.tail = F)
    }
    
    which(p.values < input$motif_significance)
    p.values[which(p.values < input$motif_significance)]
    
    ## Adjust p-values using Benjamini Hochberg
    q.values <- p.adjust(p = p.values,method = "BH")
    
    which(q.values < input$motif_significance)
    q.values[which(q.values < input$motif_significance)]
    
    ## Compute enrichments
    enrichments <- (x / k) / (m / nrow(precomputed.result))
    
    ## Final motifs, pvalues, qvalues and enrichments
    input <- list(motif_significance = 0.05, enrichment_threshold = 2 )
    
    sig.enrich.motifs <- names(which(q.values < input$motif_significance & enrichments > input$enrichment_threshold))
    sig.enrich.ids <- motif.ids[sig.enrich.motifs]
    final.q.values <- q.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    final.p.values <- p.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    final.enrichments <- enrichments[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    
    ## Determine genes for each motif
    genes.with.motif <- vector(length = length(sig.enrich.motifs))
    for (i in 1:length(sig.enrich.motifs))
    {
      print(i)
      rows.with.motif <- which(precomputed.result[,sig.enrich.motifs[i]] != 0)
      all.genes.with.motif <- rownames(precomputed.result)[rows.with.motif]
      genes.with.motif[i] <- paste(... = intersect(all.genes.with.motif,target.genes), collapse = ",")
    }
    
    ## Motifs logos
    motifs.images <- paste0("motifs_images/",sig.enrich.motifs)
    
    for (i in 1:length(motifs.images))
    {
      motifs.images[i] <- paste0("<img src='",motifs.images[i],".png', align = 'center', width = 100>")
    }
    
    
    ## Store data
    tfbs.result.table <- data.frame(sig.enrich.motifs, sig.enrich.ids, motifs.images, final.p.values, final.q.values, final.enrichments, genes.with.motif, stringsAsFactors = FALSE) 
    colnames(tfbs.result.table) <- c("DNA motifs", "Motif ID", "DNA logo", "P-values", "Q-values", "Enrichments", "Genes")
    
    ## Add links to jaspar motifs
    tfbs.result.table[["Motif ID"]] <- sapply(X=sig.enrich.ids,FUN = tfbs.link)
    
    ## Add links to genes
    for (i in 1:length(genes.with.motif))
    {
      tfbs.result.table$Genes[i] <- paste(sapply(X = strsplit(genes.with.motif[i], split = ",")[[1]],FUN = gene.link.function), collapse = ", ")
    }
    
    tfbs.result.table <- tfbs.result.table[order(final.q.values),]
    ## Output table with TFBS enrichment result
    output$output_tfbs_table <- renderDataTable({
      tfbs.result.table
    },escape = FALSE,options =list(pageLength = 10))
    
    output$download_ui_tfbs_table<- renderUI(
      tagList(downloadButton(outputId= "downloadTFBSData", "Download TFBS Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
    )
    
    output$downloadTFBSData<- downloadHandler(
      filename= function() {
        paste0(paste(c("TFBS_enrichment",selected.tfs.with.zts, input$peak, input$trough),collapse = "_"), ".tsv")
      },
      content= function(file) {
        write.table(tfbs.result.table, 
                    file=file, 
                    sep = "\t", 
                    quote = FALSE,
                    row.names = FALSE)
      })
    
    
    
    
  })

}

# Run the application 
shinyApp(ui = ui, server = server)

