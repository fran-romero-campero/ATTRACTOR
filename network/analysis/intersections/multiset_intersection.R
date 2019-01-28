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

# Date: January 2019

# Analysis of intersections among multiple sets is fundamental 
# for in-depth understanding of their complex relationships

# This script performs intersections between the targets of two transcription factors and 
# a set of genes. 

#install.packages("SuperExactTest")
library(SuperExactTest)

##Reading the two sets TF target genes
tf1 <- read.table(file = "../../../web_apps/peak_visualizer/data/targets_in_network/CCA1_ZT02_targets_in_network.txt",
                           header = FALSE, as.is = TRUE)
tf2 <- read.table(file="../../../web_apps/peak_visualizer/data/targets_in_network/LHY_targets_in_network.txt",
                          header = FALSE, as.is = TRUE)

#Reading the group of genes peaking at specific time
genes.peak.zt <- read.table(file = "../../../network/clusters/peak_ZT0.txt",
                             header = FALSE, as.is = TRUE)

#EStablishing the sets to test
sets <- c(cca1.targets, lhy.targets, genes.peak.zt)
names(sets) <- c("cca1", "lhy", "peakZT0")

#Test and visualization of intersections

#vignette("set_html")
results <- supertest(x = sets, n = 6830)
help("plot.msets")
par(mar=c(3,3,3,3))
plot(results, sort.by = "size")
plot(results, Layout = "landscape")

results.table <- summary(results)

############# Exploring the output#####


typeof(results.table)
length(results.table)
names(results.table)
results.table$Barcode
results.table$otab
typeof(results.table$otab)
final.intersection <- results.table$otab[["111"]]
results.table$etab
results.table$P.value
tail(results.table$P.value, n=1)
enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]

intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]

intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]

length(sets)

length.gene.sets <- sapply(X = sets,FUN = length)

#If you want to carry out the test introducing only numbers:
cpsets(x = 43 -1, L = length.gene.sets, n = 6830, lower.tail = FALSE)



###########################################################################################
##---------Function that returns p-value, enrichment and the---------------------------#
##---------set of genes in the intersecion (intersectSets)----------------------------------------#

##Correspondencia entre agi symbols y primary symbol
library(org.At.tair.db)
columns(org.At.tair.db)
my.key <- keys(org.At.tair.db, keytype="ENTREZID")
my.col <- c("SYMBOL", "TAIR")
alias2symbol.table <- select(org.At.tair.db, keys=my.key, columns=my.col, keytype="ENTREZID")
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 

## Correspondencia entre agi symbol y descripcion
network.data <- read.table(file="../../../web_apps/attractor_dev/data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")

description <- network.data$description
names(description) <- network.data$names
description[1:3]

intersectSets <- function(tf1,tf2,set.of.genes, alias,gene.descriptions){
  intersection.data <- list()
  sets <- list(tf1, tf2, set.of.genes)
  #names(sets) <- c("cca1", "lhy", "peakZT0")
  results <- supertest(x = sets, n = 6830)
  results.table <- summary(results)
  p.value <- tail(results.table$P.value, n=1) #Get the last p-value
  enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]
  intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]
  intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment
  
  intersection.genes.agi <- intersection.genes
  intersection.genes.primary.symbol <- alias[intersection.genes]
  names(intersection.genes.primary.symbol) <- NULL
  
  intersection.genes.description <- gene.descriptions[intersection.genes]
  names(intersection.genes.description) <- NULL
  
  intersection.data[[3]] <- data.frame(intersection.genes,intersection.genes.primary.symbol,intersection.genes.description,stringsAsFactors = F)
  
  names(intersection.data) <- c("p-value", "enrichment", "gene.table")
  return(intersection.data)
  
}

intersectSets(tf1 = tf1, tf2 = tf2, set.of.genes = genes.peak.zt, alias=alias,gene.descriptions = description)

######Exploring to obtain information (AGI, description, peak time...) #####
#-----from a list of genes.
result <- intersectSets(tf1, tf2, genes.peak.zt)
intersect.genes <- result[3][[1]]

columns(org.At.tair.db)
my.key <- keys(org.At.tair.db, keytype="ENTREZID")
my.col <- c("SYMBOL", "TAIR")
alias2symbol.table <- select(org.At.tair.db, keys=my.key, columns=my.col, keytype="ENTREZID")
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 


#Read expression data table
expression.data <- read.table(file="../../../web_apps/network_visualizer/data/athaliana_neutral_circadian_genes.txt", 
                              as.is = TRUE, header = TRUE, row.names = NULL)
head(expression.data)


subset(expression.data, genes == intersect.genes)
#AQUI ME QUEDO
#-------------------------------------------------------------------------@
#-------------------------------------------------------------------------@
#-------------------------------------------------------------------------@
#-------------------------------------------------------------------------@

#Initialize table to store information about the genes in the intersection
table.of.genes <- matrix(ncol=5, nrow = length(intersect.genes))
colnames(table.of.genes) <- c("Primary_Symbol", "AGI", "Description", "Peak", "Trough")
#Fill the table
table.of.genes[,"AGI"] <- intersect.genes
table.of.genes[,"Primary_Symbol"] <- alias[intersect.genes]



#####loop to perform all possible intersections####

tf.files <- list.files(path = "../../../web_apps/peak_visualizer/data/targets_in_network/targets_to_intersect/", pattern = "targets")
gene.files <- list.files(path = "../../../web_apps/peak_visualizer/data/clusters/by_peaks", pattern = "peak")

#Initialize matrix to store the results
intersection.table <- matrix(ncol=6)
colnames(intersection.table) <- c("TF1", "TF2", "Cluster", "P value", 
                                  "Enrichment", "Intersection Genes") 

#Initialize vector to add it as row into the matrix
current.intersection <- c()

head(intersection.table)

i <- 1
j <- 5
k <- 5
for (i in 1:length(tf.files))
{
  for (j in 1: length(tf.files))
  {
    for (k in 1:length(gene.files))
    {
      tf1 <- read.table(file=paste0("../../../web_apps/peak_visualizer/data/targets_in_network/",tf.files[i]),
                        header = TRUE, as.is = TRUE)
      tf2 <- read.table(file=paste0("../../../web_apps/peak_visualizer/data/targets_in_network/",tf.files[j]),
                        header = TRUE, as.is = TRUE)
      set.of.genes <- read.table(file=paste0("../../../web_apps/peak_visualizer/data/clusters/by_peaks/",gene.files[k]),
                                 header = FALSE, as.is = TRUE)
      
      if(tf.files[i] != tf.files[j])
      {
        print("TEST")
        result <- intersectSets(tf1,tf2,set.of.genes)
        p.value <- result[1][[1]]
        enrichment <- result[2][[1]]
        intersect.genes <- result[3][[1]]
        
        if (length(intersect.genes) !=0 & p.value < 0.000005)
        {
          print("HIT")
          current.intersection[1] <- strsplit(tf.files[i], split = "_")[[1]][1]
          current.intersection[2] <- strsplit(tf.files[j], split = "_")[[1]][1]
          current.intersection[3] <- strsplit(gene.files[k], split = ".txt")[[1]][1]
          current.intersection[4] <- p.value
          current.intersection[5] <- enrichment
          current.intersection[6] <- paste(intersect.genes, collapse= ",")
          intersection.table <- rbind(intersection.table, current.intersection)
          
        }
      }
        
      
    }
  }
  write.table(intersection.table, file="all_intersections.txt", sep="\t", row.names = FALSE)
}



######------Test of intersection between beds------#######


#Reading the bed files of the transcription factors
peaks1 <- read.table(file = "bed.files/PRR5_1_peaks.narrowPeak")
head(peaks1)

peaks2 <- read.table(file = "bed.files/PRR7_peaks.narrowPeak")
head(peaks2)

peaks3 <- read.table(file = "bed.files/PRR9_1_peaks.narrowPeak")
head(peaks3)

peaks.list <- list(peaks1, peaks2, peaks3)

length.sets <- sapply(X = peaks.list, FUN = nrow)


peaks.set1 <- peaks1
peaks.set2 <- peaks2


i <- 89

#Initialize matrix
intersection <- matrix(ncol = 3, nrow=0 )
current.intersection <- matrix(ncol = 3 )
for (i in 1:nrow(peaks.set1))
{
  #Set the current peak values of set1
  current.chr <- peaks.set1[i,1]
  current.start <- peaks.set1[i,2]
  current.end <- peaks.set1[i,3]
  #Checking if there is intersection between the current peak and any peak of set2
  option1 <- nrow(subset(peaks.set2, V1==current.chr & V2<=current.start & V3>=current.start))
  option2 <- nrow(subset(peaks.set2, V1==current.chr & V2<=current.end & V3>=current.end))
 
  print(i)
  
  if(option1+option2 > 0)
  {
    print("HIT")
    if(option1>0)
    {
      hit.peak2 <- subset(peaks.set2, V1==current.chr & V2<=current.start & V3>=current.start)
      current.intersection[1,1] <- current.chr
      current.intersection[1,2] <- current.start
      current.intersection[1,3] <- hit.peak2$V3
      
    }else
    {
      hit.peak2 <- subset(peaks.set2, V1==current.chr & V2<=current.end & V3>=current.end)
      current.intersection[1,1] <- current.chr
      current.intersection[1,2] <- hit.peak2$V2
      current.intersection[1,3] <- current.end
    }
    
    intersection <- rbind(intersection, current.intersection)
  }
}







  
#intersectBed function
intersectBed <- function(peaks.set1, peaks.set2)
{
  intersection <- matrix(ncol = 3, nrow=0 )
  current.intersection <- matrix(ncol = 3 )
  for (i in 1:nrow(peaks.set1))
  {
    #Set the current peak values of set1
    current.chr <- peaks.set1[i,1]
    current.start <- peaks.set1[i,2]
    current.end <- peaks.set1[i,3]
    #Checking if there is intersection between the current peak and any peak of set2
    option1 <- nrow(subset(peaks.set2, V1==current.chr & V2<=current.start & V3>=current.start))
    option2 <- nrow(subset(peaks.set2, V1==current.chr & V2<=current.end & V3>=current.end))
    
    print(i)
    
    if(option1+option2 > 0)
    {
      print("HIT")
      if(option1>0)
      {
        hit.peak2 <- subset(peaks.set2, V1==current.chr & V2<=current.start & V3>=current.start)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- current.start
        current.intersection[1,3] <- hit.peak2$V3
        
      }else
      {
        hit.peak2 <- subset(peaks.set2, V1==current.chr & V2<=current.end & V3>=current.end)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- hit.peak2$V2
        current.intersection[1,3] <- current.end
      }
      
      intersection <- rbind(intersection, current.intersection)
    }
  }
  return(intersection)
}

# The intersectBed function allow to get the intersections betwwen to bed files. If you want to perform the
# the intersection between three bed files, you can do it in a consecutive manner. 

first <- intersectBed(peaks1, peaks2)
second <- intersectBed(first, peaks3)
second

cpsets(x = 43 -1, L = length.gene.sets, n = 6830, lower.tail = FALSE)

