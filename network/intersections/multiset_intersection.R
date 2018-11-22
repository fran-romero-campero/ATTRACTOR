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

# Analysis of intersections among multiple sets is fundamental 
# for in-depth understanding of their complex relationships

# This script performs intersections between the targets of two transcription factors and 
# a set of genes. 

#install.packages("SuperExactTest")
library(SuperExactTest)

##Reading the two sets TF target genes
tf1 <- read.table(file = "../../data/targets_in_network/CCA1_ZT02_targets_in_network.txt",
                           header = FALSE, as.is = TRUE)
tf2 <- read.table(file="../../data/targets_in_network/LHY_targets_in_network.txt",
                          header = FALSE, as.is = TRUE)

#Reading the group of genes peaking at specific time
genes.peak.zt <- read.table(file = "../../data/clusters/peak_ZT16.txt",
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

#Si queremos hacer un test introduciendo sólo cifras. 
cpsets(x = 43 -1, L = length.gene.sets, n = 6830, lower.tail = FALSE)



###########################################################################################
##---------Function that returns p-value, enrichment and the---------------------------#
##---------set of genes in the intersecion (intersectSets)----------------------------------------#

intersectSets <- function(tf1,tf2,set.of.genes){
  intersection.data <- list()
  sets <- c(tf1, tf2, set.of.genes)
  #names(sets) <- c("cca1", "lhy", "peakZT0")
  results <- supertest(x = sets, n = 6830)
  results.table <- summary(results)
  p.value <- tail(results.table$P.value, n=1) #Get the last p-value
  enrichment <- (results.table$Table)[["FE"]][nrow(results.table$Table)]
  intersection.genes <- (results.table$Table)[["Elements"]][nrow(results.table$Table)]
  intersection.genes <- strsplit(intersection.genes, split = ", ")[[1]]
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment
  intersection.data[[3]] <- intersection.genes #hay que meter gene.table con info
  names(intersection.data) <- c("p-value", "enrichment", "gene.table")
  return(intersection.data)
  
}

intersectSets(tf1, tf2, genes.peak.zt)

#####loop to perform all possible intersections



for (i in 1:length(tfs))
{
  for (j in 1: length(tfs))
  {
    for (k in 1:length(group.of.genes))
    {
      intersectSets(i,j,k)
    }
  }
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

#Ya he creado una función para extraer las intersecciones entre dos archivos
#bed (intersectBed) Para hacerlo entre tres, primero lo hago entre las
#dos primeras y después con la tercera. 

first <- intersectBed(peaks1, peaks2)
second <- intersectBed(first, peaks3)
