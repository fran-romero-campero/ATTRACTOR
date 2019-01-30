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

library(SuperExactTest)

## Load network data
network.data <- read.table(file="../../../web_apps/attractor_dev/data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")
nrow(network.data)

## Function to extract the significance between three sets
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
  intersection.data <- list()
  intersection.data[[1]] <- p.value
  intersection.data[[2]] <- enrichment
  intersection.data[[3]] <- intersection.genes
  names(intersection.data) <- c("p-value", "enrichment", "genes")
  return(intersection.data)
}

## Loop over TFs and circadian genes
tf.names <- colnames(network.data)[6:21]

number.of.test <- 6*5*(length(tf.names) -1)*length(tf.names)/2

tfs.1 <- vector(mode="character",length=number.of.test)
tfs.2 <- vector(mode="character",length=number.of.test)
peaks <- vector(mode="character",length=number.of.test)
troughs <- vector(mode="character",length=number.of.test)
pvalues <- vector(mode="numeric",length=number.of.test)
enrichments <- vector(mode="numeric",length=number.of.test)
genes.in.intersection <- vector(mode="character",length=number.of.test)

n <- 1

for(i in 1:(length(tf.names)-1))
{
  for(j in (i+1):length(tf.names))
  {
    for(k in seq(from=0,to=20,by=4))
    {
      for(l in seq(from=0,to=20,by=4))
      {
        
        tf.i <- tf.names[i]
        tf.j <- tf.names[j]
        
        if(k != l)
        {
          targets.tf.i <- network.data$names[network.data[,tf.i] == 1]
          targets.tf.j <- network.data$names[network.data[,tf.j] == 1]
          
          circadian.gene.set <- subset(network.data, peak.zt == paste0("peak",k) & trough.zt == paste0("trough",l))$names
          
          intersection.result <- intersectSets(tf1 = targets.tf.i, tf2 = targets.tf.j, set.of.genes = circadian.gene.set)

          tfs.1[n] <- tf.i
          tfs.2[n] <- tf.j
          
          peaks[n] <- paste0("ZT",k)
          troughs[n] <- paste0("ZT",l)
                    
          pvalues[n] <- intersection.result$`p-value`
          enrichments[n] <- intersection.result$enrichment
          genes.in.intersection[n]  <- paste(intersection.result$genes,collapse = ",")
          
          n <- n + 1
        }
      }
    }
  }
}

## Due to multiple testing, it is necessary to compute FDR 
## http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R-Manual/R-Manual22.html
## https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html

fdr.values <- p.adjust(p = pvalues,method = "BH")


significance.intersection <- data.frame(tf1 = tfs.1,
                                        tf2 = tfs.2,
                                        peak.zt = peaks,
                                        trough.zt = troughs,
                                        pvalues = pvalues,
                                        fdr = fdr.values,
                                        enrichment = enrichments,
                                        genes = genes.in.intersection)

head(significance.intersection)

write.table(x = significance.intersection,file = "significance_intersection.tsv",quote = F,sep = "\t",row.names = F)

significance.results <- read.table(file="significance_intersection.tsv",header = T,sep = "\t",as.is=T)


filtered.significant.results.1 <- subset(significance.results, fdr < 0.01 & enrichment > 2)
write.table(x = filtered.significant.results.1,file = "filtered_significant_results_1.tsv",quote = F,sep = "\t",row.names = F)

filtered.significant.results.2 <- subset(significance.results, fdr < 0.001 & enrichment > 10)
write.table(x = filtered.significant.results.2,file = "filtered_significant_results_2.tsv",quote = F,sep = "\t",row.names = F)