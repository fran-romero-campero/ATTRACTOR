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
                           header = FALSE, as.is = TRUE)[[1]]
tf2 <- read.table(file="../../../web_apps/peak_visualizer/data/targets_in_network/LHY_targets_in_network.txt",
                          header = FALSE, as.is = TRUE)[[1]]

#Reading the group of genes peaking at specific time
genes.peak.zt <- read.table(file = "../../../network/clusters/peak_ZT0.txt",
                             header = FALSE, as.is = TRUE)[[1]]

#Establishing the list of sets to test
sets <- list(tf1, tf2, genes.peak.zt)


#####Test and visualization of intersections#####

#vignette("set_html")
results <- supertest(x = sets, n = 5778)
help("plot.msets")
par(mar=c(3,3,3,3))
plot(results, sort.by = "size")
plot(results, Layout = "landscape")

results.table <- summary(results)

#-- Exploring the output--#


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




###### ---- Function that returns p-value, enrichment and the set of genes in the intersecion (intersectSets)#####

##Get translation between AGI and primary symbol
library(org.At.tair.db)
columns(org.At.tair.db)
my.key <- keys(org.At.tair.db, keytype="ENTREZID")
my.col <- c("SYMBOL", "TAIR")
alias2symbol.table <- select(org.At.tair.db, keys=my.key, columns=my.col, keytype="ENTREZID")
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 


## Get description of each AGI symbol
network.data <- read.table(file="../../../web_apps/attractor_dev/data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")

description <- network.data$description
names(description) <- network.data$names
description[1:3]

#This is the function, called intersectSets
intersectSets <- function(tf1,tf2,set.of.genes, alias,gene.descriptions){
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
  intersection.data[[3]] <- gene.table
  

  intersection.genes.description <- gene.descriptions[intersection.genes]
  names(intersection.genes.description) <- NULL
  

  intersection.data[[3]] <- data.frame(intersection.genes,intersection.genes.primary.symbol,intersection.genes.description,stringsAsFactors = F)
  
  names(intersection.data) <- c("p-value", "enrichment", "gene.table")
  return(intersection.data)
  
}

intersectSets(tf1 = tf1, tf2 = tf2, set.of.genes = genes.peak.zt, alias=alias,gene.descriptions = description)


#####----Loop to perform all possible intersections between two TFs and a cluster-----####

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




### Intersection between nodes (genes) with high topological values and genes peaking at each ZT####

attractor.data <- read.table(file="../../attractor_network_representation.tsv", 
                               sep = "\t", as.is = TRUE, header = TRUE, quote = "")
head(attractor.data)
nrow(attractor.data)
gene.names <- attractor.data$names

threshold <- 0.75 #Here you can change the threshold

indegree.threshold <- quantile(attractor.data$indegree, prob=threshold)
indegree.top <- gene.names[attractor.data$indegree > indegree.threshold]

outdegree.threshold <- quantile(attractor.data$outdegree, prob=threshold)
outdegree.top <- gene.names[attractor.data$outdegree > outdegree.threshold]

attractor.degree <- attractor.data$indegree + attractor.data$outdegree
degree.threshold <- quantile(attractor.degree, prob=threshold)
degree.top <- gene.names[attractor.degree > degree.threshold]

attractor.data$transitivity[is.na(attractor.data$transitivity)] <- 0
trans.threshold <- quantile(attractor.data$transitivity, prob=threshold)
trans.top <- gene.names[attractor.data$trans > trans.threshold]

closeness.threshold <- quantile(attractor.data$closeness, prob=threshold)
closeness.top <- gene.names[attractor.data$closeness > closeness.threshold]

betweeness.threshold <- quantile(attractor.data$betweeness, prob=threshold)
betweeness.top <- gene.names[attractor.data$betweeness > betweeness.threshold]

eccentricity.threshold <- quantile(attractor.data$eccentricity, prob=threshold)
eccentricity.top <- gene.names[attractor.data$eccentricity > eccentricity.threshold]


##--Function to perform an intersection of TWO sets--##
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

intersect2sets(set1=degree.top, set2 = genes.peak.zt, alias=alias, gene.descriptions = description)

#####---Several if loops to selectize (in the app) the topological parameter and the set of genes---#####

input <- list(peak="ZT0", trough="Any", topological_parameter="Degree", threshold="0.90")

if (topological_parameter == "Degree")
{
  attractor.degree <- attractor.network$indegree + attractor.network$outdegree
  degree.threshold <- quantile(attractor.degree, prob=input$threshold)
  top.genes <- gene.names[attractor.degree > degree.threshold]
} else if (topological.parameter == "Transitivity")
{
  attractor.trans <- attractor.network$transitivity
  trans.threshold <- quantile(attractor.trans, prob= input$threshold)
  top.genes <- gene.names[attractor.trans > trans.threshold]
} else if (topological.parameter == "Closeness")
{
  attractor.closeness <- attractor.network$closeness
  closeness.threshold <- quantile(attractor.closeness, prob= input$threshold)
  top.genes <- gene.names[attractor.closeness > closeness.threshold]
} else if (topological.parameter == "Betweeness")
{
  attractor.bet <- attractor.network$betweeness
  bet.threshold <- quantile(attractor.bet, prob= input$threshold)
  top.genes <- gene.names[attractor.bet > bet.threshold]
} else if (topological.parameter == "Eccentricity")
{
  attractor.eccen <- attractor.network$eccentricity
  eccen.threshold <- quantile(attractor.eccen, prob= input$threshold)
  top.genes <- gene.names[attractor.eccen > eccen.threshold]
} 

if (input$peak == "Any")
{
  if (input$trough == "Any") 
  {
    zt.genes <- attractor.network$names
  } else
  {
    
    zt.genes <- subset(attractor.network, trough.zt == paste0("trough", substr(x = input$trough, start = 3, stop = nchar(input$trough))))$names
    
  }

} else 
{
  if (input$trough == "Any")
  {
    peak.selection <- paste0("peak", substr(x = input$peak, start = 3, stop = nchar(input$peak)))
    zt.genes <- subset(attractor.network, peak.zt == peak.selection)$names
  } else 
  {
    trough.selection <- paste0("trough", substr(x = input$trough, start = 3, stop = nchar(input$trough)))
    peak.selection <- paste0("peak", substr(x = input$peak, start = 3, stop = nchar(input$peak)))
    zt.genes <- subset(attractor.network, trough.zt == trough.selection & peak.zt == peak.selection)$names
  }
  
}

result <- intersect2sets(set1 = top.genes, set2 = zt.genes, alias = alias, gene.descriptions = description)
p.value <- result[1][[1]]
enrichment <- result[2][[1]]
intersect.genes <- result[3][[1]]$intersection.genes


##Loop to classify genes according to their peak/trough and then use them for intersections####
top.parameters <- c("Degree","Betweeness", "Closeness", "Eccentricity","Transitivity")
zts <- c("Any",paste("ZT",seq(from=0,to=20,by=4),sep=""))
possible.zts <- expand.grid(zts, zts)

for (i in 1:nrow(possible.zts))
{
  
    if (possible.zts[i, "Var1"] == "Any")
    {
      if (possible.zts[i, "Var2"] == "Any")
      {
        zt.genes <- attractor.network$names
      } else
      {

        zt.genes <- subset(attractor.network, trough.zt == paste0("trough", substr(x = possible.zts[i, "Var2"], start = 3, stop = nchar(as.character(possible.zts[i, "Var2"])))))$names

      }

    } else
    {
      if (possible.zts[i, "Var2"] == "Any")
      {
        peak.selection <- paste0("peak", substr(x = possible.zts[i, "Var1"], start = 3, stop = nchar(as.character(possible.zts[i, "Var1"]))))
        zt.genes <- subset(attractor.network, peak.zt == peak.selection)$names
      } else
      {
        trough.selection <- paste0("trough", substr(x = possible.zts[i, "Var2"], start = 3, stop = nchar(as.character(possible.zts[i, "Var2"]))))
        peak.selection <- paste0("peak", substr(x = possible.zts[i, "Var1"], start = 3, stop = nchar(as.character(possible.zts[i, "Var1"]))))
        zt.genes <- subset(attractor.network, trough.zt == trough.selection & peak.zt == peak.selection)$names
      }
    }
  
  if (length(zt.genes) == 0)
  {
    zt.genes <- NA
  }
  file.name <- paste0("peak_",possible.zts[i, "Var1"], "_trough_", possible.zts[i, "Var2"])
  write.table(zt.genes, file =paste0("clusters_ok/", file.name, ".txt"), sep= "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

##Loop to intersect the previous classified genes and high topological values genes####
clusters.files <- list.files(path = "clusters_ok/", pattern = "txt")
top.genes <- list(degree.top, trans.top, closeness.top, betweeness.top, eccentricity.top)
names(top.genes) <- c("Degree", "Transitivity", "Closeness", "Betweeness", "Eccentricity")
#Initialize matrix to store the results
intersection.table <- matrix(ncol=6, nrow = length(clusters.files))
colnames(intersection.table) <- c("peak", "through", "p-value", "fdr", "enrichment", "Intersection Genes") 
head(intersection.table)
i <- 1
j <- 2    
for (i in 1:length(top.parameters))
{
  for (j in 1:length(clusters.files))
  {
    
      current.top <- top.genes[i][[1]]
      set.of.genes <- read.table(file=paste0("clusters_ok/",clusters.files[j]),
                                 header = FALSE, as.is = TRUE)[[1]]
      
      
      print("TEST")
      result <- intersect2sets(set1 = current.top, set2 = set.of.genes, alias = alias, gene.descriptions = description)
      p.value <- result[1][[1]]
      enrichment <- result[2][[1]]
      intersect.genes <- result[3][[1]]$intersection.genes
        
      circadian.info <- strsplit(clusters.files[j], split = "peak_")[[1]][2]
      trough.info <- strsplit(circadian.info, split = "_")[[1]][3]
          
      intersection.table[j,1]<- strsplit(circadian.info, split = "_")[[1]][1]
      intersection.table[j,2] <- strsplit(trough.info, split = ".txt")[[1]][1]
      intersection.table[j,3] <- p.value
      intersection.table[j,5] <- enrichment
      intersection.table[j,6] <- paste(intersect.genes, collapse= ",")
      
        
  }
  fdr.values <- p.adjust(intersection.table[,3], method = "BH")
  intersection.table[,4] <- fdr.values
  write.table(intersection.table, 
              file=paste0("topvalues_clusters_OK/intersections_", names(top.genes[i]), as.character(threshold),".txt"), 
              sep="\t", row.names = FALSE, quote = FALSE)
}


# degree.intersections <- read.table(file="topvalues_clusters/intersections_Degree0.9.txt", header = TRUE, sep = "\t")
# head(degree.intersections)
# degree.intersections$fdr
# degree.intersections$p.value
# 
# degree.intersections <- read.table(file="topvalues_clusters/intersections_Degree0.7.txt", header = TRUE, sep = "\t")
# head(degree.intersections)
# degree.intersections$fdr
# degree.intersections$p.value




#####Intersections between binding regions in DNA (BED files)####
#Reading the bed files of the transcription factors
peaks1 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/CCA1_ZT02_peaks.narrowPeak")
head(peaks1)

peaks2 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/CCA1_peaks.narrowPeak")
head(peaks2)

peaks3 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/PRR9_1_peaks.narrowPeak")
head(peaks3)

peaks.list <- list(peaks1, peaks2, peaks3)

length.sets <- sapply(X = peaks.list, FUN = nrow)


peaks.set1 <- peaks1
peaks.set2 <- peaks2

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
        current.intersection[1,3] <- hit.peak2[,3]
        
      }else
      {
        hit.peak2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- hit.peak2[,2]
        current.intersection[1,3] <- current.end
      }
      
      intersection <- rbind(intersection, current.intersection)
    }
  }
  return(intersection)
}

# The intersectBed function allow to get the intersections between two bed files. If you want to perform the
# the intersection between three bed files, you can do it in a consecutive manner. 

first <- intersectBed(peaks1, peaks2)
nrow(first)
second <- intersectBed(first, peaks3)
nrow(second)


## Permutation of peaks.set2 (random.peaks2) and comparing with peaks.set1 ####
chromosomes.length <- read.table(file="../../../web_apps/peak_visualizer/data/bed_files/atha_chr_lengths.txt",as.is=T)[[1]]
number.randomisation <- 20 #100000

random.intersections <- vector(mode = "numeric",length=number.randomisation) #Creating vector
for(j in 1:number.randomisation)
{
  print(j)
  random.peaks2 <- matrix(nrow=nrow(peaks2),ncol=3) #Matriz con 3 columnas, una para el cromosoma, otra para el comienzo y otra para el final de la marca aleatoria.
  for(i in 1:nrow(peaks2))
  {
    current.chr <- peaks2[i,1][[1]] #Chr de la iésima marca real
    current.start <- peaks2[i,2] #Start de la iésima marca real
    current.end <- peaks2[i,3] #End de la iésima marca real
    current.length <- current.end - current.start #Longitud de la iésima marca real
    
    chr.length <- chromosomes.length[current.chr] #Length del actual cromosoma
    #Ahora genero los mismos datos para marcas aleatorias
    random.start <- floor(runif(n = 1,min = 1,max = chr.length))
    random.end <- random.start + current.length
    
    random.peaks2[i,1] <- current.chr
    random.peaks2[i,2] <- random.start
    random.peaks2[i,3] <- random.end
  }
  
  
  random.intersections[j] <- nrow(intersectBed(peaks.set1 = peaks1, peaks.set2 = random.peaks2 )) 
  
  p.value <- sum(random.intersections > nrow(first)) / number.randomisation
  p.value 

}


##Loop to check the intersection of binding regions (bed files) between all the transcription factors together and store the results in a table####
number.randomisation <- 5
bed.files <- list.files(path = "../../../web_apps/peak_visualizer/data/bed_files/", pattern = "peaks.narrowPeak")

combinations <- expand.grid(bed.files, bed.files)
bed.intersections <- matrix(ncol = 5, nrow = nrow(combinations))
colnames(bed.intersections) <- c("TF1", "TF2", "p-value", "fdr", "Genes")

total.randomisation <- number.randomisation*nrow(combinations) #just to see the progress

i <- 5
for (i in 1:nrow(combinations))
{
  peaks2 <- read.table(file = paste0("../../../web_apps/peak_visualizer/data/bed_files/", combinations[i,2]))
  peaks1 <- read.table(file = paste0("../../../web_apps/peak_visualizer/data/bed_files/", combinations[i,1]))
  real.intersection <- intersectBed(peaks.set1 = peaks1, peaks.set2 = peaks2)
  random.intersections <- vector(mode = "numeric",length=number.randomisation) #Creating vector
  for(j in 1:number.randomisation)
  {
    print(paste0(((i-1*number.randomisation)+j)/total.randomisation, " %"))
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
  
  bed.intersections[i,1] <- strsplit(x = as.character(combinations[i,1]), split = "_peaks")[[1]][1]
  bed.intersections[i,2] <- strsplit(x = as.character(combinations[i,2]), split = "_peaks")[[1]][1]
  bed.intersections[i,3] <- p.value

}

write.table(bed.intersections, file = "bed_intersections.txt", sep = "\t")

# cpsets(x = nrow(second), L = length.sets, n = sum(length.sets), 5778, lower.tail = FALSE) ##DUda, cuánto es n (population size)


