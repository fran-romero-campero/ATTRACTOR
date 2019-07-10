

## This link has a nice tutorial on the use of the hypergeometric distribution in
## enrichment analysis

## http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html

## m number of genes in the background WITH the motif
## n number of genes in the background WITHOUT the motif
## k number of genes in the gene selection
## x number of genes WITH the motif in the gene selection

## Input parameters
input <- list(promoter_length=1000,downstream_length=0,score="90",motif_significance=0.05,enrichment_threshold=2)

file.precomputed <- paste0(c("precomputed_",input$promoter_length,"_",
                           input$downstream_length, "_",
                           input$score, "_5778.tsv"),collapse="")

## Load file with precomputed results (background) and compute m and n
precomputed.result <- read.table(file=file.precomputed,header = T)
m <- colSums(precomputed.result > 0) 
n <- nrow(precomputed.result) - m

## Load file with selection, compute its size (k) and number of ocurrences (x).
target.genes <- read.table(file = "peak_ZT0_trough_ZT12.txt",as.is = T)[[1]]
target.genes <- intersect(rownames(precomputed.result),target.genes)


k <- length(target.genes)
x <- colSums(precomputed.result[target.genes,] > 0)

## Compute p-values for enrichment aocording to a hypergeometric distribution
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


## Final motifs
sig.enrich.motifs <- names(which(p.values < input$motif_significance & enrichments > input$enrichment_threshold))


## Graphical representation
library(seqLogo)

number.motifs <- 453

## Load Position Weight Matrices
## Open file connection
con <- file("../../../../web_apps/peak_visualizer/data/jaspar_motifs/pfm_plants_20180911.txt",open = "r")

## Empty list for storing PWM
motifs.pwm <- vector(mode="list",length = number.motifs)
motif.ids <- vector(mode="character",length=number.motifs)
motif.names <- vector(mode="character",length=number.motifs)

## Load 64 PWM
for(j in 1:number.motifs)
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




current.pwm <- motifs.pwm[["HAT2"]]

seqLogo(makePWM(current.pwm),xaxis = F, yaxis = F)




m["GBF3"]
x["GBF3"]

p.values["GBF3"]


## Compute significance based on randomisations


target.motif.multiplicity <- colSums(background[target.genes,])

number.randomisations <- 1000

random.multiplicities <- matrix(nrow=number.randomisations,ncol=ncol(background))
colnames(random.multiplicities) <- colnames(background)

i <- 1

for(i in 1:number.randomisations)
{
  random.selection <- sample(x = 1:nrow(background),size = length(target.genes))
  random.set <- background[random.selection,]

  random.multiplicities[i,] <- colSums(random.set)
}

head(random.multiplicities)

i <- 1
p.values <- vector(mode = "numeric", length = ncol(background)) 

for(i in 1:ncol(background))
{
  p.values[i] <- sum(target.motif.multiplicity[i] < random.multiplicities[,i]) / number.randomisations
}

q.values <- p.adjust(p = p.values,method = "BH")

colnames(background)[q.values < 0.001]



######--Apply hypergeometric distribution to perform enrichment analysis-----######
######--over several set of genes and store data in tables (Pedro) --####

## TFs and clusters intersections
tfs.intersections.data <- read.table(file = "../../intersections/filtered_significant_results_1.tsv",as.is = T, header = TRUE)
head(tfs.intersections.data)

all.sig.enrich.motifs <- vector(mode="numeric", length=nrow(tfs.intersections.data))
for (j in 1:nrow(tfs.intersections.data)) 
{
  print(j)
  target.genes <- strsplit(x = tfs.intersections.data$genes[j], split = ",")[[1]]
  k <- length(target.genes)
  x <- colSums(precomputed.result[target.genes,] > 0)
  ## Compute p-values for enrichment aocording to a hypergeometric distribution
  p.values <- vector(mode="numeric", length=length(x))
  names(p.values) <- colnames(precomputed.result)
  
  for(i in 1:length(x))
  {
    p.values[i] <- phyper(q = x[i] - 1, m = m[i], n = n[i], k = k, lower.tail = F)
  }
  
  # which(p.values < input$motif_significance)
  # p.values[which(p.values < input$motif_significance)]
  
  ## Adjust p-values using Benjamini Hochberg
  q.values <- p.adjust(p = p.values,method = "BH")
  names(q.values) <- names(p.values)
  
  which(q.values < input$motif_significance)
  
  ## Compute enrichments
  enrichments <- (x / k) / (m / nrow(precomputed.result))
  
  
  ## Final motifs
  sig.enrich.motifs <- names(which(q.values < input$motif_significance & enrichments > input$enrichment_threshold))
  
  final.q.values <- q.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
  final.p.values <- p.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
  final.enrichments <- enrichments[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
  
  ##Store data
  
  store.data <- data.frame(sig.enrich.motifs, final.p.values, final.q.values, final.enrichments) 
  write.table(store.data, sep = "\t", row.names = FALSE,
              file = paste0("TFBS_enrichment_over_intersections/TFS_cluster_genes/TFBS_enrichment_of_", tfs.intersections.data$tf1[j], "_",
                            tfs.intersections.data$tf2[j], "_",
                            tfs.intersections.data$peak.zt[j],"_", 
                            tfs.intersections.data$trough.zt[j], "_", "intersection.txt")) 
  all.sig.enrich.motifs[j] <- paste(sig.enrich.motifs, collapse = ",")
  
}
global.data <- data.frame(tfs.intersections.data$tf1, tfs.intersections.data$tf2,
                          tfs.intersections.data$peak.zt, tfs.intersections.data$trough.zt, 
                          all.sig.enrich.motifs)
head(global.data)
colnames(global.data) <- c("TF1", "TF2", "Peak ZT", "Trough ZT", "Significant enriched Motifs")
write.table(global.data, file = "TFBS_enrichment_over_intersections/TFS_cluster_genes/all_results.txt", 
            sep = "\t", row.names = FALSE)

## Bed intersections
bed.intersections.data <- read.table(file = "../../intersections/bed_intersections_filtered.txt", 
           sep = "\t", header = TRUE, as.is=TRUE)

all.sig.enrich.motifs <- vector(mode="numeric", length=nrow(bed.intersections.data))
for (j in 1:nrow(bed.intersections.data)) 
{
  print(j)
  target.genes <- strsplit(x = bed.intersections.data$Genes[j], split = ",")[[1]]
  #SInce there are genes that are not present in our network, i have to filter them:
  target.genes <- intersect(target.genes, rownames(precomputed.result)) 
  if (length(target.genes) > 0)
  {
    k <- length(target.genes)
    x <- colSums(precomputed.result[target.genes,] > 0)
    ## Compute p-values for enrichment aocording to a hypergeometric distribution
    p.values <- vector(mode="numeric", length=length(x))
    names(p.values) <- colnames(precomputed.result)
    
    for(i in 1:length(x))
    {
      p.values[i] <- phyper(q = x[i] - 1, m = m[i], n = n[i], k = k, lower.tail = F)
    }
    
    # which(p.values < input$motif_significance)
    # p.values[which(p.values < input$motif_significance)]
    
    ## Adjust p-values using Benjamini Hochberg
    q.values <- p.adjust(p = p.values,method = "BH")
    names(q.values) <- names(p.values)
    
    which(q.values < input$motif_significance)
    
    ## Compute enrichments
    enrichments <- (x / k) / (m / nrow(precomputed.result))
    
    ## Final motifs
    sig.enrich.motifs <- names(which(q.values < input$motif_significance & enrichments > input$enrichment_threshold))
    
    final.q.values <- q.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    final.p.values <- p.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    final.enrichments <- enrichments[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
    
    ##Store data
    
    store.data <- data.frame(sig.enrich.motifs, final.p.values, final.q.values, final.enrichments) 
    write.table(store.data, sep = "\t", row.names = FALSE,
                file = paste0("TFBS_enrichment_over_intersections/beds/TFBS_enrichment_of_", 
                              bed.intersections.data$TF1[j], "_and_",
                              bed.intersections.data$TF2[j], "_beds_intersection.txt"))   
    
  } else 
  {
    sig.enrich.motifs <- "The genes that are target of this intersection are not present in ATTRACTOR"
  }
  
  all.sig.enrich.motifs[j] <- paste(sig.enrich.motifs, collapse = ",")
  
}

global.data <- data.frame(bed.intersections.data$TF1, bed.intersections.data$TF2,all.sig.enrich.motifs)
head(global.data)
colnames(global.data) <- c("TF1", "TF2", "Significant enriched Motifs")
write.table(global.data, file = "TFBS_enrichment_over_intersections/beds/all_results.txt", 
            sep = "\t", row.names = FALSE)

##Top topological paramers genes and clusters intersections

degree.clusters <- read.table(file = "../../intersections/topvalues_clusters_OK/intersections_Degree0.95.txt", 
                                     sep = "\t", header = TRUE, as.is=TRUE)
degree.clusters <- subset(x = degree.clusters, subset = fdr<0.05 & enrichment>1.5)

all.sig.enrich.motifs <- vector(mode="numeric", length=nrow(degree.clusters))
for (j in 1:nrow(degree.clusters)) 
{
  print(j)
  target.genes <- strsplit(x = degree.clusters$Intersection.Genes[j], split = ",")[[1]]
  k <- length(target.genes)
  x <- colSums(precomputed.result[target.genes,] > 0)
  ## Compute p-values for enrichment aocording to a hypergeometric distribution
  p.values <- vector(mode="numeric", length=length(x))
  names(p.values) <- colnames(precomputed.result)
  
  for(i in 1:length(x))
  {
    p.values[i] <- phyper(q = x[i] - 1, m = m[i], n = n[i], k = k, lower.tail = F)
  }
  
  # which(p.values < input$motif_significance)
  # p.values[which(p.values < input$motif_significance)]
  
  ## Adjust p-values using Benjamini Hochberg
  q.values <- p.adjust(p = p.values,method = "BH")
  names(q.values) <- names(p.values)
  
  which(q.values < input$motif_significance)
  
  ## Compute enrichments
  enrichments <- (x / k) / (m / nrow(precomputed.result))
  
  
  ## Final motifs
  sig.enrich.motifs <- names(which(q.values < input$motif_significance & enrichments > input$enrichment_threshold))
  
  final.q.values <- q.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
  final.p.values <- p.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
  final.enrichments <- enrichments[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
  
  ##Store data
  
  store.data <- data.frame(sig.enrich.motifs, final.p.values, final.q.values, final.enrichments) 
  write.table(store.data, sep = "\t", row.names = FALSE,
              file = paste0("TFBS_enrichment_over_intersections/Degree095_clusters/TFBS_enrichment_of_cluster_peak", 
                            degree.clusters$peak[j], "_trough",
                            degree.clusters$through[j], "_and_degree095_genes_intersection.txt"))
  all.sig.enrich.motifs[j] <- paste(sig.enrich.motifs, collapse = ",")
  
}
global.data <- data.frame(degree.clusters$peak, degree.clusters$through,all.sig.enrich.motifs)
head(global.data)
colnames(global.data) <- c("Peak", "Trough", "Significant enriched Motifs")
write.table(global.data, file = "TFBS_enrichment_over_intersections/Degree095_clusters/all_results.txt", 
            sep = "\t", row.names = FALSE)

