

## This link has a nice tutorial on the use of the hypergeometric distribution in
## enrichment analysis

## http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html

## m number of genes in the background WITH the motif
## n number of genes in the background WITHOUT the motif
## k number of genes in the gene selection
## x number of genes WITH the motif in the gene selection

## Input parameters
input <- list(promoter_length=500,downstream_length=0,score="100",motif_significance=0.05,enrichment_threshold=2)

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


