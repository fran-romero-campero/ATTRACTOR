library(seqinr)
library(Biostrings)

sequences.file <- "background_1000_500_5778.fa"
score.value <- "95%"

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

## Read sequences file
sequences.data <- read.fasta(file = sequences.file,seqtype = "AA")

sequences <- getSequence(sequences.data)
seq.names <- getName(sequences.data)

#multiplicities <- vector(mode = "numeric",length=number.motifs)
multiplicities <- matrix(nrow=length(sequences),ncol=number.motifs)
rownames(multiplicities) <- seq.names
colnames(multiplicities) <- motif.names

head(multiplicities)

#names(multiplicities) <- motif.names

i <- 1

for(i in 1:length(seq.names))
{
  current.seq <- c2s(sequences[[i]])
  current.seq.name <- seq.names[i]
  
  
  #motifs.pwm[[j]]
  #motif.ids[[j]]
  #motif.names[[j]]
  
  for(j in 1:number.motifs)
  {
    multiplicities[current.seq.name, motif.names[[j]]] <- nrow(as.data.frame(matchPWM(pwm = motifs.pwm[[j]],subject = current.seq,min.score = score.value)))
  }
}

write.table(x = multiplicities,file = "precomputed_network.tsv",sep = "\t",quote = F)
