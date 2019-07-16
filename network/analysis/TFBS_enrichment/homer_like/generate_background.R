# This script generates a background of promoter genes to perform
# HOMER-like enrichment.

library(seqinr)
library(TxDb.Athaliana.BioMart.plantsmart28)
promoter.length <- 2000
downstream.promoter <- 500

# Extract sequences of the promoters
chr1 <- getSequence(read.fasta(file = "../../../../web_apps/peak_visualizer/data/athaliana_genome/chr1.fa",seqtype = "AA"))[[1]]
chr2 <- getSequence(read.fasta(file = "../../../../web_apps/peak_visualizer/data/athaliana_genome/chr2.fa",seqtype = "AA"))[[1]]
chr3 <- getSequence(read.fasta(file = "../../../../web_apps/peak_visualizer/data/athaliana_genome/chr3.fa",seqtype = "AA"))[[1]]
chr4 <- getSequence(read.fasta(file = "../../../../web_apps/peak_visualizer/data/athaliana_genome/chr4.fa",seqtype = "AA"))[[1]]
chr5 <- getSequence(read.fasta(file = "../../../../web_apps/peak_visualizer/data/athaliana_genome/chr5.fa",seqtype = "AA"))[[1]]


# Extract info of promoters of genes
txdb <- TxDb.Athaliana.BioMart.plantsmart28

network.data <- read.table(file="../../../../web_apps/attractor_dev/data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")
head(network.data)
network.genes <- network.data$names

#genes.data <- subset(genes(txdb,columns=c("tx_id", "tx_name","gene_id")), seqnames %in% 1:5)
genes.data <- subset(genes(txdb,columns=c("tx_id", "tx_name","gene_id")), seqnames %in% 1:5)
genes.tss <- resize(genes.data, width=1, fix='start')
genes.tss <- as.data.frame(genes.tss)
# genes.tss  <- genes.tss[network.genes,]
head(genes.tss)

genes.promoters.coordinates <- genes.tss[,c(1:3,5)]
head(genes.promoters.coordinates)

# This loop set the size of the promoters. It adds 1000 bp upstream
# of the TSS and 500 bp downstream of the TSS
for(i in 1:nrow(genes.promoters.coordinates))
{
  print(i)
  current.strand <- as.character(genes.promoters.coordinates$strand[i])
  
  if(current.strand == "+")
  {
    # genes.promoters.coordinates$start[i] <- genes.promoters.coordinates$start[i] - 1000
    # genes.promoters.coordinates$end[i] <- genes.promoters.coordinates$end[i] + 500
    genes.promoters.coordinates$start[i] <- genes.promoters.coordinates$start[i] - promoter.length
    genes.promoters.coordinates$end[i] <- genes.promoters.coordinates$end[i] + downstream.promoter
    ## Si al restarle las pares de bases del promotor da negativo (porque está al princioio del cromosoma),
    ## lo fijamos a 0. Ejemplo: i <- 8434
    if (genes.promoters.coordinates$start[i] < 0) 
    {
      genes.promoters.coordinates$start[i] <- 0
    }
  } else if(current.strand == "-")
  {
    # genes.promoters.coordinates$start[i] <- genes.promoters.coordinates$start[i] - 500
    # genes.promoters.coordinates$end[i] <- genes.promoters.coordinates$end[i] + 1000
    genes.promoters.coordinates$start[i] <- genes.promoters.coordinates$start[i] - downstream.promoter
    genes.promoters.coordinates$end[i] <- genes.promoters.coordinates$end[i] + promoter.length
    ## Si al sumarle las pares de bases del promotor da más que la longitud del cromosoma,
    ## lo fijamos en la longitud de éste. Esto lo tengo q hacer
    
  } else
  {
    print("no strand!!!!")
  }
}




background.sequences <- vector(mode = "list",length = nrow(genes.promoters.coordinates))

for(i in 1:nrow(genes.promoters.coordinates))
{
  print(i)
  current.chr <- as.numeric(genes.promoters.coordinates$seqnames[i])
  
  if(current.chr == 1)
  {
    background.sequences[[i]] <- chr1[genes.promoters.coordinates$start[i]:genes.promoters.coordinates$end[i]]
  } else if (current.chr == 2)
  {
    background.sequences[[i]] <- chr2[genes.promoters.coordinates$start[i]:genes.promoters.coordinates$end[i]]
  } else if (current.chr == 3)
  {
    background.sequences[[i]] <- chr3[genes.promoters.coordinates$start[i]:genes.promoters.coordinates$end[i]]
  } else if (current.chr == 4)
  {
    background.sequences[[i]] <- chr4[genes.promoters.coordinates$start[i]:genes.promoters.coordinates$end[i]]
  } else if (current.chr == 5)
  {
    background.sequences[[i]] <- chr5[genes.promoters.coordinates$start[i]:genes.promoters.coordinates$end[i]]
  }
}

# Write into a file
write.fasta(sequences = background.sequences,names = rownames(genes.promoters.coordinates),
            file.out = paste(paste(c("background",promoter.length,downstream.promoter,nrow(genes.promoters.coordinates)), collapse="_"),".fa",sep=""))      
