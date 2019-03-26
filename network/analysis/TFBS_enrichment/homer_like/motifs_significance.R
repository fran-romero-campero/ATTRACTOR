target.genes <- read.table(file = "peak_ZT0_trough_ZT12.txt",as.is = T)[[1]]

background <- read.table(file = "precomputed_network_1000_500_5778.tsv",header=T)
head(background)

target.genes <- intersect(rownames(background),target.genes)

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


