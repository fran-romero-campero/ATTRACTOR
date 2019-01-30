## This script performs a topological analysis over ATTRACTOR

# Authors: Francisco J. Romero-Campero
#          Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
# 
# Contact: Francisco J. Romero-Campero - fran@us.es


## Load library and graph
library(igraph)

## Load ATTRACTOR network and extract gene names
attractor.graph <- read.graph(file="../../attractor.graphml", format = "graphml")
vertex.names <- V(attractor.graph)$name
number.nodes <- length(vertex.names)

## Scale free property
in.degree <- degree(graph = attractor.graph,mode = "in")

in.degree.distribution <- table(in.degree)

x.coord <- log10(as.numeric(names(in.degree.distribution)))
y.coord <- log10(in.degree.distribution)
lm.r <- lm(y.coord ~ x.coord)
lm.res <- summary(lm.r)[[4]]

## ATTRACTOR R-squared 0.7988 0.7787

beta <- lm.res[1,1]
alpha <- lm.res[2,1]

x.coord.1 <- seq(from=0,to=2,by=0.01)
y.coord.1 <- beta + alpha*x.coord.1
x.coord.2 <- seq(from=0,to=20,by=0.01)
y.coord.2 <- 10^beta*x.coord.2^alpha

plot(x.coord,y.coord)
lines(x.coord.1,y.coord.1)

hist(in.degree,
     xlab="In Degree",ylab="Frequency",
     main="In Degree Distribution",
     cex.main=2,cex.lab=1.5,border="darkblue",lwd=2,col="lightblue")

lines(x.coord.2,y.coord.2,lwd=3,lty=4,col="darkred")
text(x = 8,y=3000,labels = expression(paste("y=", beta,"x"^alpha)),cex = 2,col="darkred",pos = 4)

text(x = 8,y=2600,labels = expression(paste(alpha," = ")),cex = 1.4,col="darkred",pos = 4)
text(x = 9,y=2650,labels = round(alpha,digits=2),cex = 1.4,col="darkred",pos = 4)
text(x = 8,y=2200,labels = expression(paste(beta, " = ")),cex = 1.4,col="darkred",pos = 4)
text(x = 9,y=2250,labels = round(10^beta,digits=2),cex = 1.4,col="darkred",pos = 4)

## Average path length
path.hist <- path.length.hist(graph = attractor.graph,directed = TRUE)$res
names(path.hist) <- 1:length(path.hist)
barplot(as.table(path.hist),ylab="Frequency",xlab="Path Length",border="cyan",col="blue",space = 0,cex.lab=1.3,lwd=2,main="Minimal Path Length Distribution",cex.main=1.5)

text(x = 3.1 ,y = 40000,
     labels = paste("Average Path Length ",round(average.path.length(graph = attractor.graph,directed = TRUE),digits = 2)),cex = 1.2)

## Scale free property
attractor.degree <- degree(graph = attractor.graph)

attractor.degree.distribution <- table(attractor.degree)

x.coord <- log10(as.numeric(names(attractor.degree.distribution)))
y.coord <- log10(attractor.degree.distribution)
lm.r <- lm(y.coord ~ x.coord)
summary(lm.r)
lm.res <- summary(lm.r)[[4]]

beta <- lm.res[1,1]
alpha <- lm.res[2,1]

x.coord.1 <- seq(from=0,to=20,by=0.01)
y.coord.1 <- beta + alpha*x.coord.1
x.coord.2 <- seq(from=0,to=20,by=0.01)
y.coord.2 <- 10^beta*x.coord.2^alpha

plot(x.coord,y.coord)
lines(x.coord.1,y.coord.1)

hist(attractor.degree,breaks=seq(from=0,to=10000,by=1),
     xlab="Degree",ylab="Frequency",
     main="Degree Distribution",
     cex.main=2,cex.lab=1.5,border="darkblue",lwd=2,col="lightblue",xlim=c(0,20))

lines(x.coord.2,y.coord.2,lwd=3,lty=4,col="darkred")
text(x = 8,y=1500,labels = "y = ax^b",cex = 2,col="darkred",pos = 4)
text(x = 8,y=1300,labels = paste0("a = ",round(10^beta,digits=2)),cex = 1.4,col="darkred",pos = 4)
text(x = 8,y=1100,labels = paste0("b = ",round(alpha,digits=2)),cex = 1.4,col="darkred",pos = 4)

## Generate random networks similar to attractor

number.randomisation <- 1000

r.square <- vector(mode = "numeric", length = number.randomisation)
p.val <- vector(mode = "numeric", length = number.randomisation)

for(j in 1:number.randomisation)
{
  random.network.adjacency <- matrix(nrow=number.nodes,ncol=number.nodes)
  
  out.degree <- degree(graph = attractor.graph,mode = "out")
  out.degree <- out.degree[out.degree != 0]
  
  node.regulators <- sample(x = 1:number.nodes,size = length(out.degree),replace = FALSE)
  
  for(i in 1:length(out.degree))
  {
    random.network.adjacency[node.regulators[i],
                             sample(x = 1:number.nodes,size = out.degree[i],replace=FALSE)] <- 1
  }
  
  random.network <- graph.adjacency(adjmatrix = random.network.adjacency,mode = "directed")
  
  in.degree <- degree(graph = random.network,mode = "in")
  in.degree.distribution <- table(in.degree)
  in.degree.distribution <- in.degree.distribution[names(in.degree.distribution) != 0]
  
  x.coord <- log10(as.numeric(names(in.degree.distribution)))
  y.coord <- log10(in.degree.distribution)
  lm.r <- lm(y.coord ~ x.coord)
  r.square[j] <- summary(lm.r)[[9]]
  p.val[j] <- summary(lm.r)$coefficients[2,4]
  print(p.val[j])
}

sum(r.square > 0.7988)
sum(r.square > 0.7787)
sum(r.square > 0.7)
sum(r.square > 0.6)

fdr.val <- p.adjust(p = p.val, method = "BH")
sum(fdr.val < 0.05)

randomisation.result <- data.frame(r.square,p.val,fdr.val)
colnames(randomisation.result) <- c("Rsquare","pvalues","FDR")

write.table(x = randomisation.result,file = "randomisation_result.tsv",quote = F,sep = "\t",row.names = F)
