## This script performs an identification of network motifs in ATTRACTOR

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

## Compute outdegree
out.degree <- degree(graph = attractor.graph,mode = "out")
out.degree <- out.degree[out.degree != 0]

max.outdegree <- max(out.degree)

## Generate random sequence for adding edges
random.seq <- rep(0,number.nodes)
points.to.add.regulation <- sample(x = (max.outdegree+1):number.nodes,
                                   size = length(out.degree),replace = FALSE)

random.seq[points.to.add.regulation] <- out.degree

random.network <- barabasi.game(n = number.nodes,out.seq = random.seq,directed = TRUE)
#attractor.graph

#diameter(random.network)

in.degree <- degree(graph = random.network,mode = "in")
#hist(in.degree)

## Network motif with a single node: autoregulation
autorregulation.in.random <- sum(diag(as.matrix(get.adjacency(random.network))))
#autorregulation.in.random

autorregulation.in.attractor <- sum(diag(as.matrix(get.adjacency(attractor.graph))))
#autorregulation.in.attractor

## Generate randomisation
number.randomisation <- 1000
autorregulation.random.graphs <- vector(length=number.randomisation, mode="numeric")

for (i in 1:number.randomisation)
{
  print(i)
  random.seq <- rep(0,number.nodes)
  points.to.add.regulation <- sample(x = (max.outdegree+1):number.nodes,
                                     size = length(out.degree),replace = FALSE)
  
  random.seq[points.to.add.regulation] <- out.degree
  
  random.network <- barabasi.game(n = number.nodes,out.seq = random.seq,directed = TRUE)
  
  autorregulation.random.graphs[i] <- sum(diag(as.matrix(get.adjacency(random.network))))
}

mean(autorregulation.random.graphs)
sd(autorregulation.random.graphs)

sum(autorregulation.random.graphs > autorregulation.in.attractor)/1000

## Network motif with three nodes

## graph.motifs inputs consists of a network and a subgraph size k
## (only k= 3 or 4 are supported). It outputs the number of occurencies
## of any subgraph of size k in the given network.

occurrency.subgraph.three.nodes.in.attractor <- graph.motifs(attractor.graph, size=3)
occurrency.subgraph.three.nodes.in.attractor
length(occurrency.subgraph.three.nodes.in.attractor)
write(occurrency.subgraph.three.nodes.in.attractor,file="occurency_subgraph_three_nodes_in_attractor.txt",ncolumns = 1)

## Subgraphs of size 3
plot.igraph(graph.isocreate(size=3, number=0))
plot.igraph(graph.isocreate(size=3, number=1))
plot.igraph(graph.isocreate(size=3, number=2))
plot.igraph(graph.isocreate(size=3, number=3))
plot.igraph(graph.isocreate(size=3, number=4))
plot.igraph(graph.isocreate(size=3, number=5))
plot.igraph(graph.isocreate(size=3, number=6))
plot.igraph(graph.isocreate(size=3, number=7))
plot.igraph(graph.isocreate(size=3, number=8))
plot.igraph(graph.isocreate(size=3, number=9))
plot.igraph(graph.isocreate(size=3, number=10))
plot.igraph(graph.isocreate(size=3, number=11))
plot.igraph(graph.isocreate(size=3, number=12))
plot.igraph(graph.isocreate(size=3, number=13))
plot.igraph(graph.isocreate(size=3, number=14))
plot.igraph(graph.isocreate(size=3, number=15))

## Generate randomisation
number.randomisation <- 10

motifs.3.random.graph <- matrix(0,nrow=number.randomisation, ncol=16)
motifs.3.random.graph[1:3,]

for (i in 1:number.randomisation)
{
  print(i)
  random.seq <- rep(0,number.nodes)
  points.to.add.regulation <- sample(x = (max.outdegree+1):number.nodes,
                                     size = length(out.degree),replace = FALSE)
  
  random.seq[points.to.add.regulation] <- out.degree
  
  random.network <- barabasi.game(n = number.nodes,out.seq = random.seq,directed = TRUE)
  
  motifs.3.random.graph[i,] <- graph.motifs(random.network, size=3)
}

write.table(x = motifs.3.random.graph,file = "motifs_three_random_graph.txt",quote = F,row.names = F,col.names = F,sep = "\t")

## Test significance for each motif
estimated.p.values <- vector(mode="numeric", length=16)
for(i in 1:16)
{
  estimated.p.values[i] <- sum(motifs.3.random.graph[,i] > occurrency.subgraph.three.nodes.in.attractor[i])/number.randomisation
}

indeces.significant.motifs.3 <- which(estimated.p.values < 0.001) - 1
write(x = indeces.significant.motifs.3,file = "indeces_significant_motifs_3.txt",ncolumns = 1)


## Motifs of size 4
occurrency.subgraph.four.nodes.in.attractor <- graph.motifs(attractor.graph, size=4)
occurrency.subgraph.four.nodes.in.attractor
write(occurrency.subgraph.four.nodes.in.attractor,file="occurency_subgraph_four_nodes_in_attractor.txt",ncolumns = 1)

## Generate randomisation
number.randomisation <- 100

motifs.4.random.graph <- matrix(0,nrow=number.randomisation, ncol=218)
motifs.4.random.graph[1:3,]

for (i in 1:number.randomisation)
{
  print(paste0("motif 4 ",i))
  random.seq <- rep(0,number.nodes)
  points.to.add.regulation <- sample(x = (max.outdegree+1):number.nodes,
                                     size = length(out.degree),replace = FALSE)
  
  random.seq[points.to.add.regulation] <- out.degree
  
  random.network <- barabasi.game(n = number.nodes,out.seq = random.seq,directed = TRUE)
  
  motifs.4.random.graph[i,] <- graph.motifs(random.network, size=4)
}

write.table(x = motifs.4.random.graph,file = "motifs_four_random_graph.txt",quote = F,row.names = F,col.names = F,sep = "\t")

## Test significance for each motif
estimated.p.values <- vector(mode="numeric", length=16)
for(i in 1:218)
{
  estimated.p.values[i] <- sum(motifs.4.random.graph[,i] > occurrency.subgraph.four.nodes.in.attractor[i])/number.randomisation
}

indeces.significant.motifs.4 <- which(estimated.p.values < 0.001) - 1
write(x = indeces.significant.motifs.4,file = "indeces_significant_motifs_4.txt",ncolumns = 1)
