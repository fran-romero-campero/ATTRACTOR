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
attractor.graph

diameter(random.network)

## Network motif with a single node: autoregulation
autorregulation.in.random <- sum(diag(as.matrix(get.adjacency(random.network))))
autorregulation.in.random

autorregulation.in.attractor <- sum(diag(as.matrix(get.adjacency(attractor.graph))))
autorregulation.in.attractor

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

## La función graph.motifs recibe como entrada una red y un tamaño de subgrafo k
## (en la actualidad sólo puede recibir tamaños 3 o 4) y devuelve el número de
## veces que se encuentra cada subgrafo con k nodos en la red. 

occurrency.subgraph.three.nodes.in.attractor <- graph.motifs(attractor.graph, size=3)
occurrency.subgraph.three.nodes.in.attractor
length(occurrency.subgraph.three.nodes.in.attractor)

## La función graph.isocreate genera todos los grafos posibles de un tamaño dado
## empezando a enumerarlos por 0. 
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

## Para ver si cada uno de los 16 posibles subgrafos es o no motivo de red
## realizamos le procedimiento anterior basado en la generación de redes aleatorias.
## En este caso el acumulador es matricial. 
motifs.3.random.graph <- matrix(0,nrow=1000, ncol=16)
motifs.3.random.graph[1:3,]

for (i in 1:1000)
{
  random.seq <- rep(0,number.nodes)
  points.to.add.regulation <- sample(x = (max.outdegree+1):number.nodes,
                                     size = length(out.degree),replace = FALSE)
  
  random.seq[points.to.add.regulation] <- out.degree
  
  random.network <- barabasi.game(n = number.nodes,out.seq = random.seq,directed = TRUE)
  
  motifs.3.random.graph[i,] <- graph.motifs(random.network, size=3)
}
