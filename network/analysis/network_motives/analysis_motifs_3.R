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

## Load three nodes motifs indeces
motifs.3.ind <- read.table(file="indeces_significant_motifs_3.txt")[[1]]

plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[1]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[2]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[3]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[4]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[5]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[6]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[7]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[8]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[9]))
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[10]))
