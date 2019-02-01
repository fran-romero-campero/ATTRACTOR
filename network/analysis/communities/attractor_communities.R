## This script performs a topological analysis over ATTRACTOR

# Authors: Francisco J. Romero-Campero
#          Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
# 
# Contact: Francisco J. Romero-Campero - fran@us.es


## Load library and graph
#install.packages("linkcomm")
library(linkcomm)

attractor.graph <- read.graph(file="../../attractor.graphml", format = "graphml")

## Convert to edge list
attractor <- as_edgelist(graph = attractor.graph,names = T)
head(attractor)
nrow(attractor)

## Extract link communities and visualize them
lc <- getLinkCommunities(network = attractor, hcmethod = "single",directed = T)
print(lc)
plot(lc, type = "graph", layout = layout.fruchterman.reingold)
plot(lc, type = "graph", layout = "spencer.circle")
plot(lc, type = "graph", layout = "spencer.circle", shownodesin = 3)
plot(lc, type = "members")
plot(lc, type = "summary")
plot(lc, type = "dend")
