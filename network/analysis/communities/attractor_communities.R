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
save(lc,file = "link_communities_directed.rda")
plot(lc, type = "graph", layout = layout.fruchterman.reingold)
plot(lc, type = "graph", layout = "spencer.circle")
plot(lc, type = "graph", layout = "spencer.circle", shownodesin = 3)
plot(lc, type = "members")
plot(lc, type = "summary")
plot(lc, type = "dend")


community.matrix <- getCommunityMatrix(lc,nodes=V(attractor.graph)$name)
dim(community.matrix)
colSums(community.matrix)

## Nested 
getAllNestedComm(lc)
plot(lc, type = "graph", clusterids = c(1,6))

cr <- getClusterRelatedness(lc, hcmethod = "ward")

cc <- getCommunityCentrality(lc)
head(sort(cc, decreasing = TRUE))

cm <- getCommunityConnectedness(lc, conn = "modularity")
plot(lc, type = "commsumm", summary = "modularity")
dev.off()
for(i in 1:9)
{
  write(x = getNodesIn(lc, clusterids = i),file = paste0("community_",i,".txt"))
}


lc.ward <- getLinkCommunities(network = attractor, hcmethod = "ward",directed = T)
save(lc.ward,file = "link_communities_directed_ward.rda")

lc.ward <- getLinkCommunities(network = attractor, hcmethod = "complete",directed = T)
save(lc.ward,file = "link_communities_directed_complete.rda")

lc.average <- getLinkCommunities(network = attractor, hcmethod = "average",directed = T)
save(lc.average,file = "link_communities_directed_average.rda")


