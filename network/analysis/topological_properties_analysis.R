## This script performs a topological analysis over ATTRACTOR

## Load library and graph
library(igraph)

## Load ATTRACTOR network and extract gene names
atha.graph <- read.graph(file="../attractor.graphml", format = "graphml")
vertex.names <- V(atha.graph)$name
length(vertex.names)

## Initialise vectors to store topological parameters
atha.neighbors <- vector(mode="character",length=length(vertex.names))
atha.indegree <- vector(mode="numeric",length=length(vertex.names))
atha.outdegree <- vector(mode="numeric",length=length(vertex.names))
atha.trans <- vector(mode="numeric",length=length(vertex.names))
atha.close <- vector(mode="numeric",length=length(vertex.names))
atha.between <- vector(mode="numeric",length=length(vertex.names))
atha.eccent <- vector(mode="numeric",length=length(vertex.names))

## Loop to compute the different topological parameters and the names of the nodes
## connect to curren node (regulators)
for (i in 1:vcount(atha.graph))
{
  print(i)
  if (length(neighbors(graph = atha.graph, v=vertex.names[i], mode="in")$name) == 0) 
  {
    atha.neighbors[i] <- NA
  } else
  {
    atha.neighbors[i] <- capture.output(cat(neighbors(graph = atha.graph, v=vertex.names[i], mode="in")$name, sep = ","))
  }
  
  atha.indegree[i] <- degree(graph = atha.graph, v=vertex.names[i],mode = "in")
  atha.outdegree[i] <- degree(graph = atha.graph, v=vertex.names[i],mode = "out")
  atha.trans[i] <- transitivity(graph = atha.graph, type = "local", vids=vertex.names[i])
  atha.close[i] <- closeness(graph = atha.graph, vids=vertex.names[i], normalized = TRUE)
  atha.between[i] <- betweenness(graph = atha.graph, v = vertex.names[i], normalized = TRUE)
  atha.eccent[i] <- eccentricity(graph = atha.graph, v = vertex.names[i])
}

## Generate data frame with the topological parameters
node.topological.parameters <- data.frame(vertex.names, atha.neighbors, atha.indegree, 
                                          atha.outdegree, atha.trans, 
                                          atha.close, atha.between, atha.eccent)
colnames(node.topological.parameters) <- c("names", "regulators", "indegree", "outdegree",
                                           "transitivity", "closeness", "betweeness", "eccentricity")

rownames(node.topological.parameters) <- vertex.names
head(node.topological.parameters)


nrow(node.topological.parameters)


## Add topological parameters to file containing the info of the network used in the web app.
attractor.network <- read.table(file="../attractor_network.tsv",header=TRUE)
head(attractor.network)
nrow(attractor.network)

head(node.topological.parameters[attractor.network$names,])

attractor.network.topological.parameters <- cbind(attractor.network,
                                                  node.topological.parameters[as.vector(attractor.network$names),2:ncol(node.topological.parameters)])

head(attractor.network.topological.parameters)

write.table(attractor.network.topological.parameters, file="../attractor_network_topological_parameters.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)

####
component_distribution(graph = atha.graph, cumulative = TRUE, mul.size = TRUE)
diameter(graph = atha.graph, directed = TRUE, unconnected = TRUE)
average.path.length(graph = atha.graph, directed = TRUE, unconnected = TRUE)

hist(atha.indegree)
hist(atha.outdegree,breaks=200)
hist(atha.outdegree+atha.indegree)
