## This script performs a topological analysis over ATTRACTOR

## Load library and graph
library(igraph)

## Load ATTRACTOR network and extract gene names
atha.graph <- read.graph(file="../attractor.graphml", format = "graphml")
vertex.names <- V(atha.graph)$name
length(vertex.names)


## Initialise vectors to store topological parameters
atha.indegree <- vector(mode="numeric",length=length(vertex.names))
atha.outdegree <- vector(mode="numeric",length=length(vertex.names))
atha.trans <- vector(mode="numeric",length=length(vertex.names))
atha.close <- vector(mode="numeric",length=length(vertex.names))
atha.between <- vector(mode="numeric",length=length(vertex.names))
atha.eccent <- vector(mode="numeric",length=length(vertex.names))

## Initialize matrix to represent the regulators for each gene
atha.neighbors <- vector(mode="character",length=length(vertex.names))

tfs <- c("CCA1","LHY", "TOC1", "PRR5", "PRR7", "PRR9", "PHYA","PHYB",
         "CRY2","FHY1","LUX","PIF3","PIF4","PIF5","ELF4","ELF3")

tfs.agi <- c("AT2G46830", "AT1G01060", "AT5G61380", "AT5G24470", "AT5G02810", "AT2G46790",
             "AT1G09570", "AT2G18790", "AT1G04400", "AT2G37678", "AT3G46640", "AT1G09530",
             "AT2G43010", "AT3G59060", "AT2G40080", "AT2G25930")

names(tfs) <- tfs.agi

regulators <- matrix(0,ncol=length(tfs),nrow=length(vertex.names))
rownames(regulators) <- vertex.names
colnames(regulators) <- tfs


## Loop to compute the different topological parameters and the names of the nodes
## connect to curren node (regulators)
for (i in 1:vcount(atha.graph))
{
  regulators[i,tfs[neighbors(graph = atha.graph, v=vertex.names[i], mode="in")$name]] <- 1

  atha.indegree[i] <- degree(graph = atha.graph, v=vertex.names[i],mode = "in")
  atha.outdegree[i] <- degree(graph = atha.graph, v=vertex.names[i],mode = "out")
  atha.trans[i] <- transitivity(graph = atha.graph, type = "local", vids=vertex.names[i])
  atha.close[i] <- closeness(graph = atha.graph, vids=vertex.names[i], normalized = TRUE)
  atha.between[i] <- betweenness(graph = atha.graph, v = vertex.names[i], normalized = TRUE)
  atha.eccent[i] <- eccentricity(graph = atha.graph, v = vertex.names[i])
}

## Generate data frame with the topological parameters
node.topological.parameters <- data.frame(vertex.names, atha.indegree, 
                                          atha.outdegree, atha.trans, 
                                          atha.close, atha.between, atha.eccent)
colnames(node.topological.parameters) <- c("names", "indegree", "outdegree",
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
                                                  regulators[as.vector(attractor.network$names),],
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

### Añadir columna de troughs en la data.frame con todos los datos para cada nodo. 
attractor.network <- read.table(file="../../network/attractor_network_topological_parameters.tsv",
                                header=TRUE, as.is = T)
head(attractor.network)
peak.troughs <- read.table(file = "../../web_apps/network_visualizer/data/athaliana_neutral_circadian_genes.txt",
                           header = TRUE, as.is = TRUE)
head(peak.troughs)
troughs <- peak.troughs$troughs
names(troughs) <- peak.troughs$genes

troughs.data <- troughs[attractor.network$names]
names(troughs.data) <- NULL



trough.zt <- c()
i <- 10
for (i in 1:length(troughs.data))
{
  trough.zt[i] <- substr(x = troughs.data[i],start = 3,stop = nchar(troughs.data[i]))
  trough.zt[i] <-paste0("trough", trough.zt[i])
}

library(tibble)

new.data <- add_column(attractor.network, trough.zt, .after = "cluster.classification")
colnames(new.data)[4] <- "peak.zt"
head(new.data)

#Rewriting the data.frame with parameters
write.table(new.data, file="../attractor_network_topological_parameters.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)



