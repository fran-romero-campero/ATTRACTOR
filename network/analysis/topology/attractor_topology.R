## This script performs a topological analysis over ATTRACTOR

# Authors: Francisco J. Romero-Campero
#          Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
# 
# Contact: Francisco J. Romero-Campero - fran@us.es


## Load library and graph
library(igraph)

## Load ATTRACTOR network and extract gene names
atha.graph <- read.graph(file="../../attractor.graphml", format = "graphml")
vertex.names <- V(atha.graph)$name
length(vertex.names)

## Scale free property
in.degree <- degree(graph = atha.graph,mode = "in")

in.degree.distribution <- table(in.degree)

x.coord <- log10(as.numeric(names(in.degree.distribution)))
y.coord <- log10(in.degree.distribution)
lm.r <- lm(y.coord ~ x.coord)
summary(lm.r)

beta <- 4.0432
alpha <- -2.8508

x.coord.1 <- seq(from=0,to=2,by=0.01)
y.coord.1 <- beta + alpha*x.coord.1
y.coord.2 <- beta*x.coord^alpha

plot(x.coord,y.coord)
lines(x.coord.1,y.coord.1)

hist(in.degree)
lines(x.coord,y.coord.2)
