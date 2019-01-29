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
length(vertex.names)

## Scale free property
in.degree <- degree(graph = attractor.graph,mode = "in")

in.degree.distribution <- table(in.degree)

x.coord <- log10(as.numeric(names(in.degree.distribution)))
y.coord <- log10(in.degree.distribution)
lm.r <- lm(y.coord ~ x.coord)
lm.res <- summary(lm.r)[[4]]

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
text(x = 8,y=3000,labels = "y = ax^b",cex = 2,col="darkred",pos = 4)
text(x = 8,y=2700,labels = paste0("a = ",round(10^beta,digits=2)),cex = 1.4,col="darkred",pos = 4)
text(x = 8,y=2500,labels = paste0("b = ",round(alpha,digits=2)),cex = 1.4,col="darkred",pos = 4)


## Average path length
path.hist <- path.length.hist(graph = train.network,directed = TRUE)$res
names(path.hist) <- 1:length(path.hist)
barplot(as.table(path.hist),ylab="Frecuencia",xlab="Longitud de camino",border="cyan",col="blue",space = 0)

text(x = 14 ,y = 550,
     labels = paste("Media =",round(average.path.length(graph = train.network,directed = FALSE),digits = 2)),cex = 1.5)


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
