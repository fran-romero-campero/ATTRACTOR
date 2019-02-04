## This script performs a topological analysis over ATTRACTOR

# Authors: Francisco J. Romero-Campero
#          Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
# 
# Contact: Francisco J. Romero-Campero - fran@us.es

## Transcription factors AGI ids and names
tfs.names <- c("CCA1","LHY", "TOC1", "PRR5", "PRR7", "PRR9", "PHYA","PHYB",
               "CRY2","FHY1","LUX","PIF3","PIF4","PIF5","ELF4","ELF3")

tf.ids <- c("AT2G46830", "AT1G01060", "AT5G61380", "AT5G24470", "AT5G02810", "AT2G46790",
            "AT1G09570", "AT2G18790", "AT1G04400", "AT2G37678", "AT3G46640", "AT1G09530",
            "AT2G43010", "AT3G59060", "AT2G40080", "AT2G25930")

names(tf.ids) <- tfs.names
names(tfs.names) <- tf.ids

## Function to extract the TF names of all the different instances associated with
## a motif
extract.different.maps <- function(motif.maps,tfs.names)
{
  motifs.instances <- list(tfs.names[sort(names(motif.maps[[1]]))])
  motif.size <- length(motifs.instances[[1]])
  
  number.different.motifs <- 2
  
  for(i in 2:length(motif.maps))
  {
    new.motif <- TRUE
    current.motif <- tfs.names[sort(names(motif.maps[[i]]))]
    print(current.motif)
    
    for(j in 1:length(motifs.instances))
    {
      if(sum(motifs.instances[[j]] == current.motif) == motif.size)
      {
        new.motif <- FALSE
        break
      }
    }
    
    if(new.motif)
    {
      motifs.instances[[number.different.motifs]] <- current.motif
      number.different.motifs <- number.different.motifs + 1
    }
  }
  
  return(motifs.instances)
}

## Load library and graph
library(igraph)

## Load ATTRACTOR network and extract gene names
attractor.graph <- read.graph(file="../../attractor.graphml", format = "graphml")

## Load three nodes motifs indeces
motifs.3.ind <- read.table(file="indeces_significant_motifs_3.txt")[[1]]

## Load three nodes motifs occurences in attractor
occurrences.3 <- read.table(file="occurency_subgraph_three_nodes_in_attractor.txt")[[1]]

## Motif number 1
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[1]))
occurrences.3[motifs.3.ind[1] + 1]
# No muy interesante

## Motif number 2
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[2]))
occurrences.3[motifs.3.ind[2] + 1]
# No muy interesante

## Motif number 3
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[3]))
occurrences.3[motifs.3.ind[3] + 1]



## Motif number 4
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[4]))
occurrences.3[motifs.3.ind[4] + 1]
## No muy interesante

## Motif number 5
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[5]))
occurrences.3[motifs.3.ind[5] + 1]
## No muy interesante

## Motif number 6
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[6]))
occurrences.3[motifs.3.ind[6] + 1]
maps.motif.6 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[6]), 
                                         attractor.graph, all.maps=TRUE)[["maps"]]




plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[7]))
occurrences.3[motifs.3.ind[7] + 1]

## Motif number 8
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[8]))
occurrences.3[motifs.3.ind[8] + 1]

maps.motif.8 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[8]), 
                                         attractor.graph, all.maps=TRUE)[["maps"]]

length(maps.motif.8)

maps.motif.8[1:3]

  first.intance <- names(maps.motif.8[[1]])

  tfs.names[c(first.intance[c(1,3)])]
  
  motifs.instances <- list(c(sort(tfs.names[c(first.intance[c(1,3)])]), first.intance[2]))
  number.different.motifs <- 2
  
  
  for(i in 2:length(maps.motif.8))
  {
    print(i)
    new.motif <- TRUE
    current.regulators <- sort(tfs.names[names(maps.motif.8[[i]])[c(1,3)]])
    current.target <- names(maps.motif.8[[i]])[2]
    #print(current.motif)

    for(j in 1:length(motifs.instances))
    {
      if(sum(motifs.instances[[j]][1:2] == current.regulators) == 2)
      {
        if(!(current.target %in% strsplit(motifs.instances[[j]][3],split=",")[[1]]))
        {
          motifs.instances[[j]][3] <- paste(motifs.instances[[j]][3],current.target,sep=",")
        }
        new.motif <- FALSE
        break
      }
    }
    
    if(new.motif)
    {
      motifs.instances[[number.different.motifs]] <- c(current.regulators,current.target)
      number.different.motifs <- number.different.motifs + 1
    }
  }
  
get.first <- function(x)  
{
  return(x[[1]])
}

get.second <- function(x)  
{
  return(x[[2]])
}

get.third <- function(x)  
{
  return(x[[3]])
}

res.motif.instances <- data.frame(sapply(motifs.instances,get.first), sapply(motifs.instances,get.second), sapply(motifs.instances,get.third))
head(res.motif.instances)

write.table(x = res.motif.instances,file = "motifs_instances_feedbackloop_multiple_output.txt",quote = F,sep = "\t",row.names = F,col.names = F)



plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[9]))
occurrences.3[motifs.3.ind[9] + 1]


## Motif number 10
plot.igraph(graph.isocreate(size=3, number=motifs.3.ind[10]))

occurrences.3[motifs.3.ind[10] + 1]

maps.motif.10 <- graph.subisomorphic.lad(graph.isocreate(size=3, number=motifs.3.ind[10]), 
                                                          attractor.graph, all.maps=TRUE)[["maps"]]
names.maps.motifs.10 <- extract.different.maps(motif.maps = maps.motif.10 ,tfs.names = tfs.names)
names.maps.motifs.10[12]
