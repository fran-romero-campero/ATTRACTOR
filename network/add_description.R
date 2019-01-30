# R script for adding a description column to the data frame representing the network

# Copyright (C) 2019  Francisco J. Romero-Campero, Pedro de los Reyes
# Ana Belén Romero Losada
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public
# License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Authors: Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
#          Francisco J. Romero-Campero
# 
# Contact: Francisco J. Romero-Campero - fran@us.es 
# Date: January 2019


## Load and extract gene description
tair10.description <- read.table(file="Athaliana_167_TAIR10.defline.txt", sep="\t",quote = "",as.is=T)
head(tair10.description)

extract.gene.name <- function(transcript.name)
{
  return(substr(x = transcript.name,start = 1,stop = 9))
}

tair10.description$V1 <- sapply(tair10.description$V1,extract.gene.name)
head(tair10.description)

tair10.description <- tair10.description[!duplicated(tair10.description),]
sum(duplicated(tair10.description))

head(tair10.description)

description <- tair10.description$V3
names(description) <- tair10.description$V1
description[1:3]

## Load previous network description and add new column description
network.df <- read.table(file="attractor_network_topological_parameters.tsv",header=T,as.is=T)

network.df <- data.frame(network.df,description[network.df$names],stringsAsFactors = F)
colnames(network.df)[ncol(network.df)] <- "description"
colnames(network.df)

#Rewriting the data.frame with description
write.table(network.df, file="attractor_network_representation.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
