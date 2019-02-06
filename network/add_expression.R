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

## Load previous network representation
network.representation <- read.table(file="../web_apps/attractor_dev/data/attractor_network_representation.tsv",header=T,sep = "\t",as.is = T,quote = "")
head(network.representation)

## Load expression profiles
expression.profiles <- read.table(file="athaliana_neutral_mean_expression.txt",header = T,as.is = T)
head(expression.profiles)
gene.names <- expression.profiles$gene

expression.profiles <- expression.profiles[!duplicated(gene.names),]
gene.names <- gene.names[!duplicated(gene.names)]

expression.profiles <- expression.profiles[,2:ncol(expression.profiles)]
rownames(expression.profiles) <- gene.names
head(expression.profiles)

## Extract network gene expression profiles
network.expression.profiles <- expression.profiles[network.representation$names,]
head(network.expression.profiles)

network.representation.with.profiles <- cbind(network.representation,network.expression.profiles)
head(network.representation.with.profiles)

rownames(network.representation.with.profiles) <- NULL
head(network.representation.with.profiles)

write.table(network.representation.with.profiles, file="attractor_network_representation.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
