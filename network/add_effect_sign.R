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
# Date: February 2019

## Load previous network representation
network.representation <- read.table(file="../web_apps/attractor_dev/data/attractor_network_representation.tsv",header=T,sep = "\t",as.is = T,quote = "")
head(network.representation)

## Extract gene expression profiles
expression.profiles <- as.matrix(network.representation[,c("ZT0","ZT4","ZT8","ZT12","ZT16","ZT20")])
rownames(expression.profiles) <- network.representation$names
head(expression.profiles)

## Extract regulatory matrix
regularoty.matrix <- network.representation[,35:53]
rownames(regularoty.matrix) <- network.representation$names
head(regularoty.matrix)

target.gene <- "AT1G22770"

gene.expression.profile <- expression.profiles[target.gene,]
