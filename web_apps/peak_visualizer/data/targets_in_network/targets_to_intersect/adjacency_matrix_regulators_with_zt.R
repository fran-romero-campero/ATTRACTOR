# R script for pre-processing 
# Copyright (C) 2018  Francisco J. Romero-Campero, Pedro de los Reyes Rodríguez
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

## Check that the set of the targets is equal to the number of genes in 
files <- list.files(path = ".")
files <- files[-1]

targets <- c()
for(i in 1:length(files))
{
  targets <- c(targets, read.table(file=files[i],header=F, as.is=T)[[1]])
}

length(unique(targets))

## Load network data
network.data <- read.table(file="../../../../attractor_dev/data/attractor_network_representation.tsv",sep="\t",header=T,quote = "",as.is=T)
head(network.data)
nrow(network.data)

adjacency.matrix <- matrix(0,ncol=19,nrow=nrow(network.data))
rownames(adjacency.matrix) <- network.data$names
colnames(adjacency.matrix) <- c("CCA1_ZT02", "CCA1_ZT14", "LHY_ZT02","TOC1_ZT15","PRR5_ZT10", "PRR7_ZT12", "PRR9_ZT??",
                                "PHYA_ZT00", "PHYB_ZT00", "CRY2_ZT08", "FHY1_ZT04", "LUX_ZT10", "LUX_ZT12", 
                                "PIF3_ZT08", "PIF4_ZT04", "PIF5_ZT04", "ELF4_ZT10", "ELF3_ZT00", "ELF3_ZT04")

## CCA1 ZT02
cca1.zt02 <- read.table(file="CCA1_ZT02_targets_in_network.txt",as.is=T,header=F)[[1]]
adjacency.matrix[cca1.zt02,"CCA1_ZT02"] <- 1

## CCA1 ZT14
cca1.zt14 <- read.table(file="CCA1_ZT14_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[cca1.zt14,"CCA1_ZT14"] <- 1

## LHY ZT02
lhy.zt02 <- read.table(file="LHY_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[lhy.zt02,"LHY_ZT02"] <- 1

## TOC1_ZT15

## PRR5_ZT10

##PRR7_ZT12

##PRR9_ZT??

## PHYA_ZT00

## PHYB_ZT00

##CRY2_ZT08

##FHY1_ZT04

##LUX_ZT10

## LUX_ZT12

## PIF3_ZT08

## PIF4_ZT04

## PIF5_ZT04

## ELF4_ZT10

## ELF3_ZT00

## ELF3_ZT04