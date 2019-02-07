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
toc1.zt15 <- read.table(file="TOC1_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[toc1.zt15,"TOC1_ZT15"] <- 1

## PRR5_ZT10
prr5.zt10 <- read.table(file="PRR5_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[prr5.zt10,"PRR5_ZT10"] <- 1

## PRR7_ZT12
prr7.zt12 <- read.table(file="PRR7_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[prr7.zt12,"PRR7_ZT12"] <- 1

## PRR9_ZT??
prr9.zt <- read.table(file="PRR9_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[prr9.zt,"PRR9_ZT??"] <- 1

## PHYA_ZT00
phya.zt00 <- read.table(file="PHYA_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[phya.zt00,"PHYA_ZT00"] <- 1

## PHYB_ZT00
phyb.zt00 <- read.table(file="PHYB_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[phyb.zt00,"PHYB_ZT00"] <- 1

## CRY2_ZT08
cry2.zt08 <- read.table(file="CRY2_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[cry2.zt08,"CRY2_ZT08"] <- 1

## FHY1_ZT04
fhy1.zt04 <- read.table(file="FHY1_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[fhy1.zt04,"FHY1_ZT04"] <- 1

## LUX_ZT10
lux.zt10 <- read.table(file="LUX_ZT10_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[lux.zt10,"LUX_ZT10"] <- 1

## LUX_ZT12
lux.zt12 <- read.table(file="LUX_ZT12_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[lux.zt12,"LUX_ZT12"] <- 1

## PIF3_ZT08
pif3.zt08 <- read.table(file="PIF3_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[pif3.zt08,"PIF3_ZT08"] <- 1

## PIF4_ZT04
pif4.zt04 <- read.table(file="PIF4_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[pif4.zt04,"PIF4_ZT04"] <- 1

## PIF5_ZT04
pif5.zt04 <- read.table(file="PIF5_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[pif5.zt04,"PIF5_ZT04"] <- 1

## ELF4_ZT10
elf4.zt10 <- read.table(file="ELF4_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[elf4.zt10,"ELF4_ZT10"] <- 1

## ELF3_ZT00
elf3.zt00 <- read.table(file="ELF3_ZT00_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[elf3.zt00,"ELF3_ZT00"] <- 1

## ELF3_ZT04
elf3.zt04 <- read.table(file="ELF3_ZT04_targets_in_network.txt", as.is=T, header=F)[[1]]
adjacency.matrix[elf3.zt04,"ELF3_ZT04"] <- 1
