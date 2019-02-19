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
expression.profiles <- as.matrix(network.representation[,c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20")])
rownames(expression.profiles) <- network.representation$names
head(expression.profiles)

## Extract regulatory matrix
regulatory.matrix <- network.representation[,35:53]
rownames(regulatory.matrix) <- network.representation$names
head(regulatory.matrix)

## Auxiliary function to determine surronding ZTs
zts.to.consider <- function(zt.point)
{
  zts <- c("ZT00","ZT04","ZT08","ZT12","ZT16","ZT20")
  zts.numeric <- seq(from=0,to=20,by=4)
  
  if(zt.point %in% zts)
  {
    return(c(zt.point, zts[which(zts == zt.point) + 1]))
  } else
  {
    current.zt.numeric <- as.numeric(substr(zt.point,start=3,stop=nchar(current.regulator.zt)))
    next.zt <- zts[which(zts.numeric >= current.zt.numeric)[1]]
    previous.zt <- zts[which(zts.numeric >= current.zt.numeric)[1] - 1]
    return(c(previous.zt, next.zt))
  }
}

gene.names <- network.representation$names

for(i in 1:length(gene.names))
{
  target.gene <- gene.names[i]
  gene.expression.profile <- expression.profiles[target.gene,]
  
  if(sum(is.na(gene.expression.profile)) == 0)
  {
    gene.regulators <- regulatory.matrix[target.gene,]
    gene.regulators.names <- names(gene.regulators)
    
    for(j in 1:length(gene.regulators))
    {
      if(gene.regulators[j] != 0)
      {
        current.regulator <- gene.regulators.names[j]
        current.regulator.zt <- strsplit(current.regulator, split="_")[[1]][2]
        zts.to.check <- zts.to.consider(zt.point = current.regulator.zt)
        
        increment <- gene.expression.profile[zts.to.check[2]] - gene.expression.profile[zts.to.check[1]]
        relative.increment <- increment / gene.expression.profile[zts.to.check[1]]
        
        if(relative.increment > 0.1)
        {
          regulatory.matrix[target.gene,current.regulator] <- 1
        } else 
        {
          regulatory.matrix[target.gene,current.regulator] <- -1
        }
      }
    }
  }
}

head(regulatory.matrix)

network.representation[,35:53] <- regulatory.matrix

write.table(network.representation, file="attractor_network_representation.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)

network.representation <- read.table(file="../web_apps/attractor_dev/data/attractor_network_representation.tsv",header=T,sep = "\t",as.is = T,quote = "")
head(network.representation)

tf.name <- "CCA1_ZT02"

percent.tf <- function(tf.name)
{
  freq.act.rep <- table(network.representation[,tf.name])
  number.rep <- freq.act.rep["-1"]
  number.act <- freq.act.rep["1"]
  total.targets <- number.act + number.rep
  return(list(activated=100*number.act/total.targets, repressed=100*number.rep/total.targets))
}

percent.tf(tf.name = "CCA1_ZT02")
percent.tf(tf.name = "CCA1_ZT14")
percent.tf(tf.name = "PHYA_ZT00")
percent.tf(tf.name = "PHYB_ZT00")
percent.tf(tf.name = "PIF5_ZT04")

table(network.representation$CCA1_ZT02)
table(network.representation$CCA1_ZT14)
table(network.representation$TOC1_ZT15)
table(network.representation$PRR5_ZT10)
table(network.representation$PHYA_ZT00)
