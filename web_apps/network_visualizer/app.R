# R script for pre-processing 
# Copyright (C) 2018  Francisco J. Romero-Campero, Pedro de los Reyes Rodríguez,
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


library(shiny)
library(org.At.tair.db)
library(igraph)

##Parameters
radius.1 <- 100 #Outer circle radius

## Read graph adjacency matrix
network.data <- read.table(file="../attractor_dev/data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")

# columns(org.At.tair.db)
## Extract gene ids
genes <- sort(network.data$name)

## Load all and circadian genes
my.key <- keys(org.At.tair.db, keytype="ENTREZID")
my.col <- c("SYMBOL", "TAIR")
alias2symbol.table <- select(org.At.tair.db, keys=my.key, columns=my.col, keytype="ENTREZID")
alias2symbol.table <- subset(alias2symbol.table, genes %in% TAIR)
alias <- alias2symbol.table$SYMBOL
names(alias) <- alias2symbol.table$TAIR
alias[is.na(alias)] <- "" 
genes <- paste(names(alias), alias, sep=" - ")

agis <-alias2symbol.table$TAIR
names(agis) <- alias2symbol.table$SYMBOL
agis[is.na(agis)] <- ""

##Functions
#Function for radian conversion
radian.conversion <- function(alpha)
{
  rad <- (alpha*pi/180)
  return(rad)
}

## Draw a circle and plot it
angle <- seq(from=0, to=2*pi, by=0.01)
x.circle.1 <- radius.1*sin(angle)
y.circle.1 <- radius.1*cos(angle)

radius.2 <- radius.1 - radius.1/12
x.circle.2 <- radius.2 * sin(angle)
y.circle.2 <- radius.2 * cos(angle)

## Read graph adjacency matrix
agi.tfs <- c("AT2G46830", "AT1G01060", "AT5G61380", "AT5G24470", "AT5G02810", "AT2G46790","AT1G09570",
             "AT2G18790", "AT1G04400", "AT2G37678", "AT3G46640", "AT1G09530", "AT2G43010", "AT3G59060",
             "AT2G40080", "AT2G25930")
name.tfs <- c("CCA1", "LHY",  "TOC1", "PRR5", "PRR7", "PRR9", "PHYA", "PHYB", "CRY2", "FHY1", "LUX", "PIF3",
              "PIF4", "PIF5", "ELF4", "ELF3")
agi.tfs.zts <- list(c("ZT02","ZT14"),
                    c("ZT02"),c("ZT15"),c("ZT10"),c("ZT12"),c("ZT04"),c("ZT00"),c("ZT00"),
                    c("ZT08"),c("ZT04"),c("ZT10","ZT12"),c("ZT08"),c("ZT04"),c("ZT04"),
                    c("ZT10"),c("ZT00","ZT04"))

tfs.selectize <- c("CCA1", "LHY",  "TOC1", "PRR5", "PRR7", "PRR9", "PHYA", "PHYB", "CRY2", "FHY1", "LUX", "PIF3",
                   "PIF4", "PIF5", "ELF4", "ELF3")

agi.tfs.zts.multiplicity <- sapply(agi.tfs.zts,length)
names(agi.tfs.zts) <- agi.tfs
names(agi.tfs.zts.multiplicity) <- agi.tfs
names(name.tfs) <- agi.tfs

adj.matrix <- as.matrix(network.data[,35:53])
rownames(adj.matrix) <- network.data$names
adj.matrix <- adj.matrix[agi.tfs,]

cca1.tf.targets <- adj.matrix[,"CCA1_ZT02"] + adj.matrix[,"CCA1_ZT14"]
lhy.tf.targets <- adj.matrix[,"LHY_ZT02"]
toc1.tf.targets <- adj.matrix[,"TOC1_ZT15"]
prr5.tf.targets <- adj.matrix[,"PRR5_ZT10"]
prr7.tf.targets <- adj.matrix[,"PRR7_ZT12"]
prr9.tf.targets <- adj.matrix[,"PRR9_ZT04"]
phya.tf.targets <- adj.matrix[,"PHYA_ZT00"]
phyb.tf.targets <- adj.matrix[,"PHYB_ZT00"]
cry2.tf.targets <- adj.matrix[,"CRY2_ZT08"]
fhy1.tf.targets <- adj.matrix[,"FHY1_ZT04"]
lux.tf.targets <- adj.matrix[,"LUX_ZT10"] + adj.matrix[,"LUX_ZT12"]
pif3.tf.targets <- adj.matrix[,"PIF3_ZT08"]
pif4.tf.targets <- adj.matrix[,"PIF4_ZT04"]
pif5.tf.targets <- adj.matrix[,"PIF5_ZT04"]
elf4.tf.targets <- adj.matrix[,"ELF4_ZT10"]
elf3.tf.targets <- adj.matrix[,"ELF3_ZT00"] + adj.matrix[,"ELF3_ZT04"]


adj.matrix <- matrix(c(cca1.tf.targets, lhy.tf.targets, toc1.tf.targets, prr5.tf.targets, prr7.tf.targets,
         prr9.tf.targets, phya.tf.targets, phyb.tf.targets, cry2.tf.targets, fhy1.tf.targets,
         lux.tf.targets, pif3.tf.targets, pif4.tf.targets, pif5.tf.targets, elf4.tf.targets,
         elf3.tf.targets),ncol=16,nrow=16)
colnames(adj.matrix) <- agi.tfs
rownames(adj.matrix) <- agi.tfs

# adj.matrix <- as.matrix(read.table(file = "data/adjacency_matrix_compressed_only_tfs.txt"))
# is.matrix(adj.matrix)
# dim(adj.matrix)
# rownames(adj.matrix) == colnames(adj.matrix)
# adj.global.matrix <- as.matrix(read.table(file = "data/adjacency_matrix_compressed.txt"))

adj.global.matrix <- as.matrix(network.data[,35:53])
rownames(adj.global.matrix) <- network.data$names

#Read expression data table
#expression.data <- network.data[,29:34]
#rownames(expression.data) <- network.data$names

# expression.data <- read.table(file="data/athaliana_neutral_circadian_genes.txt", 
#                               as.is = TRUE, header = TRUE, row.names = NULL)
#head(expression.data)

#Read mean expression data
#mean.expression <- read.table(file="data/athaliana_neutral_mean_expression.txt", header = TRUE)
#atha.genes <- as.vector(mean.expression$gene)
#mean.expression <- as.matrix(mean.expression[,2:ncol(mean.expression)])
#rownames(mean.expression) <- atha.genes

mean.expression <- as.matrix(network.data[,29:34])
rownames(mean.expression) <- network.data$names

#mean.expression <- expression.data


#Generating the network
#tfs.network <- graph.adjacency(adjmatrix = adj.matrix, mode = "directed")

## Set the angle to each transcription factor. The position (angle) in the network depends on the
## time point at which the ChIP-seq data were collected
# tfs.angles <- radian.conversion(c(2*15, 8*15, 8*15, 0, 0, 0, 4*15, 8*15, 10*15, 12*15, 
#                                   16*15, 20*15, 4*15, 10*15, 4*15, 4*15, 2*15, 14*15, 8*15, 10*15, 12*15, 
#                                   4*15, 12*15, 10*15, 15*15))
# length(tfs.angles)

tfs.names <- colnames(adj.global.matrix)
splitted.tfs.names <- strsplit(tfs.names,split="_")
tfs.angles <- vector(mode="numeric",length=length(tfs.names))
tfs.zts <- vector(mode="numeric",length=length(tfs.names))

for(i in 1:length(splitted.tfs.names))
{
  tfs.angles[i] <- radian.conversion(15*as.numeric(substr(x=splitted.tfs.names[i][[1]][2],start = 3,stop=nchar(splitted.tfs.names[i][[1]][2]))))
  tfs.zts[i] <- substr(x=splitted.tfs.names[i][[1]][2],start = 3,stop=nchar(splitted.tfs.names[i][[1]][2]))
}

zt.multiplicity <- table(tfs.zts)

#Set a radius to each TF to avoid the overlap

radius.to.multiply <- vector(mode="numeric",length=length(splitted.tfs.names))
node.labels <- vector(mode="numeric",length=length(splitted.tfs.names))
for(i in 1:length(splitted.tfs.names))
{
  node.labels[i] <- splitted.tfs.names[i][[1]][1]
  current.zt <- substr(x=splitted.tfs.names[i][[1]][2],start=3,stop=nchar(splitted.tfs.names[i][[1]][2]))
  current.multiplicity <- zt.multiplicity[current.zt]
  radius.to.multiply[i] <- (1 - (0.16*current.multiplicity))*radius.1
  zt.multiplicity[current.zt] <- zt.multiplicity[current.zt] - 1
}

# radius.to.multiply <- c(rep((radius.1*0.8),2), radius.1*0.7,radius.1*0.8,radius.1*0.7,
#                         radius.1*0.6, radius.1*0.8,radius.1*0.6,rep((radius.1*0.8),4), 
#                         rep((radius.1*0.7),2), radius.1*0.6, 
#                         radius.1*0.5, radius.1*0.7,radius.1*0.8, 
#                         radius.1*0.5, radius.1*0.6, radius.1*0.7,
#                         radius.1*0.4, radius.1*0.6, radius.1*0.5, radius.1*0.8)
# tfs.x <- (radius.2 - 3) * sin(tfs.angles)
# tfs.y <- (radius.2 - 3) * cos(tfs.angles)

#Set the x.y coordinates for the positions 
tfs.x <- radius.to.multiply * sin(tfs.angles)
tfs.y <- radius.to.multiply * cos(tfs.angles)

#Generatin a positions matrix 
matrix.pos <- matrix(data = c(tfs.x, tfs.y), nrow = length(tfs.x), ncol = 2)

# node.labels <- c("LHY1","CRY2", "PIF3", "PHYA", "PHYB", rep(x = "ELF3", times=7), "FHY1", 
#                  "ELF4", "PIF4", "PRR9", rep(x = "CCA1", times=5), "PIF5", "PRR7", "PRR5", "TOC1")

selected.colors <- c("blue4","blue","deepskyblue","gold","firebrick","gray47")
peak.times <- c("peak20","peak0","peak4","peak8","peak12","peak16")
names(selected.colors) <- peak.times

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Circadian clock network"),
   
   
   
   # Sidebar with inputs
   sidebarLayout(
      sidebarPanel(
        #Selectize input to choose the gene for represent in the network
        selectizeInput(inputId="target.gene", 
                       label="Target gene", 
                       choices=genes, 
                       selected = "AT1G22770 - GI", 
                       multiple = FALSE),
        
        checkboxGroupInput(inputId = "selected.tfs",
                           label = "Select Transcription Factors:",
                           choices = tfs.selectize,#list("LHY1","CRY2","PIF3","PHYA","PHYB","ELF3",
                                          #"FHY1","ELF4","PIF4","PRR9","CCA1","LUX","PIF5",
                                          #"PRR7","PRR5","TOC1"), 
                           inline = TRUE, 
                           width = "100%"),
        
        checkboxInput(inputId = "all",label = "Select All", width = "100%"),
        
        selectInput(inputId = "interactions", 
                    label = "Show interactions between TFs?", 
                    choices = list("Yes", "No"), 
                    selected = "Yes", multiple = FALSE),
        
        actionButton(inputId = "button", label = "GO"),
        
        
        width = 3

      ),
      
      #Plot the generated network
      
      mainPanel(
         plotOutput(outputId = "network"),
         tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
         tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),
         plotOutput(outputId = "expression"),
         
         width = 9, align = "center"
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  observeEvent(eventExpr = input$button, handlerExpr = {
    output$network <- renderPlot({
      target.agi <- strsplit(x = input$target.gene, split = " - ")[[1]][1]
      
#      if (target.agi %in% row.names(adj.global.matrix)) {
        #Plot circle
        par(mar=c(0,0,0,0))
        plot(x.circle.1,y.circle.1, type = "l", lwd=3, axes=FALSE, xlab = "", ylab="",xlim=c(-1.2 * radius.1, 1.2 * radius.1),ylim=c(-1.2 * radius.1, 1.2 * radius.1))
        lines(x.circle.2, y.circle.2, lwd=3)
        x.polygon <- c(sin(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
                       sin(seq(from=-pi, to=0, by=0.01))* radius.1)
        y.polygon <-c(cos(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
                      cos(seq(from=-pi, to=0, by=0.01))*radius.1)
        polygon(x = x.polygon, y = y.polygon, col = "black")
        for (i in 0:5)
        {
          angle.zt <- radian.conversion(alpha = 60*i)
          zt <- 4*i
          current.zt <- paste("ZT", zt,  sep = "")
          text(x = (radius.1 + radius.1/6)*sin(angle.zt), y = (radius.1 + radius.1/6)*cos(angle.zt), labels = current.zt,cex = 1.5)
          lines(x = c(radius.1 * sin(angle.zt), (radius.1 + radius.1/20)* sin(angle.zt)), 
                y = c(radius.1 * cos(angle.zt), (radius.1 + radius.1/20)* cos(angle.zt)), lwd=2)
        }
        
        if (input$all){
          sel.tfs <- c("LHY1","CRY2","PIF3","PHYA","PHYB","ELF3","FHY1","ELF4","PIF4","PRR9","CCA1","LUX","PIF5","PRR7","PRR5","TOC1")
          to.keep <- agis[sel.tfs]
          
        } else {
          selected.tfs.agi <- agis[input$selected.tfs]
          to.keep <- rep(FALSE,ncol(adj.global.matrix))
          for(i in 1:length(input$selected.tfs))
          {
            to.keep <- (to.keep | grepl(input$selected.tfs[i],colnames(adj.global.matrix)))
          }
        }
        ##First, modify the adj matrix to keep only the selected tfs marked in the app
        adj.matrix.to.represent <- adj.global.matrix[selected.tfs.agi,to.keep]
        
        new.row.names <- c()
        updated.adj.matrix.to.represent <- c()
        for(i in 1:nrow(adj.matrix.to.represent))
        {
          current.tf.name <- name.tfs[rownames(adj.matrix.to.represent)[i]]
          current.tf.zts <- agi.tfs.zts[[rownames(adj.matrix.to.represent)[i]]]
          current.tf.name.zt <- paste(current.tf.name,current.tf.zts,sep="_")
          new.row.names <- c(new.row.names,current.tf.name.zt)
          if(length(current.tf.zts) > 1)
          {
            updated.adj.matrix.to.represent <- rbind(
              rbind(updated.adj.matrix.to.represent,
                   adj.matrix.to.represent[i,]),
              adj.matrix.to.represent[i,])
          } else
          {
            updated.adj.matrix.to.represent <-rbind(updated.adj.matrix.to.represent,
                                                    adj.matrix.to.represent[i,])
          }
        }
        
        rownames(updated.adj.matrix.to.represent) <- new.row.names
        updated.adj.matrix.to.represent <- rbind(updated.adj.matrix.to.represent,adj.global.matrix[target.agi,to.keep])
        updated.adj.matrix.to.represent <- cbind(updated.adj.matrix.to.represent,rep(0,nrow(updated.adj.matrix.to.represent)))
        
        #rows.cols.to.keep <- unlist(sapply(to.keep, grep, row.names(adj.matrix)))
        #tf.adj.matrix <- adj.matrix[rows.cols.to.keep,rows.cols.to.keep]
        
        #Modify adj.matrix and matrix.pos to add the target.gene
        gene.peak.str <- subset(network.data, names == target.agi)$peak.zt 
        gene.peak <- as.numeric(substr(x=gene.peak.str,start=5,stop=nchar(gene.peak.str)))
        #gene.row.complete <- adj.global.matrix[target.agi,] 
        #gene.row <- gene.row[rows.cols.to.keep] #remove non selected tfs from the added row
        
        target.color <- selected.colors[paste0("peak",gene.peak)]
        
        if (input$interactions == "Yes")
        {
          #To show the interactions between the TFs too
          #new.matrix <- cbind(tf.adj.matrix, gene.row)
          #new.matrix <- rbind(new.matrix, rep(0,ncol(new.matrix)))
          new.matrix <- t(updated.adj.matrix.to.represent)
        } else {
          #To show only the interactions between the TFs and the selected gene.
          number.tfs <- nrow(updated.adj.matrix.to.represent) - 1
          new.matrix <- updated.adj.matrix.to.represent
          new.matrix[1:number.tfs,] <- 0
          new.matrix <- t(new.matrix)
          
          # rownames.matrix <- row.names(adj.matrix)
          # null.matrix <- matrix(data = rep(x = 0, times=length(rownames.matrix)^2),
          #                       byrow = TRUE, ncol = length(rownames.matrix), nrow = length(rownames.matrix))
          # rownames(null.matrix) <- rownames.matrix
          # colnames(null.matrix) <- rownames.matrix
          # new.matrix <- cbind(null.matrix, gene.row)
          # new.matrix <- rbind(new.matrix, rep(0,ncol(new.matrix)))
          
        }

                
        #Generating the complete network
        tfs.network <- graph.adjacency(adjmatrix = new.matrix, mode = "directed")
        
        #First, modify the angles, radius, and labels positions to keep only 
        #the selected tfs
        tfs.angles <- tfs.angles[to.keep]
        radius.to.multiply <- radius.to.multiply[to.keep]
        node.labels <- node.labels[to.keep]
        
        #Modify the angles, the radius and the positions to add the new node
        new.tfs.angles <- c(tfs.angles, radian.conversion(gene.peak*15))
        new.multiply <- c(radius.to.multiply, radius.1*0.3)
        tfs.x <- new.multiply * sin(new.tfs.angles)
        tfs.y <- new.multiply * cos(new.tfs.angles)
        
        new.matrix.pos <- matrix(data = c(tfs.x, tfs.y), nrow = nrow(new.matrix), ncol = 2)
        new.node.labels <- c(node.labels, alias[target.agi])
        
        #Plot the network
        #par(mar = c(4,4,4,4))
        plot.igraph(tfs.network, layout=new.matrix.pos, add = TRUE, rescale=FALSE, vertex.size=radius.1*13,
                    vertex.color = c(rep(x = "firebrick1",times=nrow(new.matrix)-1),target.color), vertex.label=new.node.labels, edge.arrow.size = 0.6, 
                    edge.arrow.width=1, edge.curved= TRUE, edge.width = 2, vertex.label.dist = 0,
                    vertex.label.cex=1, vertex.label.font=2,vertex.label.color="black")
        
        
      # } else {
      #   #par(mar = c(0,0,0,0))
      #   plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      #   text(x=0.5, y = 0.5, paste("The", target.agi, "gene does not present a circadian expression pattern"), cex = 0.7)
      # }
      # 
       }, height = 650, width = 650)
    
  })
  
  observeEvent(eventExpr = input$button, handlerExpr = {
    output$expression <- renderPlot({
      target.agi <- strsplit(x = input$target.gene, split = " - ")[[1]][1]
      gene.expression <- as.vector(scale(mean.expression[target.agi,]))
      gene.expression <- c(gene.expression, gene.expression[1])
      
      plot(x=seq(from=0,to=24,by=4),gene.expression,
           type="o",lwd=5,cex=1.5,
           ylim=c(-2.5,2),xlim=c(0,24),
           col="darkgreen",axes=FALSE,xlab="",ylab="", 
           main=paste(target.agi, alias[target.agi],sep=" - "))
      
      polygon(x=c(0,12,12,0),y=c(-2,-2,-2.3,-2.3),lwd=2)
      polygon(x=c(12,24,24,12),y=c(-2,-2,-2.3,-2.3),col = "black",lwd=2)
      
      axis(side = 2,at = -2:2,labels = FALSE,lwd=2)
      mtext("Normalized Gene Expression",side = 2,line = 1.3,cex = 1.3,at = 0)
      axis(side = 1,at=seq(from=0,to=24,by=2),line=-1,las=2,labels = paste("ZT",seq(from=0,to=24,by=2),sep=""),lwd=2)
      
    }, height = 600)
    
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

