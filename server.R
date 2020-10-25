## ATTRACTOR server

#library(ChIPpeakAnno)

server <- function(input, output, session) {
  
  ## Animation in main page with node size
  # rv <- reactiveValues(i = 0)
  # 
  # increase.step <- 2
  # max.steps <- 1000
  # 
  # output$networkAnimation <- renderPlot( {
  #   ggplot(network.data, aes(x.pos,y.pos)) + 
  #     theme(panel.background = element_blank(), 
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank(),
  #           axis.title = element_blank(),
  #           axis.text = element_blank(),
  #           axis.ticks.y = element_blank()) + 
  #     geom_point(fill=node.colors,size=1.65^norm.data[[(rv$i %% 48)+1]],pch=21)
  # },height = 600)
  # 
  # observeEvent(input$run, {
  #   rv$i <- 0
  #   observe({
  #     isolate({
  #       rv$i <- rv$i + increase.step
  #       print(rv$i)
  #     })
  #     
  #     if(rv$i < max.steps) {
  #       invalidateLater(5, session)
  #     }
  #   })
  # })
  
  ## Animation in main page with red gradient
  # rv <- reactiveValues(i = 0)
  # 
  # increase.step <- 2
  # increase.step.sec <- 0.2
  # max.steps <- 1000
  # 
  # output$networkAnimation <- renderPlot( {
  #   ggplot(network.data, aes(x.pos,y.pos)) + 
  #     theme(panel.background = element_blank(), 
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank(),
  #           axis.title = element_blank(),
  #           axis.text = element_blank(),
  #           axis.ticks.y = element_blank()) + 
  #     geom_point(fill=current.red.gradient[ceiling(1.7^norm.data[[(rv$i %% 48)+1]])],size=5,pch=21)
  # },height = 600, width = 600)
  # 
  # output$clockAnimation <- renderPlot({
  #   #Plot circle
  #   par(mar=c(0,0,0,0))
  #   plot(x.circle.1,y.circle.1, type = "l", lwd=3, axes=FALSE, xlab = "", ylab="",xlim=c(-1.2 * radius.1, 1.2 * radius.1),ylim=c(-1.2 * radius.1, 1.2 * radius.1))
  #   lines(x.circle.2, y.circle.2, lwd=3)
  #   x.polygon <- c(sin(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
  #                  sin(seq(from=-pi, to=0, by=0.01))* radius.1)
  #   y.polygon <-c(cos(seq(from=0, to=-pi, by=-0.01)) * radius.2, 
  #                 cos(seq(from=-pi, to=0, by=0.01))*radius.1)
  #   polygon(x = x.polygon, y = y.polygon, col = "black")
  #   for (j in 0:5)
  #   {
  #     angle.zt <- radian.conversion(alpha = 60*j)
  #     zt <- 4*j
  #     current.zt <- paste("ZT", zt,  sep = "")
  #     text(x = (radius.1 + radius.1/6)*sin(angle.zt), y = (radius.1 + radius.1/6)*cos(angle.zt), labels = current.zt,cex = 1.5,font=2)
  #     lines(x = c(radius.1 * sin(angle.zt), (radius.1 + radius.1/20)* sin(angle.zt)), 
  #           y = c(radius.1 * cos(angle.zt), (radius.1 + radius.1/20)* cos(angle.zt)), lwd=2)
  #   }
  #   
  #   radio.flecha <- 80
  #   angle.zt <- radian.conversion(alpha = 8*rv$i)
  #   arrows(x0 = 0, y0 = 0, x1 = sin(angle.zt)*radio.flecha, y1 = cos(angle.zt)*radio.flecha,lwd = 5)
  #   
  # }, height = 300, width = 300)
  # 
  # observeEvent(input$run, {
  #   rv$i <- 0
  #   observe({
  #     isolate({
  #       rv$i <- rv$i + increase.step
  #       print(rv$i)
  #     })
  #   
  #     
  #     if(rv$i < max.steps) {
  #       invalidateLater(5, session)
  #     }
  #   })
  #   
  # })
  
  ## video tutorial
  output$video_tutorial <- renderUI({
    HTML("<iframe width=\"560\" height=\"315\" src=\"https://www.youtube.com/embed/8eJN5zrMZbI\" frameborder=\"0\" allow=\"autoplay; encrypted-media\" allowfullscreen></iframe>")
  })  
  
  ## clock visualizer code
  output$clock <- renderPlot({
    
    ## Error messages for the user
    validate(
      need(input$selected.tfs, "Please select some transcription factor"),
      need(input$target.gene, "Please select a target gene")
    )
    
    ## Extracting agi ID for selected gene
    target.agi <- strsplit(x = input$target.gene, split = " - ")[[1]][1]
    
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
      text(x = (radius.1 + radius.1/6)*sin(angle.zt), y = (radius.1 + radius.1/6)*cos(angle.zt), labels = current.zt,cex = 1.5,font=2)
      lines(x = c(radius.1 * sin(angle.zt), (radius.1 + radius.1/20)* sin(angle.zt)), 
            y = c(radius.1 * cos(angle.zt), (radius.1 + radius.1/20)* cos(angle.zt)), lwd=2)
    }
    
    ## Get agis and alias of selected tfs and extracting data from adj.global.matrix 
    ## for clock visualizer
    selected.tfs.alias <- sapply(X=strsplit(input$selected.tfs,split=" "),FUN = get.first)
    selected.tfs.agi <- agis[selected.tfs.alias]
    selected.tfs.zts.pasted <- vector(mode = "character", length = length(input$selected.tfs))
    for (i in 1:length(input$selected.tfs))
    {
      selected.tfs.zts.pasted[i] <- paste(strsplit(input$selected.tfs[i], split=" ")[[1]], collapse = "_")
    }
    
    
    ## Set transcription factors to keep in the matrix
    to.keep <- rep(FALSE,ncol(adj.global.matrix))
    for(i in 1:length(input$selected.tfs))
    {
      to.keep <- (to.keep | grepl(selected.tfs.zts.pasted[i],colnames(adj.global.matrix)))
    }
    
    ## Update adjacency matrix to create the transcriptional network 
    ## represented in the clock visualizer
    updated.adj.matrix.to.represent <- adj.global.matrix[target.agi,to.keep]
    for (i in 1:length(input$selected.tfs))
    {
      updated.adj.matrix.to.represent <- rbind(rep(0,length(input$selected.tfs)),
                                               updated.adj.matrix.to.represent)
    }
    
    updated.adj.matrix.to.represent <- cbind(updated.adj.matrix.to.represent, rep(0,length(input$selected.tfs)+1))
    colnames(updated.adj.matrix.to.represent) <- c(selected.tfs.zts.pasted, alias[target.agi])
    rownames(updated.adj.matrix.to.represent) <- c(selected.tfs.zts.pasted, alias[target.agi])
    
    ## Modify adj.matrix and matrix.pos to add the target.gene
    gene.peak.str <- subset(network.data, names == target.agi)$peak.zt 
    gene.peak <- as.numeric(substr(x=gene.peak.str,start=5,stop=nchar(gene.peak.str)))
    target.color <- selected.colors[paste0("peak",gene.peak)]
    
    ## Transpose the matrix to create graph adjacency with igraph
    new.matrix <- t(updated.adj.matrix.to.represent)
    
    ## Generating the complete network
    tfs.network <- graph.adjacency(adjmatrix = new.matrix, mode = "directed",weighted = TRUE)
    edge.weights <- E(tfs.network)$weight
    
    ## Edge colors
    if (length(edge.weights) == 0)
    {
      # E(tfs.network)$color[k] <- NULL
      # print("no edge")
    }else 
    {
      for(k in 1:length(edge.weights))
      {
        if(edge.weights[k] == 1)
        {
          E(tfs.network)$color[k] <- activator.color #"darkgreen"
        } else if(edge.weights[k] == -1)
        {
          E(tfs.network)$color[k] <- repressor.color  #"darkred"
        } else if(edge.weights[k] == 2)
        {
          print(1)
          E(tfs.network)$color[k] <- neutral.color
        }
      }
      
    }
    
    
    ## Vertex colors
    node.colors <- vector(mode="character",length=nrow(new.matrix))
    for(k in 1:(length(node.colors)-1))
    {
      if(new.matrix[k,ncol(new.matrix)] == 1)
      {
        node.colors[k] <- activator.color
      } else if(new.matrix[k,ncol(new.matrix)] == -1)
      {
        node.colors[k] <- repressor.color
      } else #if (new.matrix[k,ncol(new.matrix)] == 0)
      {
        node.colors[k] <- neutral.color
      }
    }
    
    node.colors[length(node.colors)] <- target.color
    
    ## Modify the angles, radius, and labels positions to keep only the selected tfs
    tfs.angles <- tfs.angles[to.keep]
    radius.to.multiply <- radius.to.multiply[to.keep]
    node.labels <- node.labels[to.keep]
    
    ## Modify the angles, the radius and the positions to add the new node
    new.tfs.angles <- c(tfs.angles, radian.conversion(gene.peak*15))
    new.multiply <- c(radius.to.multiply, radius.1*0.3)
    tfs.x <- new.multiply * sin(new.tfs.angles)
    tfs.y <- new.multiply * cos(new.tfs.angles)
    
    new.matrix.pos <- matrix(data = c(tfs.x, tfs.y), nrow = nrow(new.matrix), ncol = 2)
    new.node.labels <- c(node.labels, alias[target.agi])
    
    ## Plot the network
    plot.igraph(tfs.network, layout=new.matrix.pos, add = TRUE, rescale=FALSE, vertex.size=radius.1*13,
                vertex.color = node.colors, vertex.label=new.node.labels, edge.arrow.size = 0.8, 
                edge.arrow.width=2, edge.curved= TRUE, edge.width = 4, vertex.label.dist = 0,
                vertex.label.cex=1, vertex.label.font=2,vertex.label.color="black",label.font=2)
    
  })
  
  ## Express visualizer code
  output$expression <- renderPlot({
    
    ## Error message for the user
    validate(
      need(input$selected.tfs, "Please select some transcription factor"),
      need(input$target.gene, "Please select a target gene")
    )
    
    
    ## Get agis and alias of selected tfs and extracting data from adj.global.matrix 
    ## for expression visualizer
    selected.tfs.alias <- sapply(X=strsplit(input$selected.tfs,split=" "),FUN = get.first)
    selected.tfs.agi <- agis[selected.tfs.alias]
    selected.tfs.zts.pasted <- vector(mode = "character", length = length(input$selected.tfs))
    for (i in 1:length(input$selected.tfs))
    {
      selected.tfs.zts.pasted[i] <- paste(strsplit(input$selected.tfs[i], split=" ")[[1]], collapse = "_")
    }
    
    ## Get expression data for the selected gene
    target.agi <- strsplit(x = input$target.gene, split = " - ")[[1]][1]
    gene.expression <- as.vector(scale(mean.expression[target.agi,]))
    gene.expression <- c(gene.expression, gene.expression[1])
    extended.gene.expression <- approx(x = seq(from=0,to=24,by=4), y = gene.expression, xout=c(0,2,4,8,10,12,14,15,16,20,24))
    extended.gene.expression.values <- extended.gene.expression$y
    names(extended.gene.expression.values) <- c("ZT00", "ZT02", "ZT04", "ZT08", "ZT10", "ZT12", "ZT14", "ZT15", "ZT16", "ZT20", "ZT24")
    line.color <- selected.colors[network.data[target.agi, "peak.zt"]]
    
    ## Plot the initial visualization of the expression profile
    plot(x=seq(from=0,to=24,by=4),gene.expression,
         type="o",lwd=5,cex=1.5,
         ylim=c(-2.5,height),xlim=c(0,24),
         col=line.color,axes=FALSE,xlab="",ylab="", 
         main=paste(target.agi, alias[target.agi],sep=" - "))
    
    ## Add TFs to expression profile
    for(i in 1:length(selected.tfs.agi))
    {
      current.tf.name <- names(selected.tfs.agi[i])
      current.zt <- strsplit(x = selected.tfs.zts.pasted[i], split="_")[[1]][2]
      current.time.point <- as.numeric(substr(x = current.zt, start = 3, stop = nchar(selected.tfs.zts.pasted[i])))
      # current.time.point <- as.numeric(substr(x = current.tf.zts[j], start = 3, stop = nchar(current.tf.zts[j])))
      current.regulation <- adj.global.matrix[target.agi,selected.tfs.zts.pasted[i]]
      
      if(current.regulation == 1)
      {
        point.color <- activator.color #"seagreen3"#"darkgreen"
        arrow.angle <- 45
        draw.tf <- TRUE
      } else if (current.regulation == -1)
      {
        point.color <- repressor.color # "firebrick1"
        arrow.angle <- 90
        draw.tf <- TRUE
      } else if (current.regulation == 2)
      {
        point.color <- neutral.color
        arrow.angle <- 45
        draw.tf <- TRUE
      } else if (current.regulation == 0)
      {
        draw.tf <- FALSE
      }
      
      if(draw.tf)
      {
        arrows(x0 = current.time.point, y0 = height.to.multiply[selected.tfs.zts.pasted[i]],
               x1 = current.time.point ,y1= extended.gene.expression.values[current.zt] + 0.3,lwd=4,angle=arrow.angle,length=0.05,col=point.color)
        points(x = current.time.point,y=height.to.multiply[selected.tfs.zts.pasted[i]],lwd=4,cex=4, col=point.color, pch = 19)  
        text(x = current.time.point,y=height.to.multiply[selected.tfs.zts.pasted[i]],labels = current.tf.name, cex=1.2, font=2 )
      }
    }
    
    polygon(x=c(0,12,12,0),y=c(-2,-2,-2.3,-2.3),lwd=2)
    polygon(x=c(12,24,24,12),y=c(-2,-2,-2.3,-2.3),col = "black",lwd=2)
    
    axis(side = 2,at = -2:2,labels = FALSE,lwd=2)
    mtext("Normalized Gene Expression",side = 2,line = 1.3,cex = 1.5,at = 0)
    axis(side = 1,at=seq(from=0,to=24,by=2),line=-1,las=2,labels = paste("ZT",seq(from=0,to=24,by=2),sep=""),lwd=2)
    
  })
  
  ## Peak visualizer code
  output$peak_plot <- renderPlot({
    
    ## Sanity checks
    validate(
      need(length(input$selected.tfs) > 0 , "Please select a set of transcription factors"),
      need(input$target.gene, "Please select a target gene")
    )
    
    ## Extract target gene annotation 
    gene.name <-  strsplit(input$target.gene,split=" - ")[[1]][1]
    
    target.gene.body <- genes.data[gene.name,]
    target.gene.chr <- as.character(target.gene.body$seqnames)
    target.gene.start <- target.gene.body$start
    target.gene.end <- target.gene.body$end
    
    target.gene.strand <- as.character(target.gene.body$strand)
    
    ## Extract cds annotation
    cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
    
    ## Extract exons annotation
    exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
    
    ## Determine the genome range to plot including promoter, gene body and 3' UTR
    ## This depends on whether the gene is on the forward or reverse strand
    range.to.plot <- target.gene.body
    
    if(target.gene.strand == "+")
    {
      range.to.plot$start <- range.to.plot$start - input$promoter.length
      range.to.plot$end <- range.to.plot$end + input$threeprime.length
    } else if (target.gene.strand == "-")
    {
      range.to.plot$end <- range.to.plot$end + input$promoter.length
      range.to.plot$start <- range.to.plot$start - input$threeprime.length
    }
    
    ## Compute the length of the genome range to represent
    current.length <- range.to.plot$end - range.to.plot$start
    
    ## Determine upper limit of the graph
    number.tfs <- length(input$selected.tfs)
    upper.lim <- 25 * length(input$selected.tfs)
    
    ## Draw DNA strand
    gene.height <- -25
    cord.x <- 1:current.length
    
    plot(cord.x, rep(gene.height,length(cord.x)),type="l",col="black",lwd=3,ylab="",
         cex.lab=2,axes=FALSE,xlab="",main="",cex.main=2,
         ylim=c(-30,25 * length(input$selected.tfs)),
         xlim=c(-3000,max(cord.x)))
    
    ## Extract exons for target gene
    exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
    
    ## Transform exon coordinates to current range
    min.pos <- min(exons.data.target.gene$start)
    
    if(target.gene.strand == "+")
    {
      exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$promoter.length
      exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$promoter.length
    } else if(target.gene.strand == "-")
    {
      exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$threeprime.length
      exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$threeprime.length
    }
    
    ## Represent exons
    exon.width <- 2
    for(i in 1:nrow(exons.data.target.gene))
    {
      # Determine start/end for each exon
      current.exon.start <- exons.data.target.gene$start[i]
      current.exon.end <- exons.data.target.gene$end[i]
      
      ## Determine coordinates for each exon polygon and represent it
      exon.x <- c(current.exon.start,current.exon.end,current.exon.end,current.exon.start)
      exon.y <- c(gene.height + exon.width, gene.height + exon.width, gene.height - exon.width, gene.height - exon.width)
      
      polygon(x = exon.x, y = exon.y, col = "blue",border = "blue")
    }
    
    ## Extract cds for target gene
    cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
    
    ## Transform cds coordinates to current range
    if(target.gene.strand == "+")
    {
      cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$promoter.length
      cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$promoter.length
    } else if (target.gene.strand == "-")
    {
      cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$threeprime.length
      cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$threeprime.length
    }
    
    cds.width <- 3
    for(i in 1:nrow(cds.data.target.gene))
    {
      # Determine current cds start/end
      current.cds.start <- cds.data.target.gene$start[i]
      current.cds.end <- cds.data.target.gene$end[i]
      
      # Determine curret cds coordinates for the polygon and represent it
      cds.x <- c(current.cds.start,current.cds.end,current.cds.end,current.cds.start)
      cds.y <- c(gene.height + cds.width, gene.height + cds.width, gene.height - cds.width, gene.height - cds.width)
      
      polygon(x = cds.x, y = cds.y, col = "blue",border = "blue")
    }
    
    ## Draw arrow to represent transcription direction 
    if(target.gene.strand == "+")
    {
      lines(c(input$promoter.length,input$promoter.length,input$promoter.length+100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
      lines(c(input$promoter.length+50,input$promoter.length+100),y=c(gene.height+6,gene.height+5),lwd=3)
      lines(c(input$promoter.length+50,input$promoter.length+100),y=c(gene.height+4,gene.height+5),lwd=3)
    } else if (target.gene.strand == "-")
    {
      lines(c(current.length - input$promoter.length, current.length - input$promoter.length, current.length - input$promoter.length-100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
      lines(c(current.length - input$promoter.length-50, current.length - input$promoter.length - 100),y=c(gene.height + 6, gene.height + 5),lwd=3)
      lines(c(current.length - input$promoter.length-50, current.length - input$promoter.length - 100),y=c(gene.height + 4, gene.height + 5),lwd=3)
    }
    
    ## Draw promoter range
    if(target.gene.strand == "+")
    {
      axis(side = 1,labels = c(- input$promoter.length, - input$promoter.length / 2,"TSS"),at = c(1,input$promoter.length/2,input$promoter.length),lwd=2,cex=1.5,las=2,cex=2)
    } else if(target.gene.strand == "-")
    {
      axis(side = 1,labels = c("TSS",- input$promoter.length / 2,- input$promoter.length),at = c(current.length-input$promoter.length,current.length-input$promoter.length/2, current.length),lwd=2,cex=1.5,las=2,cex=2)
    }
    
    selected.bigwig.files <- bigwig.files[input$selected.tfs]
    selected.bed.files <- bed.files[input$selected.tfs]
    
    ## Since ChIPpeakAnno needs more than one region to plot our region
    ## is duplicated 
    regions.plot <- GRanges(rbind(range.to.plot,range.to.plot))
    
    ## Import signal from the bigwig files
    cvglists <- sapply(selected.bigwig.files, import, 
                       format="BigWig", 
                       which=regions.plot, 
                       as="RleList")
    
    names(cvglists) <- input$selected.tfs
    
    ## Compute signal in the region to plot
    chip.signal <- ChIPpeakAnno::featureAlignedSignal(cvglists, regions.plot, 
                                                      upstream=ceiling(current.length/2), 
                                                      downstream=ceiling(current.length/2),
                                                      n.tile=current.length) 
    
    ## Compute mean signal 
    chip.signal.means <- matrix(nrow=length(input$selected.tfs), ncol=ncol(chip.signal[[1]]))
    
    for(i in 1:length(input$selected.tfs))
    {
      if(target.gene.strand == "+")
      {
        chip.signal.means[i, ] <- colMeans(chip.signal[[i]],na.rm = TRUE)
      } else if (target.gene.strand == "-")
      {
        chip.signal.means[i, ] <- rev(colMeans(chip.signal[[i]],na.rm = TRUE))
      }
    }
    
    ## Draw peak regions for each TF and determing TF binding sequences
    
    ## Determine TFBS motifs to search for and Selecting an
    ## example motif if the user does not select any of them
    if(input$all.motifs)
    {
      selected.motifs.pwm <- motifs.pwm
    } else if(length(input$selected.motifs)==0)
    {
      selected.motifs.pwm <- motifs.pwm["EE"]
    }else
    {
      selected.motifs.pwm <- motifs.pwm[input$selected.motifs]
    }
    
    selected.motif.names <- names(selected.motifs.pwm)
    selected.motif.ids <- motif.ids[selected.motif.names]
    
    ## Initialize data frame containing TF binding sequences in the peak regions
    df.hits <- data.frame(0,0,"","","")
    colnames(df.hits) <- c("tf_number","position","id","name","seq")
    
    ## Width of the rectangle representing the peak region
    peak.width <- 1
    for(i in 1:length(input$selected.tfs))
    {
      ## Extract bed file name 1 and read it
      current.bed.file <- selected.bed.files[i]
      current.peaks <- read.table(file=current.bed.file,header = F, as.is = T)
      peak.coordinates <- subset(current.peaks, V1 == range.to.plot$seqnames & V2 >= range.to.plot$start & V3 <= range.to.plot$end) 
      current.peaks.to.plot <- peak.coordinates[,2:3]
      
      ## Transform coordinates 
      current.peaks.to.plot <- current.peaks.to.plot - range.to.plot$start
      
      ## Check if there are peaks for the target gene
      if(nrow(current.peaks.to.plot) > 0)
      {
        ## Normalization
        chip.signal.means[i, ] <- 10 * chip.signal.means[i, ] / max(chip.signal.means[i, ])
        
        #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
        for(j in 1:nrow(current.peaks.to.plot))
        {
          ## Extract start and end point of each peak region
          current.peak.start <- current.peaks.to.plot[j,1]
          current.peak.end <- current.peaks.to.plot[j,2]
          
          ## Computer coordinates for polygon and draw it
          peak.x <- c(current.peak.start,current.peak.end,
                      current.peak.end,current.peak.start)
          peak.y <- c(25*(i - 1) - 5 + peak.width, 25*(i - 1) - 5 + peak.width, 
                      25*(i - 1) - 5 - peak.width, 25*(i - 1) - 5 - peak.width)  
          
          polygon(x = peak.x, y = peak.y, col = area.colors[i], border = line.colors[i],lwd=2)
          
          ## Identify TF binding DNA motifs 
          peak.chr <- peak.coordinates[j, 1]
          peak.start <- peak.coordinates[j, 2]
          peak.end <- peak.coordinates[j, 3]
          
          ## Extract peak sequence
          if(peak.chr == "1")
          {
            peak.sequence <- toString(getSeq(x = chr1)[[1]][peak.start:peak.end]) #c2s(chr1[peak.start:peak.end])
          } else if(peak.chr == "2")
          {
            peak.sequence <- toString(getSeq(x = chr2)[[1]][peak.start:peak.end]) #c2s(chr2[peak.start:peak.end])
          } else if(peak.chr == "3")
          {
            peak.sequence <- toString(getSeq(x = chr3)[[1]][peak.start:peak.end]) #c2s(chr3[peak.start:peak.end])
          } else if(peak.chr == "4")
          {
            peak.sequence <- toString(getSeq(x = chr4)[[1]][peak.start:peak.end]) #c2s(chr4[peak.start:peak.end])
          } else if(peak.chr == "5")
          {
            peak.sequence <- toString(getSeq(x = chr5)[[1]][peak.start:peak.end]) #c2s(chr5[peak.start:peak.end])
          }
          
          peak.rev.comp.sequence <- reverse.complement(peak.sequence)
          
          for(k in 1:length(selected.motifs.pwm))
          {
            motif.pwm <- selected.motifs.pwm[[k]]
            
            hits.fw <- matchPWM(motif.pwm, peak.sequence, 
                                min.score = paste0(input$min.score.pwm,"%"))
            hits.fw.seqs <- as.data.frame(hits.fw)[[1]]
            hits.fw <- as(hits.fw, "IRanges")
            hits.fw.start <- start(hits.fw)
            hits.fw.end <- end(hits.fw)
            
            if(length(hits.fw.start) > 0)
            {
              df.hits.fw <- data.frame(rep(i,length(hits.fw.start)),
                                       ((hits.fw.start+hits.fw.end)/2) + current.peak.start,
                                       rep(selected.motif.ids[k],length(hits.fw.start)),
                                       rep(selected.motif.names[k],length(hits.fw.start)),
                                       hits.fw.seqs)
              colnames(df.hits.fw)  <- c("tf_number","position","id","name","seq")
              df.hits <- rbind(df.hits,df.hits.fw)
            }
            
            hits.rev <- matchPWM(motif.pwm, peak.rev.comp.sequence, 
                                 min.score = paste0(input$min.score.pwm,"%"))
            hits.rev.seqs <- as.data.frame(hits.rev)[[1]]
            hits.rev.seqs <- sapply(hits.rev.seqs,reverse.complement)
            names(hits.rev.seqs) <- NULL
            
            hits.rev <- as(hits.rev, "IRanges")
            hits.rev.start <- nchar(peak.sequence) - end(hits.rev) + 1
            hits.rev.end <- nchar(peak.sequence) - start(hits.rev) + 1
            
            if(length(hits.rev.start) > 0)
            {
              df.hits.rev <- data.frame(rep(i,length(hits.rev.start)),
                                        ((hits.rev.start+hits.rev.end)/2) + current.peak.start,
                                        rep(selected.motif.ids[k],length(hits.rev.start)),
                                        rep(selected.motif.names[k],length(hits.rev.start)),
                                        hits.rev.seqs)
              colnames(df.hits.rev)  <- c("tf_number","position","id","name","seq")
              df.hits <- rbind(df.hits,df.hits.rev)
            }
            
          }
          
        }
      }
    }
    
    ## Remove first line of the data frame added just for technical reason
    df.hits <- df.hits[-1,]
    
    ## Draw TF binding sites
    detected.tfbs <- unique(as.vector(df.hits$name))
    
    number.of.shapes <- ceiling(length(detected.tfbs) / length(symbol.color))
    
    necessary.shapes <- rep(symbol.shapes[1:number.of.shapes],each = length(detected.tfbs)/number.of.shapes)
    necessary.colors <- rep(symbol.color,number.of.shapes)
    
    if(length(detected.tfbs) > 0)
    {
      for(i in 1:length(detected.tfbs))
      {
        current.tfbs <- detected.tfbs[i]
        current.shape <- necessary.shapes[i]
        current.color <- necessary.colors[i]
        
        positions <- subset(df.hits, name == current.tfbs)
        
        for(j in 1:nrow(positions))
        {
          tf.to.draw <- positions$tf_number[j]
          pos.to.draw <- positions$position[j]
          
          points(x = pos.to.draw, y = 25*(tf.to.draw - 1) - 5 - 5*peak.width,
                 pch = current.shape, col = current.color, cex = 1)
        }
      }
      
      ## Add legend for TFBS
      legend.step <- 10
      for(i in 1:length(detected.tfbs))
      {
        points(x = -3000, y = upper.lim - (i-1)*legend.step, 
               pch=necessary.shapes[i], col = necessary.colors[i],cex = 1)
        
        
        current.seq <- as.character(subset(df.hits,name == detected.tfbs[i])[["seq"]][[1]])
        current.label <- paste(c(detected.tfbs[i], "  -  ", current.seq ),collapse="")
        
        text(x = -2900, y = upper.lim - (i-1)*legend.step, labels = current.label,
             adj = 0,cex = 0.7)
      }
    }
    
    ## Draw profiles for TF binding
    for(i in 1:length(input$selected.tfs))
    {
      ## Compute base line for current TF
      current.base.line <- 25 * (i - 1)
      
      ## Represent signal from the current TF
      lines(chip.signal.means[i,]+current.base.line,type="l",col=line.colors[i],lwd=3)
      
      ## Determine polygon coordinates and represent it
      cord.y <- c(current.base.line,chip.signal.means[i,]+current.base.line,current.base.line)
      cord.x <- 1:length(cord.y)
      
      polygon(cord.x,cord.y,col=area.colors[i])
      
      text(x = -50,y = 25*(i-1) + 12,labels = input$selected.tfs[i],adj = 1,col = line.colors[i],font = 2)
    }
    
  })
  
  ## Multiple transcription factor code
  
  ## Initial/default visualization of ATTRACTOR
  output$networkPlot <- renderPlot({
    ggplot(network.data, aes(x.pos,y.pos)) + 
      theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks.y = element_blank()) + 
      geom_point(color=node.colors,size=1)
  },height = 700)
  
  ## Determine common targets and perform analysis when button is clicked
  observeEvent(input$go_multiple, {
    
    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    if(input$peak != "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, (peak.zt == input$peak & trough.zt == input$trough))
    } else if(input$peak == "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, trough.zt == input$trough)
    } else if(input$peak != "any" && input$trough == "any")
    {
      selected.genes.df <- subset(selected.genes.df, peak.zt == input$peak)
    }
    
    if(nrow(selected.genes.df) == 0)
    {
      output$empty_overlap_message_2 <- renderText(expr = "The intersection between the target genes of the
                                                   selected transcription factors is empty. 
                                                   The selected TFs do no have any common target genes. No
                                                   further analysis can be performed.")
      
      network.representation <- ggplot(network.data, aes(x.pos,y.pos)) + 
        theme(panel.background = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) + 
        geom_point(color=node.colors,size=1)
      
      output$overlap.significance.text <- renderText(expr = {""})
      output$venn.diagram.plot <- renderPlot(expr = {})
      
    } else
    {
      
      output$empty_overlap_message_2 <- renderText(expr = "")
      
      ## Node colors for representation
      selected.nodes.colors <- selected.colors[selected.genes.df$peak.zt]
      
      ## Determine significance of overlap
      ## Number of sets to overlap
      if(input$peak != "any" || input$trough != "any")
      {
        number.of.sets <- length(input$selected.multiple.tfs) + 1
      } else
      {
        number.of.sets <- length(input$selected.multiple.tfs)
      }
      
      list.of.sets <- vector(mode = "list",length=number.of.sets)
      
      ## TFs targets
      for(i in 1:length(selected.tfs.with.zts))
      {
        list.of.sets[[i]] <- rownames(network.data)[which(network.data[,selected.tfs.with.zts[[i]]] != 0)]
      }
      
      ## Gene cluster with specified expression pattern
      if(input$peak != "any" && input$trough != "any")
      {
        list.of.sets[[length(list.of.sets)]] <- rownames(subset(network.data, (peak.zt == input$peak & trough.zt == input$trough)))
        
        expression.gene.set.name <- paste0(paste0("peak at ZT",paste(substr(start = 5,stop=nchar(input$peak),x = input$peak))),
                                           paste0(" and trough at ZT",paste(substr(start = 5,stop=nchar(input$trough),x = input$trough))))
        name.of.sets <- c(input$selected.multiple.tfs, expression.gene.set.name)
      } else if(input$peak == "any" && input$trough != "any")
      {
        list.of.sets[[length(list.of.sets)]] <- rownames(subset(network.data, trough.zt == input$trough))
        
        expression.gene.set.name <- paste0(paste0("trough at ZT",paste(substr(start = 5,stop=nchar(input$trough),x = input$trough))))
        name.of.sets <- c(input$selected.multiple.tfs, expression.gene.set.name)
      } else if(input$peak != "any" && input$trough == "any")
      {
        list.of.sets[[length(list.of.sets)]] <- rownames(subset(network.data, peak.zt == input$peak))
        
        expression.gene.set.name <- paste0(paste0("peak at ZT",paste(substr(start = 5,stop=nchar(input$peak),x = input$peak))))
        name.of.sets <- c(input$selected.multiple.tfs, expression.gene.set.name)
      } else
      {
        name.of.sets <- input$selected.multiple.tfs
      }
      
      names(list.of.sets) <- name.of.sets
      
      ## Compute overlap p value and enrichment
      overlap.output <- supertest(x = list.of.sets, n = nrow(network.data))
      overlap.results <- summary(overlap.output)
      overlap.p.value <- tail(overlap.results$P.value, n=1)
      overlap.p.value <- round(overlap.p.value,digits = -log10(overlap.p.value)+2)
      overlap.enrichment <- round((overlap.results$Table)[["FE"]][nrow(overlap.results$Table)],digits=2)
      
      ## Output text for the overlap
      overlap.message <- "The overlap between the targets of the transcription factors"
      text.first.tfs <- paste(input$selected.multiple.tfs[1:(length(input$selected.multiple.tfs)-1)],collapse=",")
      text.all.tfs <- paste(c(text.first.tfs, "and", input$selected.multiple.tfs[length(input$selected.multiple.tfs)]),collapse=" ")
      overlap.message <- paste(overlap.message,text.all.tfs,sep=" ")
      
      if(input$peak != "any" && input$trough != "any")
      {
        overlap.message <- paste(overlap.message, paste(c("and the genes with expression peaks at ZT", substr(x=input$peak, start = 5,stop = nchar(input$peak)), 
                                                          "and expression troughs at", substr(x=input$trough, start = 5,stop = nchar(input$trough))), collapse=" "),sep= " ")
      } else if(input$peak == "any" && input$trough != "any")
      {
        overlap.message <- paste(overlap.message, paste(c("and the genes with expression troughs at", input$trough), collapse=" "), sep = " ")
      } else if(input$peak != "any" && input$trough == "any")
      {
        overlap.message <- paste(overlap.message, paste(c("and the genes with expression peaks at ZT", substr(x=input$peak, start = 5,stop = nchar(input$peak))), collapse=" "), sep = " ")
      }
      
      if(overlap.p.value < 0.05)
      {
        overlap.message <- paste(overlap.message, paste0(paste(c("is significant with a p-value of",overlap.p.value,"and an enrichment of", overlap.enrichment),collapse=" "),"."), sep=" ")
      } else
      {
        overlap.message <- paste(overlap.message, paste0(paste(c("is NOT significant according to a p-value of",overlap.p.value),collapse=" "),"."), sep = " ")
      }
      
      
      output$overlap.significance.text <- renderText(expr = {
        overlap.message
      })
      
      ## Draw venn diagram
      if(number.of.sets == 1)
      {
        output$overlap.significance.text <- renderText(expr = {
          "A single transcription factor or set of genes with specific expression pattern was selected. The analysis of 
        overlap significance is not applicable."
        })
      }
      if(number.of.sets == 2)
      {
        output$venn.diagram.plot <- renderPlot(expr = {
          grid.newpage()
          draw.pairwise.venn(area1 = length(list.of.sets[[1]]),
                             area2 = length(list.of.sets[[2]]),
                             cross.area = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                             category = name.of.sets,
                             scaled = TRUE,lwd = 2,col = "black",
                             fill = c("blue","red"),alpha = 0.7,
                             cex = 2,cat.cex = 1.5,cat.pos = c(0,0))
        },width = 600,height = 600)
      } else if (number.of.sets == 3)
      {
        output$venn.diagram.plot <- renderPlot(expr = {
          grid.newpage()
          draw.triple.venn(area1 = length(list.of.sets[[1]]),
                           area2 = length(list.of.sets[[2]]),
                           area3 = length(list.of.sets[[3]]),
                           n12 = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                           n23 = length(intersect(list.of.sets[[2]],list.of.sets[[3]])),
                           n13 = length(intersect(list.of.sets[[1]],list.of.sets[[3]])), 
                           n123 =  length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[3]])),
                           category = name.of.sets,
                           scaled = TRUE,lwd = 2,col = "black",
                           fill = c("blue","red","darkgreen"),alpha = 0.7,
                           cex = 1.5,cat.cex = 1.5,cat.pos = c(-30,30,180))
        }, width = 600, height = 600)
      } else if (number.of.sets == 4)
      {
        output$venn.diagram.plot <- renderPlot(expr = {
          grid.newpage()
          draw.quad.venn(area1 = length(list.of.sets[[1]]),
                         area2 = length(list.of.sets[[2]]),
                         area3 = length(list.of.sets[[3]]),
                         area4 = length(list.of.sets[[4]]),
                         n12 = length(intersect(list.of.sets[[1]],list.of.sets[[2]])),
                         n13 = length(intersect(list.of.sets[[1]],list.of.sets[[3]])), 
                         n14 = length(intersect(list.of.sets[[1]],list.of.sets[[4]])),
                         n23 = length(intersect(list.of.sets[[2]],list.of.sets[[3]])),
                         n24 = length(intersect(list.of.sets[[2]],list.of.sets[[4]])),
                         n34 = length(intersect(list.of.sets[[3]],list.of.sets[[4]])),
                         n123 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[3]])), 
                         n124 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]),list.of.sets[[4]])),
                         n134 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[3]]),list.of.sets[[4]])),
                         n234 = length(intersect(intersect(list.of.sets[[2]],list.of.sets[[3]]),list.of.sets[[4]])),
                         n1234 = length(intersect(intersect(list.of.sets[[1]],list.of.sets[[2]]), intersect(list.of.sets[[3]],list.of.sets[[4]]))),
                         category = name.of.sets,
                         scaled = TRUE,lwd = 2,col = "black",
                         fill = c("blue","red","darkgreen","orange"),alpha = 0.7,
                         cex = 1.5,cat.cex = 1.5,cat.pos = c(0,0,0,0))
        }, width = 600, height = 600)
      } else if (number.of.sets > 4)
      {
        output$venn.diagram.plot <- renderPlot(expr = {
          plot(overlap.output, Layout = "landscape")
        })
      }
      ## Target gene representation on the network
      network.representation <- ggplot(network.data, aes(x.pos,y.pos)) + 
        theme(panel.background = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks.y = element_blank()) + 
        geom_point(color=node.colors,size=1) +
        geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
      
      ## Add edges to the network when selected
      if(input$edges)
      {
        for(i in 1:length(input$selected.multiple.tfs))
        {
          tf.xpos <- subset(network.data, names == tf.ids[selected.only.tfs[i]])[["x.pos"]]
          tf.ypos <- subset(network.data, names == tf.ids[selected.only.tfs[i]])[["y.pos"]]
          network.representation <- network.representation +
            annotate("segment",
                     x=rep(tf.xpos,nrow(selected.genes.df)),
                     y=rep(tf.ypos,nrow(selected.genes.df)),
                     xend=selected.genes.df$x.pos,
                     yend=selected.genes.df$y.pos, 
                     color="grey", arrow=arrow(type="closed",length=unit(0.1, "cm")))
        }
        
        network.representation <- network.representation + 
          geom_point(data = selected.genes.df,aes(x.pos,y.pos), size=4, fill=selected.nodes.colors,colour="black",pch=21)
      }
      
    }
    
    ## Update network representation on the app
    output$networkPlot <- renderPlot({
      network.representation
    },height = 700)
    
    
    ## Output message when empty information
    if(nrow(selected.genes.df) == 0)
    {
      output$empty_overlap_message_1 <- renderText(expr = "The intersection between the target genes of the
                                                   selected transcription factors is empty. 
                                                   The selected TFs do no have any common target genes. No
                                                   further analysis can be performed.")
      output$outputTable <- renderDataTable({})
      output$download_ui_for_table<- renderUI("")
    } else
    {
      ## Output table with gene info
      output$outputTable <- renderDataTable({
        create.output.table(input.gene.df=selected.genes.df,alias,tfs.names)
      },escape=FALSE)
      
      ## Generate UI to download table and creating a downlodable table
      output$download_ui_for_table<- renderUI(
        tagList(downloadButton(outputId= "downloadData", "Download Selected Genes"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
      )
      
      output$downloadData<- downloadHandler(
        filename= function() {
          paste0(paste(input$selected.tfs,collapse = "_"), ".tsv")
        },
        content= function(file) {
          write.table(create.downloadable.output.table(input.gene.df=selected.genes.df,alias,tfs.names), 
                      file=file, 
                      sep = "\t", 
                      quote = FALSE,
                      row.names = FALSE)
        })
    }
    
  })
  
  ##Perform GO terms enrichment analysis when button is clicked
  observeEvent(input$goterm,{
    
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )
    
    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    ## Check whether or not the selected TFs share some common targets
    if(nrow(selected.genes.df) == 0)
    {
      output$empty_overlap_message_3 <- renderText(expr = "The intersection between the target genes of the
                                                   selected transcription factors is empty. 
                                                   The selected TFs do no have any common target genes. No
                                                   further analysis can be performed.")
    } else
    {
      ## GO enrichment analysis
      
      output$empty_overlap_message_3 <- renderText(expr = "")
      
      ## Set the background to perform the GO terms enrichment analysis depending on the user selection
      if (input$go.background == "allgenome")
      {
        go.universe <- ath.universe
      } else
      {
        go.universe <- network.data$name
      }
      
      
      ## Show element when GO term enrichment analysis starts
      shinyjs::showElement(id = 'loading.div')
      
      enrich.go <- enrichGO(gene          = selected.genes.df$name,
                            universe      = go.universe,
                            OrgDb         = org.At.tair.db,
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE,
                            keyType = "TAIR")
      
      ## Hide loading element when GO term enrichment is done
      shinyjs::hideElement(id = 'loading.div')
      
      ## Generate output table
      enrich.go.result <- as.data.frame(enrich.go)
      
      if(nrow(enrich.go.result) > 0)
      {
        ## Set to empty the message stating no go enrichment found
        output$no_go_results <- renderText({""})
        
        ## GO term Description P-value Q-value Enrichment (SetRatio, BgRatio) Genes
        go.term.enrichments <- compute.enrichments(gene.ratios = enrich.go.result$GeneRatio,
                                                   bg.ratios = enrich.go.result$BgRatio)
        
        go.result.table <- data.frame(enrich.go.result$ID, enrich.go.result$Description,
                                      enrich.go.result$pvalue, enrich.go.result$qvalue,
                                      go.term.enrichments, 
                                      gsub(pattern = "/",replacement = " ",x = enrich.go.result$geneID),
                                      stringsAsFactors = FALSE)
        
        colnames(go.result.table) <- c("GO ID", "Description", "p-value", "q-value",
                                       "Enrichment (Target Ratio; BG Ration)","Genes")
        
        go.result.table.with.links <- go.result.table
        ## Add links to the genes
        genes.in.go.enrichment <- go.result.table$Genes
        
        ## Add link to genes
        for(i in 1:length(genes.in.go.enrichment))
        {
          go.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(genes.in.go.enrichment[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }
        
        ## Add links to GO ids
        go.result.table.with.links[["GO ID"]] <- sapply(X = go.result.table.with.links[["GO ID"]], FUN = go.link)
        
        ## Introductory text for GO enrichment table
        go.table.text <- "The table below summarizes the result of the GO term
      enrichment analysis. Each row represents a GO term significantly enriched in the selected
      gene set. The first column represents the GO term
      identifier. The second column contains a human readable description. For more details on the
      corresponding GO term, click on the identifier in the first column. The third and fourth
      column presents the p-value and q-value (adjusted p-value or FDR) capturing the level
      of significance. The fifth column displays the corresponding enrichment value E (m/n; M/N) where
      n is the number of genes with annotation from the target set, N is the number of genes with
      annotation from the gene universe, m is the number of genes from the target set annotated with the
      corresponding GO term and M is the number of genes from the gene universe annotated with
      the GO term associated with the corresponding row. The enrichment is then computed as
      E = (m/n) / (M/N). Finally, the last column, contains the genes from the target set
      annotated with the GO term represented in the corresponding row."
        
        output$textGOTable <- renderText(expr = go.table.text)
        
        ## Output table with GO enrichment result
        output$output_go_table <- renderDataTable({
          ## Error message for the user
          validate(
            need(input$selected.multiple.tfs, "Please select some transcription factor")
          )
          go.result.table.with.links #go.result.table
        },escape=FALSE,options =list(pageLength = 5)) 
        
        output$download_ui_for_go_table<- renderUI(
          tagList(downloadButton(outputId= "downloadGOData", "Download GO Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
        )
        
        output$downloadGOData<- downloadHandler(
          filename= function() {
            paste0(paste(c("GO_enrichment",input$selected.tfs),collapse = "_"), ".tsv")
          },
          content= function(file) {
            write.table(go.result.table, 
                        file=file, 
                        sep = "\t", 
                        quote = FALSE,
                        row.names = FALSE)
          })
        
        ## Download result
        output$downloadData<- downloadHandler(
          filename= function() {
            paste("godata" , ".csv", sep="")
          },
          content= function(file) {
            write.csv(go.result.table,
                      file,
                      row.names=TRUE
            )
          })
        
        ## Link to REVIGO 
        revigo.data <- paste(revigo.data <- apply(go.result.table[,c("GO ID", "q-value")], 1, paste, collapse = " "), collapse="\n")
        
        url1 <- tags$a("here", href="#", onclick="document.revigoForm.submit();")
        url2 <- tags$form(
          name="revigoForm", action="http://revigo.irb.hr/", method="post", target="_blank",
          tags$textarea(name="inputGoList", rows="1", cols="8", class="revigoText",
                        style="visibility: hidden", revigo.data)
        )
        
        output$revigo<- renderUI(
          tagList("The enriched GO terms above may be redundant. Visualize these results in REViGO in order to remove redundancy. Click", url1, url2)
        )
        
        
        output$barplot_text <- renderText("In the following barplot each bar represents a significantly enriched 
GO term. The length of the bar corresponds to the number of genes in the
                                        target set annotated with the given GO term. The bar color captures the level
                                        of significance from blue, less significant, to red, more significant.")
        
        
        ## GO map
        output$gomap_text <- renderText("The following figure corresponds to a subgraph
                                      induced by most significant GO terms from 
                                      Biological Process subcategory. Enriched terms 
                                      are colored and the color depends on the 
                                      adjusted p-value according to the Benjamini & Hochberg 
                                      method, increasing the p-value from purple to red")
        
        output$gomap <- renderPlot(
          width     = 1040,
          height    = 1000,
          res       = 120,
          expr = {
            ## Error message for the user
            validate(
              need(input$selected.multiple.tfs, "Please select some transcription factor")
            )
            goplot(enrich.go,showCategory = 10)
          })
        
        ## Barplot
        output$bar.plot <- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            ## Error message for the user
            validate(
              need(input$selected.multiple.tfs, "Please select some transcription factor")
            )
            barplot(enrich.go,drop=TRUE,showCategory = 10)
          })
        
        ##EMAP plot
        output$emapplot_text <- renderText("The following figure consists of an enrichment map where nodes represent enriched GO terms. The
        size of a node is proportional to the number of genes annotated with the corresponding GO term in the target set.
The node colors represent the level of significance from less signficant in blue to more significant in red. Edges are drawn
between two nodes when the corresponding GO terms are semantically related.")
        
        output$emap.plot <- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            emapplot(enrich.go)
          })
        
        ## CNET plot
        output$cnetplot_text <- renderText("The following figure corresponds to a gene-concept network. The beige
nodes represents GO terms and the grey nodes genes. An edge is drawn from a gene to a GO term when the gene is annotated
with the corresponding gene. The size of nodes representing GO terms is proportional to the number of genes annotated
with the corresponding GO term.")
        
        output$cnet.plot <- renderPlot(
          width     = 870,
          height    = 600,
          res       = 120,
          expr = {
            cnetplot(enrich.go)
          })
      } else
      {
        output$no_go_results <- renderText({"No GO term enrichment was found."})
        output$textGOTable <- renderText("")
        output$output_go_table <- renderDataTable({}) 
        output$download_ui_for_go_table<- renderUI("")
        #output$downloadGOData<- downloadHandler("")
        #output$downloadData<- downloadHandler("")
        output$revigo<- renderUI("")
        output$barplot_text <- renderText("")
        output$gomap_text <- renderText("")
        output$gomap <- renderPlot("")
        output$bar.plot <- renderPlot("")
        output$emapplot_text <- renderText("")
        output$emap.plot <- renderPlot("")
        output$cnetplot_text <- renderText("")
        output$cnet.plot <- renderPlot("")
      }
    }
  })
  
  ##Perform KEGG pathway enrichment analysis when button is clicked
  observeEvent(input$pathway_button,{
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )
    
    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]    
    
    if(nrow(selected.genes.df) == 0)
    {
      output$empty_overlap_message_4 <- renderText(expr = "The intersection between the target genes of the
                                                   selected transcription factors is empty. 
                                                   The selected TFs do no have any common target genes. No
                                                   further analysis can be performed.")
      output$output_pathway_table <- renderDataTable({""})
      output$download_ui_for_kegg_table<- renderUI("")
      output$kegg_selectize <- renderUI({""})
      output$output_module_table <- renderDataTable({""}) 
      output$kegg_image <- renderImage({filename = "blank.png"})
    } else
    {
      output$empty_overlap_message_4 <- renderText(expr = "")
      
      ## Set the background to perform pathway enrichment analysis depending on the user selection
      if (input$pathway_background == "allgenome")
      {
        pathway.universe <- ath.universe
      } else
      {
        pathway.universe <- network.data$name
      }
      
      ## Show element when kegg pathway enrichment starts
      shinyjs::showElement(id = 'loading.div.kegg')
      
      ## Compute KEGG pathway enrichment
      pathway.enrichment <- enrichKEGG(gene = selected.genes.df$name, 
                                       organism = "ath", 
                                       keyType = "kegg",
                                       universe = pathway.universe,
                                       qvalueCutoff = 0.05)
      pathway.enrichment.result <- as.data.frame(pathway.enrichment)
      
      ## Hide loading element when KEGG pathway enrichment is done
      shinyjs::hideElement(id = 'loading.div.kegg')
      
      ## Generate output table
      if(nrow(pathway.enrichment.result) > 0)
      {
        pathways.enrichment <- compute.enrichments(gene.ratios = pathway.enrichment.result$GeneRatio,
                                                   bg.ratios = pathway.enrichment.result$BgRatio)
        
        ## Separate genes with blank spaces
        kegg.enriched.genes <- pathway.enrichment.result$geneID
        for(i in 1:length(kegg.enriched.genes))
        {
          kegg.enriched.genes[i] <- paste(strsplit(kegg.enriched.genes[i],split="/")[[1]],collapse=" ")
        }
        
        ## Generate data frame with output table
        pathways.result.table <- data.frame(pathway.enrichment.result$ID, pathway.enrichment.result$Description,
                                            pathway.enrichment.result$pvalue, pathway.enrichment.result$qvalue,
                                            pathways.enrichment, 
                                            kegg.enriched.genes,
                                            stringsAsFactors = FALSE)
        
        colnames(pathways.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                             "Enrichment (Target Ratio; BG Ration)","Genes")
        
        kegg.result.table.with.links <- pathways.result.table
        
        ## Add links to genes
        for(i in 1:length(kegg.enriched.genes))
        {
          kegg.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(kegg.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }
        
        ## Add links to kegg pathways
        kegg.result.table.with.links[["KEGG ID"]] <- sapply(X=kegg.result.table.with.links[["KEGG ID"]],FUN = kegg.pathway.link)
        
        output$output_pathway_table <- renderDataTable({
          kegg.result.table.with.links
        },escape=FALSE,options =list(pageLength = 5)) 
        
        
        output$download_ui_for_kegg_table<- renderUI(
          tagList(downloadButton(outputId= "downloadKEGGData", "Download KEGG Pathway Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
        )
        
        output$downloadKEGGData<- downloadHandler(
          filename= function() {
            paste0(paste(c("KEGG_enrichment",input$selected.tfs),collapse = "_"), ".tsv")
          },
          content= function(file) {
            write.table(pathways.result.table, 
                        file=file, 
                        sep = "\t", 
                        quote = FALSE,
                        row.names = FALSE)
          })
        
        ## Download result
        output$downloadData<- downloadHandler(
          filename= function() {
            paste("keggdata" , ".csv", sep="")
          },
          content= function(file) {
            write.csv(pathways.result.table,
                      file,
                      row.names=TRUE
            )
          })
        
        
        ## Visualization of specific enriched pathways
        genes.pathway <- rep(0, length(pathway.universe))
        names(genes.pathway) <- pathway.universe
        
        genes.pathway[selected.genes.df$name] <- 1
        
        pathways.for.select <- paste(pathways.result.table[["KEGG ID"]], pathways.result.table[["Description"]], sep=" - ")
        
        output$kegg_selectize <- renderUI({
          selectInput(inputId = "kegg_pathway", 
                      label="Choose Pathway for Representation",
                      multiple = FALSE,
                      selected = pathways.for.select[1],
                      choices=pathways.for.select)
        })
        
        ## Enriched pathway image
        output$kegg_image <- renderImage({
          pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
                   pathway.id = strsplit(input$kegg_pathway,split=" - ")[[1]][1],
                   species = "ath",
                   limit = list(gene=max(abs(genes.pathway)), cpd=1),
                   gene.idtype ="kegg")
          
          list(src = paste(c(strsplit(input$kegg_pathway,split=" - ")[[1]][1],"pathview","png"), collapse="."),
               contentType="image/png",width=1200,height=900)
        },deleteFile = T)
        
      } else
      {
        output$no_kegg_enrichment <- renderText(expr = "No enriched KEGG pathway was detected in the selected genes.")
        output$no_pathway_visualization <- renderText(expr = "No enriched KEGG pathway was detected in the selected genes.")
        output$output_pathway_table <- renderDataTable({""})
        output$download_ui_for_kegg_table<- renderUI("")
        output$kegg_selectize <- renderUI("")
        print("llego aqui perfe")
        output$kegg_image <- renderImage({filename = "blank.png"})
      }
      
      ## KEGG module enrichment analysis
      modules.enrichment <- enrichMKEGG(gene = selected.genes.df$name, 
                                        universe = pathway.universe, 
                                        organism = "ath", 
                                        keyType = "kegg",
                                        minGSSize = 4)
      
      modules.enrichment.result <- as.data.frame(modules.enrichment)
      if(nrow(modules.enrichment.result) > 0)
      {
        modules.enrichment <- compute.enrichments(gene.ratios = modules.enrichment.result$GeneRatio,
                                                  bg.ratios = modules.enrichment.result$BgRatio)
        
        modules.enriched.genes <- modules.enrichment.result$geneID
        for(i in 1:length(modules.enriched.genes))
        {
          modules.enriched.genes[i] <- paste(strsplit(modules.enriched.genes[i],split="/")[[1]],collapse=" ")
        }
        
        modules.result.table <- data.frame(modules.enrichment.result$ID, modules.enrichment.result$Description,
                                           modules.enrichment.result$pvalue, modules.enrichment.result$qvalue,
                                           modules.enrichment, 
                                           modules.enriched.genes,
                                           stringsAsFactors = FALSE)
        
        colnames(modules.result.table) <- c("KEGG ID", "Description", "p-value", "q-value",
                                            "Enrichment (Target Ratio; BG Ration)","Genes")
        
        modules.result.table.with.links <- modules.result.table
        
        ## Add links to genes
        for(i in 1:length(modules.enriched.genes))
        {
          modules.result.table.with.links$Genes[i] <- paste(sapply(X = strsplit(modules.enriched.genes[i],split=" ")[[1]],FUN = gene.link.function),collapse=" ")
        }
        
        ## Add links to kegg pathways
        modules.result.table.with.links[["KEGG ID"]] <- sapply(X=modules.result.table.with.links[["KEGG ID"]],FUN = kegg.module.link)
        
        ## Generate output table
        output$output_module_table <- renderDataTable({
          modules.result.table.with.links
        },escape=FALSE,options =list(pageLength = 5)) 
      } else
      {
        output$text_module_kegg <- renderText(expr = "No enriched KEGG module was detected in the selected genes.")
        output$output_module_table <- renderDataTable({""}) 
      }
      
      
    }
    
    
  })
  
  ##Perform TFBS enrichment analysis when button is clicked
  observeEvent(input$tfbs_button,{
    ## Sanity checks
    validate(
      need(input$selected.multiple.tfs, "Please select a set of transcription factors")
    )
    
    ## Determine targets of selected TFs
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    if(length(selected.tfs.with.zts) > 1)
    {
      gene.selection <- rowSums(selected.tfs.adj) == length(selected.tfs.with.zts)
    } else if (length(selected.tfs.with.zts) == 1)
    {
      gene.selection <- as.vector(selected.tfs.adj)
    }
    
    ## Determine targets with the specified expression profile
    selected.genes.df <- network.data[gene.selection,]  
    
    if(input$peak != "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, (peak.zt == input$peak & trough.zt == input$trough))
    } else if(input$peak == "any" && input$trough != "any")
    {
      selected.genes.df <- subset(selected.genes.df, trough.zt == input$trough)
    } else if(input$peak != "any" && input$trough == "any")
    {
      selected.genes.df <- subset(selected.genes.df, peak.zt == input$peak)
    }
    
    ## Set the background to perform TFBS enrichment analysis depending on the user selection
    if (input$tfbs_background == "allgenome")
    {
      tfbs.universe <- 33323
    } else
    {
      tfbs.universe <- 5778
    }
    
    # ## Show element when TFBS enrichment starts
    # shinyjs::showElement(id = 'loading.div.tfbs')
    
    file.precomputed <- paste0(c("precomputed_",input$up_promoter,"_",
                                 input$down_promoter, "_",
                                 input$score, "_",tfbs.universe,".tsv"),collapse="")
    
    ## Load file with precomputed results (background) and compute m and n
    precomputed.result <- read.table(file=paste0("data/precomputed_results_tfbs/",file.precomputed),header = T)
    m <- colSums(precomputed.result > 0) 
    n <- nrow(precomputed.result) - m
    
    ## Compute selection size (k) and number of ocurrences (x).
    # target.genes <- read.table(file = "peak_ZT0_trough_ZT12.txt",as.is = T)[[1]]
    target.genes <- intersect(rownames(precomputed.result),selected.genes.df$names)
    
    if(length(target.genes) == 0)
    {
      output$output_tfbs_table <- renderDataTable("")
      output$download_ui_tfbs_table<- renderUI("")
      output$empty_overlap_message_5 <- renderText(expr = "The intersection between the target genes of the
                                                   selected transcription factors is empty. 
                                                   The selected TFs do no have any common target genes. No
                                                   further analysis can be performed.")
    } else
    {
      output$empty_overlap_message_5 <- renderText(expr = "")
      k <- length(target.genes)
      x <- colSums(precomputed.result[target.genes,] > 0)
      
      ## Compute p-values for enrichment according to a hypergeometric distribution
      p.values <- vector(mode="numeric", length=length(x))
      names(p.values) <- colnames(precomputed.result)
      
      for(i in 1:length(x))
      {
        p.values[i] <- phyper(q = x[i] - 1, m = m[i], n = n[i], k = k, lower.tail = F)
      }
      
      which(p.values < input$motif_significance)
      p.values[which(p.values < input$motif_significance)]
      
      ## Adjust p-values using Benjamini Hochberg
      q.values <- p.adjust(p = p.values,method = "BH")
      
      which(q.values < input$motif_significance)
      q.values[which(q.values < input$motif_significance)]
      
      ## Compute enrichments
      enrichments <- (x / k) / (m / nrow(precomputed.result))
      
      ## Final motifs, pvalues, qvalues and enrichments
      input <- list(motif_significance = 0.05, enrichment_threshold = 2 )
      
      sig.enrich.motifs <- names(which(q.values < input$motif_significance & enrichments > input$enrichment_threshold))
      sig.enrich.ids <- motif.ids[sig.enrich.motifs]
      final.q.values <- q.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
      final.p.values <- p.values[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
      final.enrichments <- enrichments[which(q.values < input$motif_significance & enrichments > input$enrichment_threshold)]
      
      if(length(sig.enrich.motifs) > 0)
      { 
        ## Set message of no TFBS enrichment to the empty string
        output$no_tfbs_enrichment <- renderText("")
        ## Determine genes for each motif
        genes.with.motif <- vector(length = length(sig.enrich.motifs))
        for (i in 1:length(sig.enrich.motifs))
        {
          #print(i)
          rows.with.motif <- which(precomputed.result[,sig.enrich.motifs[i]] != 0)
          all.genes.with.motif <- rownames(precomputed.result)[rows.with.motif]
          genes.with.motif[i] <- paste(intersect(all.genes.with.motif,target.genes), collapse = ",")
        }
        
        ## Motifs logos
        motifs.images <- paste0("motifs_images/",sig.enrich.motifs)
        
        for (i in 1:length(motifs.images))
        {
          motifs.images[i] <- paste0("<img src='",motifs.images[i],".png', align = 'center', width = 100>")
        }
        
        ## Store data
        tfbs.result.table <- data.frame(sig.enrich.motifs, sig.enrich.ids, motifs.images, final.p.values, final.q.values, final.enrichments, genes.with.motif, stringsAsFactors = FALSE,row.names = NULL) 
        colnames(tfbs.result.table) <- c("DNA motifs", "Motif ID", "DNA logo", "P-values", "Q-values", "Enrichments", "Genes")
        
        ## Add links to jaspar motifs
        tfbs.result.table[["Motif ID"]] <- sapply(X=sig.enrich.ids,FUN = tfbs.link)
        
        ## Add links to genes
        for (i in 1:length(genes.with.motif))
        {
          tfbs.result.table$Genes[i] <- paste(sapply(X = strsplit(genes.with.motif[i], split = ",")[[1]],FUN = gene.link.function), collapse = ", ")
        }
        
        tfbs.result.table <- tfbs.result.table[order(final.q.values),]
        ## Output table with TFBS enrichment result
        output$output_tfbs_table <- renderDataTable({
          tfbs.result.table
        },escape = FALSE,options =list(pageLength = 10))
        
        output$download_ui_tfbs_table<- renderUI(
          tagList(downloadButton(outputId= "downloadTFBSData", "Download TFBS Enrichment"),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br())
        )
        
        output$downloadTFBSData<- downloadHandler(
          filename= function() {
            paste0(paste(c("TFBS_enrichment",selected.tfs.with.zts, input$peak, input$trough),collapse = "_"), ".tsv")
          },
          content= function(file) {
            write.table(tfbs.result.table, 
                        file=file, 
                        sep = "\t", 
                        quote = FALSE,
                        row.names = FALSE)
          })
      } else
      {
        output$no_tfbs_enrichment <- renderText("No TFBS enrichment found on the selected set of genes.")
        output$output_tfbs_table <- renderDataTable({""})
        output$download_ui_tfbs_table<- renderUI("")
      }
      
    }
    
    
    
  })
  
  ##Retrieve bulk data when button is clicked
  observeEvent(input$retrieve_button, {
    selected.tfs.with.zts <- str_replace(string = input$selected.multiple.tfs,pattern = " ",replacement = "_")
    selected.only.tfs <- sapply(X = strsplit(x = input$selected.multiple.tfs,split = " "), FUN = get.first)
    selected.tfs.adj <- (network.data[,selected.tfs.with.zts] != 0)
    
    print("selected.tfs.with.zts:")
    print(selected.tfs.with.zts)
    print("selected.only.tfs:")
    print(selected.only.tfs)
    
    print("input$select_data_retrieval")
    print(input$select_data_retrieval)
    
    
    
    if("genomic_locations" %in% input$select_data_retrieval)
    {
      download.bed.files <-  bed.files[input$selected.multiple.tfs]
      output$bed_files_download <- renderUI(expr = {
        tagList(
          tags$ol(
            tags$ul( downloadLink(outputId = "first_file",label = download.bed.files[1])  )
          )
        )
      })
      
    } 
    
    if("genomic_signal"%in% input$select_data_retrieval)
    {
      download.bigwig.files <- bigwig.files[input$selected.multiple.tfs]
    }
    
    
    
    
  })
  
}
