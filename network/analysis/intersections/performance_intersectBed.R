# Authors: Pedro de los Reyes Rodríguez
#          Ana Belén Romero-Losada
#          Francisco J. Romero-Campero
# 
# Contact: Francisco J. Romero-Campero - fran@us.es

# Date: February 2019

#####Intersections between binding regions in DNA (BED files)####
#Reading the bed files of the transcription factors
peaks1 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/CCA1_ZT02_peaks.narrowPeak")
head(peaks1)

peaks2 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/CCA1_peaks.narrowPeak")
head(peaks2)

peaks3 <- read.table(file = "../../../web_apps/peak_visualizer/data/bed_files/PRR9_1_peaks.narrowPeak")
head(peaks3)

peaks.list <- list(peaks1, peaks2, peaks3)

length.sets <- sapply(X = peaks.list, FUN = nrow)


peaks.set1 <- peaks1
peaks.set2 <- peaks2
peaks.set2 <- random.peaks2

#intersectBed function 
intersectBed <- function(peaks.set1, peaks.set2)
{
  intersection <- matrix(ncol = 3, nrow=0 )
  current.intersection <- matrix(ncol = 3 )
  for (i in 1:nrow(peaks.set1))
  {
    #Set the current peak values of set1
    current.chr <- as.numeric(peaks.set1[i,1])
    current.start <- peaks.set1[i,2]
    current.end <- peaks.set1[i,3]
    #Checking if there is intersection between the current peak and any peak of set2
    option1 <- nrow(subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.start & peaks.set2[,3]>=current.start))
    option2 <- nrow(subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end))
    
    # print(i)
    
    if(option1+option2 > 0)
    {
      # print("HIT")
      if(option1>0)
      {
        hit.peak2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.start & peaks.set2[,3]>=current.start)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- current.start
        current.intersection[1,3] <- hit.peak2[1,3]
        
      }else
      {
        hit.peak2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end)
        current.intersection[1,1] <- current.chr
        current.intersection[1,2] <- hit.peak2[1,2]
        current.intersection[1,3] <- current.end
      }
      
      intersection <- rbind(intersection, current.intersection)
    }
  }
  return(intersection)
}

intersectBed.fran <- function(peaks.set1, peaks.set2)
{
  intersection <- matrix(ncol = 3, nrow=min(nrow(peaks.set1),nrow(peaks.set2)) )
  j <- 0
  for (i in 1:nrow(peaks.set1))
  {
    #Set the current peak values of set1
    current.chr <- as.numeric(peaks.set1[i,1])
    current.start <- peaks.set1[i,2]
    current.end <- peaks.set1[i,3]
    #Checking if there is intersection between the current peak and any peak of set2
    hit.1 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.start & peaks.set2[,3]>=current.start)
    option1 <- nrow(hit.1)
    hit.2 <- subset(peaks.set2, peaks.set2[,1]==current.chr & peaks.set2[,2]<=current.end & peaks.set2[,3]>=current.end)    
    option2 <- nrow(hit.2)
    
    # print(i)
    if( option1 + option2 > 0)    
    {
      # print("HIT")
      if(option1>0)
      {
        j <- j + 1
        intersection[j,1] <- current.chr
        intersection[j,2] <- current.start
        intersection[j,3] <- hit.1[1,3]
      } else 
      {
        j <- j + 1
        intersection[j,1] <- current.chr
        intersection[j,2] <- hit.2[1,2]
        intersection[j,3] <- current.end
      }
    }
      
  }
  return(intersection[1:j,])
}



# The intersectBed function allow to get the intersections between two bed files. If you want to perform the
# the intersection between three bed files, you can do it in a consecutive manner. 

first <- intersectBed(peaks1, peaks2)
nrow(first)

first <- intersectBed.fran(peaks.set1 = peaks1, peaks.set2 = peaks2)
nrow(first)

start.time <- Sys.time()
first <- intersectBed(peaks1, peaks2)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
first <- intersectBed.fran(peaks1, peaks2)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken




second <- intersectBed(first, peaks3)
nrow(second)

convert.df.to.granges <- function(input.df)
{
  output.granges <- input.df[,1:3]
  colnames(output.granges) <- c("seqnames", "start", "end")
  return(toGRanges(data = output.granges))
}


library(ChIPpeakAnno)

granges.peaks1 <- convert.df.to.granges(input.df = peaks1)
granges.peaks2 <- convert.df.to.granges(input.df = peaks2)

#intersect.1.2 <- 
system.time(intersectBed(peaks1, peaks2))
system.time(findOverlapsOfPeaks(granges.peaks1,granges.peaks2))

nrow(intersect.1.2)

help("findOverlapsOfPeaks")

start.time <- Sys.time()
res <- intersectBed(peaks1, peaks2)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
res <- findOverlapsOfPeaks(granges.peaks1,granges.peaks2)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



## Permutation of peaks.set2 (random.peaks2) and comparing with peaks.set1 ####
chromosomes.length <- read.table(file="../../../web_apps/peak_visualizer/data/bed_files/atha_chr_lengths.txt",as.is=T)[[1]]
number.randomisation <- 20 #100000

random.intersections <- vector(mode = "numeric",length=number.randomisation) #Creating vector
for(j in 1:number.randomisation)
{
  print(j)
  random.peaks2 <- matrix(nrow=nrow(peaks2),ncol=3) #Matriz con 3 columnas, una para el cromosoma, otra para el comienzo y otra para el final de la marca aleatoria.
  for(i in 1:nrow(peaks2))
  {
    current.chr <- peaks2[i,1][[1]] #Chr de la iésima marca real
    current.start <- peaks2[i,2] #Start de la iésima marca real
    current.end <- peaks2[i,3] #End de la iésima marca real
    current.length <- current.end - current.start #Longitud de la iésima marca real
    
    chr.length <- chromosomes.length[current.chr] #Length del actual cromosoma
    #Ahora genero los mismos datos para marcas aleatorias
    random.start <- floor(runif(n = 1,min = 1,max = chr.length))
    random.end <- random.start + current.length
    
    random.peaks2[i,1] <- current.chr
    random.peaks2[i,2] <- random.start
    random.peaks2[i,3] <- random.end
  }
  
  
  random.intersections[j] <- nrow(intersectBed(peaks.set1 = peaks1, peaks.set2 = random.peaks2 )) 
  
  p.value <- sum(random.intersections > nrow(first)) / number.randomisation
  p.value 
  
}


##Loop to check the intersection of binding regions (bed files) between all the transcription factors together and store the results in a table####
chromosomes.length <- read.table(file="../../../web_apps/peak_visualizer/data/bed_files/atha_chr_lengths.txt",as.is=T)[[1]]
number.randomisation <- 1000
bed.files <- list.files(path = "../../../web_apps/peak_visualizer/data/bed_files/", pattern = "peaks.narrowPeak")

combinations <- expand.grid(bed.files, bed.files)
bed.intersections <- matrix(ncol = 6, nrow = nrow(combinations))
colnames(bed.intersections) <- c("TF1", "TF2", "p-value", "fdr", "number of intersections", "Genes" )


txdb <- TxDb.Athaliana.BioMart.plantsmart28
i <- 26
total.tests <- nrow(combinations)
# total.tests <- 50

# Start the clock!
ptm <- proc.time()

for (i in 1:total.tests)
  # for (i in 1:nrow(combinations))
{
  # print(paste0("test number ", i, " of ", nrow(combinations)))
  print(paste0((i/total.tests)*100, " %"))
  peaks1 <- read.table(file = paste0("../../../web_apps/peak_visualizer/data/bed_files/", combinations[i,1]))
  peaks2 <- read.table(file = paste0("../../../web_apps/peak_visualizer/data/bed_files/", combinations[i,2]))
  real.intersection <- intersectBed(peaks.set1 = peaks1, peaks.set2 = peaks2)
  if (nrow(real.intersection) > 0)
  {
    random.intersections <- vector(mode = "numeric",length=number.randomisation) #Creating vector
    for(j in 1:number.randomisation)
    {
      print(j)
      random.peaks2 <- matrix(nrow=nrow(peaks2),ncol=3) #Matriz con 3 columnas, una para el cromosoma, otra para el comienzo y otra para el final de la región aleatoria.
      for(k in 1:nrow(peaks2))
      {
        current.chr <- peaks2[k,1][[1]] #Chr de la iésima marca real
        current.start <- peaks2[k,2] #Start de la iésima marca real
        current.end <- peaks2[k,3] #End de la iésima marca real
        current.length <- current.end - current.start #Longitud de la iésima marca real
        
        chr.length <- chromosomes.length[current.chr] #Length del actual cromosoma
        #Ahora genero los mismos datos para regiones aleatorias
        random.start <- floor(runif(n = 1,min = 1,max = chr.length))
        random.end <- random.start + current.length
        
        random.peaks2[k,1] <- current.chr
        random.peaks2[k,2] <- random.start
        random.peaks2[k,3] <- random.end
      }
      
      
      random.intersections[j] <- nrow(intersectBed(peaks.set1 = peaks1, peaks.set2 = random.peaks2 )) 
      
      
    }
    
    p.value <- sum(random.intersections > nrow(real.intersection)) / number.randomisation
    if( p.value == 0)
    {
      p.value <- 1/number.randomisation
    }
    
    colnames(real.intersection) <- c("chromosome", "start", "end")
    granges.intersection <- makeGRangesFromDataFrame(real.intersection,
                                                     keep.extra.columns=FALSE,
                                                     ignore.strand=FALSE,
                                                     seqinfo=NULL,
                                                     seqnames.field="chromosome",
                                                     start.field="start",
                                                     end.field="end",
                                                     starts.in.df.are.0based=FALSE)
    
    
    
    peakAnno <- annotatePeak(granges.intersection, tssRegion=c(-2000, 2000),
                             TxDb=txdb, annoDb="org.At.tair.db")
    
    annot.peaks <- as.data.frame(peakAnno)
    target.genes <- subset(annot.peaks, distanceToTSS >= 2000 | distanceToTSS <= -2000)$geneId
    target.genes <- paste(target.genes, collapse = ",")
    
    bed.intersections[i,1] <- strsplit(x = as.character(combinations[i,1]), split = "_peaks")[[1]][1]
    bed.intersections[i,2] <- strsplit(x = as.character(combinations[i,2]), split = "_peaks")[[1]][1]
    bed.intersections[i,3] <- p.value
    bed.intersections[i,5] <- nrow(real.intersection)
    bed.intersections[i,6] <- target.genes
    
  } else 
  {
    bed.intersections[i,1] <- strsplit(x = as.character(combinations[i,1]), split = "_peaks")[[1]][1]
    bed.intersections[i,2] <- strsplit(x = as.character(combinations[i,2]), split = "_peaks")[[1]][1]
    bed.intersections[i,3] <- NA
    bed.intersections[i,5] <- "No intersection"
    bed.intersections[i,6] <- NA
  }
  
  
  
}

write.table(bed.intersections, file = "bed_intersections.txt", sep = "\t", row.names = FALSE)
# Stop the clock
proc.time() - ptm

library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)
library(ChIPseeker)

txdb <- TxDb.Athaliana.BioMart.plantsmart28
colnames(real.intersection) <- c("chromosome", "start", "end")
granges.intersection <- makeGRangesFromDataFrame(real.intersection,
                                                 keep.extra.columns=FALSE,
                                                 ignore.strand=FALSE,
                                                 seqinfo=NULL,
                                                 seqnames.field="chromosome",
                                                 start.field="start",
                                                 end.field="end",
                                                 starts.in.df.are.0based=FALSE)



peakAnno <- annotatePeak(granges.intersection, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.At.tair.db")

annot.peaks <- as.data.frame(peakAnno)
target.genes <- subset(annot.peaks, distanceToTSS >= 2000 | distanceToTSS <= -2000)$geneId
target.genes <- paste(target.genes, collapse = ",")
