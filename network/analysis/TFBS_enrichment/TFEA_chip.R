#Load package
library(TFEA.ChIP)


######---Step1: Filter peaks from source and store them as a GRanges object---###########

#Read the metadata table
# A Metadata table (storing at least, Accession ID, name of the file, 
# and TF tested in the ChIP-Seq experiment). The metadata table included 
# with this package has the following fields: “Name”, “Accession”, “Cell”, 
# “Cell Type”, “Treatment”, “Antibody”, and “TF”.
myMetadata <- read.table(file = "metadata.txt", header = T)
head(myMetadata)
bed <- read.table(file = "bed_files/CCA1_peaks.narrowPeak")
head(bed)

#Convert to narrowpeak format (ten columns)
# extra.columns <- matrix(nrow = nrow(bed), ncol = 7)
# new.bed <- cbind(bed, extra.columns)
# write.table(new.bed, file="bed_files/CCA1_new_peaks.narrowPeak", sep = "\t", 
#             col.names = FALSE, row.names = FALSE)


#txt2GR is the function to convert a file.txt with peaks to GRanges object
txt2GR(fileTable = bed, format = "narrowpeak", fileMetaData = myMetadata)

folder<-"bed_files/"
File.list<-dir(folder)
format<-"narrowpeak"

#This function read the bed file, select the info from myMetaData and apply the txt2GR function
txt2GR.fun <- function( File.list, myMetaData, format ){
  
  tmp<-read.table( file = paste0(folder,File.list[i]), stringsAsFactors = FALSE )
  
  file.metadata <- myMetadata[ myMetadata$Name == File.list[i], ]
  
  ChIP.dataset.gr<-txt2GR(tmp, format, file.metadata)
  
  return(ChIP.dataset.gr)
}

#For loop
gr.list <- list()
for (i in 1:length(File.list))
{
  gr.list[i] <- txt2GR.fun(File.list = File.list, myMetaData = myMetadata, format="narrowpeak")

}

#######---Step2: Assign TFBS peaks from ChIP dataset to specific genes---#####

library(GSEABase)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
atha.genes <- genes(txdb)

TF.gene.binding.db <- GR2tfbs_db(atha.genes, gr.list, distanceMargin = 0) #asignation of peaks to genes
str(TF.gene.binding.db)

#######---Step3: Generation of the TFBS database ---#####

# The function makeTFBSmatrix generates a binary matrix with a row for each gene 
# and a column for each independent ChIPseq dataset. The cell of the matrix contain
# a value of 1 if the gene is bound by the TF and 0 otherwise. 
# This matrix and the metadata table can be used instead of the default TFBS database.

gen.list <- genes(txdb)$gene_id # selecting all the genes in knownGene

myTFBSmatrix <- makeTFBSmatrix(gen.list,TF.gene.binding.db)
myTFBSmatrix[2530:2531,] # The gene AT1G23080 has TFBS for this five ChIP-Seq datasets
