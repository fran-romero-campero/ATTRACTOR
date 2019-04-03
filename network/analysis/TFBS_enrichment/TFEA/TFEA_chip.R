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
bed <- read.table(file = "CCA1_peaks.narrowPeak")
head(bed)

#Convert to narrowpeak format (ten columns)
# extra.columns <- matrix(nrow = nrow(bed), ncol = 7)
# new.bed <- cbind(bed, extra.columns)
# write.table(new.bed, file="bed_files/CCA1_new_peaks.narrowPeak", sep = "\t", 
#             col.names = FALSE, row.names = FALSE)


#txt2GR is the function to convert a file.txt with peaks to GRanges object
txt2GR(fileTable = bed, format = "narrowpeak", fileMetaData = myMetadata)

folder<-"../bed_files/"
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

#######---Step4: Sustitute the default database by a custom generated table ---#####

# At the beginning of a session, use the function set_user_data 
# to use your TFBS binary matrix and metadata table with the rest of the package.

set_user_data(binary_matrix = myTFBSmatrix, metadata = myMetadata)
help("set_user_data")

#########################################################################
###########Analysis of the association of TFBS and any condition#########
#########################################################################

# We must provide a list of genes are considered differentially induced (ANY CONDITION, for example,
# genes peaking at ZT0) and a list of control genes whose expression is not altered 
# in the analyzed experiment (The control genes would be the entire network)

# Crea una tabla de contingencia para cada TF de nuestra base de datos. Como ésta:

#*************************
#         | TFbound_yes |  TFbound_no
# -------------------------------------
# my_list |   number    |   number
# -------------------------------------
# control|    number    |   number
#--------------------------------------

# Entonces aplica el test de FIsher a cada tabla de contingencia para testear si la hipótesis
# nula de que la unión de TF y la condición (por ejemplo que pique a ZT0) es independiente. 
# Devuelve los p-valores y los FDRs.


zt0.peak <- read.table("../../../clusters/peak_ZT0_trough_ZT12.txt", as.is = T)[[1]]
length(zt0.peak)
network.data <- read.table(file="../../../../web_apps/attractor_dev/data/attractor_network_representation.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "")
head(network.data)

all.genes <- network.data$names
length(all.genes)

my.list <- intersect(zt0.peak, all.genes)
length(my.list)

control.list <- setdiff(all.genes, my.list)
length(control.list) 

CM.list <- contingency_matrix(test_list = my.list, 
                              control_list = control.list, chip_index = myMetadata)


pvalues <- getCMstats(CM.list)
head(pvalues)

plot_CM(pvalues)  #plot p-values against ORs (odd-ratios)

##Hay que buscar log10(adj.pval) altos y odds distintos de uno. 
rows.to.keep <- subset(x = pvalues, subset = adj.p.value < 0.001)


######Bucle para testear todos los conjuntos de genes (clusters)####

clusters.folder<-"../../../clusters/"
cluster.list<-dir(clusters.folder)

for (i in 1:length(cluster.list))
{
  gene.list <- read.table(paste0("../../../clusters/", cluster.list[i]), as.is = T)[[1]]
  my.list <- intersect(gene.list, all.genes)
  control.list <- setdiff(all.genes, my.list)
  CM.list <- contingency_matrix(test_list = my.list, 
                                control_list = control.list, chip_index = myMetadata)
  
  
  pvalues <- getCMstats(CM.list)
  rows.to.keep <- subset(x = pvalues, subset = adj.p.value < 0.001)
  
  split.name <- strsplit(x = cluster.list[i], split = ".txt")[[1]]
  write.table(rows.to.keep, 
              file=paste0("TFEA_tests/TFBS_enrichment_",split.name),
              sep="\t", row.names = FALSE)
  
  
}



# An odds ratio is a relative measure of effect, which allows the 
# comparison of the intervention group of a study relative to the 
# comparison or placebo group. So if the outcome is the same in both 
# groups the ratio will be 1, which implies there is no difference 
# between the two arms of the study.

#' #NOTE: Since the contingency matrix function does not work for me, I edit it to work
#' # Basically, i replace the lapply for a loop and remove a function
#' 
#' test_list <- my.list
#' control_list <- control.list
#' chip_index <- myMetadata
#' Mat01 <- myTFBSmatrix
#' 
#' 
#' 
#' contingency_matrix_pedro <- function(test_list, control_list,
#'                                chip_index = get_chip_index()) {
#'   
#'   #' @title Computes 2x2 contingency matrices
#'   #' @description Function to compute contingency 2x2 matrix by the partition
#'   #' of the two gene ID lists according to the presence or absence of the
#'   #' terms in these list in a ChIP-Seq binding database.
#'   #' @param test_list List of gene Entrez IDs
#'   #' @param control_list If not provided, all human genes not present in
#'   #' test_list will be used as control.
#'   #' @param chip_index Output of the function “get_chip_index”, a data frame
#'   #' containing accession IDs of ChIPs on the database and the TF each one
#'   #' tests. If not provided, the whole internal database will be used
#'   #' @return List of contingency matrices, one CM per element in chip_index
#'   #' (i.e. per ChIP-seq dataset).
#'   #' @export contingency_matrix
#'   #' @examples
#'   #' data('Genes.Upreg',package = 'TFEA.ChIP')
#'   #' CM_list_UP <- contingency_matrix(Genes.Upreg)
#'   
#'   if (!exists("Mat01")) {
#'     Mat01 <- NULL
#'     data("Mat01", package = "TFEA.ChIP", envir = environment())
#'   }
#'   if (missing(control_list)) {
#'     # Generating control gene list in case is not provided.
#'     control_list <- rownames(Mat01)
#'   }
#'   
#'   control_list <- control_list[!(control_list %in% test_list)]
#'   
#'   Matrix1 <- Mat01[rownames(Mat01) %in% test_list, colnames(Mat01) %in%
#'                      chip_index$Accession]
#'   Matrix2 <- Mat01[rownames(Mat01) %in% control_list, colnames(Mat01) %in%
#'                      chip_index$Accession]
#'   
#'   contMatrix_list <- lapply(1:nrow(myMetadata), matrix, data= NA, nrow=2, ncol=2) #Create an empty list of matrices
#'   # contMatrix_list <- vector("list", nrow(myMetadata))
#'   for (i in 1:(nrow(myMetadata)-1))
#'   {
#'     
#'       chip.vector1 <- Matrix1[, chip_index$Accession[i] ]
#'       chip.vector2 <- Matrix2[, chip_index$Accession[i] ]
#'       
#'       pos1 <- sum( chip.vector1 == 1 )
#'       pos2 <- sum(chip.vector2 == 1 )
#'       neg1 <- sum( chip.vector1 == 0 )
#'       neg2 <- sum( chip.vector2 == 0 )
#'       
#'       contMatrix <- cbind(c(pos1, pos2), c(neg1, neg2))
#'       rownames(contMatrix) <- c("Test", "Control")
#'       colnames(contMatrix) <- c("Positive", "Negative")
#'       contMatrix_list[[i]] <- contMatrix
#'   }
#' 
#'   names(contMatrix_list) <- as.character(chip_index$Accession)
#'   return(contMatrix_list)
#' }
#' 
#' # El bucle lo hace bien pero falla en el TF número 20 (TOC1) ¿Por qué?? CORREGIIIR
#' 
#' CM.list <- contingency_matrix_pedro(test_list = my.list, 
#'                               control_list = control.list, chip_index = myMetadata)
#' 
#' pvalues <- getCMstats(CM.list)
#'     
#'     
#' 
#'   
