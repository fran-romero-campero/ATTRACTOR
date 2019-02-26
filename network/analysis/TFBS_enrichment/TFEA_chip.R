#Load package
library(TFEA.ChIP)

#Read the metadata table
# A Metadata table (storing at least, Accession ID, name of the file, 
# and TF tested in the ChIP-Seq experiment). The metadata table included 
# with this package has the following fields: “Name”, “Accession”, “Cell”, 
# “Cell Type”, “Treatment”, “Antibody”, and “TF”.
myMetadata <- read.table(file = "metadata.txt", header = T)
head(myMetadata)
# bed <- read.table(file = "bed_files/PIF3_peaks.narrowPeak")
# head(bed)


txt2GR(fileTable = bed, format = "narrowpeak", fileMetaData = myMetadata)

folder<-"bed_files/"
File.list<-dir(folder)
format<-"narrowpeak"

gr.list <- lapply(
  seq_along(File.list),
  function( File.list, myMetaData, format ){
    
    tmp<-read.table( file = paste0(folder,File.list[i]), stringsAsFactors = FALSE )
    
    file.metadata <- myMetadata[ myMetadata$Name == File.list[i], ]
    
    ChIP.dataset.gr<-txt2GR(tmp, format, file.metadata)
    
    return(ChIP.dataset.gr)
  },
  File.list = File.list,
  myMetadata = myMetadata,
  format = format
)
