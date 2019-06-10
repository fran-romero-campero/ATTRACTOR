library(clusterProfiler)
library(org.At.tair.db)

columns(org.At.tair.db)
select(org.At.tair.db,columns = c("SYMBOL"),keys=keys(org.At.tair.db,keytype = "TAIR"))[["SYMBOL"]]

compute.enrichments <- function(gene.ratios, bg.ratios)
{
  gene.ratios.eval <- sapply(parse(text=gene.ratios),FUN = eval)
  bg.ratios.eval <- sapply(parse(text=bg.ratios),FUN = eval)
  enrichments <- round(x=gene.ratios.eval/bg.ratios.eval,digits = 2)
  enrichments.text <- paste(enrichments, " (", gene.ratios, "; ", bg.ratios, ")",sep="")
  
  return(enrichments.text)  
}

#genes <- read.table(file="genes_CCA1_PRR5_PIF5.txt",as.is=T)[[1]]
genes <- read.table(file="genes_CCA1_PRR5_PIF5.txt",as.is=T)[[1]]

length(genes)

atha.universe <- unique(AnnotationDbi::select(org.At.tair.db,columns = c("GO"),keys=keys(org.At.tair.db,keytype = "TAIR"))[["TAIR"]])
length(atha.universe)

ego <- enrichGO(gene          = genes,
                universe      = atha.universe,
                OrgDb         = org.At.tair.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = FALSE,
                keyType = "TAIR")


barplot(ego,drop=TRUE,showCategory = 20)
goplot(x=ego,showCategory = 10)
dotplot(ego,showCategory=20)
emapplot(ego,showCategory=20)
cnetplot(ego,showCategory=20)

## Generate ouput table
enrich.go.result <- as.data.frame(ego)
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

head(go.result.table)  
write.table(x = go.result.table,file = "hubs_go.txt",row.names = F,quote = F,sep = "\t")

pathway.enrichment <- as.data.frame(enrichKEGG(gene = genes, 
                                               organism = "ath", universe = atha.universe,
                                               qvalueCutoff = 0.05))


write.table(x = pathway.enrichment,file = "hubs_pathways.txt",row.names = F,col.names=F,quote = F,sep = "\t")

genes.pathway <- rep(0,length(atha.universe))
names(genes.pathway) <- atha.universe

genes.pathway[genes] <- 1

pathview(gene.data = sort(genes.pathway,decreasing = TRUE), pathway.id =pathway.enrichment$ID[6], species = "ath",
         limit = list(gene=max(abs(genes.pathway)), cpd=1),gene.idtype ="kegg")
head(pathway.enrichment)
