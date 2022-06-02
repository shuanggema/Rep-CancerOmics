options(stringsAsFactors = FALSE)

######## input ########
Final_gene_data <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_BRCA/Gene/Final_gene_data.csv",
                            row.names = 1)

######## output ########
# File - "EnsemblID_EntrezID_Pathway.csv": column names: Entrez_ID, Ensembl_ID, Pathway_ID, Description, GeneRatio, BgRatio
# File - "gene_list_full.csv": gene list with a linked pathway


########################## Get gene related pathways #######################
# related library
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")

library("AnnotationDbi")
library("org.Hs.eg.db")
library(clusterProfiler)

#columns(org.Hs.eg.db) # returns list of available keytypes

############## Get Ensembl ID linked Entrez ID ###########
EnsemblID_TO_EntrezID_list = mapIds(org.Hs.eg.db,
                                    keys=as.matrix(Final_gene_data$Ensembl_ID), #Column containing Ensembl gene ids
                                    column="ENTREZID",
                                    keytype="ENSEMBL",
                                    multiVals="list")

EnsemblID_TO_EntrezID_df <- data.frame(Ensembl_ID = names(unlist(EnsemblID_TO_EntrezID_list)),
                                       Entrez_ID = unlist(EnsemblID_TO_EntrezID_list))

EnsemblID_TO_EntrezID_df <- EnsemblID_TO_EntrezID_df[!is.na(EnsemblID_TO_EntrezID_df$Entrez_ID),]

############## Get Entrez ID linked pathways ###########
# get the linked matrix for gene-pathway
# will get the newest version from https://www.genome.jp/kegg/mapper.html

KEGG_EntrezID <- enrichKEGG(
  gene = unique(EnsemblID_TO_EntrezID_df$Entrez_ID),
  keyType = "kegg",
  organism = "hsa",
  pvalueCutoff = 1,
  pAdjustMethod = "none",
  qvalueCutoff=1
)

KEGG_EntrezID <- as.data.frame(KEGG_EntrezID)

###
new <- NA
for (i in 1:nrow(KEGG_EntrezID)){
  old <- unlist(strsplit(KEGG_EntrezID[i,"geneID"], "/"))
  for (id in old) {
    new <- rbind(new, cbind(id,KEGG_EntrezID[i,c("ID","Description","GeneRatio","BgRatio")]))
  }
}
new <- new[-1,]

colnames(new) <- c("Entrez_ID","Pathway_ID","Description","GeneRatio","BgRatio")
EnsemblID_EntrezID_Pathway <- merge(EnsemblID_TO_EntrezID_df, new, by = "Entrez_ID")
write.csv(EnsemblID_EntrezID_Pathway,file = "/Users/name/Desktop/TCGA-2021/Data/EnsemblID_EntrezID_Pathway.csv")

############## only retain genes with linked pathways ###########
gene_list_full <- Final_gene_data$Ensembl_ID[Final_gene_data$Ensembl_ID %in% EnsemblID_EntrezID_Pathway$Ensembl_ID]
write.csv(gene_list_full,file = "/Users/name/Desktop/TCGA-2021/Data/gene_list_full.csv")

