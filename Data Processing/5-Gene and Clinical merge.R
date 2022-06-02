options(stringsAsFactors = TRUE)

######## input ########
Final_gene_data <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Gene/Final_gene_data.csv",
                            row.names = 1)
Final_demo <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Clinic/Final_demo.csv",
                       row.names = 1)
variable_list <- c("age_at_diagnosis","race","gender","prior_malignancy","ajcc_pathologic_stage",
                   "radiation_therapy","pharmaceutical_therapy","cigarettes_per_day")

setwd("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Merged")

######## output ########
# File - "Final_data.csv": gene expression [unit: log(TPM+1)] and demo and survival dataframe, 
#                          colname: barcode (delete normal sample, duplicated sample/portion), demo variables, and outcomes
#                          first col: EnsemblID without version number,
#                          second col: GeneID,
#                          third col: Gene type, protein encoding only

# File - "x_demo.csv": demo variables matrix, based on the input variable list order
# File - "x_gene.csv": gene expression [unit: log(TPM+1)] matrix
# File - "x_survival.csv": survival matrix

########
library(fastDummies)

############# data merge - data frame ################
### extract transpose gene matrix
# retain un-duplicated Ensembl ID
Final_gene_data <- Final_gene_data[!duplicated(Final_gene_data$Ensembl_ID),]
gene_expression <- Final_gene_data[,4:ncol(Final_gene_data)]
gene_list <- Final_gene_data$Ensembl_ID
rownames(gene_expression) <- gene_list
Final_gene_data <- data.frame(t(as.matrix(gene_expression)))
Final_gene_data$subject_id <- substr(rownames(Final_gene_data), start = 9, stop = 12)

### redefine subjecte id for demo data
Final_demo$subject_id <- substr(Final_demo$case_submitter_id, start = 9, stop = 12)
Final_demo <- Final_demo[,c("subject_id",variable_list,"time","status")]

### merge gene and demo data
Final_data <- merge(Final_gene_data, Final_demo, by = "subject_id")

write.csv(Final_data, "./Final_data.csv")

############# data merge - matrix ################
x_demo <- NA
for (var in variable_list){
  if(class(Final_data[,var]) == "numeric"){
    x_demo <- cbind(x_demo,Final_data[,var])
  }
  
  if(class(Final_data[,var]) == "factor"){
    if(length(levels(Final_data[,var])) == 2){
      x_demo <- cbind(x_demo,as.numeric(Final_data[,var])-1)
    } else{
      x_demo <- cbind(x_demo,dummy_cols(Final_data[,var], remove_first_dummy = TRUE)[-1])
    }
  }
}
x_demo <- as.matrix(x_demo[,-1])
x_gene <- as.matrix(Final_data[,gene_list])
x_survival <- as.matrix(Final_data[,c("time","status")])

write.csv(x_demo, "./x_demo.csv")
write.csv(x_gene, "./x_gene.csv")
write.csv(x_survival, "./x_survival.csv")



