options(stringsAsFactors = FALSE)

######## input ########
data_before_sample_filter <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUSC/Gene/GeneID_EnsebmleID_barcode_log(TPM+1)_data_protein.csv",
                                      row.names = 1)

######## output ########
# File - "Final_gene_data.csv": log(TPM+1) dataframe, colname: barcode (delete normal sample, duplicated sample/portion), 
#                               first row: EnsemblID without version number,
#                               second row: GeneID,
#                               third row: Gene type, protein encoding only

############# sample filtering ################
# only retain primary solid tumor sample: sample 1
retain_barcode_list <- colnames(data_before_sample_filter[,-c(1:3)])
retain_barcode_list <- retain_barcode_list[substr(retain_barcode_list, start = 14, stop = 15) %in%
                                             c("01")]
# "01","02","03","04","05","06","07","08","09","10"

# only retain first (A>B>C...) sample of the sequence samples if one subject has more than one sample: vial A
subject_list <- substr(retain_barcode_list, start = 9, stop = 12)
# whether one subject has more than one sample
table(subject_list)[table(subject_list) > 1]
subject_list_more_than_one_sample <- names(table(subject_list)[table(subject_list) > 1])

delete_barcode_list <- NA
for (subject in subject_list_more_than_one_sample){
  delete_barcode_list <- c(delete_barcode_list, retain_barcode_list[substr(retain_barcode_list, start = 9, stop = 12) == subject])
}
delete_barcode_list <- delete_barcode_list[-1]

# only retain the first portion of the sample: portion 1
delete_barcode_list <- c("TCGA.21.1076.01A.02R.0692.07")

retain_barcode_list <- retain_barcode_list[!retain_barcode_list %in% delete_barcode_list]

############# finally we have 501 subjects ################
Final_gene_data <- data.frame(data_before_sample_filter[,1:3],data_before_sample_filter[,retain_barcode_list])
write.csv(Final_gene_data, "/Users/name/Desktop/TCGA-2021/Data/TCGA_LUSC/Gene/Final_gene_data.csv")




