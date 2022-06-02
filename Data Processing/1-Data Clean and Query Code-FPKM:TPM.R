options(stringsAsFactors = FALSE)

######## input ########
main_path <- "/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Gene"
gene_path <- "/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Gene/gdc_download_20210918_131706.519786"
metadata_file <- "/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Gene/metadata.cart.2021-09-18.json"
GeneID_EensemblID_file <- "/Users/name/Desktop/TCGA-2021/Data/gencode.v38.annotation.gtf"

######## output ########
# Folder - "SampleFiles": store upzip TCGA raw gene data
# File - "EnsebmleID_barcode_FPKM_data.csv": FPKM dataframe, colname: barcode, 
#                                            first row: EnsemblID without version number
# File - "Ensembl_ID_TO_Genename.csv": store mapping from EnsemblID to Genename
# File - "GeneID_EnsebmleID_barcode_log(TPM+1)_data_protein.csv": log(TPM+1) dataframe, colname: barcode, 
#                                            first row: EnsemblID without version number,
#                                            second row: GeneID,
#                                            third row: Gene type, protein encoding only

########
library(R.utils)
library(rjson)

################### Unzip raw gene data #######################
# move all files into same folder "SampleFiles"
setwd(main_path)
dir.create("SampleFiles")

# set file path to downloaded data (cart)
filepath <- dir(path = gene_path, full.names = TRUE)

for (wd in filepath) {
  files <- dir(path = wd, pattern = "gz$")
  fromfilepath <- paste(wd, "/", files, sep = "")
  tofilepath <- paste("./SampleFiles/", files, sep = "")
  file.copy(fromfilepath,tofilepath)
}

# unzip all files and delete original files
setwd("./SampleFiles")
countsFiles <- dir(path = "./", pattern = "gz$")

sapply(countsFiles, gunzip)

################### Barcode & Case ID #######################
# deal with json files
metadata_json_File <- fromJSON(file = metadata_file)
json_File_Info <- data.frame(filesName = c(), TCGA_Barcode = c())
for (i in 1:length(metadata_json_File)) {
  TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name <- metadata_json_File[[i]][["file_name"]]
  json_File_Info <- rbind(json_File_Info,data.frame(filesName = file_name,TCGA_Barcode = TCGA_Barcode))
}

rownames(json_File_Info) <- json_File_Info[,1]
filesName_To_TCGA_BarcodeFile <- json_File_Info[-1]


########################## Gene data matrix #######################
countsFileNames <- dir(pattern="FPKM.txt$") # OR function: list.files()

allSampleRawCounts <- data.frame()
for(txtFile in countsFileNames){
  # each loop = one file
  SampleCounts <- read.table(txtFile, header = FALSE)
  rownames(SampleCounts) <- as.matrix(SampleCounts[1])
  SampleCounts <- SampleCounts[-1]
  # according to caseID in filesName_To_TCGA_BarcodeFile, find column name (barcode)
  colnames(SampleCounts) <- filesName_To_TCGA_BarcodeFile[paste(txtFile,".gz",sep = ""),]
  if (dim(allSampleRawCounts)[1]== 0){
    allSampleRawCounts <- SampleCounts
  }
  else
  {allSampleRawCounts<- cbind(allSampleRawCounts,SampleCounts)}
}


Ensembl_ID <- substr(row.names(allSampleRawCounts),1,15) # EnsemblID without version number
EnsebmleID_barcode_FPKM_data <- cbind(Ensembl_ID,allSampleRawCounts)
rownames(EnsebmleID_barcode_FPKM_data) <- 1:nrow(EnsebmleID_barcode_FPKM_data) # to avoid duplicated row names
write.csv(EnsebmleID_barcode_FPKM_data, "../EnsebmleID_barcode_FPKM_data.csv")

########################## Ensemble ID & Gene ID ##########################
# function: by gtf file, get mapping relationship between EnsemblID and gene names
get_map = function(input) {
  input = data.table::fread(input, header = FALSE)

  input = input[input[[3]] == "gene", ]

  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  pattern_type = ".*gene_type \"([^;]+)\";.*"

  gene_id = sub(pattern_id, "\\1", input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  gene_type = sub(pattern_type, "\\1", input[[9]])

  Ensembl_ID_TO_Genename <- data.frame(gene_id = gene_id,
                                       gene_name = gene_name,
                                       gene_type = gene_type,
                                       stringsAsFactors = FALSE)
  return(Ensembl_ID_TO_Genename)
}

Ensembl_ID_TO_Genename <- get_map(GeneID_EensemblID_file)

gtf_Ensembl_ID <- substr(Ensembl_ID_TO_Genename[,1],1,15) # delete version number of EnsemblID
Ensembl_ID_TO_Genename <- data.frame(Ensembl_ID = gtf_Ensembl_ID,
                                     gene_id = Ensembl_ID_TO_Genename[,2],
                                     gene_type = Ensembl_ID_TO_Genename[,3])

write.csv(Ensembl_ID_TO_Genename,file = "../Ensembl_ID_TO_Genename.csv")

########################## merge data ##########################
Ensembl_ID_TO_Genename_protein <- Ensembl_ID_TO_Genename[Ensembl_ID_TO_Genename$gene_type == "protein_coding",]
GeneID_EnsebmleID_barcode_FPKM_data_protein <- merge(Ensembl_ID_TO_Genename_protein,EnsebmleID_barcode_FPKM_data,by="Ensembl_ID")

########################## final data ##########################
# log(TPM+1), protein_coding gene information
GeneID_EnsebmleID_barcode_FPKM_data_protein[,-c(1:3)] <- log(GeneID_EnsebmleID_barcode_FPKM_data_protein[,-c(1:3)]*1000000/
  colSums(GeneID_EnsebmleID_barcode_FPKM_data_protein[,-c(1:3)])+1)

write.csv(GeneID_EnsebmleID_barcode_FPKM_data_protein,file = "../GeneID_EnsebmleID_barcode_log(TPM+1)_data_protein.csv")







