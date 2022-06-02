options(stringsAsFactors = TRUE)

######## input ########
Final_data <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Merged/Final_data.csv",
                            row.names = 1)
gene_list_full <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/gene_list_full.csv",row.names = 1))

# write.csv(sample(1000:9999,100), "/Users/name/Desktop/TCGA-2021/Data/seed_list.csv")
seed_list <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/seed_list.csv",row.names = 1))

setwd("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Significant Gene")

variable_list <- c("age_at_diagnosis","race","gender","prior_malignancy","ajcc_pathologic_stage",
                   "radiation_therapy","pharmaceutical_therapy","cigarettes_per_day")

survival_formula <- "age_at_diagnosis + race + gender + prior_malignancy + ajcc_pathologic_stage + radiation_therapy + pharmaceutical_therapy + cigarettes_per_day"

number_var <- 9 # the dimensions of all covariates

######## output ########
# File - "Final_data_selected.csv": gene expression [unit: log(TPM+1)] and demo and survival dataframe, 
#                          only gene that has variance > 1st quantile,
#                          colname: barcode (delete normal sample, duplicated sample/portion), demo variables, and outcomes
#                          first col: EnsemblID without version number,
#                          second col: GeneID,
#                          third col: Gene type, protein encoding only

# File - "gene_list_subset.csv": delete genes with variance <= 1st Quantile
# File - "Selected_gene_first_half.csv"
# File - "Selected_gene_last_half.csv"

########
library(clusterSim)
library(survival)

############# delete genes with variance <= 1st Quantile ################
# only retain gene that has different Ensembl ID (delete repeated gene expression for the same Ensembl ID)
table(gene_list_full)[table(gene_list_full) >= 2]
gene_list_full <- unique(gene_list_full)

gene_variance <- apply(as.matrix(Final_data[,gene_list_full]), 2, var)
gene_list_subset <- gene_list_full[!gene_variance <= quantile(gene_variance, 0.25)]
Final_data_selected <- Final_data[,c(gene_list_subset,variable_list,"time","status")]

write.csv(gene_list_subset,"./gene_list_subset.csv")
write.csv(Final_data_selected,"./Final_data_selected.csv")


############# conduct survival analysis for the remaining gene ################
# standardize all the gene expression
Final_data_selected[,gene_list_subset] <- data.Normalization(Final_data_selected[,gene_list_subset],type="n1",normalization="column")

number_obs <- nrow(Final_data_selected)
ind <- 1:number_obs

Selected_gene_first_half <- list()
Selected_gene_last_half <- list()

i <- 1
for (seed in seed_list){
  set.seed(seed)
  ind_half <- sample(ind, round(number_obs/2,0), replace = FALSE)
  Final_data_selected_first_half <- Final_data_selected[ind_half,]
  Final_data_selected_last_half <- Final_data_selected[-ind_half,]
 
  ### for the first half of samples 
  p_value <- NA
  j <- 1
  Final_data_selected_first_half <- rbind(Final_data_selected_first_half,Final_data_selected_first_half)
  for(gene in gene_list_subset){
    model <- coxph(as.formula(paste("Surv(time, status) ~",gene,"+",survival_formula)), 
                   Final_data_selected_first_half)
    p_value[j] <- summary(model)$coefficients[4*(number_var+1)+1]
    j <- j+1
  }
  # p_value_adjust <- p.adjust(p_value, method = "bonferroni", n = length(p_value)) 
  p_value_adjust <- p.adjust(p_value, method = "fdr", n = length(p_value)) 
  Selected_gene_first_half[[i]] <- gene_list_subset[p_value_adjust < 0.1]
  
  ### for the last half of samples 
  p_value <- NA
  j <- 1
  Final_data_selected_last_half <- rbind(Final_data_selected_last_half,Final_data_selected_last_half)
  for(gene in gene_list_subset){
    model <- coxph(as.formula(paste("Surv(time, status) ~",gene,"+",survival_formula)), 
                   Final_data_selected_last_half)
    p_value[j] <- summary(model)$coefficients[4*(number_var+1)+1]
    j <- j+1
  }
  # p_value_adjust <- p.adjust(p_value, method = "bonferroni", n = length(p_value))
  p_value_adjust <- p.adjust(p_value, method = "fdr", n = length(p_value))
  Selected_gene_last_half[[i]] <- gene_list_subset[p_value_adjust < 0.1]
  
  i <- i+1
}

mean(c(lengths(Selected_gene_first_half),lengths(Selected_gene_last_half)))
sd(c(lengths(Selected_gene_first_half),lengths(Selected_gene_last_half)))

############# store the list as dataframe  ################
Selected_gene_first_half_df <- data.frame(ind = 1:100, gene = rep(NA,100))
Selected_gene_last_half_df <- data.frame(ind = 1:100, gene = rep(NA,100))

for (i in 1:100){
  temp <- NA
  for (j in Selected_gene_first_half[[i]]){
    if (is.na(temp)){
      temp <- j
    }else{
      temp <- paste(temp, "/", j, sep = "")
      }
  }
  Selected_gene_first_half_df[i,"gene"] <- temp 
}

for (i in 1:100){
  temp <- NA
  for (j in Selected_gene_last_half[[i]]){
    if (is.na(temp)){
      temp <- j
    }else{
      temp <- paste(temp, "/", j, sep = "")
    }
  }
  Selected_gene_last_half_df[i,"gene"] <- temp 
}

write.csv(Selected_gene_first_half_df,"./Selected_gene_first_half.csv")
write.csv(Selected_gene_last_half_df,"./Selected_gene_last_half.csv")
