options(stringsAsFactors = FALSE)

######## input ########
Final_data <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Merged/Final_data.csv",
                       row.names = 1)
gene_list_subset <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Significant Gene/gene_list_subset.csv",row.names = 1))

Selected_gene_first_half <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Significant Gene/Selected_gene_first_half.csv",
                                     row.names = 1)
Selected_gene_last_half <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Significant Gene/Selected_gene_last_half.csv",
                                    row.names = 1)

EnsemblID_EntrezID_Pathway <- read.csv("/Users/name/Desktop/TCGA-2021/Data/EnsemblID_EntrezID_Pathway.csv", row.names = 1)

setwd("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Significant Gene/")

######## output ########
# File - "Overlap_measure.csv": colname: no_genes_first, no_genes_last; 
#                               measure1, measure1_expect, measure1_baseline, measure1_025, measure1_950, measure1_975; 
#                               measure2, measure2_baseline, measure2_025, measure2_950, measure2_975; 
#                               no_path_first, no_path_last; 
#                               measure3, measure3_baseline, measure3_025, measure3_950, measure3_975; 


############# FIRST HALF VS LAST HALF ################
### How many significant genes are selected

# restore data frame back into list
Selected_gene_first_half_list <- list()
Selected_gene_last_half_list <- list()

for (i in 1:100){
  Selected_gene_first_half_list[[i]] <- strsplit(Selected_gene_first_half2[i,"gene"], "/")[[1]]
  Selected_gene_last_half_list[[i]] <- strsplit(Selected_gene_last_half2[i,"gene"], "/")[[1]]
}


####################################################################
################### gene name measure ##############################
####################################################################
# initialize vectors
OM1 <- NA
OM1_expect <- NA
OM1_baseline <- NA
OM1_baseline_025 <- NA
OM1_baseline_950 <- NA
OM1_baseline_975 <- NA

for (i in 1:100){
  first <- Selected_gene_first_half_list[[i]]
  last <- Selected_gene_last_half_list[[i]]
  
  if (sum(is.na(first),is.na(last)) > 0){
    OM1[i] <- 0
    OM1_expect[i] <- 0
    OM1_baseline[i] <- 0
    OM1_baseline_025[i] <- 0
    OM1_baseline_950[i] <- 0
    OM1_baseline_975[i] <- 0
  } else{
    OM1[i] <- sum(first %in% last)/min(length(first),length(last))
  
    # calculate baseline (expected)
    OM1_expect[i] <- max(length(first),length(last))/length(gene_list_subset)
  
    # calculate baseline (from data)
    OM1_baseline_boot <- rep(NA,1000)
    for (t in 1:1000){
      gene_list_temp1 <- sample(gene_list_subset, length(first), replace = FALSE)
      gene_list_temp2 <- sample(gene_list_subset, length(last), replace = FALSE)
    
      OM1_baseline_boot[t] <- sum(gene_list_temp1 %in% gene_list_temp2)/min(length(first),length(last))
    }
  
    OM1_baseline[i] <- mean(OM1_baseline_boot)
    OM1_baseline_025[i] <- quantile(OM1_baseline_boot, 0.025)
    OM1_baseline_950[i] <- quantile(OM1_baseline_boot, 0.95)
    OM1_baseline_975[i] <- quantile(OM1_baseline_boot, 0.975)
  }
}

####################################################################
################### gene correlation measure #######################
####################################################################
# initialize vectors
OM2 <- NA
OM2_baseline <- NA
OM2_baseline_025 <- NA
OM2_baseline_950 <- NA
OM2_baseline_975 <- NA

for (i in 1:100){
  first <- Selected_gene_first_half_list[[i]]
  last <- Selected_gene_last_half_list[[i]]
  
  if(sum(is.na(first),is.na(last)) > 0){
    OM2[i] <- 0
    OM2_baseline[i] <- 0
    OM2_baseline_025[i] <- 0
    OM2_baseline_950[i] <- 0
    OM2_baseline_975[i] <- 0
  } else{
    gene_cor_comparison <- abs(cor(Final_data[first],Final_data[last]))
    OM2[i] <- (sum(apply(gene_cor_comparison,1,max)) + sum(apply(gene_cor_comparison,2,max)))/(length(first)+length(last))

    # calculate baseline (from data)
    OM2_baseline_boot <- rep(NA,1000)
    for (t in 1:1000){
      gene_list_temp1 <- sample(gene_list_subset, length(first), replace = FALSE)
      gene_list_temp2 <- sample(gene_list_subset, length(last), replace = FALSE)
    
      gene_cor_comparison_temp <- abs(cor(Final_data[gene_list_temp1],Final_data[gene_list_temp2]))

      OM2_baseline_boot[t] <- (sum(apply(gene_cor_comparison_temp,1,max)) +
                                          sum(apply(gene_cor_comparison_temp,2,max)))/(length(first) + length(last))
    }
  
    OM2_baseline[i] <- mean(OM2_baseline_boot)
    OM2_baseline_025[i] <- quantile(OM2_baseline_boot, 0.025)
    OM2_baseline_950[i] <- quantile(OM2_baseline_boot, 0.95)
    OM2_baseline_975[i] <- quantile(OM2_baseline_boot, 0.975)
  }
}


#################################################################
################### gene pathway measure ########################
#################################################################

############# FIRST HALF VS LAST HALF ################
### How many significant genes related pathways are selected

# get gene linked pathway
Selected_pathway_first_half_list <- list()
Selected_pathway_last_half_list <- list()

# store the number of pathways selected
number_pathway_first_half <- NA
number_pathway_last_half <- NA

for (i in 1:100){
  Selected_pathway_first_half_list[[i]] <- EnsemblID_EntrezID_Pathway[EnsemblID_EntrezID_Pathway$Ensembl_ID %in% Selected_gene_first_half_list[[i]],]
  number_pathway_first_half[i] <- length(unique(Selected_pathway_first_half_list[[i]]$Pathway_ID))
  
  Selected_pathway_last_half_list[[i]] <- EnsemblID_EntrezID_Pathway[EnsemblID_EntrezID_Pathway$Ensembl_ID %in% Selected_gene_last_half_list[[i]],]
  number_pathway_last_half[i] <- length(unique(Selected_pathway_last_half_list[[i]]$Pathway_ID))
}


############# pathway measure ################

# initialize vectors
OM3 <- NA
OM3_baseline <- NA
OM3_baseline_025 <- NA
OM3_baseline_950 <- NA
OM3_baseline_975 <- NA

for (i in 1:100){
  first_gene <- Selected_gene_first_half_list[[i]]
  last_gene <- Selected_gene_last_half_list[[i]]
  
  if(sum(is.na(first_gene),is.na(last_gene)) > 0){
    OM3[i] <- 0
    OM3_baseline[i] <- 0
    OM3_baseline_025[i] <- 0
    OM3_baseline_950[i] <- 0
    OM3_baseline_975[i] <- 0
  } else{
    first_path <- Selected_pathway_first_half_list[[i]]
    last_path <- Selected_pathway_last_half_list[[i]]
  
    ####   ####   ####   ####   ####    ####   ####   ####   ####   ####    ####   ####   ####   ####   ####
    unique_pathway_first_half <- table(first_path$Pathway_ID)
    unique_pathway_first_half <- data.frame(Pathway_ID = names(unique_pathway_first_half), count1 = as.numeric(unique_pathway_first_half))
    unique_pathway_first_half$count2 <- NA
    
    j <- 1
    for(id in unique_pathway_first_half$Pathway_ID){
      unique_pathway_first_half$count2[j] <- first_path$GeneRatio[first_path$Pathway_ID == id][1]
      j <- j + 1
    }
  
    unique_pathway_first_half$count2 <- as.numeric(sub("/.*", "", unique_pathway_first_half$count2))
    unique_pathway_first_half$weight <- unique_pathway_first_half$count1/unique_pathway_first_half$count2
  
    ####   ####   ####   ####   ####    ####   ####   ####   ####   ####    ####   ####   ####   ####   ####
    unique_pathway_last_half <- table(last_path$Pathway_ID)
    unique_pathway_last_half <- data.frame(Pathway_ID = names(unique_pathway_last_half), count1 = as.numeric(unique_pathway_last_half))
    unique_pathway_last_half$count2 <- NA
  
    j <- 1
    for(id in unique_pathway_last_half$Pathway_ID){
      unique_pathway_last_half$count2[j] <- last_path$GeneRatio[last_path$Pathway_ID == id][1]
      j <- j + 1
    }
  
    unique_pathway_last_half$count2 <- as.numeric(sub("/.*", "", unique_pathway_last_half$count2))
    unique_pathway_last_half$weight <- unique_pathway_last_half$count1/unique_pathway_last_half$count2
  
    ####   ####   ####   ####   ####   
    ####   ####   ####   ####   ####
    denominator <- 0
    nominator <- 0
    for(q in 1:nrow(unique_pathway_first_half)){
      for(p in 1:nrow(unique_pathway_last_half)){
        if(unique_pathway_last_half$Pathway_ID[p] == unique_pathway_first_half$Pathway_ID[q]){
          nominator <- nominator + unique_pathway_last_half$weight[p]*unique_pathway_first_half$weight[q]
        }
      }
    }
  
    gene_number1 <- nrow(unique_pathway_first_half)
    gene_number2 <- nrow(unique_pathway_last_half)
    if (gene_number1 <= gene_number2){
      denominator <- sort(unique_pathway_first_half$weight, decreasing = T) %*%
        sort(unique_pathway_last_half$weight, decreasing = T)[1:gene_number1]
    } else {
      denominator <- sort(unique_pathway_last_half$weight, decreasing = T) %*%
        sort(unique_pathway_first_half$weight, decreasing = T)[1:gene_number2]
    }
  
    OM3[i] <- nominator/denominator
    
    ####   ####   ####   ####   ####   
    ####   ####   ####   ####   ####  
    # calculate baseline (from data)
    OM3_baseline_boot <- rep(NA,1000)
  
    for (t in 1:1000){
      gene_list_temp1 <- sample(gene_list_subset, length(first_gene), replace = FALSE)
      gene_list_temp2 <- sample(gene_list_subset, length(last_gene), replace = FALSE)
    
      sig_gene_pathway_temp1 <- EnsemblID_EntrezID_Pathway[EnsemblID_EntrezID_Pathway$Ensembl_ID %in% gene_list_temp1,]
      sig_gene_pathway_temp2 <- EnsemblID_EntrezID_Pathway[EnsemblID_EntrezID_Pathway$Ensembl_ID %in% gene_list_temp2,]
    
      ###
      unique_pathway_temp1 <- table(sig_gene_pathway_temp1$Pathway_ID)
      unique_pathway_temp1 <- data.frame(Pathway_ID = names(unique_pathway_temp1), count1 = as.numeric(unique_pathway_temp1))
    
      j <- 1
      for(id in unique_pathway_temp1$Pathway_ID){
        unique_pathway_temp1$count2[j] <- sig_gene_pathway_temp1$GeneRatio[sig_gene_pathway_temp1$Pathway_ID == id][1]
        j <- j + 1
      }
    
      unique_pathway_temp1$count2 <- as.numeric(sub("/.*", "", unique_pathway_temp1$count2))
      unique_pathway_temp1$weight <- unique_pathway_temp1$count1/unique_pathway_temp1$count2
    
      ###
      unique_pathway_temp2 <- table(sig_gene_pathway_temp2$Pathway_ID)
      unique_pathway_temp2 <- data.frame(Pathway_ID = names(unique_pathway_temp2), count1 = as.numeric(unique_pathway_temp2))
    
      j <- 1
      for(id in unique_pathway_temp2$Pathway_ID){
        unique_pathway_temp2$count2[j] <- sig_gene_pathway_temp2$GeneRatio[sig_gene_pathway_temp2$Pathway_ID == id][1]
        j <- j + 1
      }
    
      unique_pathway_temp2$count2 <- as.numeric(sub("/.*", "", unique_pathway_temp2$count2))
      unique_pathway_temp2$weight <- unique_pathway_temp2$count1/unique_pathway_temp2$count2
    
      denominator <- 0
      nominator <- 0
      for(p in 1:nrow(unique_pathway_temp1)){
        for(q in 1:nrow(unique_pathway_temp2)){
          if(unique_pathway_temp2$Pathway_ID[q] == unique_pathway_temp1$Pathway_ID[p]){
            nominator <- nominator + unique_pathway_temp2$weight[q]*unique_pathway_temp1$weight[p]
          }
        }
      }
    
      gene_number1 <- nrow(unique_pathway_temp1)
      gene_number2 <- nrow(unique_pathway_temp2)
      if (gene_number1 <= gene_number2){
        denominator <- sort(unique_pathway_temp1$weight, decreasing = T) %*%
          sort(unique_pathway_temp2$weight, decreasing = T)[1:gene_number1]
      } else {
        denominator <- sort(unique_pathway_temp2$weight, decreasing = T) %*%
          sort(unique_pathway_temp1$weight, decreasing = T)[1:gene_number2]
      }
    
      OM3_baseline_boot[t] <- nominator/denominator
    }
  
    OM3_baseline[i] <- mean(OM3_baseline_boot)
    OM3_baseline_025[i] <- quantile(OM3_baseline_boot, 0.025)
    OM3_baseline_950[i] <- quantile(OM3_baseline_boot, 0.95)
    OM3_baseline_975[i] <- quantile(OM3_baseline_boot, 0.975)
  }
}


no_genes_first <- rep(0,100)
no_genes_first[!is.na(Selected_gene_first_half$gene)] <- lengths(Selected_gene_first_half_list)[!is.na(Selected_gene_first_half$gene)]

no_genes_last <- rep(0,100)
no_genes_last[!is.na(Selected_gene_last_half$gene)] <- lengths(Selected_gene_last_half_list)[!is.na(Selected_gene_last_half$gene)]


Overlap_measure <- data.frame(no_genes_first = no_genes_first, no_genes_last = no_genes_last,
                              measure1 = OM1, measure1_expect = OM1_expect, measure1_baseline = OM1_baseline, 
                              measure1_025 = OM1_baseline_025, measure1_950 = OM1_baseline_950, measure1_975 = OM1_baseline_975,
                              measure2 = OM2, measure2_baseline = OM2_baseline, 
                              measure2_025 = OM2_baseline_025, measure2_950 = OM2_baseline_950, measure2_975 = OM2_baseline_975,
                              no_path_first = number_pathway_first_half, no_path_last = number_pathway_last_half,
                              measure3 = OM3, measure3_baseline = OM3_baseline, 
                              measure3_025 = OM3_baseline_025, measure3_950 = OM3_baseline_950, measure3_975 = OM3_baseline_975)

write.csv(Overlap_measure,"./Overlap_measure.csv")


