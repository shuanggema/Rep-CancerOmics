options(stringsAsFactors = TRUE)

######## input ########
x_demo <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Merged/x_demo.csv",
                   row.names = 1))
x_gene <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Merged/x_gene.csv",
                   row.names = 1))
x_survival <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Merged/x_survival.csv",
                   row.names = 1))

gene_list_subset <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Significant Gene/gene_list_subset.csv",row.names = 1))

seed_list <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/seed_list.csv",row.names = 1))

setwd("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Significant Gene")

######## output ########
# File - "Selected_gene_first_half.csv"
# File - "Selected_gene_last_half.csv"

########
library(glmnet)

############# delete genes with variance <= 1st Quantile ################
x_gene <- x_gene[,gene_list_subset]
x_gene_num <- matrix(as.numeric(x_gene), ncol = ncol(x_gene))

############# conduct survival analysis for the remaining gene ################
# Setting alpha = 0.5 implements half lasso and half ridge regression
# Penalize both genes and demos


# through cross-validation select best lambda 
# (gives lowest cross validation partial-likelihood for the Cox model )

###

number_obs <- nrow(x_gene)
ind <- 1:number_obs

Selected_gene_first_half <- list()
Selected_gene_last_half <- list()

i <- 1
for (seed in seed_list){
  set.seed(seed)
  ind_half <- sample(ind, round(number_obs/2,0), replace = FALSE)
  gene_selected_first_half <- x_gene_num[ind_half,]
  gene_selected_last_half <- x_gene_num[-ind_half,]
  
  demo_selected_first_half <- x_demo[ind_half,]
  demo_selected_last_half <- x_demo[-ind_half,]
  
  survival_selected_first_half <- x_survival[ind_half,]
  survival_selected_last_half <- x_survival[-ind_half,]
  
  gene_selected_first_half <- rbind(gene_selected_first_half,gene_selected_first_half)
  gene_selected_last_half <- rbind(gene_selected_last_half,gene_selected_last_half)
  
  demo_selected_first_half <- rbind(demo_selected_first_half,demo_selected_first_half)
  demo_selected_last_half <- rbind(demo_selected_last_half,demo_selected_last_half)
  
  survival_selected_first_half <- rbind(survival_selected_first_half,survival_selected_first_half)
  survival_selected_last_half <- rbind(survival_selected_last_half,survival_selected_last_half)
  
  ### for the first half of samples 
  EN_cv <- cv.glmnet(cbind(gene_selected_first_half,demo_selected_first_half), 
                     Surv(survival_selected_first_half[,"time"] + runif(1,0,0.01), 
                          survival_selected_first_half[,"status"], type = 'right'),
                     family = "cox", alpha = 0.5, standardize = TRUE, nfolds = 3)
  
  Selected_gene_first_half[[i]] <- gene_list_subset[as.vector(!coef(EN_cv) == 0)[1:ncol(x_gene)]]

  ### for the last half of samples 
  EN_cv <- cv.glmnet(cbind(gene_selected_last_half,demo_selected_last_half), 
                     Surv(survival_selected_last_half[,"time"] + runif(1,0,0.01), 
                          survival_selected_last_half[,"status"], type = 'right'),
                     family = "cox", alpha = 0.5, standardize = TRUE, nfolds = 3)
  
  Selected_gene_last_half[[i]] <- gene_list_subset[as.vector(!coef(EN_cv) == 0)[1:ncol(x_gene)]]
  
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















