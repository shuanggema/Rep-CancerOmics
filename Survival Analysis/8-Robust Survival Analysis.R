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

setwd("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD")

number_var <- 9 # for marginal analysis, dimensions of all covariates

target_number_joint <- 71 # target number of genes being selected

lambda <- c(0.3, 0.2, 0.1, 0.05, 0.01) # potential lambda for joint analysis

######## output ########
# File - "Selected_gene_first_half.csv"
# File - "Selected_gene_last_half.csv"

########
library(clusterSim)
library(quantreg)

############# delete genes with variance <= 1st Quantile ################
x_gene <- x_gene[,gene_list_subset]
x_gene_num <- matrix(as.numeric(x_gene), ncol = ncol(x_gene))

# standardize all the gene expression
x_gene_num <- data.Normalization(x_gene_num,type="n1",normalization="column")

############# conduct survival analysis for the remaining gene ################

### calculate weight to account censoring ###
AFT_weight <- function(survtime, status){
  n <- length(survtime)
  weight <- rep(NA,n)
  
  order <- order(survtime)
  time_order <- survtime[order]
  status_order <- status[order]
  
  temp_vec <- ((n-1:n)/(n-1:n+1))^status_order
  
  weight[order[1]] <- status_order[1]/n
  for (i in 2:n){
    weight[order[i]] <- status_order[i]/(n-i+1) * prod(temp_vec[1:(i-1)])
  }
  return(weight)
}

###
number_obs <- nrow(x_gene)
ind <- 1:number_obs

Selected_gene_first_half_marginal <- list()
Selected_gene_last_half_marginal <- list()
Selected_gene_first_half_joint <- list()
Selected_gene_last_half_joint <- list()

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
  
  ##################### for the first half of samples #####################
  weight <- AFT_weight(survival_selected_first_half[,1], survival_selected_first_half[,2])
  X <- cbind(gene_selected_first_half,demo_selected_first_half) * weight
  X <- X[weight != 0,]
  Y <- survival_selected_first_half[,1] * weight
  Y <- Y[weight != 0]
  
  ### marginal analysis
  p_value <- NA

  for(j in 1:length(gene_list_subset)){
    model <- rq(Y ~ X[,c(j,(dim(X)[2]-number_var + 1):dim(X)[2])], tau=.5, method="br", model = TRUE)
    p_value[j] <- summary(model, se = "boot")$coefficients[3*(number_var+2)+2]
  }
  
  p_value_adjust <- p.adjust(p_value, method = "bonferroni", n = length(p_value))
  Selected_gene_first_half_marginal[[i]] <- gene_list_subset[p_value_adjust < 0.1]
  
  ### joint analysis
  gene_k <- list()
  for (k in 1:length(lambda)){
    model <- rq.fit.lasso(X, Y, tau = 0.5, eps = 1e-06, lambda = lambda[k])
    gene_k[[k]] <- gene_list_subset[abs(model$coefficients[1:length(gene_list_subset)]) >= 
                                    max(abs(model$coefficients[1:length(gene_list_subset)]))*0.01]
  }
  Selected_gene_first_half_joint[[i]] <- gene_k[[which.min(abs(lengths(gene_k)-target_number_joint))]]

  
  ##################### for the last half of samples #####################
  weight <- AFT_weight(survival_selected_last_half[,1], survival_selected_last_half[,2])
  X <- cbind(gene_selected_last_half,demo_selected_last_half) * weight
  X <- X[weight != 0,]
  Y <- survival_selected_last_half[,1] * weight
  Y <- Y[weight != 0]
  
  ### marginal analysis
  p_value <- NA
  
  for(j in 1:length(gene_list_subset)){
    model <- rq(Y ~ X[,c(j,(dim(X)[2]-number_var + 1):dim(X)[2])], tau=.5, method="br", model = TRUE)
    p_value[j] <- summary(model, se = "boot")$coefficients[3*(number_var+2)+2]
  }
  
  p_value_adjust <- p.adjust(p_value, method = "bonferroni", n = length(p_value))
  Selected_gene_last_half_marginal[[i]] <- gene_list_subset[p_value_adjust < 0.1]
  
  ### joint analysis
  gene_k <- list()
  for (k in 1:length(lambda)){
    model <- rq.fit.lasso(X, Y, tau = 0.5, eps = 1e-06, lambda = lambda[k])
    gene_k[[k]] <- gene_list_subset[abs(model$coefficients[1:length(gene_list_subset)]) >= 
                                      max(abs(model$coefficients[1:length(gene_list_subset)]))*0.01]
  }
  Selected_gene_last_half_joint[[i]] <- gene_k[[which.min(abs(lengths(gene_k)-target_number_joint))]]
  i <- i+1
}


############# store the list as dataframe  ################
# select to store either the marginal analysis results or the joint analysis results
Selected_gene_first_half <- Selected_gene_first_half_marginal
Selected_gene_last_half <- Selected_gene_last_half_marginal

Selected_gene_first_half <- Selected_gene_first_half_joint
Selected_gene_last_half <- Selected_gene_last_half_joint

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

# select to store either the marginal analysis results or the joint analysis results
write.csv(Selected_gene_first_half_df,"./Significant Gene_Robust/Selected_gene_first_half.csv") # store marginal
write.csv(Selected_gene_last_half_df,"./Significant Gene_Robust/Selected_gene_last_half.csv") # store marginal

write.csv(Selected_gene_first_half_df,"./Significant Gene_Robust2/Selected_gene_first_half.csv") # store joint
write.csv(Selected_gene_last_half_df,"./Significant Gene_Robust2/Selected_gene_last_half.csv") # store joint
