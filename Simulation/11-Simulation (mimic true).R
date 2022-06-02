options(stringsAsFactors = FALSE)

######## input ########
Final_data <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Merged/Final_data.csv", row.names = 1)
gene_list_full <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/gene_list_full.csv",row.names = 1))
EnsemblID_EntrezID_Pathway <- read.csv("/Users/name/Desktop/TCGA-2021/Data/EnsemblID_EntrezID_Pathway.csv", row.names = 1)

# use the same seeds for data splitting in the data analysis
seed_list <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/seed_list.csv",row.names = 1))

setwd("/Users/name/Desktop/TCGA-2021/Data/Simulation-LUAD")

######## output ########
# File - "Overlap_measure.csv": marginal or joint
# File - "Overlap_measure_baseline.csv": marginal or joint


########
library(clusterSim)
library(ggfortify)
library(survival)
library(glmnet)

################ Pre-process #################
# only retain gene that has different Ensembl ID (delete repeated gene expression for the same Ensembl ID)
table(gene_list_full)[table(gene_list_full) >= 2]
gene_list_full <- unique(gene_list_full)

# only retain gene expression information (variance > 0.25) and discard original covariates
gene_variance <- apply(as.matrix(Final_data[,gene_list_full]), 2, var)
gene_list_subset <- gene_list_full[!gene_variance <= quantile(gene_variance, 0.25)]
Simu_data <- Final_data[,gene_list_subset]
Simu_data <- data.Normalization(Simu_data,type="n1",normalization="column")

############# True data part (get observed beta) #################
# True survival information
Surv_info <- Final_data[,c("time","status")]

# plot the survival curve
fit <- survfit(Surv(time, status) ~ 1, data = Surv_info)
autoplot(fit)

beta_p_val_from_data <- data.frame(gene = gene_list_subset, beta = NA, p_value = NA)
j <- 1
for(gene in gene_list_subset){
  model <- coxph(as.formula(paste("Surv(time, status) ~",gene)), cbind(Simu_data, Surv_info))
  beta_p_val_from_data[j, "beta"] <- model$coefficients
  beta_p_val_from_data[j, "p_value"] <- summary(model)$coefficients[5]
  j <- j+1
}

par(mfrow = c(2,1))
plot(1:5799, sort(beta_p_val_from_data$beta,decreasing = TRUE), main = "Distribution of true beta (ordered)", xlab = "", ylab = "true beta")
plot(1:5799, sort(abs(beta_p_val_from_data$beta),decreasing = TRUE), main = "Distribution of absolute true beta (ordered)", xlab = "", ylab = "abs true beta")

p_value_adjust <- p.adjust(beta_p_val_from_data$p_value, method = "bonferroni", n = length(beta_p_val_from_data$p_value))
Selected_gene_from_data <- gene_list_subset[p_value_adjust < 0.1]

############# Select gene and get true beta (mimic observed beta) #################
# store selected genes
select_gene <- sample(gene_list_subset, 30, replace = FALSE)
# select_gene <- as.matrix(read.csv("./selected_gene.csv", row.names = 1))

# set parameters
lambda <- 10^(-3)
beta <- rep(0,ncol(Simu_data))
beta[which(gene_list_subset %in% select_gene)] <- 0.13

# generate duplicated data and get estimated true beta
Simu_data_pop_ind <- sample(nrow(Simu_data), nrow(Simu_data)*100, replace = TRUE)
Simu_data_pop <- Simu_data[Simu_data_pop_ind,]

# show the distribution of beta*X
summary(rowSums(0.4*Simu_data_pop[,select_gene]))

# create simulated population data (no censoring)
Simu_data_pop$time <-  -log(runif(nrow(Simu_data_pop)))/(lambda*exp(as.matrix(Simu_data_pop) %*% as.matrix(beta)))
Simu_data_pop$status <-  1

# plot the survival curve
fit <- survfit(Surv(time, status) ~ 1, data = Simu_data_pop[Simu_data_pop$time < 10000,])
autoplot(fit)

fit <- survfit(Surv(time, status) ~ 1, data = Simu_data_temp)
autoplot(fit)

beta_from_simu <- data.frame(gene = gene_list_subset, beta = NA)
j <- 1
for(gene in gene_list_subset){
  model <- coxph(as.formula(paste("Surv(time, status) ~",gene)), Simu_data_pop)
  beta_from_simu[j,"beta"] <- model$coefficients
  j <- j+1
}
par(mfrow = c(1,1))
plot(1:5799, sort(beta_from_simu[,2],decreasing = TRUE), main = "Distribution of true beta (ordered)", xlab = "", ylab = "true beta")
plot(1:5799, sort(abs(beta_from_simu[,2]),decreasing = TRUE), main = "Distribution of absolute true beta (ordered)", xlab = "", ylab = "abs true beta")

par(mfrow = c(2,1))
plot(1:5799, sort(abs(beta_p_val_from_data$beta),decreasing = TRUE), main = "Distribution of absolute true beta (ordered)", xlab = "from data", ylab = "absolute true beta")
plot(1:5799, sort(abs(beta_from_simu[,2]),decreasing = TRUE), xlab = "from simulation", ylab = "absolute true beta", main = "")

# write.csv(select_gene,"./selected_gene.csv")


###################### Another WAY ##################################
###################### Another WAY ##################################
###################### Another WAY ##################################
Selected_gene_from_simu_list_marginal1 <- list()
Selected_gene_from_simu_list_marginal2 <- list()
Selected_gene_from_simu_list_en1 <- list()
Selected_gene_from_simu_list_en2 <- list()

i <- 1
for (seed in seed_list){
  Simu_data_temp <- Simu_data
  Simu_data_temp$time1 <-  -log(runif(nrow(Simu_data_temp)))/(lambda*exp(as.matrix(Simu_data_temp) %*% as.matrix(beta)))
  Simu_data_temp$time2 <- runif(n = nrow(Simu_data_temp),  min = 0, max = 2000)
  Simu_data_temp$status <-  0
  Simu_data_temp$status[Simu_data_temp$time1 <= Simu_data_temp$time2] <- 1
  Simu_data_temp$time <- Simu_data_temp$time1
  Simu_data_temp$time[Simu_data_temp$status == 0] <- Simu_data_temp$time2[Simu_data_temp$status == 0]
  
  set.seed(seed)
  ind_half <- sample(1:nrow(Simu_data_temp), round(nrow(Simu_data_temp)/2,0), replace = FALSE)
  Simu_data_temp1 <- Simu_data_temp[ind_half,]
  Simu_data_temp2 <- Simu_data_temp[-ind_half,]
  
  Simu_data_temp1 <- rbind(Simu_data_temp1,Simu_data_temp1)
  Simu_data_temp2 <- rbind(Simu_data_temp2,Simu_data_temp2)
  ############# Marginal: conduct survival analysis for the remaining gene ################
  p_value <- NA
  j <- 1
  for(gene in gene_list_subset){
    model <- coxph(as.formula(paste("Surv(time, status) ~",gene)), Simu_data_temp1)
    p_value[j] <- summary(model)$coefficients[5]
    j <- j+1
  }
  p_value_adjust <- p.adjust(p_value, method = "bonferroni", n = length(p_value))
  
  Selected_gene_from_simu_list_marginal1[[i]] <- gene_list_subset[p_value_adjust < 0.1]
  
  p_value <- NA
  j <- 1
  for(gene in gene_list_subset){
    model <- coxph(as.formula(paste("Surv(time, status) ~",gene)), Simu_data_temp2)
    p_value[j] <- summary(model)$coefficients[5]
    j <- j+1
  }
  p_value_adjust <- p.adjust(p_value, method = "bonferroni", n = length(p_value))
  
  Selected_gene_from_simu_list_marginal2[[i]] <- gene_list_subset[p_value_adjust < 0.1]
  
  ############# EN: conduct survival analysis for the remaining gene ################ 
  EN_cv <- cv.glmnet(as.matrix(Simu_data_temp1[,1:ncol(Simu_data)]), 
                     Surv(as.matrix(Simu_data_temp1$time), as.matrix(Simu_data_temp1$status), type = 'right'),
                     family = "cox", alpha = 0.5, standardize = TRUE, nfolds = 3)
  
  Selected_gene_from_simu_list_en1[[i]] <- gene_list_subset[as.vector(!coef(EN_cv) == 0)]
  
  EN_cv <- cv.glmnet(as.matrix(Simu_data_temp2[,1:ncol(Simu_data)]), 
                     Surv(as.matrix(Simu_data_temp2$time), as.matrix(Simu_data_temp2$status), type = 'right'),
                     family = "cox", alpha = 0.5, standardize = TRUE, nfolds = 3)
  
  Selected_gene_from_simu_list_en2[[i]] <- gene_list_subset[as.vector(!coef(EN_cv) == 0)]
  
  i <- i+1
}

# choose either the marginal analysis results or the joint analysis results
Selected_gene_first_half_list <- Selected_gene_from_simu_list_en1
Selected_gene_last_half_list <- Selected_gene_from_simu_list_en2


####################################################################
################### gene name overlap ##############################
####################################################################
# initialize vectors
OM1 <- NA
OM1_expect <- NA
OM1_baseline <- NA
OM1_baseline_025 <- NA
OM1_baseline_975 <- NA

for (i in 1:100){
  first <- Selected_gene_first_half_list[[i]]
  last <- Selected_gene_last_half_list[[i]]
  
  if (min(length(first),length(last)) == 0){
    OM1[i] <- 0
    OM1_expect[i] <- 0
    OM1_baseline[i] <- 0
    OM1_baseline_025[i] <- 0
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
    OM1_baseline_975[i] <- quantile(OM1_baseline_boot, 0.975)
  }
}

####################################################################
################### gene correlation overlap #######################
####################################################################
# initialize vectors
OM2 <- NA
OM2_baseline <- NA
OM2_baseline_025 <- NA
OM2_baseline_975 <- NA

for (i in 1:100){
  first <- Selected_gene_first_half_list[[i]]
  last <- Selected_gene_last_half_list[[i]]
  
  if(min(length(first),length(last)) == 0){
    OM2[i] <- 0
    OM2_baseline[i] <- 0
    OM2_baseline_025[i] <- 0
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
    OM2_baseline_975[i] <- quantile(OM2_baseline_boot, 0.975)
  }
}


#################################################################
################### gene pathway overlap ########################
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


############# pathway overlap ################

# initialize vectors
OM3 <- NA
OM3_baseline <- NA
OM3_baseline_025 <- NA
OM3_baseline_975 <- NA

for (i in 1:100){
  first_gene <- Selected_gene_first_half_list[[i]]
  last_gene <- Selected_gene_last_half_list[[i]]
  
  if(min(length(first_gene),length(last_gene)) == 0){
    OM3[i] <- 0
    OM3_baseline[i] <- 0
    OM3_baseline_025[i] <- 0
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
    OM3_baseline_975[i] <- quantile(OM3_baseline_boot, 0.975)
  }
}

no_genes_first <- lengths(Selected_gene_first_half_list)
no_genes_last <- lengths(Selected_gene_last_half_list)


Overlap_measure <- data.frame(no_genes_first = no_genes_first, no_genes_last = no_genes_last,
                              measure1 = OM1, 
                              measure2 = OM2, 
                              no_path_first = number_pathway_first_half, 
                              no_path_last = number_pathway_last_half,
                              measure3 = OM3)

Overlap_measure_baseline <- data.frame(measure1_expect = OM1_expect, measure1_baseline = OM1_baseline, 
                                       measure1_025 = OM1_baseline_025, measure1_975 = OM1_baseline_975,
                                       measure2_baseline = OM2_baseline, 
                                       measure2_025 = OM2_baseline_025, measure2_975 = OM2_baseline_975,
                                       measure3_baseline = OM3_baseline, 
                                       measure3_025 = OM3_baseline_025, measure3_975 = OM3_baseline_975)
###
write.csv(Overlap_measure,"./Overlap_measure.csv")
write.csv(Overlap_measure_baseline,"./Overlap_measure_baseline.csv")






