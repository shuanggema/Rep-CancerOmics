options(stringsAsFactors = TRUE)

######## input ########
Final_data <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Merged/Final_data.csv", row.names = 1)
gene_list_full <- as.matrix(read.csv("/Users/name/Desktop/TCGA-2021/Data/gene_list_full.csv",row.names = 1))

EnsemblID_EntrezID_Pathway <- read.csv("/Users/name/Desktop/TCGA-2021/Data/EnsemblID_EntrezID_Pathway.csv", row.names = 1)

setwd("/Users/nannie/Desktop/TCGA-2021/Data/Simulation/beta0.15")
######## output ########
# File - "./Marginal/Selected_gene_from_data_x.csv": x = 1,2,3,4,5
# File - "./EN/Selected_gene_from_data_x.csv": x = 1,2,3,4,5
# File - "./Marginal/Overlap_measure_vs100_x.csv": x = 1,2,3,4,5
# File - "./Marginal/Overlap_measure_baseline_vs100_x.csv": x = 1,2,3,4,5
# File - "./EN/Overlap_measure_vs100_x.csv": x = 1,2,3,4,5
# File - "./EN/Overlap_measure_baseline_vs100_x.csv": x = 1,2,3,4,5

########
library(clusterSim)
library(ggfortify)
library(survival)
library(glmnet)
library(ggplot2)
library(dplyr)
library(cowplot)

############# Select gene #################
# only retain gene that has different Ensembl ID (delete repeated gene expression for the same Ensembl ID)
table(gene_list_full)[table(gene_list_full) >= 2]
gene_list_full <- unique(gene_list_full)

# only retain gene expression information (variance > 0.25) and discard original covariates and survival
gene_variance <- apply(as.matrix(Final_data[,gene_list_full]), 2, var)
gene_list_subset <- gene_list_full[!gene_variance <= quantile(gene_variance, 0.25)]
Simu_data <- Final_data[,gene_list_subset]
Simu_data <- data.Normalization(Simu_data,type="n1",normalization="column")


############# Set parameter and find true parameter #################
# set parameters
lambda <- 10^(-3)
beta <- rep(0,ncol(Simu_data))
beta[which(gene_list_subset %in% select_gene)] <- 0.1

# generate duplicated data and get estimated true beta
Simu_data_pop_ind <- sample(nrow(Simu_data), nrow(Simu_data)*100, replace = TRUE)
Simu_data_pop <- Simu_data[Simu_data_pop_ind,]

# create simulated population data (no censoring)
Simu_data_pop$time <-  -log(runif(nrow(Simu_data_pop)))/(lambda*exp(as.matrix(Simu_data_pop) %*% as.matrix(beta)))
Simu_data_pop$status <-  1

# plot the survival curve
fit <- survfit(Surv(time, status) ~ 1, data = Simu_data_pop[Simu_data_pop$time < 10000,])
autoplot(fit)

# calculate true beta
true_beta <- data.frame(gene = gene_list_subset, beta = NA)

j <- 1
for(gene in gene_list_subset){
  model <- coxph(as.formula(paste("Surv(time, status) ~",gene)), Simu_data_pop)
  true_beta[j,"beta"] <- model$coefficients
  j <- j+1
}
par(mfrow = c(2,1))
plot(1:5799, sort(true_beta[,2],decreasing = TRUE), main = "Distribution of true beta (ordered)", xlab = "", ylab = "true beta")
plot(1:5799, sort(abs(true_beta[,2]),decreasing = TRUE), main = "Distribution of absolute true beta (ordered)", xlab = "", ylab = "abs true beta")

# write.csv(true_beta,"./true_beta.csv")

############# Create simulated data for analysis #################
lambda <- 10^(-3)
beta <- rep(0,ncol(Simu_data))
beta[which(gene_list_subset %in% select_gene)] <- 0.1

Simu_data_temp <- Simu_data
Simu_data_temp$time1 <-  -log(runif(nrow(Simu_data_temp)))/(lambda*exp(as.matrix(Simu_data_temp) %*% as.matrix(beta)))
Simu_data_temp$time2 <- runif(n = nrow(Simu_data_temp),  min = 0, max = 2000)
Simu_data_temp$status <-  0
Simu_data_temp$status[Simu_data_temp$time1 <= Simu_data_temp$time2] <- 1
Simu_data_temp$time <- Simu_data_temp$time1
Simu_data_temp$time[Simu_data_temp$status == 0] <- Simu_data_temp$time2[Simu_data_temp$status == 0]

# plot the survival curve
fit <- survfit(Surv(time, status) ~ 1, data = Simu_data_temp)
autoplot(fit)


############# Marginal: conduct survival analysis for the remaining gene ################
number_genes_compare <- 100

Selected_gene_first_half <- list()
Selected_gene_last_half <- list()

p_value <- NA
j <- 1
for(gene in gene_list_subset){
  model <- coxph(as.formula(paste("Surv(time, status) ~",gene)), Simu_data_temp)
  p_value[j] <- summary(model)$coefficients[5]
  j <- j+1
}
p_value_adjust <- p.adjust(p_value, method = "bonferroni", n = length(p_value))
Selected_gene_from_data <- gene_list_subset[p_value_adjust < 0.1]

write.csv(Selected_gene_from_data,"./Marginal/Selected_gene_from_data_1.csv")


gene_true_beta_ordered_list <- true_beta[order(abs(true_beta[,2]),decreasing = TRUE),"gene"]
for (i in 1:(length(gene_list_subset)-number_genes_compare+1)){
  Selected_gene_first_half[[i]] <- Selected_gene_from_data
  Selected_gene_last_half[[i]] <- as.character(gene_true_beta_ordered_list[i:(i+number_genes_compare-1)])
}


############# EN: conduct survival analysis for the remaining gene ################
number_genes_compare <- 100

Selected_gene_first_half <- list()
Selected_gene_last_half <- list()

EN_cv <- cv.glmnet(as.matrix(Simu_data), 
                     Surv(Simu_data_temp$time, Simu_data_temp$status, type = 'right'),
                     family = "cox", alpha = 0.5, standardize = TRUE, nfolds = 3)
  
Selected_gene_from_data <- gene_list_subset[as.vector(!coef(EN_cv) == 0)]

write.csv(Selected_gene_from_data,"./EN/Selected_gene_from_data_1.csv")

gene_true_beta_ordered_list <- true_beta[order(abs(true_beta[,2]),decreasing = TRUE),"gene"]
for (i in 1:(length(gene_list_subset)-number_genes_compare+1)){
  Selected_gene_first_half[[i]] <- Selected_gene_from_data
  Selected_gene_last_half[[i]] <- as.character(gene_true_beta_ordered_list[i:(i+number_genes_compare-1)])
}

####################################################################
################### gene name overlap ##############################
####################################################################
# initialize vectors
OM1 <- NA
OM1_expect <- NA
OM1_baseline <- NA
OM1_baseline_025 <- NA
OM1_baseline_975 <- NA

for (i in 1:(length(gene_list_subset)-number_genes_compare+1)){
  first <- Selected_gene_first_half[[i]]
  last <- Selected_gene_last_half[[i]]
  
  OM1[i] <- sum(first %in% last)/min(length(first),length(last))
}   

# calculate baseline (expected)
OM1_expect <- max(length(Selected_gene_from_data),number_genes_compare)/length(gene_list_subset)
    
# calculate baseline (from data)
OM1_baseline_boot <- rep(NA,1000)
for (t in 1:1000){
  gene_list_temp1 <- sample(gene_list_subset, length(Selected_gene_from_data), replace = FALSE)
  gene_list_temp2 <- sample(gene_list_subset, number_genes_compare, replace = FALSE)
      
  OM1_baseline_boot[t] <- sum(gene_list_temp1 %in% gene_list_temp2)/min(length(Selected_gene_from_data),number_genes_compare)
}
    
OM1_baseline <- mean(OM1_baseline_boot)
OM1_baseline_025 <- quantile(OM1_baseline_boot, 0.025)
OM1_baseline_975 <- quantile(OM1_baseline_boot, 0.975)


####################################################################
################### gene correlation overlap #######################
####################################################################
# initialize vectors
OM2 <- NA
OM2_baseline <- NA
OM2_baseline_025 <- NA
OM2_baseline_975 <- NA

for (i in 1:(length(gene_list_subset)-number_genes_compare+1)){
  first <- Selected_gene_first_half[[i]]
  last <- Selected_gene_last_half[[i]]

  gene_cor_comparison <- abs(cor(Simu_data_temp[first],Simu_data_temp[last]))
  OM2[i] <- (sum(apply(gene_cor_comparison,1,max)) + sum(apply(gene_cor_comparison,2,max)))/(length(first)+length(last))
}    

# calculate baseline (from data)
OM2_baseline_boot <- rep(NA,1000)
for (t in 1:1000){
  gene_list_temp1 <- sample(gene_list_subset, length(Selected_gene_from_data), replace = FALSE)
  gene_list_temp2 <- sample(gene_list_subset, number_genes_compare, replace = FALSE)
      
  gene_cor_comparison_temp <- abs(cor(Simu_data_temp[gene_list_temp1],Simu_data_temp[gene_list_temp2]))
      
  OM2_baseline_boot[t] <- (sum(apply(gene_cor_comparison_temp,1,max)) +
                                 sum(apply(gene_cor_comparison_temp,2,max)))/(length(Selected_gene_from_data) + number_genes_compare)
}
    
OM2_baseline <- mean(OM2_baseline_boot)
OM2_baseline_025 <- quantile(OM2_baseline_boot, 0.025)
OM2_baseline_975 <- quantile(OM2_baseline_boot, 0.975)


#################################################################
################### gene pathway overlap ########################
#################################################################

############# FIRST HALF VS LAST HALF ################
### How many significant genes related pathways are selected

# get gene linked pathway
Selected_pathway_last_half_list <- list()

# store the number of pathways selected
number_pathway_first_half <- NA
number_pathway_last_half <- NA

Selected_pathway_first_half <- EnsemblID_EntrezID_Pathway[EnsemblID_EntrezID_Pathway$Ensembl_ID %in% Selected_gene_first_half[[1]],]
number_pathway_first_half_number <- length(unique(Selected_pathway_first_half$Pathway_ID))

for (i in 1:(length(gene_list_subset)-number_genes_compare+1)){
  number_pathway_first_half[i] <- number_pathway_first_half_number
  
  Selected_pathway_last_half_list[[i]] <- EnsemblID_EntrezID_Pathway[EnsemblID_EntrezID_Pathway$Ensembl_ID %in% Selected_gene_last_half[[i]],]
  number_pathway_last_half[i] <- length(unique(Selected_pathway_last_half_list[[i]]$Pathway_ID))
}


############# pathway overlap ################

# initialize vectors
OM3 <- NA
OM3_baseline <- NA
OM3_baseline_025 <- NA
OM3_baseline_975 <- NA

##### for first half
first_gene <- Selected_gene_first_half[[1]]
first_path <- Selected_pathway_first_half

unique_pathway_first_half <- table(as.character(first_path$Pathway_ID))
unique_pathway_first_half <- data.frame(Pathway_ID = names(unique_pathway_first_half), count1 = as.numeric(unique_pathway_first_half))
unique_pathway_first_half$count2 <- NA

j <- 1
for(id in unique_pathway_first_half$Pathway_ID){
  unique_pathway_first_half$count2[j] <- as.character(first_path$GeneRatio)[first_path$Pathway_ID == id][1]
  j <- j + 1
}

unique_pathway_first_half$count2 <- as.numeric(sub("/.*", "", unique_pathway_first_half$count2))
unique_pathway_first_half$weight <- unique_pathway_first_half$count1/unique_pathway_first_half$count2


##### for last half
for (i in 1:(length(gene_list_subset)-number_genes_compare+1)){
  last_gene <- Selected_gene_last_half[[i]]
  last_path <- Selected_pathway_last_half_list[[i]]
    
  ####   ####   ####   ####   ####    ####   ####   ####   ####   ####    ####   ####   ####   ####   ####
  unique_pathway_last_half <- table(as.character(last_path$Pathway_ID))
  unique_pathway_last_half <- data.frame(Pathway_ID = names(unique_pathway_last_half), count1 = as.numeric(unique_pathway_last_half))
  unique_pathway_last_half$count2 <- NA
    
  j <- 1
  for(id in unique_pathway_last_half$Pathway_ID){
    unique_pathway_last_half$count2[j] <- as.character(last_path$GeneRatio)[last_path$Pathway_ID == id][1]
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
}   

####   ####   ####   ####   ####   
####   ####   ####   ####   ####  
# calculate baseline (from data)
OM3_baseline_boot <- rep(NA,1000)
    
for (t in 1:1000){
  gene_list_temp1 <- sample(gene_list_subset, length(Selected_gene_from_data), replace = FALSE)
  gene_list_temp2 <- sample(gene_list_subset, number_genes_compare, replace = FALSE)
      
  sig_gene_pathway_temp1 <- EnsemblID_EntrezID_Pathway[EnsemblID_EntrezID_Pathway$Ensembl_ID %in% gene_list_temp1,]
  sig_gene_pathway_temp2 <- EnsemblID_EntrezID_Pathway[EnsemblID_EntrezID_Pathway$Ensembl_ID %in% gene_list_temp2,]
      
  ###
  unique_pathway_temp1 <- table(as.character(sig_gene_pathway_temp1$Pathway_ID))
  unique_pathway_temp1 <- data.frame(Pathway_ID = names(unique_pathway_temp1), count1 = as.numeric(unique_pathway_temp1))
      
  j <- 1
  for(id in unique_pathway_temp1$Pathway_ID){
    unique_pathway_temp1$count2[j] <- as.character(sig_gene_pathway_temp1$GeneRatio)[sig_gene_pathway_temp1$Pathway_ID == id][1]
    j <- j + 1
  }
      
  unique_pathway_temp1$count2 <- as.numeric(sub("/.*", "", unique_pathway_temp1$count2))
  unique_pathway_temp1$weight <- unique_pathway_temp1$count1/unique_pathway_temp1$count2
      
  ###
  unique_pathway_temp2 <- table(as.character(sig_gene_pathway_temp2$Pathway_ID))
  unique_pathway_temp2 <- data.frame(Pathway_ID = names(unique_pathway_temp2), count1 = as.numeric(unique_pathway_temp2))
      
  j <- 1
  for(id in unique_pathway_temp2$Pathway_ID){
    unique_pathway_temp2$count2[j] <- as.character(sig_gene_pathway_temp2$GeneRatio)[sig_gene_pathway_temp2$Pathway_ID == id][1]
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
    
OM3_baseline <- mean(OM3_baseline_boot)
OM3_baseline_025 <- quantile(OM3_baseline_boot, 0.025)
OM3_baseline_975 <- quantile(OM3_baseline_boot, 0.975)


no_genes_first <- lengths(Selected_gene_first_half)
no_genes_last <- lengths(Selected_gene_last_half)


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

write.csv(Overlap_measure,"./EN/Overlap_measure_vs100_1.csv") # or change to Marginal/...
write.csv(Overlap_measure_baseline,"./EN/Overlap_measure_baseline_vs100_1.csv") # or change to Marginal/...


################################ Draw plots #####################################
Overlap_measure1 <- read.csv("./Marginal/Overlap_measure_vs100_1.csv", row.names = 1)
Overlap_measure2 <- read.csv("./Marginal/Overlap_measure_vs100_2.csv", row.names = 1)
Overlap_measure3 <- read.csv("./Marginal/Overlap_measure_vs100_3.csv", row.names = 1)
Overlap_measure4 <- read.csv("./Marginal/Overlap_measure_vs100_4.csv", row.names = 1)
Overlap_measure5 <- read.csv("./Marginal/Overlap_measure_vs100_5.csv", row.names = 1)
Overlap_measure_baseline1 <- read.csv("./Marginal/Overlap_measure_baseline_vs100_1.csv", row.names = 1)
Overlap_measure_baseline2 <- read.csv("./Marginal/Overlap_measure_baseline_vs100_2.csv", row.names = 1)
Overlap_measure_baseline3 <- read.csv("./Marginal/Overlap_measure_baseline_vs100_3.csv", row.names = 1)
Overlap_measure_baseline4 <- read.csv("./Marginal/Overlap_measure_baseline_vs100_4.csv", row.names = 1)
Overlap_measure_baseline5 <- read.csv("./Marginal/Overlap_measure_baseline_vs100_5.csv", row.names = 1)

Overlap_measure1 <- read.csv("./EN/Overlap_measure_vs100_1.csv", row.names = 1)
Overlap_measure2 <- read.csv("./EN/Overlap_measure_vs100_2.csv", row.names = 1)
Overlap_measure3 <- read.csv("./EN/Overlap_measure_vs100_3.csv", row.names = 1)
Overlap_measure4 <- read.csv("./EN/Overlap_measure_vs100_4.csv", row.names = 1)
Overlap_measure5 <- read.csv("./EN/Overlap_measure_vs100_5.csv", row.names = 1)
Overlap_measure_baseline1 <- read.csv("./EN/Overlap_measure_baseline_vs100_1.csv", row.names = 1)
Overlap_measure_baseline2 <- read.csv("./EN/Overlap_measure_baseline_vs100_2.csv", row.names = 1)
Overlap_measure_baseline3 <- read.csv("./EN/Overlap_measure_baseline_vs100_3.csv", row.names = 1)
Overlap_measure_baseline4 <- read.csv("./EN/Overlap_measure_baseline_vs100_4.csv", row.names = 1)
Overlap_measure_baseline5 <- read.csv("./EN/Overlap_measure_baseline_vs100_5.csv", row.names = 1)

################################ Draw plots #####################################

############# true beta ##########
data.frame(beta = sort(abs(true_beta[,2]),decreasing = TRUE)) %>% 
  ggplot(aes(x = 1:5799, y = beta)) + geom_point(aes(y = beta)) + 
  ggtitle("Ordered abs true beta")

############# OM 1 ################
p1 <- ggplot(Overlap_measure1, aes(1:(5799-100+1))) +  
  geom_line(aes(y = measure1), color = "pink3", alpha = 0.5, size = 2.5) +
  geom_line(aes(y = Overlap_measure_baseline1$measure1_baseline), color = "pink3") +
  geom_line(aes(y = Overlap_measure_baseline1$measure1_025), color = "pink3") +
  geom_line(aes(y = Overlap_measure_baseline1$measure1_975), color = "pink3") + 
  geom_line(aes(y = Overlap_measure2$measure1), color = "69b3a2", alpha = 0.5, size = 2) +
  geom_line(aes(y = Overlap_measure_baseline2$measure1_baseline), color = "69b3a2") +
  geom_line(aes(y = Overlap_measure_baseline2$measure1_025), color = "69b3a2") +
  geom_line(aes(y = Overlap_measure_baseline2$measure1_975), color = "69b3a2", size = 1) + 
  geom_line(aes(y = Overlap_measure3$measure1), color = "seagreen3", alpha = 0.5, size = 1.5) +
  geom_line(aes(y = Overlap_measure_baseline3$measure1_baseline), color = "seagreen3") +
  geom_line(aes(y = Overlap_measure_baseline3$measure1_025), color = "seagreen3") +
  geom_line(aes(y = Overlap_measure_baseline3$measure1_975), color = "seagreen3") + 
  geom_line(aes(y = Overlap_measure4$measure1), color = "chocolate3", alpha = 0.5, size = 1) +
  geom_line(aes(y = Overlap_measure_baseline4$measure1_baseline), color = "chocolate3") +
  geom_line(aes(y = Overlap_measure_baseline4$measure1_025), color = "chocolate3") +
  geom_line(aes(y = Overlap_measure_baseline4$measure1_975), color = "chocolate3") + 
  geom_line(aes(y = Overlap_measure5$measure1), color = "royalblue4", alpha = 0.3, size = 0.5) +
  geom_line(aes(y = Overlap_measure_baseline5$measure1_baseline), color = "royalblue4") +
  geom_line(aes(y = Overlap_measure_baseline5$measure1_025), color = "royalblue4") +
  geom_line(aes(y = Overlap_measure_baseline5$measure1_975), color = "royalblue4") + 
  xlab("Gene sets ordered by abs true coefficient") + ylab("S1")

############# OM 2 ################
p2 <- ggplot(Overlap_measure1, aes(1:(5799-100+1))) +  
  geom_line(aes(y = measure2), color = "pink3", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline1$measure2_baseline), color = "pink3") +
  geom_line(aes(y = Overlap_measure_baseline1$measure2_025), color = "pink3") +
  geom_line(aes(y = Overlap_measure_baseline1$measure2_975), color = "pink3") + 
  geom_line(aes(y = Overlap_measure2$measure2), color = "69b3a2", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline2$measure2_baseline), color = "69b3a2") +
  geom_line(aes(y = Overlap_measure_baseline2$measure2_025), color = "69b3a2") +
  geom_line(aes(y = Overlap_measure_baseline2$measure2_975), color = "69b3a2") + 
  geom_line(aes(y = Overlap_measure3$measure2), color = "seagreen3", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline3$measure2_baseline), color = "seagreen3") +
  geom_line(aes(y = Overlap_measure_baseline3$measure2_025), color = "seagreen3") +
  geom_line(aes(y = Overlap_measure_baseline3$measure2_975), color = "seagreen3") + 
  geom_line(aes(y = Overlap_measure4$measure2), color = "chocolate3", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline4$measure2_baseline), color = "chocolate3") +
  geom_line(aes(y = Overlap_measure_baseline4$measure2_025), color = "chocolate3") +
  geom_line(aes(y = Overlap_measure_baseline4$measure2_975), color = "chocolate3") + 
  geom_line(aes(y = Overlap_measure5$measure2), color = "royalblue4", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline5$measure2_baseline), color = "royalblue4") +
  geom_line(aes(y = Overlap_measure_baseline5$measure2_025), color = "royalblue4") +
  geom_line(aes(y = Overlap_measure_baseline5$measure2_975), color = "royalblue4") + 
  xlab("Gene sets ordered by abs true coefficient") + ylab("S2")


############# OM 3 ################
p3 <- ggplot(Overlap_measure1, aes(1:(5799-100+1))) +  
  geom_line(aes(y = measure3), color = "pink3", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline1$measure3_baseline), color = "pink3") +
  geom_line(aes(y = Overlap_measure_baseline1$measure3_025), color = "pink3") +
  geom_line(aes(y = Overlap_measure_baseline1$measure3_975), color = "pink3") + 
  geom_line(aes(y = Overlap_measure2$measure3), color = "69b3a2", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline2$measure3_baseline), color = "69b3a2") +
  geom_line(aes(y = Overlap_measure_baseline2$measure3_025), color = "69b3a2") +
  geom_line(aes(y = Overlap_measure_baseline2$measure3_975), color = "69b3a2") + 
  geom_line(aes(y = Overlap_measure3$measure3), color = "seagreen3", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline3$measure3_baseline), color = "seagreen3") +
  geom_line(aes(y = Overlap_measure_baseline3$measure3_025), color = "seagreen3") +
  geom_line(aes(y = Overlap_measure_baseline3$measure3_975), color = "seagreen3") + 
  geom_line(aes(y = Overlap_measure4$measure3), color = "chocolate3", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline4$measure3_baseline), color = "chocolate3") +
  geom_line(aes(y = Overlap_measure_baseline4$measure3_025), color = "chocolate3") +
  geom_line(aes(y = Overlap_measure_baseline4$measure3_975), color = "chocolate3") + 
  geom_line(aes(y = Overlap_measure5$measure3), color = "royalblue4", alpha = 0.7, size = 0.75) +
  geom_line(aes(y = Overlap_measure_baseline5$measure3_baseline), color = "royalblue4") +
  geom_line(aes(y = Overlap_measure_baseline5$measure3_025), color = "royalblue4") +
  geom_line(aes(y = Overlap_measure_baseline5$measure3_975), color = "royalblue4") + 
  xlab("Gene sets ordered by abs true coefficient") + ylab("S3")

plot_grid(p1, p2, p3, nrow = 3)


