options(stringsAsFactors = FALSE)

######## input ########
setwd("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Clinic/clinical.cart.2021-09-18")
clinical <- as.data.frame(fread("clinical.tsv"))
exposure <- as.data.frame(fread("exposure.tsv"))
family_history <- as.data.frame(fread("family_history.tsv"))


######## output ########
# File - "Final_demo.csv": contains all related clinical, exposure, family history information; 
#                          data frame format

########
library(data.table)

########################## clean clinical data #############################
colnames(clinical)

# all meanings of variables are searched at https://docs.gdc.cancer.gov/Data_Dictionary/gdcmvs/

# select retained variable
identity_list <- c("case_submitter_id")

demographic_list <- c("age_at_diagnosis","gender","race")

outcome_list <- c("days_to_death","vital_status","days_to_last_follow_up")

med_history_list <- c("prior_malignancy","prior_treatment")

med_current_list <- c("ajcc_pathologic_stage")

treatment_list <- c("treatment_or_therapy","treatment_type")

clinical_filtered <- clinical[,c(identity_list,demographic_list,outcome_list,med_history_list,med_current_list)]
clinical_filtered <- clinical_filtered[!duplicated(clinical_filtered),]
clinical_filtered$radiation_therapy <- clinical$treatment_or_therapy[clinical$treatment_type == "Radiation Therapy, NOS"]
clinical_filtered$pharmaceutical_therapy <- clinical$treatment_or_therapy[clinical$treatment_type == "Pharmaceutical Therapy, NOS"]


Final_demo <- clinical_filtered

########################## process variables one by one #############################

# age: impute missing using median value, as numeric, divided by 365.25, standardize
Final_demo$age_at_diagnosis <- as.numeric(Final_demo$age_at_diagnosis)
Final_demo$age_at_diagnosis[is.na(Final_demo$age_at_diagnosis)] <- median(Final_demo$age_at_diagnosis[!is.na(Final_demo$age_at_diagnosis)])
Final_demo$age_at_diagnosis <- Final_demo$age_at_diagnosis/365.25
Final_demo$age_at_diagnosis <- (Final_demo$age_at_diagnosis - mean(Final_demo$age_at_diagnosis))/
  sd(Final_demo$age_at_diagnosis)

# race: white and non-white
Final_demo$race <- ifelse(Final_demo$race == "white",1,0)
Final_demo$race <- factor(Final_demo$race, levels = c(0,1), labels = c("non-white","white"))

# gender: female and male (only female)
Final_demo$gender <- ifelse(Final_demo$gender == "female",1,0)
Final_demo$gender <- factor(Final_demo$gender, levels = c(0,1), labels = c("male","female"))

# prior_malignancy: yes and no/not reported
Final_demo$prior_malignancy <- ifelse(Final_demo$prior_malignancy == "yes",1,0)
Final_demo$prior_malignancy <- factor(Final_demo$prior_malignancy, levels = c(0,1), labels = c("no or not reported","yes"))

# prior_treatment: yes and no/not reported
Final_demo$prior_treatment <- ifelse(Final_demo$prior_treatment == "Yes",1,0)
Final_demo$prior_treatment <- factor(Final_demo$prior_treatment, levels = c(0,1), labels = c("no or not reported","yes"))
Final_demo <- Final_demo[Final_demo$prior_treatment == "no or not reported",]

# cancer stage: I, II, III or IV
ajcc_pathologic_stage <- Final_demo$ajcc_pathologic_stage
Final_demo$ajcc_pathologic_stage <- 0
Final_demo$ajcc_pathologic_stage[ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB")] <- 1
Final_demo$ajcc_pathologic_stage[ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB")] <- 2
Final_demo$ajcc_pathologic_stage[ajcc_pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB","Stage IIIC","Stage IV")] <- 3
Final_demo <- Final_demo[Final_demo$ajcc_pathologic_stage != 0,]
Final_demo$ajcc_pathologic_stage <- factor(Final_demo$ajcc_pathologic_stage,
                                                                levels = c(1,2,3),
                                                                labels = c("Stage I,IA,IB","Stage II,IIA,IIB","Stage III,IIIA,IIIB,IIIC,IV"))

# therapy: yes, no/not reported
Final_demo$radiation_therapy <- ifelse(Final_demo$radiation_therapy == "yes",1,0)
Final_demo$radiation_therapy <- factor(Final_demo$radiation_therapy, levels = c(0,1), labels = c("no or not reported","yes"))

Final_demo$pharmaceutical_therapy <- ifelse(Final_demo$pharmaceutical_therapy == "yes",1,0)
Final_demo$pharmaceutical_therapy <- factor(Final_demo$pharmaceutical_therapy, levels = c(0,1), labels = c("no or not reported","yes"))

# Survival outcome
status <- rep(0,nrow(Final_demo)) # those alive
status[Final_demo$vital_status == "Dead"] <- 1 # those died
time <- Final_demo$days_to_death
time[status == 0] <- Final_demo$days_to_last_follow_up[status == 0]
data_survival <- data.frame(time, status)
Final_demo <- cbind(Final_demo,data_survival)

Final_demo <- Final_demo[!Final_demo$time == "'--",] # delete 1 subjects missing death time or loss to follow-up time
Final_demo$time <- as.numeric(Final_demo$time)
Final_demo <- Final_demo[Final_demo$time >= 0,] # delete 1 subjects has survival time < 0

########################## clean exposure data #############################
exposure_filtered <- exposure[,c("case_submitter_id","cigarettes_per_day")]
exposure_filtered$cigarettes_per_day[exposure_filtered$cigarettes_per_day == "'--"] <- 0
exposure_filtered$cigarettes_per_day <- as.numeric(exposure_filtered$cigarettes_per_day)

rownames(exposure_filtered) <- exposure_filtered$case_submitter_id
Final_demo$cigarettes_per_day <- exposure_filtered[Final_demo$case_submitter_id,"cigarettes_per_day"]

#########################
Final_demo <- Final_demo[,c("case_submitter_id","age_at_diagnosis","race","gender","prior_malignancy","ajcc_pathologic_stage",
                            "radiation_therapy","pharmaceutical_therapy","cigarettes_per_day",
                            "time","status")]

write.csv(Final_demo,"../Final_demo.csv")





