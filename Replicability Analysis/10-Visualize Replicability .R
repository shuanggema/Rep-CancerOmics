######## input ########
# if simulation
setwd("/Users/name/Desktop/TCGA-2021/Data/Simulation-LUAD")
OM1 <- read.csv("./Overlap_measure_marginal_within_samples.csv", row.names = 1)
OM2 <- read.csv("./Overlap_measure_baseline_marginal_within_samples.csv", row.names = 1)
OM <- cbind(OM1,OM2)

# if data analysis
OM <- read.csv("/Users/name/Desktop/TCGA-2021/Data/TCGA_LUAD/Significant Gene/Overlap_measure.csv", row.names = 1)

########
library(ggplot2)
library(dplyr)
library(cowplot)

############# Summary stat of number of genes identified ################
# account for 0s
round(mean(c(OM$no_genes_first, OM$no_genes_last)))
round(sd(c(OM$no_genes_first, OM$no_genes_last)))

round(mean(c(OM$no_path_first, OM$no_path_last)))
round(sd(c(OM$no_path_first, OM$no_path_last)))

# delete all 0s
round(mean(c(OM$no_genes_first, OM$no_genes_last)[c(OM$no_genes_first, OM$no_genes_last) !=0]))
round(sd(c(OM$no_genes_first, OM$no_genes_last)[c(OM$no_genes_first, OM$no_genes_last) !=0]))

round(mean(c(OM$no_path_first, OM$no_path_last)[c(OM$no_genes_first, OM$no_genes_last) !=0]))
round(sd(c(OM$no_path_first, OM$no_path_last)[c(OM$no_genes_first, OM$no_genes_last) !=0]))

############# Summary stat of OMs ################
OM <- OM[OM$no_genes_first > 0 & OM$no_genes_last > 0,]
round(mean(OM$measure1),2); round(sd(OM$measure1),2)
round(mean(OM$measure2),2); round(sd(OM$measure2),2)
round(mean(OM$measure3),2); round(sd(OM$measure3),2)

round(mean(OM$measure1_025),2); round(sd(OM$measure1_025),2)
round(mean(OM$measure1_baseline),2); round(sd(OM$measure1_baseline),2)
round(mean(OM$measure1_950),2); round(sd(OM$measure1_950),2)
round(mean(OM$measure1_975),2); round(sd(OM$measure1_975),2)

round(mean(OM$measure2_025),2); round(sd(OM$measure2_025),2)
round(mean(OM$measure2_baseline),2); round(sd(OM$measure2_baseline),2)
round(mean(OM$measure2_950),2); round(sd(OM$measure2_950),2)
round(mean(OM$measure2_975),2); round(sd(OM$measure2_975),2)

round(mean(OM$measure3_025),2); round(sd(OM$measure3_025),2)
round(mean(OM$measure3_baseline),2); round(sd(OM$measure3_baseline),2)
round(mean(OM$measure3_950),2); round(sd(OM$measure3_950),2)
round(mean(OM$measure3_975),2); round(sd(OM$measure3_975),2)

sum(OM$measure1 < OM$measure1_950)/nrow(OM)
sum(OM$measure1 >= OM$measure1_950)/nrow(OM)

sum(OM$measure1 < OM$measure1_025)/nrow(OM)
sum(OM$measure1 >= OM$measure1_025 & OM$measure1 <= OM$measure1_975)/nrow(OM)
sum(OM$measure1 > OM$measure1_975)/nrow(OM)

sum(OM$measure2 < OM$measure2_950)/nrow(OM)
sum(OM$measure2 >= OM$measure2_950)/nrow(OM)

sum(OM$measure2 < OM$measure2_025)/nrow(OM)
sum(OM$measure2 >= OM$measure2_025 & OM$measure2 <= OM$measure2_975)/nrow(OM)
sum(OM$measure2 > OM$measure2_975)/nrow(OM)

sum(OM$measure3 < OM$measure3_950)/nrow(OM)
sum(OM$measure3 >= OM$measure3_950)/nrow(OM)

sum(OM$measure3 < OM$measure3_025)/nrow(OM)
sum(OM$measure3 >= OM$measure3_025 & OM$measure3 <= OM$measure3_975)/nrow(OM)
sum(OM$measure3 > OM$measure3_975)/nrow(OM)

############# No. genes histogram ################
data.frame(Number = c(OM$no_genes_first, OM$no_genes_last), 
           class = c(rep("First Half",100),rep("Last Half",100))) %>%
  ggplot(aes(x=Number, fill=class)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="") + ggtitle("How many significant genes are selected")

############# OM 1 ################
highlight_df_1 <- OM %>% filter(measure1 >= measure1_025 & measure1 <= measure1_975)
highlight_ind_1 <- c(1:nrow(OM))[OM$measure1 >= OM$measure1_025 & OM$measure1 <= OM$measure1_975]

p1 <- ggplot(OM, aes(1:nrow(OM))) +  
  geom_point(aes(y = measure1, color = "Obs above 97.5%"), size = 1) + 
  geom_point(data = highlight_df_1, aes(x = highlight_ind_1, y = measure1, color="Obs within 95%"), size = 1) + 
  # geom_line(aes(y = measure1_expect, color = "Expected Base")) +
  geom_line(aes(y = measure1_baseline, color = "Mean")) +
  geom_line(aes(y = measure1_025, color = "2.5%")) +
  geom_line(aes(y = measure1_950, color = "95%")) +
  geom_line(aes(y = measure1_975, color = "97.5%")) + 
  geom_ribbon(aes(ymin=measure1_025, ymax=measure1_975), fill="pink2", alpha=0.2) +
  coord_flip() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("") + ylab("S1") + ylim(c(0,1)) +
  scale_color_manual(name = "", 
                     values = c("Obs above 97.5%" = "black", 
                                "Obs within 95%" = "red",
                                # "Expected Base" = "#69b3a2", 
                                "Mean" = "navy", 
                                "2.5%" = "pink3",
                                "95%" = "green3",
                                "97.5%" = "pink3"))


############# OM 2 ################
highlight_df_21 <- OM %>% filter(measure2 >= measure2_025 & measure2 <= measure2_975)
highlight_ind_21 <- c(1:nrow(OM))[OM$measure2 >= OM$measure2_025 & OM$measure2 <= OM$measure2_975]

highlight_df_22 <- OM %>% filter(measure2 < measure2_025)
highlight_ind_22 <- c(1:nrow(OM))[OM$measure2 < OM$measure2_025]

p2 <- ggplot(OM, aes(1:nrow(OM))) +  
  geom_point(aes(y = measure2, color = "Obs above 97.5%"), size = 1) +
  geom_point(data = highlight_df_21, aes(x = highlight_ind_21, y = measure2, color="Obs within 95%"), size = 1) + 
  geom_point(data = highlight_df_22, aes(x = highlight_ind_22, y = measure2, color="Obs below 2.5%"), size = 1) + 
  geom_line(aes(y = measure2_baseline, color = "Mean")) +
  geom_line(aes(y = measure2_025, color = "2.5%")) +
  geom_line(aes(y = measure2_950, color = "95%")) +
  geom_line(aes(y = measure2_975, color = "97.5%")) + 
  geom_ribbon(aes(ymin=measure2_025, ymax=measure2_975), fill="pink2", alpha=0.2) +
  coord_flip() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("") + ylab("S2") + ylim(c(0,1)) +
  scale_color_manual(name = "", 
                     values = c("Obs above 97.5%" = "black", 
                                "Obs within 95%" = "red",
                                "Obs below 2.5%" = "green4",
                                "Mean" = "navy", 
                                "2.5%" = "pink3", 
                                "95%" = "green3",
                                "97.5%" = "pink3"))

############# No. path histogram ################
data.frame(Number = c(OM$no_path_first, OM$no_path_last), 
           class = c(rep("First Half",100),rep("Last Half",100))) %>%
  ggplot(aes(x=Number, fill=class)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() + 
  labs(fill="") + ggtitle("How many significant genes related unique pathways are selected")

############# OM 3 ################
highlight_df_31 <- OM %>% filter(measure3 >= measure3_025 & measure3 <= measure3_975)
highlight_ind_31 <- c(1:nrow(OM))[OM$measure3 >= OM$measure3_025 & OM$measure3 <= OM$measure3_975]

highlight_df_32 <- OM %>% filter(measure3 < measure3_025)
highlight_ind_32 <- c(1:nrow(OM))[OM$measure3 < OM$measure3_025]


p3 <- ggplot(OM, aes(1:nrow(OM))) +  
  geom_point(aes(y = measure3, color = "Obs above 97.5%"), size = 1) +
  geom_point(data = highlight_df_31, aes(x = highlight_ind_31, y = measure3, color = "Obs within 95%"), size = 1) + 
  geom_point(data = highlight_df_32, aes(x = highlight_ind_32, y = measure3, color = "Obs below 2.5%"), size = 1) + 
  geom_line(aes(y = measure3_baseline, color = "Mean")) +
  geom_line(aes(y = measure3_025, color = "2.5%")) +
  geom_line(aes(y = measure3_950, color = "95%")) +
  geom_line(aes(y = measure3_975, color = "97.5%")) + 
  geom_ribbon(aes(ymin = measure3_025, ymax = measure3_975), fill="pink2", alpha=0.2) +
  coord_flip() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("") + ylab("S3") + ylim(c(0,1)) +
  scale_color_manual(name = "", 
                     values = c("Obs above 97.5%" = "black", 
                                "Obs within 95%" = "red",
                                "Obs below 2.5%" = "green4",
                                "Mean" = "navy", 
                                "2.5%" = "pink3", 
                                "95%" = "green3",
                                "97.5%" = "pink3"))


plot_grid(p1 + theme(legend.key = element_rect(fill = "white"), 
                     legend.text = element_text(color = "white"), 
                     legend.title = element_text(color = "white")) + 
                     guides(color = guide_legend(override.aes = list(color = NA))), 
          p2, 
          p3 + theme(legend.key = element_rect(fill = "white"), 
                     legend.text = element_text(color = "white"), 
                     legend.title = element_text(color = "white")) + 
            guides(color = guide_legend(override.aes = list(color = NA))), nrow = 3)


