rm(list=ls())

# * R packages
library(data.table)
library(stringr)
library(tidyr)
library(dplyr)
library(utils)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggstance)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(Hmisc)
library(boot)
library(modelbased)
library(MESS)

# *** Reading and aggregating site simulation data from the first script ----
################################################################################
list <- list.files("./") # location of files
data <- c()
for (i in list) {
  temp = fread(paste0("./",i))
  data <- rbind(data,temp)
}
data <- data %>% mutate_at(c(18:21), as.numeric)
colnames(data)[4]<-"data_type"
rm(temp)

# *** Correction of energy losses by the night duration to get a loss over 24h
################################################################################
data2 <- fread("./data_2_5_t1.csv") # reading the global training dataset
data2 <- subset(data2, !data2$materiel_tot %in% c("Batcorder_-30","SM4_6","SM4_24","SM4_12","SM3_12")) # removing rare material combinations
temp <-  data2 %>% # computing the night proportion over 24h
  dplyr::select(Code_eolienne_meteo,id_site_date,sunrise, sunset) %>%
  group_by(Code_eolienne_meteo,id_site_date,sunrise, sunset)  %>%
  summarise()
temp$night_length <- as.numeric(temp$sunrise - temp$sunset)
temp$percentage <- temp$night_length/24
rm(data2)
temp2 <-  temp %>%
  dplyr::select(Code_eolienne_meteo,percentage) %>%
  group_by(Code_eolienne_meteo)  %>%
  summarise(mean_percentage = mean(percentage))
for (i in unique(data$site)) { # correcting energy losses 
  print(i)
  data[site==i,c(9,12)] <- data[site==i,c(9,12)]*temp2$mean_percentage[temp2$Code_eolienne_meteo==i]
  data[site==i,c(11,14)] <- data[site==i,c(11,14)]/temp2$mean_percentage[temp2$Code_eolienne_meteo==i]
}

# The three following parts starting with "***", named XX, XX and XX, allow to 
# compute the efficiency of an algorithm-based curtailment. All three parts 
# contains the same code, but are applied to a different part of the dataset 
# (either based on rotation losses at all sites, or based on rotation losses at 
# curtailed sites only, or based on energy losses for sites with such data). 
# As a consequence, only the first part and first guild (LRE) is annotated.

# *** Curtailment efficiency based on rotation lost at all sites ----
################################################################################
# LRE guild ----
################################################################################
data2 <- data[which(data$sumOcc_LRE_predict > 0 & data$guild == "LRE" & data$period == "all"),] # subset on simulation made for all periods and containing LRE occurrences
data3 <- data[which(data$sumOcc_LRE_predict > 0 & data$guild == "LRE" & data$period != "all"),] # subset on simulation made by season and containing LRE occurrences
nsites <- length(unique(data2$site)) # storing the number of sites involved
data2 <- data2[-which(data2$proba==0 & data2$protection<1),] # removing artefacts (incomplete protection despite the lowest probability threshold)
data3 <- data3[-which(data3$proba==0 & data3$protection<1),] # the same
# National threshold method ----
################################################################################
# vector of probability thresholds on which the efficiency will be computed
proba_boot <- c(unique(data2$proba[which(data2$proba <= 0.025)]),
                unique(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)][seq(1, length(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)]), 2)]),
                unique(data2$proba[which(data2$proba > 0.05 & data2$proba <=0.10)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 4)]),
                unique(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 8)]),
                unique(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)][seq(1, length(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)]), 16)]),
                unique(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)][seq(1, length(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)]), 32)]),
                unique(data2$proba[which(data2$proba > 0.40)][seq(1, length(data2$proba[which(data2$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]
# vectors for storing protection and loss values
temp_protection <- c()
temp_loss <- c()
cov <- c()
for (i in c("inclusion","none")) { # computing efficiency for both data inclusion methods (inclusion of site-specific data in model training or not)
  print(i)
  for (j in proba_boot) { # for each proba threshold above which wind turbines are stopped
    print(j)
    sub = data2[which(data2$data_type==i & data2$proba==j),] # subset of the right threshold and method
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE))) # boostrapped mean protection and 95% CI on sites with 1,000 iterations
    temp_loss <- as.data.frame(rbind(temp_loss,smean.cl.boot(sub$rotation_loss, conf.int = .95, B = 1000, na.rm = TRUE)))  # boostrapped mean loss and 95% CI on sites with 1,000 iterations
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="LRE"))) # storing related information (method, threshold level, guild)
  }
}
# aggregating vectors and renaming columns
result_national <- cbind(temp_protection,temp_loss,cov,method="national")
colnames(result_national)[1]<-"mean_protection"
colnames(result_national)[2]<-"low_protection"
colnames(result_national)[3]<-"upp_protection"
colnames(result_national)[4]<-"mean_loss"
colnames(result_national)[5]<-"low_loss"
colnames(result_national)[6]<-"upp_loss"
# set max CI to 1 for graphical representation
for(i in 1:nrow(result_national)){
  if(result_national$upp_protection[i]>1){result_national$upp_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# Site-specific threshold method ----
################################################################################
# round of rotation losses on which the average of the protection will be done
data2$rotation_loss_round <- round(data2$rotation_loss,2)
# vectors for storing efficiency
temp_protection <- c()
mean_loss <- c()
# vectors for storing involved threshold probability by site for a given efficiency level
proba_mean <- c()
proba_sd <- c()
proba_min <- c()
proba_max <- c()
cov <- c()
for (i in c("inclusion","none")) {  # computing efficiency for both data inclusion methods (inclusion of site-specific data in model training or not)
  print(i)
  for (j in unique(data2$rotation_loss_round)) { # for each rotation loss level
    print(j)
    sub = data2[which(data2$data_type==i & data2$rotation_loss_round==j),] # subset of the right loss level and method
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE))) # boostrapped mean protection and 95% CI on sites with 1,000 iterations
    mean_loss <- c(mean_loss,j) # storing the rotation loss j
    proba_mean <- c(proba_mean,mean(sub$proba)) # storing the mean probability threshold required for a given site
    proba_sd <- c(proba_sd,sd(sub$proba)) # storing the associated standard deviation 
    proba_min <- c(proba_min,min(sub$proba)) # storing the minimum probability threshold required among site
    proba_max <- c(proba_max,max(sub$proba)) # storing the maximum probability threshold required among site
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="LRE"))) # storing related information (method, threshold level, guild)
  }
}
# aggregating vectors and renaming columns
result_site <- cbind(temp_protection,mean_loss,proba_mean,proba_sd,proba_min,proba_max,cov,method="site")
colnames(result_site)[1]<-"mean_protection"
colnames(result_site)[2]<-"low_protection"
colnames(result_site)[3]<-"upp_protection"
# removing NAs and ordering by loss level
result_site <- result_site %>% drop_na()
result_site <- result_site[order(result_site$data_type,result_site$mean_loss, decreasing = T), ]
# smooting the protection value across rotation loss levels for graphical representation
result_site$mean_protection[which(result_site$data_type=="none")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="none")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="none")] <- smoothing(result_site$low_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$mean_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$low_protection[which(result_site$data_type=="inclusion")], method = "smooth")
# set max CI to 1 for graphical representation
for(i in 1:nrow(result_site)){
  if(result_site$upp_protection[i]>1){result_site$upp_protection[i]=1}
}

# Season-specific threshold method  ----
################################################################################
# vector of probability thresholds on which the efficiency will be computed
proba_boot <- c(unique(data3$proba[which(data3$proba <= 0.025)]),
                unique(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)][seq(1, length(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)]), 2)]),
                unique(data3$proba[which(data3$proba > 0.05 & data3$proba <=0.10)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 4)]),
                unique(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 8)]),
                unique(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)][seq(1, length(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)]), 16)]),
                unique(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)][seq(1, length(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)]), 32)]),
                unique(data3$proba[which(data3$proba > 0.40)][seq(1, length(data3$proba[which(data3$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]
# dataframe for storing results
result_national_season <- c()
for(s in unique(data3$period)){ # for each season
  print(s)
  # vectors for storing the % of protection and loss, and the number of rotation lost and the number of bat passes protected
  temp_protection <- c()
  temp_loss <- c()
  temp_nb_loss <- c()
  temp_nb_passes <- c()
  cov <- c()
  for (i in c("inclusion","none")) { # computing efficiency for both data inclusion methods (inclusion of site-specific data in model training or not)
    print(i)
    for (j in proba_boot) { # for each proba threshold above which wind turbines are stopped
      sub = data3[which(data3$data_type==i & data3$proba==j & data3$period==s),] # subset if the right period, method, and proba threshold
      # temporary vectors in which protection and loss are computed 1,000 times based on randomly selected sites
      boot_protection <- c()
      boot_loss <- c()
      boot_nb_loss <- c()
      boot_nb_loss_tot <- c()
      boot_nb_passes <- c()
      boot_nb_passes_tot <- c()
      for(l in 1:1000){ # computing 1,000 times protection ad loss based on randomly selected sites
              ind <- sample(1:nrow(sub), replace = TRUE) # random selection of sites
              boot_protection <- c(boot_protection,sum(sub$nb_protection[ind], na.rm = T)/sum(sub$nb_passes[ind], na.rm = T)) # computing protection on randomly selected sites
              boot_loss <- c(boot_loss,sum(sub$nb_rotation_loss[ind], na.rm = T)/sum(sub$nb_rotation[ind], na.rm = T)) # computing loss on randomly selected sites
              boot_nb_loss <- c(boot_nb_loss,sum(sub$nb_rotation_loss[ind], na.rm = T))  # computing the number of rotation lost on randomly selected sites
              boot_nb_loss_tot <- c(boot_nb_loss_tot,sum(sub$nb_rotation[ind], na.rm = T)) # computing the total number of rotations on randomly selected sites
              boot_nb_passes <- c(boot_nb_passes,sum(sub$nb_protection[ind], na.rm = T)) # computing the number of bat detections protected on randomly selected sites
              boot_nb_passes_tot <- c(boot_nb_passes_tot,sum(sub$nb_passes[ind], na.rm = T)) # computing the total number of bat detections on randomly selected sites
      }
      # computing the mean and 95% CI on the 1,000 previous iterations
      temp_protection <- rbind(temp_protection,cbind(mean(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)-1.96*sd(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)+1.96*sd(boot_protection, na.rm = T)))
      temp_loss <- rbind(temp_loss,cbind(mean(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)-1.96*sd(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)+1.96*sd(boot_loss, na.rm = T)))
      temp_nb_loss <- rbind(temp_nb_loss,cbind(mean(boot_nb_loss, na.rm = T),mean(boot_nb_loss, na.rm = T)-1.96*sd(boot_nb_loss, na.rm = T),mean(boot_nb_loss, na.rm = T)+1.96*sd(boot_nb_loss, na.rm = T),
                                               mean(boot_nb_loss_tot, na.rm = T),mean(boot_nb_loss_tot, na.rm = T)-1.96*sd(boot_nb_loss_tot, na.rm = T),mean(boot_nb_loss_tot, na.rm = T)+1.96*sd(boot_nb_loss_tot, na.rm = T)))
      temp_nb_passes <- rbind(temp_nb_passes,cbind(mean(boot_nb_passes, na.rm = T),mean(boot_nb_passes, na.rm = T)-1.96*sd(boot_nb_passes, na.rm = T),mean(boot_nb_passes, na.rm = T)+1.96*sd(boot_nb_passes, na.rm = T),
                                               mean(boot_nb_passes_tot, na.rm = T),mean(boot_nb_passes_tot, na.rm = T)-1.96*sd(boot_nb_passes_tot, na.rm = T),mean(boot_nb_passes_tot, na.rm = T)+1.96*sd(boot_nb_passes_tot, na.rm = T)))
      cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="LRE",method="national_season",season=s))) # storing associated information (methods, threshold, guild, season)
    }
  }
  # aggregating in a sigle dataframe
  result_national_season <- rbind(result_national_season,cbind(temp_protection,temp_loss,temp_nb_loss,temp_nb_passes,cov))
}
# renaming columns
colnames(result_national_season)[1]<-"mean_protection"
colnames(result_national_season)[2]<-"low_protection"
colnames(result_national_season)[3]<-"upp_protection"
colnames(result_national_season)[4]<-"mean_loss"
colnames(result_national_season)[5]<-"low_loss"
colnames(result_national_season)[6]<-"upp_loss"
colnames(result_national_season)[7]<-"mean_nb_loss"
colnames(result_national_season)[8]<-"low_nb_loss"
colnames(result_national_season)[9]<-"upp_nb_loss"
colnames(result_national_season)[10]<-"mean_nb_loss_tot"
colnames(result_national_season)[11]<-"low_nb_loss_tot"
colnames(result_national_season)[12]<-"upp_nb_loss_tot"
colnames(result_national_season)[13]<-"mean_nb_passes"
colnames(result_national_season)[14]<-"low_nb_passes"
colnames(result_national_season)[15]<-"upp_nb_passes"
colnames(result_national_season)[16]<-"mean_nb_passes_tot"
colnames(result_national_season)[17]<-"low_nb_passes_tot"
colnames(result_national_season)[18]<-"upp_nb_passes_tot"
# removing -Inf values
result_national_season <- result_national_season %>% 
  filter_all(all_vars(!is.infinite(.)))
# creating vectors for storing probability thresholds required for each season, method, and rotation loss level
result_national_season$rotation_loss_round <- round(result_national_season$mean_loss,2)
mean_protection <- c()
low_protection <- c()
upp_protection <- c()
round_loss <- c()
proba_mean_winter <- c()
proba_mean_spring <- c()
proba_mean_summer <- c()
proba_mean_autumn <- c()
proba_sd_winter <- c()
proba_sd_spring <- c()
proba_sd_summer <- c()
proba_sd_autumn <- c()
proba_min_winter <- c()
proba_min_spring <- c()
proba_min_summer <- c()
proba_min_autumn <- c()
proba_max_winter <- c()
proba_max_spring <- c()
proba_max_summer <- c()
proba_max_autumn <- c()
cov <- c()
for (i in c("inclusion","none")) { # storing probability thresholds for both data inclusion methods (inclusion of site-specific data in model training or not)
  print(i)
  for (j in unique(result_national_season$rotation_loss_round)) { # for each rotation loss level
    sub = result_national_season[which(result_national_season$data_type==i & result_national_season$rotation_loss_round==j),] # selecting the right season and rotation loss level
    mean_protection <- c(mean_protection,mean(sub$mean_protection)) # computing the mean protection associated with the rotation loss level j
    low_protection <- c(low_protection,mean(sub$low_protection)) # computing the lower 95% CI protection associated with the rotation loss level j
    upp_protection <- c(upp_protection,mean(sub$upp_protection)) # computing the upper 95% CI protection associated with the rotation loss level j
    round_loss <- c(round_loss,j) # storing the rotation loss level
    # computing mean, standard deviation, min, and max of the required probability threshold to reach the rotation 
    # loss level j and for each season separately
    proba_mean_winter <- c(proba_mean_winter,mean(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_mean_spring <- c(proba_mean_spring,mean(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_mean_summer <- c(proba_mean_summer,mean(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_mean_autumn <- c(proba_mean_autumn,mean(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_sd_winter <- c(proba_sd_winter,sd(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_sd_spring <- c(proba_sd_spring,sd(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_sd_summer <- c(proba_sd_summer,sd(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_sd_autumn <- c(proba_sd_autumn,sd(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_min_winter <- c(proba_min_winter,min(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_min_spring <- c(proba_min_spring,min(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_min_summer <- c(proba_min_summer,min(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_min_autumn <- c(proba_min_autumn,min(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_max_winter <- c(proba_max_winter,max(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_max_spring <- c(proba_max_spring,max(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_max_summer <- c(proba_max_summer,max(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_max_autumn <- c(proba_max_autumn,max(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="LRE")))
  }
}
# aggregating in a single datafame
result_national_season_avg <- cbind(mean_protection,low_protection,upp_protection,round_loss,proba_mean_winter,proba_mean_spring,proba_mean_summer,proba_mean_autumn,
                                    proba_sd_winter,proba_sd_spring,proba_sd_summer,proba_sd_autumn,proba_min_winter,proba_min_spring,proba_min_summer,proba_min_autumn,
                                    proba_max_winter,proba_max_spring,proba_max_summer,proba_max_autumn,cov,method="national_season")
# removing NAs and ordering by loss level
result_national_season_avg <- result_national_season_avg %>% drop_na(mean_protection)
result_national_season_avg <- result_national_season_avg[order(result_national_season_avg$data_type,result_national_season_avg$round_loss, decreasing = T), ]
# smooting the protection value across rotation loss levels for graphical representation
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
# set min and max CI to 0 and 1 for graphical representation
  for(i in 1:nrow(result_national_season_avg)){
    if(result_national_season_avg$upp_protection[i]>1){result_national_season_avg$upp_protection[i]=1}
    if(result_national_season_avg$low_protection[i]<0){result_national_season_avg$low_protection[i]=0}
    if(result_national_season_avg$mean_protection[i]>1){result_national_season_avg$mean_protection[i]=1}
  }
  rm(temp_loss,temp_protection,sub,cov)
  
# merging results from the 3 methods in a single dataframe ----
################################################################################
result_site_merge <- result_site[,-c(6:8)]
colnames(result_site_merge)[5] <- "proba"
result_national_merge <- result_national[,-c(5:6)]
result_national_season_avg_merge <- result_national_season_avg[,-c(6:20)]
colnames(result_national_season_avg_merge)[4] <- "mean_loss"
colnames(result_national_season_avg_merge)[5] <- "proba"
result_total_LRE <- rbind(result_national_merge,result_site_merge,result_national_season_avg_merge)
result_site_LRE <- result_site
result_national_LRE <- result_national
result_national_season_LRE <- result_national_season
result_national_season_avg_LRE <- result_national_season_avg
# plots ----
################################################################################
LRE_inclusion <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="inclusion"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(C) LRE - Inclusion of site-specific data to training")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

LRE_none <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(D) LRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"),legend.title = element_text(size = 9),legend.text = element_text(size = 8))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

################################################################################
# MRE guild ----
################################################################################
data2 <- data[which(data$sumOcc_MRE_predict > 0 & data$guild == "MRE" & data$period == "all"),]
data3 <- data[which(data$sumOcc_MRE_predict > 0 & data$guild == "MRE" & data$period %in% c("spring","summer","autumn")),]
nsites <- length(unique(data2$site))
data2 <- data2[-which(data2$proba==0 & data2$protection<1),]
data3 <- data3[-which(data3$proba==0 & data3$protection<1),]
# National threshold method ----
################################################################################
proba_boot <- c(unique(data2$proba[which(data2$proba <= 0.025)]),
                unique(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)][seq(1, length(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)]), 2)]),
                unique(data2$proba[which(data2$proba > 0.05 & data2$proba <=0.10)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 4)]),
                unique(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 8)]),
                unique(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)][seq(1, length(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)]), 16)]),
                unique(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)][seq(1, length(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)]), 32)]),
                unique(data2$proba[which(data2$proba > 0.40)][seq(1, length(data2$proba[which(data2$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]
temp_protection <- c()
temp_loss <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in proba_boot) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$proba==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    temp_loss <- as.data.frame(rbind(temp_loss,smean.cl.boot(sub$rotation_loss, conf.int = .95, B = 1000, na.rm = TRUE)))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="MRE")))
  }
}
result_national <- cbind(temp_protection,temp_loss,cov,method="national")
colnames(result_national)[1]<-"mean_protection"
colnames(result_national)[2]<-"low_protection"
colnames(result_national)[3]<-"upp_protection"
colnames(result_national)[4]<-"mean_loss"
colnames(result_national)[5]<-"low_loss"
colnames(result_national)[6]<-"upp_loss"

for(i in 1:nrow(result_national)){
  if(result_national$upp_protection[i]>1){result_national$upp_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# Site-specific threshold method ----
################################################################################
data2$rotation_loss_round <- round(data2$rotation_loss,2)
temp_protection <- c()
mean_loss <- c()
proba_mean <- c()
proba_sd <- c()
proba_min <- c()
proba_max <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(data2$rotation_loss_round)) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$rotation_loss_round==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    mean_loss <- c(mean_loss,j)
    proba_mean <- c(proba_mean,mean(sub$proba))
    proba_sd <- c(proba_sd,sd(sub$proba))
    proba_min <- c(proba_min,min(sub$proba))
    proba_max <- c(proba_max,max(sub$proba))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="MRE")))
  }
}
result_site <- cbind(temp_protection,mean_loss,proba_mean,proba_sd,proba_min,proba_max,cov,method="site")
colnames(result_site)[1]<-"mean_protection"
colnames(result_site)[2]<-"low_protection"
colnames(result_site)[3]<-"upp_protection"
result_site <- result_site %>% drop_na()
result_site <- result_site[order(result_site$data_type,result_site$mean_loss, decreasing = T), ]
result_site$mean_protection[which(result_site$data_type=="none")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="none")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="none")] <- smoothing(result_site$low_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$mean_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$low_protection[which(result_site$data_type=="inclusion")], method = "smooth")

for(i in 1:nrow(result_site)){
  if(result_site$upp_protection[i]>1){result_site$upp_protection[i]=1}
}

# Season-specific threshold method ----
################################################################################
proba_boot <- c(unique(data3$proba[which(data3$proba <= 0.025)]),
                unique(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)][seq(1, length(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)]), 2)]),
                unique(data3$proba[which(data3$proba > 0.05 & data3$proba <=0.10)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 4)]),
                unique(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 8)]),
                unique(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)][seq(1, length(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)]), 16)]),
                unique(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)][seq(1, length(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)]), 32)]),
                unique(data3$proba[which(data3$proba > 0.40)][seq(1, length(data3$proba[which(data3$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]

result_national_season <- c()
for(s in unique(data3$period)){
  print(s)
  temp_protection <- c()
  temp_loss <- c()
  temp_nb_loss <- c()
  temp_nb_passes <- c()
  cov <- c()
  for (i in c("inclusion","none")) {
    print(i)
    for (j in proba_boot) {
      
      sub = data3[which(data3$data_type==i & data3$proba==j & data3$period==s),]
      boot_protection <- c()
      boot_loss <- c()
      boot_nb_loss <- c()
      boot_nb_loss_tot <- c()
      boot_nb_passes <- c()
      boot_nb_passes_tot <- c()
      
      for(l in 1:1000){
        
        ind <- sample(1:nrow(sub), replace = TRUE)
        boot_protection <- c(boot_protection,sum(sub$nb_protection[ind], na.rm = T)/sum(sub$nb_passes[ind], na.rm = T))
        boot_loss <- c(boot_loss,sum(sub$nb_rotation_loss[ind], na.rm = T)/sum(sub$nb_rotation[ind], na.rm = T))
        boot_nb_loss <- c(boot_nb_loss,sum(sub$nb_rotation_loss[ind], na.rm = T))
        boot_nb_loss_tot <- c(boot_nb_loss_tot,sum(sub$nb_rotation[ind], na.rm = T))
        boot_nb_passes <- c(boot_nb_passes,sum(sub$nb_protection[ind], na.rm = T))
        boot_nb_passes_tot <- c(boot_nb_passes_tot,sum(sub$nb_passes[ind], na.rm = T))
        
      }
      temp_protection <- rbind(temp_protection,cbind(mean(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)-1.96*sd(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)+1.96*sd(boot_protection, na.rm = T)))
      temp_loss <- rbind(temp_loss,cbind(mean(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)-1.96*sd(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)+1.96*sd(boot_loss, na.rm = T)))
      temp_nb_loss <- rbind(temp_nb_loss,cbind(mean(boot_nb_loss, na.rm = T),mean(boot_nb_loss, na.rm = T)-1.96*sd(boot_nb_loss, na.rm = T),mean(boot_nb_loss, na.rm = T)+1.96*sd(boot_nb_loss, na.rm = T),
                                               mean(boot_nb_loss_tot, na.rm = T),mean(boot_nb_loss_tot, na.rm = T)-1.96*sd(boot_nb_loss_tot, na.rm = T),mean(boot_nb_loss_tot, na.rm = T)+1.96*sd(boot_nb_loss_tot, na.rm = T)))
      temp_nb_passes <- rbind(temp_nb_passes,cbind(mean(boot_nb_passes, na.rm = T),mean(boot_nb_passes, na.rm = T)-1.96*sd(boot_nb_passes, na.rm = T),mean(boot_nb_passes, na.rm = T)+1.96*sd(boot_nb_passes, na.rm = T),
                                                   mean(boot_nb_passes_tot, na.rm = T),mean(boot_nb_passes_tot, na.rm = T)-1.96*sd(boot_nb_passes_tot, na.rm = T),mean(boot_nb_passes_tot, na.rm = T)+1.96*sd(boot_nb_passes_tot, na.rm = T)))
            cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="MRE",method="national_season",season=s)))
    }
  }
  result_national_season <- rbind(result_national_season,cbind(temp_protection,temp_loss,temp_nb_loss,temp_nb_passes,cov))
}
colnames(result_national_season)[1]<-"mean_protection"
colnames(result_national_season)[2]<-"low_protection"
colnames(result_national_season)[3]<-"upp_protection"
colnames(result_national_season)[4]<-"mean_loss"
colnames(result_national_season)[5]<-"low_loss"
colnames(result_national_season)[6]<-"upp_loss"
colnames(result_national_season)[7]<-"mean_nb_loss"
colnames(result_national_season)[8]<-"low_nb_loss"
colnames(result_national_season)[9]<-"upp_nb_loss"
colnames(result_national_season)[10]<-"mean_nb_loss_tot"
colnames(result_national_season)[11]<-"low_nb_loss_tot"
colnames(result_national_season)[12]<-"upp_nb_loss_tot"
colnames(result_national_season)[13]<-"mean_nb_passes"
colnames(result_national_season)[14]<-"low_nb_passes"
colnames(result_national_season)[15]<-"upp_nb_passes"
colnames(result_national_season)[16]<-"mean_nb_passes_tot"
colnames(result_national_season)[17]<-"low_nb_passes_tot"
colnames(result_national_season)[18]<-"upp_nb_passes_tot"
result_national_season <- result_national_season %>% 
  filter_all(all_vars(!is.infinite(.)))

result_national_season$rotation_loss_round <- round(result_national_season$mean_loss,2)
mean_protection <- c()
low_protection <- c()
upp_protection <- c()
round_loss <- c()
# proba_mean_winter <- c()
proba_mean_spring <- c()
proba_mean_summer <- c()
proba_mean_autumn <- c()
# proba_sd_winter <- c()
proba_sd_spring <- c()
proba_sd_summer <- c()
proba_sd_autumn <- c()
# proba_min_winter <- c()
proba_min_spring <- c()
proba_min_summer <- c()
proba_min_autumn <- c()
# proba_max_winter <- c()
proba_max_spring <- c()
proba_max_summer <- c()
proba_max_autumn <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(result_national_season$rotation_loss_round)) {
    sub = result_national_season[which(result_national_season$data_type==i & result_national_season$rotation_loss_round==j),]
    mean_protection <- c(mean_protection,mean(sub$mean_protection))
    low_protection <- c(low_protection,mean(sub$low_protection))
    upp_protection <- c(upp_protection,mean(sub$upp_protection))
    round_loss <- c(round_loss,j)
    # proba_mean_winter <- c(proba_mean_winter,mean(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_mean_spring <- c(proba_mean_spring,mean(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_mean_summer <- c(proba_mean_summer,mean(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_mean_autumn <- c(proba_mean_autumn,mean(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    # proba_sd_winter <- c(proba_sd_winter,sd(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_sd_spring <- c(proba_sd_spring,sd(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_sd_summer <- c(proba_sd_summer,sd(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_sd_autumn <- c(proba_sd_autumn,sd(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    # proba_min_winter <- c(proba_min_winter,min(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_min_spring <- c(proba_min_spring,min(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_min_summer <- c(proba_min_summer,min(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_min_autumn <- c(proba_min_autumn,min(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    # proba_max_winter <- c(proba_max_winter,max(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_max_spring <- c(proba_max_spring,max(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_max_summer <- c(proba_max_summer,max(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_max_autumn <- c(proba_max_autumn,max(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="MRE")))
  }
}
result_national_season_avg <- cbind(mean_protection,low_protection,upp_protection,round_loss,proba_mean_spring,proba_mean_summer,proba_mean_autumn,
                                    proba_sd_spring,proba_sd_summer,proba_sd_autumn,proba_min_spring,proba_min_summer,proba_min_autumn,
                                    proba_max_spring,proba_max_summer,proba_max_autumn,cov,method="national_season")

result_national_season_avg <- result_national_season_avg %>% drop_na(mean_protection)
result_national_season_avg <- result_national_season_avg[order(result_national_season_avg$data_type,result_national_season_avg$round_loss, decreasing = T), ]
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")

for(i in 1:nrow(result_national_season_avg)){
  if(result_national_season_avg$upp_protection[i]>1){result_national_season_avg$upp_protection[i]=1}
  if(result_national_season_avg$low_protection[i]<0){result_national_season_avg$low_protection[i]=0}
  if(result_national_season_avg$mean_protection[i]>1){result_national_season_avg$mean_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# merging results from the 3 methods in a single dataframe ----
################################################################################
result_site_merge <- result_site[,-c(6:8)]
colnames(result_site_merge)[5] <- "proba"
result_national_merge <- result_national[,-c(5:6)]
result_national_season_avg_merge <- result_national_season_avg[,-c(6:16)]
colnames(result_national_season_avg_merge)[4] <- "mean_loss"
colnames(result_national_season_avg_merge)[5] <- "proba"
result_total_MRE <- rbind(result_national_merge,result_site_merge,result_national_season_avg_merge)
result_site_MRE <- result_site
result_national_MRE <- result_national
result_national_season_MRE <- result_national_season
result_national_season_avg_MRE <- result_national_season_avg
# plots ----
################################################################################
MRE_inclusion <- ggplot(data = result_total_MRE[which(result_total_MRE$data_type=="inclusion"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(A) MRE - Inclusion of site-specific data to training")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

MRE_none <- ggplot(data = result_total_MRE[which(result_total_MRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(B) MRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

# Saving figures and results ----
################################################################################
legend <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(D) LRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "right",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

df_leg <- get_legend(legend)
tiff("./Efficiency_rotations_bootstrap1000_smooth.tiff", width = 1300, height = 1400, units = "px",res = 150)
grid.arrange(arrangeGrob(MRE_inclusion,MRE_none,LRE_inclusion,LRE_none, nrow = 2),
             arrangeGrob(nullGrob(), df_leg, nullGrob(), nrow = 1), ncol = 1, heights = c(4,1),
             top = textGrob("", gp = gpar(fotsize = 12, font = 2)))
dev.off()

write_csv(rbind(result_national_LRE,result_national_MRE),"./ALL_result_national.csv")
write_csv(rbind(result_national_season_avg_LRE,cbind(result_national_season_avg_MRE,proba_mean_winter=NA,proba_sd_winter=NA,proba_min_winter=NA,proba_max_winter=NA)),"./ALL_result_national_season_avg.csv")
write_csv(rbind(result_national_season_LRE,result_national_season_MRE),"./ALL_result_national_season_nb passes_rotations.csv")
write_csv(rbind(result_site_LRE,result_site_MRE),"./ALL_result_site.csv")


# *** Curtailment efficiency based on rotation lost at non-curtailed sites ----
################################################################################
# LRE guild ----
################################################################################
data2 <- data[which(data$sumOcc_LRE_predict > 0 & data$guild == "LRE" & data$period == "all" & data$bridage_presence == "non"),]
data3 <- data[which(data$sumOcc_LRE_predict > 0 & data$guild == "LRE" & data$period != "all" & data$bridage_presence == "non"),]
nsites <- length(unique(data2$site))
data2 <- data2[-which(data2$proba==0 & data2$protection<1),]
data3 <- data3[-which(data3$proba==0 & data3$protection<1),]
# National threshold method ----
################################################################################
proba_boot <- c(unique(data2$proba[which(data2$proba <= 0.025)]),
                unique(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)][seq(1, length(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)]), 2)]),
                unique(data2$proba[which(data2$proba > 0.05 & data2$proba <=0.10)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 4)]),
                unique(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 8)]),
                unique(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)][seq(1, length(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)]), 16)]),
                unique(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)][seq(1, length(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)]), 32)]),
                unique(data2$proba[which(data2$proba > 0.40)][seq(1, length(data2$proba[which(data2$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]
temp_protection <- c()
temp_loss <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in proba_boot) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$proba==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    temp_loss <- as.data.frame(rbind(temp_loss,smean.cl.boot(sub$rotation_loss, conf.int = .95, B = 1000, na.rm = TRUE)))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="LRE")))
  }
}
result_national <- cbind(temp_protection,temp_loss,cov,method="national")
colnames(result_national)[1]<-"mean_protection"
colnames(result_national)[2]<-"low_protection"
colnames(result_national)[3]<-"upp_protection"
colnames(result_national)[4]<-"mean_loss"
colnames(result_national)[5]<-"low_loss"
colnames(result_national)[6]<-"upp_loss"

for(i in 1:nrow(result_national)){
  if(result_national$upp_protection[i]>1){result_national$upp_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# Site-specific threshold method ----
################################################################################
data2$rotation_loss_round <- round(data2$rotation_loss,2)
temp_protection <- c()
mean_loss <- c()
proba_mean <- c()
proba_sd <- c()
proba_min <- c()
proba_max <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(data2$rotation_loss_round)) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$rotation_loss_round==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    mean_loss <- c(mean_loss,j)
    proba_mean <- c(proba_mean,mean(sub$proba))
    proba_sd <- c(proba_sd,sd(sub$proba))
    proba_min <- c(proba_min,min(sub$proba))
    proba_max <- c(proba_max,max(sub$proba))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="LRE")))
  }
}
result_site <- cbind(temp_protection,mean_loss,proba_mean,proba_sd,proba_min,proba_max,cov,method="site")
colnames(result_site)[1]<-"mean_protection"
colnames(result_site)[2]<-"low_protection"
colnames(result_site)[3]<-"upp_protection"
result_site <- result_site %>% drop_na()
result_site <- result_site[order(result_site$data_type,result_site$mean_loss, decreasing = T), ]
result_site$mean_protection[which(result_site$data_type=="none")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="none")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="none")] <- smoothing(result_site$low_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$mean_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$low_protection[which(result_site$data_type=="inclusion")], method = "smooth")

# result_site$mean_protection[which(result_site$data_type=="none")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="none")], method = "loess", strength = 0.1)
# result_site$upp_protection[which(result_site$data_type=="none")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="none")], method = "loess", strength = 0.1)
# result_site$low_protection[which(result_site$data_type=="none")] <- smoothing(result_site$low_protection[which(result_site$data_type=="none")], method = "loess", strength = 0.1)
# result_site$mean_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="inclusion")], method = "loess", strength = 0.1)
# result_site$upp_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="inclusion")], method = "loess", strength = 0.1)
# result_site$low_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$low_protection[which(result_site$data_type=="inclusion")], method = "loess", strength = 0.1)

for(i in 1:nrow(result_site)){
  if(result_site$upp_protection[i]>1){result_site$upp_protection[i]=1}
}

# Season-specific threshold method  ----
################################################################################
proba_boot <- c(unique(data3$proba[which(data3$proba <= 0.025)]),
                unique(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)][seq(1, length(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)]), 2)]),
                unique(data3$proba[which(data3$proba > 0.05 & data3$proba <=0.10)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 4)]),
                unique(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 8)]),
                unique(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)][seq(1, length(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)]), 16)]),
                unique(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)][seq(1, length(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)]), 32)]),
                unique(data3$proba[which(data3$proba > 0.40)][seq(1, length(data3$proba[which(data3$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]

result_national_season <- c()
for(s in unique(data3$period)){
  print(s)
  temp_protection <- c()
  temp_loss <- c()
  cov <- c()
  for (i in c("inclusion","none")) {
    print(i)
    for (j in proba_boot) {
      
      sub = data3[which(data3$data_type==i & data3$proba==j & data3$period==s),]
      boot_protection <- c()
      boot_loss <- c()
      
      for(l in 1:1000){
        
        ind <- sample(1:nrow(sub), replace = TRUE)
        boot_protection <- c(boot_protection,sum(sub$nb_protection[ind], na.rm = T)/sum(sub$nb_passes[ind], na.rm = T))
        boot_loss <- c(boot_loss,sum(sub$nb_rotation_loss[ind], na.rm = T)/sum(sub$nb_rotation[ind], na.rm = T))
        
      }
      temp_protection <- rbind(temp_protection,cbind(mean(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)-1.96*sd(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)+1.96*sd(boot_protection, na.rm = T)))
      temp_loss <- rbind(temp_loss,cbind(mean(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)-1.96*sd(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)+1.96*sd(boot_loss, na.rm = T)))
      
      cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="LRE",method="national_season",season=s)))
    }
  }
  result_national_season <- rbind(result_national_season,cbind(temp_protection,temp_loss,cov))
}
colnames(result_national_season)[1]<-"mean_protection"
colnames(result_national_season)[2]<-"low_protection"
colnames(result_national_season)[3]<-"upp_protection"
colnames(result_national_season)[4]<-"mean_loss"
colnames(result_national_season)[5]<-"low_loss"
colnames(result_national_season)[6]<-"upp_loss"
result_national_season <- result_national_season %>% 
  filter_all(all_vars(!is.infinite(.)))

result_national_season$rotation_loss_round <- round(result_national_season$mean_loss,2)
mean_protection <- c()
low_protection <- c()
upp_protection <- c()
round_loss <- c()
proba_mean_winter <- c()
proba_mean_spring <- c()
proba_mean_summer <- c()
proba_mean_autumn <- c()
proba_sd_winter <- c()
proba_sd_spring <- c()
proba_sd_summer <- c()
proba_sd_autumn <- c()
proba_min_winter <- c()
proba_min_spring <- c()
proba_min_summer <- c()
proba_min_autumn <- c()
proba_max_winter <- c()
proba_max_spring <- c()
proba_max_summer <- c()
proba_max_autumn <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(result_national_season$rotation_loss_round)) {
    sub = result_national_season[which(result_national_season$data_type==i & result_national_season$rotation_loss_round==j),]
    mean_protection <- c(mean_protection,mean(sub$mean_protection))
    low_protection <- c(low_protection,mean(sub$low_protection))
    upp_protection <- c(upp_protection,mean(sub$upp_protection))
    round_loss <- c(round_loss,j)
    proba_mean_winter <- c(proba_mean_winter,mean(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_mean_spring <- c(proba_mean_spring,mean(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_mean_summer <- c(proba_mean_summer,mean(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_mean_autumn <- c(proba_mean_autumn,mean(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_sd_winter <- c(proba_sd_winter,sd(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_sd_spring <- c(proba_sd_spring,sd(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_sd_summer <- c(proba_sd_summer,sd(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_sd_autumn <- c(proba_sd_autumn,sd(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_min_winter <- c(proba_min_winter,min(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_min_spring <- c(proba_min_spring,min(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_min_summer <- c(proba_min_summer,min(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_min_autumn <- c(proba_min_autumn,min(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_max_winter <- c(proba_max_winter,max(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_max_spring <- c(proba_max_spring,max(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_max_summer <- c(proba_max_summer,max(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_max_autumn <- c(proba_max_autumn,max(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="LRE")))
  }
}
result_national_season_avg <- cbind(mean_protection,low_protection,upp_protection,round_loss,proba_mean_winter,proba_mean_spring,proba_mean_summer,proba_mean_autumn,
                                    proba_sd_winter,proba_sd_spring,proba_sd_summer,proba_sd_autumn,proba_min_winter,proba_min_spring,proba_min_summer,proba_min_autumn,
                                    proba_max_winter,proba_max_spring,proba_max_summer,proba_max_autumn,cov,method="national_season")

result_national_season_avg <- result_national_season_avg %>% drop_na(mean_protection)
result_national_season_avg <- result_national_season_avg[order(result_national_season_avg$data_type,result_national_season_avg$round_loss, decreasing = T), ]
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")

for(i in 1:nrow(result_national_season_avg)){
  if(result_national_season_avg$upp_protection[i]>1){result_national_season_avg$upp_protection[i]=1}
  if(result_national_season_avg$low_protection[i]<0){result_national_season_avg$low_protection[i]=0}
  if(result_national_season_avg$mean_protection[i]>1){result_national_season_avg$mean_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# merging results from the 3 methods in a single dataframe ----
################################################################################
result_site_merge <- result_site[,-c(6:8)]
colnames(result_site_merge)[5] <- "proba"
result_national_merge <- result_national[,-c(5:6)]
result_national_season_avg_merge <- result_national_season_avg[,-c(6:20)]
colnames(result_national_season_avg_merge)[4] <- "mean_loss"
colnames(result_national_season_avg_merge)[5] <- "proba"
result_total_LRE <- rbind(result_national_merge,result_site_merge,result_national_season_avg_merge)
result_site_LRE <- result_site
result_national_LRE <- result_national
result_national_season_LRE <- result_national_season
result_national_season_avg_LRE <- result_national_season_avg
# plots ----
################################################################################
LRE_inclusion <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="inclusion"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(C) LRE - Inclusion of site-specific data to training")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

LRE_none <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(D) LRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"),legend.title = element_text(size = 9),legend.text = element_text(size = 8))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

################################################################################
# MRE guild ----
################################################################################
data2 <- data[which(data$sumOcc_MRE_predict > 0 & data$guild == "MRE" & data$period == "all" & data$bridage_presence == "non"),]
data3 <- data[which(data$sumOcc_MRE_predict > 0 & data$guild == "MRE" & data$period %in% c("spring","summer","autumn") & data$bridage_presence == "non"),]
nsites <- length(unique(data2$site))
data2 <- data2[-which(data2$proba==0 & data2$protection<1),]
data3 <- data3[-which(data3$proba==0 & data3$protection<1),]
# National threshold method ----
################################################################################
proba_boot <- c(unique(data2$proba[which(data2$proba <= 0.025)]),
                unique(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)][seq(1, length(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)]), 2)]),
                unique(data2$proba[which(data2$proba > 0.05 & data2$proba <=0.10)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 4)]),
                unique(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 8)]),
                unique(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)][seq(1, length(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)]), 16)]),
                unique(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)][seq(1, length(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)]), 32)]),
                unique(data2$proba[which(data2$proba > 0.40)][seq(1, length(data2$proba[which(data2$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]
temp_protection <- c()
temp_loss <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in proba_boot) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$proba==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    temp_loss <- as.data.frame(rbind(temp_loss,smean.cl.boot(sub$rotation_loss, conf.int = .95, B = 1000, na.rm = TRUE)))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="MRE")))
  }
}
result_national <- cbind(temp_protection,temp_loss,cov,method="national")
colnames(result_national)[1]<-"mean_protection"
colnames(result_national)[2]<-"low_protection"
colnames(result_national)[3]<-"upp_protection"
colnames(result_national)[4]<-"mean_loss"
colnames(result_national)[5]<-"low_loss"
colnames(result_national)[6]<-"upp_loss"

for(i in 1:nrow(result_national)){
  if(result_national$upp_protection[i]>1){result_national$upp_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# Site-specific threshold method  ----
################################################################################
data2$rotation_loss_round <- round(data2$rotation_loss,2)
temp_protection <- c()
mean_loss <- c()
proba_mean <- c()
proba_sd <- c()
proba_min <- c()
proba_max <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(data2$rotation_loss_round)) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$rotation_loss_round==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    mean_loss <- c(mean_loss,j)
    proba_mean <- c(proba_mean,mean(sub$proba))
    proba_sd <- c(proba_sd,sd(sub$proba))
    proba_min <- c(proba_min,min(sub$proba))
    proba_max <- c(proba_max,max(sub$proba))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="MRE")))
  }
}
result_site <- cbind(temp_protection,mean_loss,proba_mean,proba_sd,proba_min,proba_max,cov,method="site")
colnames(result_site)[1]<-"mean_protection"
colnames(result_site)[2]<-"low_protection"
colnames(result_site)[3]<-"upp_protection"
result_site <- result_site %>% drop_na()
result_site <- result_site[order(result_site$data_type,result_site$mean_loss, decreasing = T), ]
result_site$mean_protection[which(result_site$data_type=="none")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="none")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="none")] <- smoothing(result_site$low_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$mean_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$low_protection[which(result_site$data_type=="inclusion")], method = "smooth")

for(i in 1:nrow(result_site)){
  if(result_site$upp_protection[i]>1){result_site$upp_protection[i]=1}
}

# Season-specific threshold method  ----
################################################################################
proba_boot <- c(unique(data3$proba[which(data3$proba <= 0.025)]),
                unique(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)][seq(1, length(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)]), 2)]),
                unique(data3$proba[which(data3$proba > 0.05 & data3$proba <=0.10)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 4)]),
                unique(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 8)]),
                unique(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)][seq(1, length(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)]), 16)]),
                unique(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)][seq(1, length(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)]), 32)]),
                unique(data3$proba[which(data3$proba > 0.40)][seq(1, length(data3$proba[which(data3$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]

result_national_season <- c()
for(s in unique(data3$period)){
  print(s)
  temp_protection <- c()
  temp_loss <- c()
  cov <- c()
  for (i in c("inclusion","none")) {
    print(i)
    for (j in proba_boot) {
      
      sub = data3[which(data3$data_type==i & data3$proba==j & data3$period==s),]
      boot_protection <- c()
      boot_loss <- c()
      
      for(l in 1:1000){
        
        ind <- sample(1:nrow(sub), replace = TRUE)
        boot_protection <- c(boot_protection,sum(sub$nb_protection[ind], na.rm = T)/sum(sub$nb_passes[ind], na.rm = T))
        boot_loss <- c(boot_loss,sum(sub$nb_rotation_loss[ind], na.rm = T)/sum(sub$nb_rotation[ind], na.rm = T))
        
      }
      temp_protection <- rbind(temp_protection,cbind(mean(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)-1.96*sd(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)+1.96*sd(boot_protection, na.rm = T)))
      temp_loss <- rbind(temp_loss,cbind(mean(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)-1.96*sd(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)+1.96*sd(boot_loss, na.rm = T)))
      
      cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="MRE",method="national_season",season=s)))
    }
  }
  result_national_season <- rbind(result_national_season,cbind(temp_protection,temp_loss,cov))
}
colnames(result_national_season)[1]<-"mean_protection"
colnames(result_national_season)[2]<-"low_protection"
colnames(result_national_season)[3]<-"upp_protection"
colnames(result_national_season)[4]<-"mean_loss"
colnames(result_national_season)[5]<-"low_loss"
colnames(result_national_season)[6]<-"upp_loss"
result_national_season <- result_national_season %>% 
  filter_all(all_vars(!is.infinite(.)))

result_national_season$rotation_loss_round <- round(result_national_season$mean_loss,2)
mean_protection <- c()
low_protection <- c()
upp_protection <- c()
round_loss <- c()
proba_mean_spring <- c()
proba_mean_summer <- c()
proba_mean_autumn <- c()
proba_sd_spring <- c()
proba_sd_summer <- c()
proba_sd_autumn <- c()
proba_min_spring <- c()
proba_min_summer <- c()
proba_min_autumn <- c()
proba_max_spring <- c()
proba_max_summer <- c()
proba_max_autumn <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(result_national_season$rotation_loss_round)) {
    sub = result_national_season[which(result_national_season$data_type==i & result_national_season$rotation_loss_round==j),]
    mean_protection <- c(mean_protection,mean(sub$mean_protection))
    low_protection <- c(low_protection,mean(sub$low_protection))
    upp_protection <- c(upp_protection,mean(sub$upp_protection))
    round_loss <- c(round_loss,j)
    proba_mean_spring <- c(proba_mean_spring,mean(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_mean_summer <- c(proba_mean_summer,mean(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_mean_autumn <- c(proba_mean_autumn,mean(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_sd_spring <- c(proba_sd_spring,sd(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_sd_summer <- c(proba_sd_summer,sd(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_sd_autumn <- c(proba_sd_autumn,sd(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_min_spring <- c(proba_min_spring,min(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_min_summer <- c(proba_min_summer,min(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_min_autumn <- c(proba_min_autumn,min(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_max_spring <- c(proba_max_spring,max(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_max_summer <- c(proba_max_summer,max(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_max_autumn <- c(proba_max_autumn,max(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="MRE")))
  }
}
result_national_season_avg <- cbind(mean_protection,low_protection,upp_protection,round_loss,proba_mean_spring,proba_mean_summer,proba_mean_autumn,
                                    proba_sd_spring,proba_sd_summer,proba_sd_autumn,proba_min_spring,proba_min_summer,proba_min_autumn,
                                    proba_max_spring,proba_max_summer,proba_max_autumn,cov,method="national_season")

result_national_season_avg <- result_national_season_avg %>% drop_na(mean_protection)
result_national_season_avg <- result_national_season_avg[order(result_national_season_avg$data_type,result_national_season_avg$round_loss, decreasing = T), ]
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")

for(i in 1:nrow(result_national_season_avg)){
  if(result_national_season_avg$upp_protection[i]>1){result_national_season_avg$upp_protection[i]=1}
  if(result_national_season_avg$low_protection[i]<0){result_national_season_avg$low_protection[i]=0}
  if(result_national_season_avg$mean_protection[i]>1){result_national_season_avg$mean_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# merging results from the 3 methods in a single dataframe  ----
################################################################################
result_site_merge <- result_site[,-c(6:8)]
colnames(result_site_merge)[5] <- "proba"
result_national_merge <- result_national[,-c(5:6)]
result_national_season_avg_merge <- result_national_season_avg[,-c(6:16)]
colnames(result_national_season_avg_merge)[4] <- "mean_loss"
colnames(result_national_season_avg_merge)[5] <- "proba"
result_total_MRE <- rbind(result_national_merge,result_site_merge,result_national_season_avg_merge)
result_site_MRE <- result_site
result_national_MRE <- result_national
result_national_season_MRE <- result_national_season
result_national_season_avg_MRE <- result_national_season_avg
# plots ----
################################################################################
MRE_inclusion <- ggplot(data = result_total_MRE[which(result_total_MRE$data_type=="inclusion"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(A) MRE - Inclusion of site-specific data to training")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

MRE_none <- ggplot(data = result_total_MRE[which(result_total_MRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(B) MRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

# Saving figures and results ----
################################################################################
legend <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost blade rotations",
    y = "Proportion of protected bat activity")+
  ggtitle("(D) LRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "right",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

df_leg <- get_legend(legend)
tiff("./Efficiency_rotations_nocurt_bootstrap1000_smooth.tiff", width = 1300, height = 1400, units = "px",res = 150)
grid.arrange(arrangeGrob(MRE_inclusion,MRE_none,LRE_inclusion,LRE_none, nrow = 2),
             arrangeGrob(nullGrob(), df_leg, nullGrob(), nrow = 1), ncol = 1, heights = c(4,1),
             top = textGrob("", gp = gpar(fotsize = 12, font = 2)))
dev.off()

write_csv(rbind(result_national_LRE,result_national_MRE),"./NOCURT_result_national.csv")
write_csv(rbind(result_national_season_avg_LRE,cbind(result_national_season_avg_MRE,proba_mean_winter=NA,proba_sd_winter=NA,proba_min_winter=NA,proba_max_winter=NA)),"./NOCURT_result_national_season_avg.csv")
write_csv(rbind(result_national_season_LRE,result_national_season_MRE),"./NOCURT_result_national_season.csv")
write_csv(rbind(result_site_LRE,result_site_MRE),"./NOCURT_result_site.csv")


# *** Curtailment efficiency based on energy losses at all sites ----
################################################################################
# LRE guild ----
################################################################################
data2 <- data[which(data$sumOcc_LRE_predict > 0 & data$guild == "LRE" & data$period == "all" & data$energy_loss != "NA"),]
data3 <- data[which(data$sumOcc_LRE_predict > 0 & data$guild == "LRE" & data$period != "all" & data$energy_loss != "NA"),]
nsites <- length(unique(data2$site[which(data2$energy_loss != "NA")]))
data2 <- data2[-which(data2$proba==0 & data2$protection<1),]
data3 <- data3[-which(data3$proba==0 & data3$protection<1),]
# National threshold method ----
################################################################################
proba_boot <- c(unique(data2$proba[which(data2$proba <= 0.025)]),
                unique(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)][seq(1, length(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)]), 2)]),
                unique(data2$proba[which(data2$proba > 0.05 & data2$proba <=0.10)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 4)]),
                unique(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 8)]),
                unique(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)][seq(1, length(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)]), 16)]),
                unique(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)][seq(1, length(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)]), 32)]),
                unique(data2$proba[which(data2$proba > 0.40)][seq(1, length(data2$proba[which(data2$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]
temp_protection <- c()
temp_loss <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in proba_boot) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$proba==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    temp_loss <- as.data.frame(rbind(temp_loss,smean.cl.boot(sub$energy_loss, conf.int = .95, B = 1000, na.rm = TRUE)))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="LRE")))
  }
}
result_national <- cbind(temp_protection,temp_loss,cov,method="national")
colnames(result_national)[1]<-"mean_protection"
colnames(result_national)[2]<-"low_protection"
colnames(result_national)[3]<-"upp_protection"
colnames(result_national)[4]<-"mean_loss"
colnames(result_national)[5]<-"low_loss"
colnames(result_national)[6]<-"upp_loss"

for(i in 1:nrow(result_national)){
  if(result_national$upp_protection[i]>1){result_national$upp_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# Site-specific threshold method ----
################################################################################
data2$rotation_loss_round <- round(data2$energy_loss,2)
temp_protection <- c()
mean_loss <- c()
proba_mean <- c()
proba_sd <- c()
proba_min <- c()
proba_max <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(data2$rotation_loss_round)) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$rotation_loss_round==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    mean_loss <- c(mean_loss,j)
    proba_mean <- c(proba_mean,mean(sub$proba))
    proba_sd <- c(proba_sd,sd(sub$proba))
    proba_min <- c(proba_min,min(sub$proba))
    proba_max <- c(proba_max,max(sub$proba))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="LRE")))
  }
}
result_site <- cbind(temp_protection,mean_loss,proba_mean,proba_sd,proba_min,proba_max,cov,method="site")
colnames(result_site)[1]<-"mean_protection"
colnames(result_site)[2]<-"low_protection"
colnames(result_site)[3]<-"upp_protection"
result_site <- result_site %>% drop_na()
result_site <- result_site[order(result_site$data_type,result_site$mean_loss, decreasing = T), ]
result_site$mean_protection[which(result_site$data_type=="none")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="none")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="none")] <- smoothing(result_site$low_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$mean_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$low_protection[which(result_site$data_type=="inclusion")], method = "smooth")

for(i in 1:nrow(result_site)){
  if(result_site$upp_protection[i]>1){result_site$upp_protection[i]=1}
}

# Season-specific threshold method  ----
################################################################################
proba_boot <- c(unique(data3$proba[which(data3$proba <= 0.025)]),
                unique(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)][seq(1, length(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)]), 2)]),
                unique(data3$proba[which(data3$proba > 0.05 & data3$proba <=0.10)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 4)]),
                unique(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 8)]),
                unique(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)][seq(1, length(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)]), 16)]),
                unique(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)][seq(1, length(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)]), 32)]),
                unique(data3$proba[which(data3$proba > 0.40)][seq(1, length(data3$proba[which(data3$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]

result_national_season <- c()
for(s in unique(data3$period)){
  print(s)
  temp_protection <- c()
  temp_loss <- c()
  cov <- c()
  for (i in c("inclusion","none")) {
    print(i)
    for (j in proba_boot) {
      
      sub = data3[which(data3$data_type==i & data3$proba==j & data3$period==s),]
      boot_protection <- c()
      boot_loss <- c()
      
      for(l in 1:1000){
        
        ind <- sample(1:nrow(sub), replace = TRUE)
        boot_protection <- c(boot_protection,sum(sub$nb_protection[ind], na.rm = T)/sum(sub$nb_passes[ind], na.rm = T))
        boot_loss <- c(boot_loss,sum(sub$nb_energy_loss[ind], na.rm = T)/sum(sub$nb_energy[ind], na.rm = T))
        
      }
      temp_protection <- rbind(temp_protection,cbind(mean(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)-1.96*sd(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)+1.96*sd(boot_protection, na.rm = T)))
      temp_loss <- rbind(temp_loss,cbind(mean(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)-1.96*sd(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)+1.96*sd(boot_loss, na.rm = T)))
      
      cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="LRE",method="national_season",season=s)))
    }
  }
  result_national_season <- rbind(result_national_season,cbind(temp_protection,temp_loss,cov))
}
colnames(result_national_season)[1]<-"mean_protection"
colnames(result_national_season)[2]<-"low_protection"
colnames(result_national_season)[3]<-"upp_protection"
colnames(result_national_season)[4]<-"mean_loss"
colnames(result_national_season)[5]<-"low_loss"
colnames(result_national_season)[6]<-"upp_loss"
result_national_season <- result_national_season %>% 
  filter_all(all_vars(!is.infinite(.)))

result_national_season$rotation_loss_round <- round(result_national_season$mean_loss,2)
mean_protection <- c()
low_protection <- c()
upp_protection <- c()
round_loss <- c()
proba_mean_winter <- c()
proba_mean_spring <- c()
proba_mean_summer <- c()
proba_mean_autumn <- c()
proba_sd_winter <- c()
proba_sd_spring <- c()
proba_sd_summer <- c()
proba_sd_autumn <- c()
proba_min_winter <- c()
proba_min_spring <- c()
proba_min_summer <- c()
proba_min_autumn <- c()
proba_max_winter <- c()
proba_max_spring <- c()
proba_max_summer <- c()
proba_max_autumn <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(result_national_season$rotation_loss_round)) {
    sub = result_national_season[which(result_national_season$data_type==i & result_national_season$rotation_loss_round==j),]
    mean_protection <- c(mean_protection,mean(sub$mean_protection))
    low_protection <- c(low_protection,mean(sub$low_protection))
    upp_protection <- c(upp_protection,mean(sub$upp_protection))
    round_loss <- c(round_loss,j)
    proba_mean_winter <- c(proba_mean_winter,mean(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_mean_spring <- c(proba_mean_spring,mean(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_mean_summer <- c(proba_mean_summer,mean(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_mean_autumn <- c(proba_mean_autumn,mean(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_sd_winter <- c(proba_sd_winter,sd(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_sd_spring <- c(proba_sd_spring,sd(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_sd_summer <- c(proba_sd_summer,sd(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_sd_autumn <- c(proba_sd_autumn,sd(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_min_winter <- c(proba_min_winter,min(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_min_spring <- c(proba_min_spring,min(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_min_summer <- c(proba_min_summer,min(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_min_autumn <- c(proba_min_autumn,min(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_max_winter <- c(proba_max_winter,max(as.numeric(sub$proba[which(sub$season=="winter")]), na.rm=T))
    proba_max_spring <- c(proba_max_spring,max(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_max_summer <- c(proba_max_summer,max(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_max_autumn <- c(proba_max_autumn,max(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="LRE")))
  }
}
result_national_season_avg <- cbind(mean_protection,low_protection,upp_protection,round_loss,proba_mean_winter,proba_mean_spring,proba_mean_summer,proba_mean_autumn,
                                    proba_sd_winter,proba_sd_spring,proba_sd_summer,proba_sd_autumn,proba_min_winter,proba_min_spring,proba_min_summer,proba_min_autumn,
                                    proba_max_winter,proba_max_spring,proba_max_summer,proba_max_autumn,cov,method="national_season")

result_national_season_avg <- result_national_season_avg %>% drop_na(mean_protection)
result_national_season_avg <- result_national_season_avg[order(result_national_season_avg$data_type,result_national_season_avg$round_loss, decreasing = T), ]
# result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
# result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
# result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
# result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
# result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
# result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")], method = "loess")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")], method = "loess")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")], method = "loess")
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")], method = "loess")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")], method = "loess")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")], method = "loess")

for(i in 1:nrow(result_national_season_avg)){
  if(result_national_season_avg$upp_protection[i]>1){result_national_season_avg$upp_protection[i]=1}
  if(result_national_season_avg$low_protection[i]<0){result_national_season_avg$low_protection[i]=0}
  if(result_national_season_avg$mean_protection[i]>1){result_national_season_avg$mean_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# merging results from the 3 methods in a single dataframe ----
################################################################################
result_site_merge <- result_site[,-c(6:8)]
colnames(result_site_merge)[5] <- "proba"
result_national_merge <- result_national[,-c(5:6)]
result_national_season_avg_merge <- result_national_season_avg[,-c(6:20)]
colnames(result_national_season_avg_merge)[4] <- "mean_loss"
colnames(result_national_season_avg_merge)[5] <- "proba"
result_total_LRE <- rbind(result_national_merge,result_site_merge,result_national_season_avg_merge)
result_site_LRE <- result_site
result_national_LRE <- result_national
result_national_season_LRE <- result_national_season
result_national_season_avg_LRE <- result_national_season_avg
# plots ----
################################################################################
LRE_inclusion <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="inclusion"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost energy",
    y = "Proportion of protected bat activity")+
  ggtitle("(C) LRE - Inclusion of site-specific data to training")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

LRE_none <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost energy",
    y = "Proportion of protected bat activity")+
  ggtitle("(D) LRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"),legend.title = element_text(size = 9),legend.text = element_text(size = 8))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

################################################################################
# MRE guild ----
################################################################################
data2 <- data[which(data$sumOcc_MRE_predict > 0 & data$guild == "MRE" & data$period == "all" & data$energy_loss != "NA"),]
data3 <- data[which(data$sumOcc_MRE_predict > 0 & data$guild == "MRE" & data$period %in% c("spring","summer","autumn") & data$energy_loss != "NA"),]
nsites <- length(unique(data2$site[which(data2$energy_loss != "NA")]))
data2 <- data2[-which(data2$proba==0 & data2$protection<1),]
data3 <- data3[-which(data3$proba==0 & data3$protection<1),]
# National threshold method ----
################################################################################
proba_boot <- c(unique(data2$proba[which(data2$proba <= 0.025)]),
                unique(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)][seq(1, length(data2$proba[which(data2$proba > 0.025 & data2$proba <=0.05)]), 2)]),
                unique(data2$proba[which(data2$proba > 0.05 & data2$proba <=0.10)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 4)]),
                unique(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)][seq(1, length(data2$proba[which(data2$proba > 0.10 & data2$proba <=0.20)]), 8)]),
                unique(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)][seq(1, length(data2$proba[which(data2$proba > 0.20 & data2$proba <=0.30)]), 16)]),
                unique(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)][seq(1, length(data2$proba[which(data2$proba > 0.30 & data2$proba <=0.40)]), 32)]),
                unique(data2$proba[which(data2$proba > 0.40)][seq(1, length(data2$proba[which(data2$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]
temp_protection <- c()
temp_loss <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in proba_boot) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$proba==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    temp_loss <- as.data.frame(rbind(temp_loss,smean.cl.boot(sub$energy_loss, conf.int = .95, B = 1000, na.rm = TRUE)))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="MRE")))
  }
}
result_national <- cbind(temp_protection,temp_loss,cov,method="national")
colnames(result_national)[1]<-"mean_protection"
colnames(result_national)[2]<-"low_protection"
colnames(result_national)[3]<-"upp_protection"
colnames(result_national)[4]<-"mean_loss"
colnames(result_national)[5]<-"low_loss"
colnames(result_national)[6]<-"upp_loss"

for(i in 1:nrow(result_national)){
  if(result_national$upp_protection[i]>1){result_national$upp_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# Site-specific threshold method ----
################################################################################
data2$rotation_loss_round <- round(data2$energy_loss,2)
temp_protection <- c()
mean_loss <- c()
proba_mean <- c()
proba_sd <- c()
proba_min <- c()
proba_max <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(data2$rotation_loss_round)) {
    print(j)
    sub = data2[which(data2$data_type==i & data2$rotation_loss_round==j),]
    temp_protection <- as.data.frame(rbind(temp_protection,smean.cl.boot(sub$protection, conf.int = .95, B = 1000, na.rm = TRUE)))
    mean_loss <- c(mean_loss,j)
    proba_mean <- c(proba_mean,mean(sub$proba))
    proba_sd <- c(proba_sd,sd(sub$proba))
    proba_min <- c(proba_min,min(sub$proba))
    proba_max <- c(proba_max,max(sub$proba))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="MRE")))
  }
}
result_site <- cbind(temp_protection,mean_loss,proba_mean,proba_sd,proba_min,proba_max,cov,method="site")
colnames(result_site)[1]<-"mean_protection"
colnames(result_site)[2]<-"low_protection"
colnames(result_site)[3]<-"upp_protection"
result_site <- result_site %>% drop_na()
result_site <- result_site[order(result_site$data_type,result_site$mean_loss, decreasing = T), ]
result_site$mean_protection[which(result_site$data_type=="none")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="none")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="none")] <- smoothing(result_site$low_protection[which(result_site$data_type=="none")], method = "smooth")
result_site$mean_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$mean_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$upp_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$upp_protection[which(result_site$data_type=="inclusion")], method = "smooth")
result_site$low_protection[which(result_site$data_type=="inclusion")] <- smoothing(result_site$low_protection[which(result_site$data_type=="inclusion")], method = "smooth")

for(i in 1:nrow(result_site)){
  if(result_site$upp_protection[i]>1){result_site$upp_protection[i]=1}
}

# Season-specific threshold method ----
################################################################################
proba_boot <- c(unique(data3$proba[which(data3$proba <= 0.025)]),
                unique(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)][seq(1, length(data3$proba[which(data3$proba > 0.025 & data3$proba <=0.05)]), 2)]),
                unique(data3$proba[which(data3$proba > 0.05 & data3$proba <=0.10)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 4)]),
                unique(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)][seq(1, length(data3$proba[which(data3$proba > 0.10 & data3$proba <=0.20)]), 8)]),
                unique(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)][seq(1, length(data3$proba[which(data3$proba > 0.20 & data3$proba <=0.30)]), 16)]),
                unique(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)][seq(1, length(data3$proba[which(data3$proba > 0.30 & data3$proba <=0.40)]), 32)]),
                unique(data3$proba[which(data3$proba > 0.40)][seq(1, length(data3$proba[which(data3$proba > 0.40)]), 64)]))
proba_boot <- proba_boot[!is.na(proba_boot)]

result_national_season <- c()
for(s in unique(data3$period)){
  print(s)
  temp_protection <- c()
  temp_loss <- c()
  cov <- c()
  for (i in c("inclusion","none")) {
    print(i)
    for (j in proba_boot) {
      
      sub = data3[which(data3$data_type==i & data3$proba==j & data3$period==s),]
      boot_protection <- c()
      boot_loss <- c()
      
      for(l in 1:1000){
        
        ind <- sample(1:nrow(sub), replace = TRUE)
        boot_protection <- c(boot_protection,sum(sub$nb_protection[ind], na.rm = T)/sum(sub$nb_passes[ind], na.rm = T))
        boot_loss <- c(boot_loss,sum(sub$nb_energy_loss[ind], na.rm = T)/sum(sub$nb_energy[ind], na.rm = T))
        
      }
      temp_protection <- rbind(temp_protection,cbind(mean(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)-1.96*sd(boot_protection, na.rm = T),mean(boot_protection, na.rm = T)+1.96*sd(boot_protection, na.rm = T)))
      temp_loss <- rbind(temp_loss,cbind(mean(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)-1.96*sd(boot_loss, na.rm = T),mean(boot_loss, na.rm = T)+1.96*sd(boot_loss, na.rm = T)))
      
      cov <- as.data.frame(rbind(cov,cbind(data_type=i,proba=j,guild="MRE",method="national_season",season=s)))
    }
  }
  result_national_season <- rbind(result_national_season,cbind(temp_protection,temp_loss,cov))
}
colnames(result_national_season)[1]<-"mean_protection"
colnames(result_national_season)[2]<-"low_protection"
colnames(result_national_season)[3]<-"upp_protection"
colnames(result_national_season)[4]<-"mean_loss"
colnames(result_national_season)[5]<-"low_loss"
colnames(result_national_season)[6]<-"upp_loss"
result_national_season <- result_national_season %>% 
  filter_all(all_vars(!is.infinite(.)))

result_national_season$rotation_loss_round <- round(result_national_season$mean_loss,2)
mean_protection <- c()
low_protection <- c()
upp_protection <- c()
round_loss <- c()
proba_mean_spring <- c()
proba_mean_summer <- c()
proba_mean_autumn <- c()
proba_sd_spring <- c()
proba_sd_summer <- c()
proba_sd_autumn <- c()
proba_min_spring <- c()
proba_min_summer <- c()
proba_min_autumn <- c()
proba_max_spring <- c()
proba_max_summer <- c()
proba_max_autumn <- c()
cov <- c()
for (i in c("inclusion","none")) {
  print(i)
  for (j in unique(result_national_season$rotation_loss_round)) {
    sub = result_national_season[which(result_national_season$data_type==i & result_national_season$rotation_loss_round==j),]
    mean_protection <- c(mean_protection,mean(sub$mean_protection))
    low_protection <- c(low_protection,mean(sub$low_protection))
    upp_protection <- c(upp_protection,mean(sub$upp_protection))
    round_loss <- c(round_loss,j)
    proba_mean_spring <- c(proba_mean_spring,mean(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_mean_summer <- c(proba_mean_summer,mean(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_mean_autumn <- c(proba_mean_autumn,mean(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_sd_spring <- c(proba_sd_spring,sd(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_sd_summer <- c(proba_sd_summer,sd(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_sd_autumn <- c(proba_sd_autumn,sd(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_min_spring <- c(proba_min_spring,min(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_min_summer <- c(proba_min_summer,min(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_min_autumn <- c(proba_min_autumn,min(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    proba_max_spring <- c(proba_max_spring,max(as.numeric(sub$proba[which(sub$season=="spring")]), na.rm=T))
    proba_max_summer <- c(proba_max_summer,max(as.numeric(sub$proba[which(sub$season=="summer")]), na.rm=T))
    proba_max_autumn <- c(proba_max_autumn,max(as.numeric(sub$proba[which(sub$season=="autumn")]), na.rm=T))
    cov <- as.data.frame(rbind(cov,cbind(data_type=i,guild="MRE")))
  }
}
result_national_season_avg <- cbind(mean_protection,low_protection,upp_protection,round_loss,proba_mean_spring,proba_mean_summer,proba_mean_autumn,
                                    proba_sd_spring,proba_sd_summer,proba_sd_autumn,proba_min_spring,proba_min_summer,proba_min_autumn,
                                    proba_max_spring,proba_max_summer,proba_max_autumn,cov,method="national_season")

result_national_season_avg <- result_national_season_avg %>% drop_na(mean_protection)
result_national_season_avg <- result_national_season_avg[order(result_national_season_avg$data_type,result_national_season_avg$round_loss, decreasing = T), ]
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="none")], method = "smooth")
result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$mean_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$upp_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")
result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")] <- smoothing(result_national_season_avg$low_protection[which(result_national_season_avg$data_type=="inclusion")], method = "smooth")

for(i in 1:nrow(result_national_season_avg)){
  if(result_national_season_avg$upp_protection[i]>1){result_national_season_avg$upp_protection[i]=1}
  if(result_national_season_avg$low_protection[i]<0){result_national_season_avg$low_protection[i]=0}
  if(result_national_season_avg$mean_protection[i]>1){result_national_season_avg$mean_protection[i]=1}
}
rm(temp_loss,temp_protection,sub,cov)

# merging results from the 3 methods in a single dataframe ----
################################################################################
result_site_merge <- result_site[,-c(6:8)]
colnames(result_site_merge)[5] <- "proba"
result_national_merge <- result_national[,-c(5:6)]
result_national_season_avg_merge <- result_national_season_avg[,-c(6:16)]
colnames(result_national_season_avg_merge)[4] <- "mean_loss"
colnames(result_national_season_avg_merge)[5] <- "proba"
result_total_MRE <- rbind(result_national_merge,result_site_merge,result_national_season_avg_merge)
result_site_MRE <- result_site
result_national_MRE <- result_national
result_national_season_MRE <- result_national_season
result_national_season_avg_MRE <- result_national_season_avg

# plots ----
################################################################################
MRE_inclusion <- ggplot(data = result_total_MRE[which(result_total_MRE$data_type=="inclusion"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost energy",
    y = "Proportion of protected bat activity")+
  ggtitle("(A) MRE - Inclusion of site-specific data to training")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

MRE_none <- ggplot(data = result_total_MRE[which(result_total_MRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost energy",
    y = "Proportion of protected bat activity")+
  ggtitle("(B) MRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "none",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

# Saving figures and results  ----
################################################################################
legend <- ggplot(data = result_total_LRE[which(result_total_LRE$data_type=="none"),]) +
  geom_line(aes(x = mean_loss, y = mean_protection,colour = method, fill = method))+
  geom_ribbon(aes(x = mean_loss, ymin=low_protection,ymax=upp_protection, fill = method), alpha=0.2 ) +
  theme_bw()+
  scale_fill_discrete(labels = c("Fixed probability across sites and seasons",
                                 "Fixed probability across sites adapted by seasons",
                                 "Variable probability adapted to sites"))+
  scale_colour_discrete(labels = c("Fixed probability across sites and seasons",
                                   "Fixed probability across sites adapted by seasons",
                                   "Variable probability adapted to sites"))+
  labs(
    fill="Curtailment thresholding method",
    color="Curtailment thresholding method",
    x = "Proportion of lost energy",
    y = "Proportion of protected bat activity")+
  ggtitle("(D) LRE - No inclusion")+
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  theme(legend.position = "right",plot.title = element_text(size = 11, face = "bold"))+
  annotate(geom = "text", x = 0.9, y = 0.05, label = paste0("N=",nsites," sites"))

df_leg <- get_legend(legend)
tiff("./Efficiency_energy_bootstrap1000.tiff", width = 1300, height = 1400, units = "px",res = 150)
grid.arrange(arrangeGrob(MRE_inclusion,MRE_none,LRE_inclusion,LRE_none, nrow = 2),
             arrangeGrob(nullGrob(), df_leg, nullGrob(), nrow = 1), ncol = 1, heights = c(4,1),
             top = textGrob("", gp = gpar(fotsize = 12, font = 2)))
dev.off()

write_csv(rbind(result_national_LRE,result_national_MRE),"./results/ENERGY_result_national.csv")
write_csv(rbind(result_national_season_avg_LRE,cbind(result_national_season_avg_MRE,proba_mean_winter=NA,proba_sd_winter=NA,proba_min_winter=NA,proba_max_winter=NA)),"./ENERGY_result_national_season_avg.csv")
write_csv(rbind(result_national_season_LRE,result_national_season_MRE),"./ENERGY_result_national_season.csv")
write_csv(rbind(result_site_LRE,result_site_MRE),"./ENERGY_result_site.csv")

# *** Checking possible extrapolation of energy results to all sites ----
################################################################################

list <- list.files("./") # location of curtailement simulation files produced in the first script
data <- c()
for (i in list) { # aggregating files
  temp = fread(paste0("./",i))
  data <- rbind(data,temp)
}
data <- data %>% mutate_at(c(18:21), as.numeric)
colnames(data)[4]<-"data_type"
rm(temp)
# list of sites with energy data vs others
list_energy <- unique(data$site[which((data$sumOcc_LRE_predict > 0 | data$sumOcc_MRE_predict > 0) & data$energy_loss != "NA")])
list_others <- unique(data$site[which((data$sumOcc_LRE_predict > 0 | data$sumOcc_MRE_predict > 0) & data$site != list_energy)])
rm(data)
# aggregating prediction dataset for sites with energy data and others separately
data_energy <- c()
for (i in list_energy) {
  temp = fread(paste0("./predict_50_",i,".csv"))
  data_energy <- rbind(data_energy,temp)
}
data_other <- c()
for (i in list_others) {
  temp = fread(paste0("./predict_50_",i,".csv"))
  data_other <- rbind(data_other,temp)
}
# initializing vectors for three bat protection levels and associated required predicted probability for each guild
protection <- c(0.8,0.9,0.95,0.8,0.9,0.95,0.8,0.9,0.95,0.8,0.9,0.95)
seuil <- c(0.0028,0.0014,0.0007,0.0016,0.0009,0.0005,0.0028,0.0014,0.0007,0.0016,0.0009,0.0005)
guild <- c("MRE","MRE","MRE","LRE","LRE","LRE")
data <- c("energy","energy","energy","energy","energy","energy","others","others","others","others","others","others")
mean_speed_loss <- c(mean(data_energy$vitesse_10min[which(data_energy$MRE_pred_none >= 0.0028)], na.rm = T), # computing average blade speed for sites with energy data
                     mean(data_energy$vitesse_10min[which(data_energy$MRE_pred_none >= 0.0014)], na.rm = T),
                     mean(data_energy$vitesse_10min[which(data_energy$MRE_pred_none >= 0.0007)], na.rm = T),
                     mean(data_energy$vitesse_10min[which(data_energy$LRE_pred_none >= 0.0016)], na.rm = T),
                     mean(data_energy$vitesse_10min[which(data_energy$LRE_pred_none >= 0.0009)], na.rm = T),
                     mean(data_energy$vitesse_10min[which(data_energy$LRE_pred_none >= 0.0005)], na.rm = T),
                     mean(data_other$vitesse_10min[which(data_other$MRE_pred_none >= 0.0028)], na.rm = T), # computing average blade speed for sites without energy data
                     mean(data_other$vitesse_10min[which(data_other$MRE_pred_none >= 0.0014)], na.rm = T),
                     mean(data_other$vitesse_10min[which(data_other$MRE_pred_none >= 0.0007)], na.rm = T),
                     mean(data_other$vitesse_10min[which(data_other$LRE_pred_none >= 0.0016)], na.rm = T),
                     mean(data_other$vitesse_10min[which(data_other$LRE_pred_none >= 0.0009)], na.rm = T),
                     mean(data_other$vitesse_10min[which(data_other$LRE_pred_none >= 0.0005)], na.rm = T))
sd_speed_loss <- c(sd(data_energy$vitesse_10min[which(data_energy$MRE_pred_none >= 0.0028)], na.rm = T), # computing standard deviation of blade speed for sites with energy data
                     sd(data_energy$vitesse_10min[which(data_energy$MRE_pred_none >= 0.0014)], na.rm = T),
                     sd(data_energy$vitesse_10min[which(data_energy$MRE_pred_none >= 0.0007)], na.rm = T),
                     sd(data_energy$vitesse_10min[which(data_energy$LRE_pred_none >= 0.0016)], na.rm = T),
                     sd(data_energy$vitesse_10min[which(data_energy$LRE_pred_none >= 0.0009)], na.rm = T),
                     sd(data_energy$vitesse_10min[which(data_energy$LRE_pred_none >= 0.0005)], na.rm = T),
                    sd(data_other$vitesse_10min[which(data_other$MRE_pred_none >= 0.0028)], na.rm = T), # computing standard deviation of blade speed for sites without energy data
                    sd(data_other$vitesse_10min[which(data_other$MRE_pred_none >= 0.0014)], na.rm = T),
                    sd(data_other$vitesse_10min[which(data_other$MRE_pred_none >= 0.0007)], na.rm = T),
                    sd(data_other$vitesse_10min[which(data_other$LRE_pred_none >= 0.0016)], na.rm = T),
                    sd(data_other$vitesse_10min[which(data_other$LRE_pred_none >= 0.0009)], na.rm = T),
                    sd(data_other$vitesse_10min[which(data_other$LRE_pred_none >= 0.0005)], na.rm = T))
compil <- as.data.frame(cbind(guild,data,protection,seuil,mean_speed_loss,sd_speed_loss))
compil <- compil %>% mutate_at(c(5:6), as.numeric)

# graphical representation of blade speed range for sites with and without energy data
blade_speed <- ggplot(compil, aes(protection, mean_speed_loss)) +
  geom_pointrange(
    aes(ymin = mean_speed_loss-sd_speed_loss, ymax = mean_speed_loss+sd_speed_loss, color = data),
    position = position_dodge(0.3)
  )+
  labs(color="Dataset",y = "Blade speed lost with curtailment (km/h)", x = "Proportion of protected bat activity") + 
  scale_color_manual(labels = c("Energy loss (n=13)",
                                "Rotations loss (n=93)"), values = c("#00AFBB", "#E7B800"))+
  facet_grid(. ~ guild)+
  ggtitle("(A)")
rm(blade_speed,compil,data_energy,data_other,temp)

# limiting energy losses between 0 and 1
data2 <- data[which(data$period == "all" & data$energy_loss != "NA"),]
for(i in 1:nrow(data2)){
  if(data2$energy_loss[i]>1){data2$energy_loss[i]=1}
  if(data2$energy_loss[i]<0){data2$energy_loss[i]=0}
}

# graphical representation of the relationship between the blade rotation loss and energy loss
energy_rotation <- ggplot(data2, aes(rotation_loss, energy_loss, colour = site)) + 
  geom_point(alpha = 0.15) +
  labs(color="Sites", y = "Proportion of lost energy", x = "Proportion of lost blade rotations") +
  geom_smooth(aes(group=1))+
  ggtitle("(B)")+
  facet_grid(. ~ guild)

# saving
tiff("./Energy_rotations_generalization.tiff", width = 1000, height = 1200, units = "px",res = 150)
grid.arrange(arrangeGrob(blade_speed,energy_rotation, nrow = 2))
dev.off()

# *** Key numbers extraction ----
################################################################################

# the function below allows to compute for each subset analysis (all sites, 
# non-curtailed only, or sites with energy data), and each protection level
# and threshold definition method (national, site_specific, season-specific),
# the predicted probability above which wind turbines have to be curtailed, 
# and related AUC values. The function saves all these information in a file.
extract_summary=function(folder,
                         protection_levels){
  
  for (i in c("ALL","NOCURT","ENERGY")) {
    
    guild <- c()
    method <- c()
    protection <- c()
    type <- c()
    proba_avg <- c()
    proba_min <- c()
    proba_max <- c()
    proba_avg_winter <- c()
    proba_min_winter <- c()
    proba_max_winter <- c()
    proba_avg_spring <- c()
    proba_min_spring <- c()
    proba_max_spring <- c()
    proba_avg_summer <- c()
    proba_min_summer <- c()
    proba_max_summer <- c()
    proba_avg_autumn <- c()
    proba_min_autumn <- c()
    proba_max_autumn <- c()
    loss_avg <- c()
    loss_min <- c()
    loss_max <- c()
    auc_avg <- c()
    auc_min <- c()
    auc_max <- c()
    
    for (s in protection_levels) {
      
      national <- fread(paste0(folder,i,"_result_national.csv"))
      national <- national[order(national$guild,national$data_type,national$mean_loss, decreasing = T), ]
      national_season <- fread(paste0(folder,i,"_result_national_season.csv"))
      national_season_avg <- fread(paste0(folder,i,"_result_national_season_avg.csv"))
      site <- fread(paste0(folder,i,"_result_site.csv"))
      national$diff <- NA
      national_season$diff <- NA
      national_season_avg$diff <- NA
      site$diff <- NA
      national$diff_low <- NA
      national_season$diff_low <- NA
      national_season_avg$diff_low <- NA
      site$diff_low <- NA
      national$diff_upp <- NA
      national_season$diff_upp <- NA
      national_season_avg$diff_upp <- NA
      site$diff_upp <- NA
      national$diff = national$mean_protection-s
      national_season$diff = national_season$mean_protection-s
      national_season_avg$diff = national_season_avg$mean_protection-s
      site$diff = site$mean_protection-s
      national$diff_low = national$low_protection-s
      national_season$diff_low = national_season$low_protection-s
      national_season_avg$diff_low = national_season_avg$low_protection-s
      site$diff_low = site$low_protection-s
      national$diff_upp = national$upp_protection-s
      national_season$diff_upp = national_season$upp_protection-s
      national_season_avg$diff_upp = national_season_avg$upp_protection-s
      site$diff_upp = site$upp_protection-s
      
      for (g in unique(national$guild)) {
        
        for (m in c("national","national_season","site")) {
          
          for (t in unique(national$data_type)) {
            
            guild <- c(guild, g)
            method <- c(method, m)
            protection <- c(protection, s)
            type <- c(type, t)
            
            national_sub = national[which(national$guild==g & national$data_type==t),]
            national_season_avg_sub= national_season_avg[which(national_season_avg$guild==g & national_season_avg$data_type==t),]
            site_sub = site[which(site$guild==g & site$data_type==t),]
            
            if(m=="national"){
              proba_avg <- c(proba_avg,max(national_sub$proba[which(abs(national_sub$diff)==min(abs(national_sub$diff)))]))
              proba_min <- c(proba_min,max(national_sub$proba[which(abs(national_sub$diff_low)==min(abs(national_sub$diff_low)))]))
              proba_max <- c(proba_max,max(national_sub$proba[which(abs(national_sub$diff_upp)==min(abs(national_sub$diff_upp)))]))
              proba_avg_winter <- c(proba_avg_winter,NA)
              proba_min_winter <- c(proba_min_winter,NA)
              proba_max_winter <- c(proba_max_winter,NA)
              proba_avg_spring <- c(proba_avg_spring,NA)
              proba_min_spring <- c(proba_min_spring,NA)
              proba_max_spring <- c(proba_max_spring,NA)
              proba_avg_summer <- c(proba_avg_summer,NA)
              proba_min_summer <- c(proba_min_summer,NA)
              proba_max_summer <- c(proba_max_summer,NA)
              proba_avg_autumn <- c(proba_avg_autumn,NA)
              proba_min_autumn <- c(proba_min_autumn,NA)
              proba_max_autumn <- c(proba_max_autumn,NA)
              loss_avg <- c(loss_avg,min(national_sub$mean_loss[which(abs(national_sub$diff)==min(abs(national_sub$diff)))]))
              loss_min <- c(loss_min,min(national_sub$mean_loss[which(abs(national_sub$diff_upp)==min(abs(national_sub$diff_upp)))]))
              loss_max <- c(loss_max,min(national_sub$mean_loss[which(abs(national_sub$diff_low)==min(abs(national_sub$diff_low)))]))
              auc_avg <- c(auc_avg,auc(national_sub$mean_loss/max(national_sub$mean_loss),national_sub$mean_protection))
              auc_min <- c(auc_min,auc(national_sub$low_loss/max(national_sub$low_loss),national_sub$low_protection))
              auc_max <- c(auc_max,auc(national_sub$upp_loss/max(national_sub$upp_loss),national_sub$upp_protection))
            }
            
            if(m=="national_season"){
              proba_avg <- c(proba_avg,NA)
              proba_min <- c(proba_min,NA)
              proba_max <- c(proba_max,NA)
              winter = national_season[which(national_season$guild==g & national_season$data_type==t & national_season$season=="winter"),]
              spring = national_season[which(national_season$guild==g & national_season$data_type==t & national_season$season=="spring"),]
              summer = national_season[which(national_season$guild==g & national_season$data_type==t & national_season$season=="summer"),]
              autumn = national_season[which(national_season$guild==g & national_season$data_type==t & national_season$season=="autumn"),]
              proba_avg_winter <- c(proba_avg_winter,max(winter$proba[which(abs(winter$diff)==min(abs(winter$diff)))]))
              proba_min_winter <- c(proba_min_winter,max(winter$proba[which(abs(winter$diff_low)==min(abs(winter$diff_low)))]))
              proba_max_winter <- c(proba_max_winter,max(winter$proba[which(abs(winter$diff_upp)==min(abs(winter$diff_upp)))]))
              proba_avg_spring <- c(proba_avg_spring,max(spring$proba[which(abs(spring$diff)==min(abs(spring$diff)))]))
              proba_min_spring <- c(proba_min_spring,max(spring$proba[which(abs(spring$diff_low)==min(abs(spring$diff_low)))]))
              proba_max_spring <- c(proba_max_spring,max(spring$proba[which(abs(spring$diff_upp)==min(abs(spring$diff_upp)))]))
              proba_avg_summer <- c(proba_avg_summer,max(summer$proba[which(abs(summer$diff)==min(abs(summer$diff)))]))
              proba_min_summer <- c(proba_min_summer,max(summer$proba[which(abs(summer$diff_low)==min(abs(summer$diff_low)))]))
              proba_max_summer <- c(proba_max_summer,max(summer$proba[which(abs(summer$diff_upp)==min(abs(summer$diff_upp)))]))
              proba_avg_autumn <- c(proba_avg_autumn,max(autumn$proba[which(abs(autumn$diff)==min(abs(autumn$diff)))]))
              proba_min_autumn <- c(proba_min_autumn,max(autumn$proba[which(abs(autumn$diff_low)==min(abs(autumn$diff_low)))]))
              proba_max_autumn <- c(proba_max_autumn,max(autumn$proba[which(abs(autumn$diff_upp)==min(abs(autumn$diff_upp)))]))
              loss_avg <- c(loss_avg,min(national_season_avg_sub$round_loss[which(abs(national_season_avg_sub$diff)==min(abs(national_season_avg_sub$diff)))]))
              loss_min <- c(loss_min,min(national_season_avg_sub$round_loss[which(abs(national_season_avg_sub$diff_upp)==min(abs(national_season_avg_sub$diff_upp)))]))
              loss_max <- c(loss_max,min(national_season_avg_sub$round_loss[which(abs(national_season_avg_sub$diff_low)==min(abs(national_season_avg_sub$diff_low)))]))
              auc_avg <- c(auc_avg,auc(national_season_avg_sub$round_loss/max(national_season_avg_sub$round_loss),national_season_avg_sub$mean_protection))
              auc_min <- c(auc_min,auc(national_season_avg_sub$round_loss/max(national_season_avg_sub$round_loss),national_season_avg_sub$low_protection))
              auc_max <- c(auc_max,auc(national_season_avg_sub$round_loss/max(national_season_avg_sub$round_loss),national_season_avg_sub$upp_protection))
            }
            
            if(m=="site"){
              proba_avg <- c(proba_avg,max(site_sub$proba_mean[which(abs(site_sub$diff)==min(abs(site_sub$diff)))]))
              proba_min <- c(proba_min,max(site_sub$proba_mean[which(abs(site_sub$diff_low)==min(abs(site_sub$diff_low)))]))
              proba_max <- c(proba_max,max(site_sub$proba_mean[which(abs(site_sub$diff_upp)==min(abs(site_sub$diff_upp)))]))
              proba_avg_winter <- c(proba_avg_winter,NA)
              proba_min_winter <- c(proba_min_winter,NA)
              proba_max_winter <- c(proba_max_winter,NA)
              proba_avg_spring <- c(proba_avg_spring,NA)
              proba_min_spring <- c(proba_min_spring,NA)
              proba_max_spring <- c(proba_max_spring,NA)
              proba_avg_summer <- c(proba_avg_summer,NA)
              proba_min_summer <- c(proba_min_summer,NA)
              proba_max_summer <- c(proba_max_summer,NA)
              proba_avg_autumn <- c(proba_avg_autumn,NA)
              proba_min_autumn <- c(proba_min_autumn,NA)
              proba_max_autumn <- c(proba_max_autumn,NA)
              loss_avg <- c(loss_avg,min(site_sub$mean_loss[which(abs(site_sub$diff)==min(abs(site_sub$diff)))]))
              loss_min <- c(loss_min,min(site_sub$mean_loss[which(abs(site_sub$diff_upp)==min(abs(site_sub$diff_upp)))]))
              loss_max <- c(loss_max,min(site_sub$mean_loss[which(abs(site_sub$diff_low)==min(abs(site_sub$diff_low)))]))
              auc_avg <- c(auc_avg,auc(site_sub$mean_loss/max(site_sub$mean_loss),site_sub$mean_protection))
              auc_min <- c(auc_min,auc(site_sub$mean_loss/max(site_sub$mean_loss),site_sub$low_protection))
              auc_max <- c(auc_max,auc(site_sub$mean_loss/max(site_sub$mean_loss),site_sub$upp_protection))
            }
          }
        }
      }
    }
    write.csv(cbind(guild,method,protection,type,proba_avg,proba_min,proba_max,
                    proba_avg_winter,proba_min_winter,proba_max_winter,
                    proba_avg_spring,proba_min_spring,proba_max_spring,
                    proba_avg_summer,proba_min_summer,proba_max_summer,
                    proba_avg_autumn,proba_min_autumn,proba_max_autumn,
                    loss_avg,loss_min,loss_max,auc_avg,auc_min,auc_max),paste0(folder,i,"_summary.csv"))
  }
}

extract_summary(folder = "./",
                protection_levels = c(0.8,0.9,0.95))



# *** Computing correction coefficients to extrapolate  energy losses to an ----
#     annual energy loss 
################################################################################
# computing the average number of nights monitored in the winter period
temp <- data[data$jday < 85 | data$jday>325,]
temp2 <-  temp %>%
  dplyr::select(Code_eolienne_meteo,jday) %>%
  group_by(Code_eolienne_meteo)  %>%
  summarise(n_jday = length(unique(jday)))
summary(temp2$n_jday)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    4.50   17.00   26.87   30.50  125.00 

# selecting data with energy production within and outside winter
data_energy_winter <- data[which(data$production != "NA" & (data$jday < 85 | data$jday>325)),]
data_energy_other <- data[which(data$production != "NA" & (data$jday > 85 & data$jday<325)),]

# computing the total energy production in non-monitored winter days for LRE and MRE
temp_w <-  data_energy_winter %>%
  dplyr::select(id_site_date, pourcentage_nuit, production,sunrise, sunset) %>%
  group_by(id_site_date,sunrise, sunset)  %>%
  summarise(sum = sum(production))
temp_w$night_length <- as.numeric(temp_w$sunrise - temp_w$sunset)
temp_w$percentage <- temp_w$night_length/24
temp_w$prod_corrected <- temp_w$sum/temp_w$percentage
mean(temp_w$prod_corrected)*(125-26.87)
# 16373344 MW for LRE
mean(temp_w$prod_corrected)*(125)
# 20856700 MW for MRE

# computing the total energy production during spring, summer, and autumn
temp_o <-  data_energy_other %>%
  dplyr::select(id_site_date, pourcentage_nuit, production,sunrise, sunset) %>%
  group_by(id_site_date,sunrise, sunset)  %>%
  summarise(sum = sum(production))
temp_o$night_length <- as.numeric(temp_o$sunrise - temp_o$sunset)
temp_o$percentage <- temp_o$night_length/24
temp_o$prod_corrected <- temp_o$sum/temp_o$percentage
mean(temp_o$prod_corrected)*240
# 24037002 MW

# computing the proportion of the annual energy production covered by the monitoring
24037002*100/(24037002+16373344)
# 59.4823 % for LRE (energy losses over the study period should be multiplied by this value to get an estimated)
24037002*100/(24037002+20856700)
# 53.54204 % for MRE (energy losses over the study period should be multiplied by this value to get an estimated)

