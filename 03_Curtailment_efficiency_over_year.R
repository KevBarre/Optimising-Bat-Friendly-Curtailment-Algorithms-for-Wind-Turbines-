# R pacakges
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(modelbased)

list <- list.files("./") # location of predict files from the script 1

##########################################################################
# *** Graphical representation of algorithm efficiency over the year ----
##########################################################################
# vectors for computing losses and protection by 5 days intervals
cut <- seq(0,365,5)[-1]
cut_jday <- c()
site <- c()
prop_yprotection_LRE <- c()
prop_yloss_LRE <- c()
prop_yprotection_MRE <- c()
prop_yloss_MRE <- c()

for (i in list) { # site by site
  
  print(i)
  
  temp <- fread(paste0("/Users/kbarre/Documents/08_WinWind/analyses/Efficiency tests/V4_no spatial_proba_t_1/predict_inclusion/",i)) # reading prediction dataset of the site i (produced in the script 1)
  
  for (j in cut) { # for each julian day interval
    
    sub <- temp[which(temp$jday <= j & temp$jday > j-5),] # selecting the last 5 days only
    
    site <- c(site, unique(temp$Code_eolienne_meteo)) # storing site identifier
    cut_jday <- c(cut_jday, j) # storing the julian day upper limit
    
    # computing the proportion of annual bat protection and annual rotation loss on the 5-days interval 
    if(nrow(sub) > 0){
      prop_yprotection_LRE <- c(prop_yprotection_LRE, sum(sub$LRE[which(sub$LRE_pred_none>=0.0009 & sub$rpm_10min_next > 1 & sub$rpm_10min_next != "NA")+1], na.rm = TRUE)/sum(temp$LRE[which(sub$LRE_pred_none>=0.0009 & temp$rpm_10min > 1 & temp$rpm_10min_next != "NA")], na.rm = TRUE))
      prop_yloss_LRE <- c(prop_yloss_LRE, sum(sub$nb_rotations[which(sub$LRE_pred_none>=0.0009 & sub$rpm_10min_next > 1 & sub$rpm_10min_next != "NA")+1], na.rm = TRUE)/sum(temp$nb_rotations[which(sub$LRE_pred_none>=0.0009 & temp$rpm_10min > 1 & temp$rpm_10min_next != "NA")], na.rm = TRUE))
      prop_yprotection_MRE <- c(prop_yprotection_MRE, sum(sub$MRE[which(sub$MRE_pred_none>=0.0014 & sub$rpm_10min_next > 1 & sub$rpm_10min_next != "NA")+1], na.rm = TRUE)/sum(temp$MRE[which(sub$MRE_pred_none>=0.0014 & temp$rpm_10min > 1 & temp$rpm_10min_next != "NA")], na.rm = TRUE))
      prop_yloss_MRE <- c(prop_yloss_MRE, sum(sub$nb_rotations[which(sub$MRE_pred_none>=0.0014 & sub$rpm_10min_next > 1 & sub$rpm_10min_next != "NA")+1], na.rm = TRUE)/sum(temp$nb_rotations[which(sub$MRE_pred_none>=0.0014 & temp$rpm_10min > 1 & temp$rpm_10min_next != "NA")], na.rm = TRUE))
    }else{
      prop_yprotection_LRE <- c(prop_yprotection_LRE, NA)
      prop_yloss_LRE <- c(prop_yloss_LRE, NA)
      prop_yprotection_MRE <- c(prop_yprotection_MRE, NA)
      prop_yloss_MRE <- c(prop_yloss_MRE, NA)
    }
    if(prop_yprotection_LRE[length(prop_yprotection_LRE)] %in% "NaN"){prop_yprotection_LRE[length(prop_yprotection_LRE)] = 0}
    if(prop_yloss_LRE[length(prop_yloss_LRE)] %in% "NaN"){prop_yloss_LRE[length(prop_yloss_LRE)] = 0}
    if(prop_yprotection_MRE[length(prop_yprotection_MRE)] %in% "NaN"){prop_yprotection_MRE[length(prop_yprotection_MRE)] = 0}
    if(prop_yloss_MRE[length(prop_yloss_MRE)] %in% "NaN"){prop_yloss_MRE[length(prop_yloss_MRE)] = 0}
    
  }
  
}

# merging results
result_final <- as.data.frame(cbind(cut_jday,site,prop_protection_LRE,prop_loss_LRE,prop_protection_MRE,prop_loss_MRE,prop_yprotection_LRE,prop_yloss_LRE,prop_yprotection_MRE,prop_yloss_MRE))
result_final <- result_final %>% 
  mutate_at(c(1,3:6), as.numeric)

# data formatting
data_plot = result_final %>%
  pivot_longer(names_to = "type", cols=prop_protection_LRE:prop_yloss_MRE)
data_plot$value <- as.numeric(data_plot$value)

# plots
plot_LRE <- ggplot(data_plot[which(data_plot$type %in% c("prop_yprotection_LRE","prop_yloss_LRE")),], aes(cut_jday, value, colour=type)) +
  theme_bw()+
  stat_summary(fun.data = "mean_cl_normal", size=0.2)+
  scale_color_discrete(labels = c("Lost rotations",
                                  "Protected bat activity"),
                       name="")+
  labs(y = "Proportion of annual lost rotations and protected bat activity \nby 5-days intervals", x = "Julian day")+
  ggtitle("(A) Long-range echolocators")+
  theme(legend.position = c(0.20,0.85))+
  scale_x_continuous(breaks=c(50,100,150,200,250,300,350), minor_breaks = c(50,100,150,200,250,300,350), expand = c(0.05,0))

plot_MRE <- ggplot(data_plot[which(data_plot$type %in% c("prop_yprotection_MRE","prop_yloss_MRE") & data_plot$cut_jday != 165),], aes(as.factor(cut_jday), value, colour=type)) +
  theme_bw()+
  stat_summary(fun.data = "mean_cl_normal", size=0.2)+
  scale_color_discrete(labels = c("Lost rotations",
                                  "Protected bat activity"),
                       name="")+
  labs(y = "Proportion of annual lost rotations and protected bat activity \nby 5-days intervals", x = "Julian day")+
  ggtitle("(B) Mid-range echolocators")+
  theme(legend.position = c(0.20,0.85))+
  scale_x_discrete(breaks=c(0,50,100,150,200,250,300,350), expand = c(0.05,0))

# saving
tiff("/Users/kbarre/Documents/08_WinWind/MS/5days_prop_loss_prot_ofannual_NONE_NATIONAL_P90.tiff", width = 1800, height = 850, units = "px",res = 150)
cowplot::plot_grid(plot_LRE, plot_MRE,
                   nrow = 1,
                   align = "h")
dev.off()

#################################################################################
# *** % of annual energy losses that occur when bat protection is negligible ----
#################################################################################

# computing the average energy loss on each 5-days intervals
jday <- c()
loss_boot <- c()
type <- c()
for (i in c("prop_yloss_MRE","prop_yloss_LRE")) {
  sub <- data_plot[which(data_plot$type==i),]
  for (j in unique(data_plot$cut_jday)) {
    loss_boot <- rbind(loss_boot,smean.cl.boot(sub$value[which(sub$cut_jday==j)], conf.int = .95, B = 1000, na.rm = TRUE))
    jday <- c(jday,j)
    type <- c(type,i)
  }
}
sum_eff_loss <- as.data.frame(cbind(type, jday,loss_boot))
sum_eff_loss$jday <- as.numeric(sum_eff_loss$jday)
sum_eff_loss$Mean <- as.numeric(sum_eff_loss$Mean)

# Computing the % of annual energy losses that occur between mid-April and mid-November, when bat protection is zero or negligible
sum(sum_eff_loss$Mean[which(sum_eff_loss$jday<105 | sum_eff_loss$jday>305 & sum_eff_loss$type == "prop_yloss_MRE")], na.rm = T)
# 0.3200896
sum(sum_eff_loss$Mean[which(sum_eff_loss$jday<105 | sum_eff_loss$jday>305 & sum_eff_loss$type == "prop_yloss_LRE")], na.rm = T)
# 0.3529375

