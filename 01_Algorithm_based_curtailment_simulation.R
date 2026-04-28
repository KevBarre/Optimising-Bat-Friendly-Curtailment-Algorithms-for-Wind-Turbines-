rm(list=ls())

# * R packages
library(data.table)
library(stringr)
library(tidyr)
library(dplyr)
library(utils)
library(glmmTMB)

# * model formula for MRE and LRE guilds ----
f_MRE_2_5 <- paste0("MRE_binom_t1 ~ scale(jday) + scale(I(jday^2)) + 
                                       scale(pourcentage_nuit) + scale(I(pourcentage_nuit^2)) + 
                                       scale(rpm_10min) + scale(I(rpm_10min^2)) + 
                                       scale(vent_10min) + scale(I(vent_10min^2)) + 
                                       scale(temperature_10min) + scale(I(temperature_10min^2)) + 
                                       scale(rr_10min) + scale(I(rr_10min^2)) + 
                             
                                       scale(jday) : scale(I(pourcentage_nuit^2)) + 
                                       scale(I(jday^2)) : scale(pourcentage_nuit) + 
                                       scale(I(jday^2)) : scale(I(pourcentage_nuit^2)) + 
                      
                                       scale(jday) : scale(temperature_10min) + 
                                       scale(jday) : scale(I(temperature_10min^2)) + 
                                       scale(I(jday^2)) : scale(temperature_10min) + 
                                       scale(I(jday^2)) : scale(I(temperature_10min^2)) + 
                                     
                                       scale(I(jday^2)) : scale(vent_10min) + 
                                       scale(I(jday^2)) : scale(I(vent_10min^2)) + 
                                       
                                       scale(pourcentage_nuit) : scale(temperature_10min) + 
                                       scale(pourcentage_nuit) : scale(I(temperature_10min^2)) + 
                                       scale(I(pourcentage_nuit^2)) : scale(temperature_10min) +
                                       
                                       (1|annee) +
                                       (1|Code_eolienne_meteo) +
                                       (1|materiel_tot)")

f_LRE_2_5 <- paste0("LRE_binom_t1 ~ scale(jday) + scale(I(jday^2)) + 
                                          scale(pourcentage_nuit) + scale(I(pourcentage_nuit^2)) + 
                                          scale(rpm_10min) + scale(I(rpm_10min^2)) + 
                                          scale(vent_10min) + scale(I(vent_10min^2)) + 
                                          scale(temperature_10min) + scale(I(temperature_10min^2)) + 
                                          scale(rr_10min) + scale(I(rr_10min^2)) + 
                             
                                          scale(jday) : scale(pourcentage_nuit) + 
                                          scale(jday) : scale(I(pourcentage_nuit^2)) + 
                                          scale(I(jday^2)) : scale(pourcentage_nuit) + 
                                          scale(I(jday^2)) : scale(I(pourcentage_nuit^2)) + 
                                          scale(temperature_10min) : scale(I(vent_10min^2)) + 
                                          scale(jday) : scale(vent_10min) + 
                                          scale(jday) : scale(temperature_10min) + 
                                          scale(I(jday^2)) : scale(temperature_10min) + 
                                          scale(I(jday^2)) : scale(I(temperature_10min^2)) + 
                                          scale(pourcentage_nuit) : scale(vent_10min) + 
                                          scale(I(pourcentage_nuit^2)) : scale(vent_10min) + 
                                          scale(pourcentage_nuit) : scale(temperature_10min) + 
                                          scale(I(pourcentage_nuit^2)) : scale(temperature_10min) + 
                              
                                          (1|annee) +
                                          (1|Code_eolienne_meteo) +
                                          (1|materiel_tot)")


# ** algorithm-based curtailment comparisons ----
setwd() # working directory
# data preparation
P2_5_ER90 <- fread("./data_2_5_t1.csv") # reading dataset 
assign("data", P2_5_ER90)
rm(P2_5_ER90)
data <- subset(data, !data$materiel_tot %in% c("Batcorder_-30","SM4_6","SM4_24","SM4_12","SM3_12")) # removing rare material combinations

# selection of sites with at least 10 occurrences of one the guilds
site <- c()
sumOcc_MRE <- c()
sumOcc_LRE <- c()
data$annee <- as.numeric(data$annee)
for(i in unique(data$Code_eolienne_meteo)){
  s = data[which(data$Code_eolienne_meteo == i),]
  site = c(site, i)
  sumOcc_MRE <- c(sumOcc_MRE, sum(s$MRE_binom))
  sumOcc_LRE <- c(sumOcc_LRE, sum(s$LRE_binom))
}
concat <- as.data.frame(cbind(site,sumOcc_MRE,sumOcc_LRE))
list_site <- concat$site[which((as.numeric(concat$sumOcc_LRE) >= 10 | as.numeric(concat$sumOcc_MRE) >=10))] 

# starting the simulation
start.time.all <- Sys.time()

for (i in list_site[1:93]) { # site by site
  
  # i = "62572_5"
  
  start.time <- Sys.time()
  
  print(paste0(i," - ",which(list_site[1:93] == i),"/93 - start : ", Sys.time()))
  
  # prediction subset on which simulate curtailment efficiency
  data_predict <- data[which(data$Code_eolienne_meteo == i & 
                             data$date_debut_UTC %in% sample(unique(data$date_debut_UTC[which(data$Code_eolienne_meteo == i)]), 
                                             round(length(unique(data$date_debut_UTC[which(data$Code_eolienne_meteo == i)]))/2,0))),] # random selection of 50% of dates on the targeted site for curtailment efficiency test   
  data_predict$rpm_10min_next <- NA
  data_predict <- data_predict[order(data_predict$date_heure_UTC),]
  for (m in 1:nrow(data_predict)) {
    if(m < nrow(data_predict) & data_predict$pourcentage_nuit[m+1]-2.5 == data_predict$pourcentage_nuit[m]){data_predict$rpm_10min_next[m] = data_predict$rpm_10min[m+1]}
  } # computing blade RPM of the next time window to only curtail on t+1 time windows if RPM > 1 (WT in operation)

  # training dataset
  data_train_none <- data[-which(data$Code_eolienne_meteo == i),]
  data_train_none <- data_train_none[-which(data_train_none$date_debut_UTC %in% unique(data_predict$date_debut_UTC)),] # all data without the targeted site neither dates included in the prediction subset
  data_train_inclusion <- data[-which(data$date_debut_UTC %in% unique(data_predict$date_debut_UTC)),] # all data without the prediction subset (so including the monitoring dates on the site not included in the prediction subset L101)    

  # training models without or without site-specific data in training data
  MRE_train_none <- glmmTMB(as.formula(f_MRE_2_5), family = binomial(link = "cloglog"), data=data_train_none)
  LRE_train_none <- glmmTMB(as.formula(f_LRE_2_5), family = binomial(link = "cloglog"), data=data_train_none)
  print(paste0("Training without inclusion done - ", format(round(difftime(Sys.time(),start.time),2))))
  MRE_train_inclusion <- glmmTMB(as.formula(f_MRE_2_5), family = binomial(link = "cloglog"), data=data_train_inclusion)
  LRE_train_inclusion <- glmmTMB(as.formula(f_LRE_2_5), family = binomial(link = "cloglog"), data=data_train_inclusion)
  print(paste0("Training with inclusion done - ", format(round(difftime(Sys.time(),start.time),2))))
  
  # prediction of trained models on prediction subset
  data_predict$annee <- NA # fixing year to NA (to not favor algorithm)
  data_predict_none <- data_predict
  data_predict_none$Code_eolienne_meteo <- NA # fixing site to NA (because of no inclusion of site-specific data in training)
  data_predict_none$materiel_tot <- NA # fixing material to NA (because of no inclusion of site-specific data in training)
  data_predict_inclusion <- data_predict
  data_predict$MRE_pred_none <- predict(MRE_train_none, data_predict_none, type = "response", allow.new.levels=TRUE)
  data_predict$LRE_pred_none <- predict(LRE_train_none, data_predict_none, type = "response", allow.new.levels=TRUE)
  data_predict$MRE_pred_inclusion <- predict(MRE_train_inclusion, data_predict_inclusion, type = "response", allow.new.levels=TRUE)
  data_predict$LRE_pred_inclusion <- predict(LRE_train_inclusion, data_predict_inclusion, type = "response", allow.new.levels=TRUE)
  
  write.csv(data_predict[,c(1:2,20:22,27,33,40:58,78:84,88:94)], paste0("./predict_",i,".csv"), row.names=F) # saving

  print(paste0("Prediction done - ",  format(round(difftime(Sys.time(),start.time),2))))
  
  # dataframe creation to compute curtailment efficiency
  thresholds <- as.data.frame( 
    cbind(site = i, # site
          bridage_presence = unique(data$bridage[which(data$Code_eolienne_meteo == i)]), # already under curtailment or not
          guild = c(rep("MRE",5001),rep("LRE",5001),rep("MRE",5001),rep("LRE",5001)), # guild
          data = c(rep("none",5001),rep("none",5001),rep("inclusion",5001),rep("inclusion",5001)), # training/prediction method
          proba = c(seq(0,0.5,0.0001),seq(0,0.5,0.0001),seq(0,0.5,0.0001),seq(0,0.5,0.0001)), # predicted probability from models
          protection = NA, # % of protected bat detection
          nb_protection = NA, # number of protected bat detection
          nb_passes = NA, # total number of bat detections
          rotation_loss = NA, # % of lost blade rotations under RPM > 1
          nb_rotation_loss = NA,  # number of lost blade rotations under RPM > 1
          nb_rotation = NA, # total number of blade rotations with under > 1
          energy_loss = NA, # % of lost energy production under RPM > 1
          nb_energy_loss = NA, # amount of energy production under RPM > 1
          nb_energy = NA, # total amount of energy production under RPM > 1
          sumOcc_LRE_predict = sum(data_predict$LRE_binom), # sum of LRE occurrences in the prediction dataset
          sumOcc_MRE_predict = sum(data_predict$MRE_binom),  # sum of MRE occurrences in the prediction dataset
          sum_LRE_predict = sum(data_predict$LRE),# sum of LRE detections in the prediction dataset
          sum_MRE_predict = sum(data_predict$MRE) # sum of MRE detections in the prediction dataset
          )
  )
  # Duplicating the dataframe for all periods and each season 
  thresholds <- rbind(cbind(thresholds, period = "all"),cbind(thresholds, period = "winter"),cbind(thresholds, period = "spring"),cbind(thresholds, period = "summer"),cbind(thresholds, period = "autumn"))
  thresholds <- thresholds %>% mutate_at(c(5,18:21), as.numeric)

  # curtailment simulation
  pb <- txtProgressBar(1, nrow(thresholds), style = 3) # progression bar
  for (s in 1:nrow(thresholds)) { # for each predicted probability threshold
    # when the predicted probability exceeds the threshold s, computing the number of bat detections preserved, the total number of bat detections,
    # the number/amount of rotations/energy lost, and the total number/amount of rotations/energy, that are contained in the dataset
    # lines +1 (because algorithms are trained to predicted the expected presence probability at t+1 based on t0 weather and time information, 
    # the curtailment being made at t+1 based on t0 conditions) with more than one blade rotation per minute.
    setTxtProgressBar(pb, s)
    # for all periods when no site-specific data was included in model training
    if(thresholds$period[s] == "all" & thresholds$data[s] == "none" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
    }
    if(thresholds$period[s] == "all" & thresholds$data[s] == "none" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
    }
    # for all periods when site-specific data was included in model training
    if(thresholds$period[s] == "all" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
    }
    if(thresholds$period[s] == "all" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA")+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA")], na.rm = TRUE)
    }
    # for winter when no site-specific data was included in model training
    if(thresholds$period[s] == "winter" & thresholds$data[s] == "none" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
    }
    if(thresholds$period[s] == "winter" & thresholds$data[s] == "none" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
    }
    # for winter when site-specific data was included in model training
    if(thresholds$period[s] == "winter" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
    }
    if(thresholds$period[s] == "winter" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday < 85 | data_predict$jday>325))], na.rm = TRUE)
    }
    # for spring when no site-specific data was included in model training
    if(thresholds$period[s] == "spring" & thresholds$data[s] == "none" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
    }
    if(thresholds$period[s] == "spring" & thresholds$data[s] == "none" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
    }
    # for spring when site-specific data was included in model training
    if(thresholds$period[s] == "spring" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
    }
    if(thresholds$period[s] == "spring" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 85 & data_predict$jday<179))], na.rm = TRUE)
    }
    # for summer when no site-specific data was included in model training
    if(thresholds$period[s] == "summer" & thresholds$data[s] == "none" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
    }
    if(thresholds$period[s] == "summer" & thresholds$data[s] == "none" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
    }
    # for summer when site-specific data was included in model training
    if(thresholds$period[s] == "summer" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
    }
    if(thresholds$period[s] == "summer" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 179 & data_predict$jday<265))], na.rm = TRUE)
    }
    # for autumn when no site-specific data was included in model training
    if(thresholds$period[s] == "autumn" & thresholds$data[s] == "none" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
    }
    if(thresholds$period[s] == "autumn" & thresholds$data[s] == "none" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_none >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
    }
    # for autumn when site-specific data was included in model training
    if(thresholds$period[s] == "autumn" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "MRE"){
      thresholds$nb_protection[s] = sum(data_predict$MRE[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$MRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$MRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
    }
    if(thresholds$period[s] == "autumn" & thresholds$data[s] == "inclusion" & thresholds$guild[s] == "LRE"){
      thresholds$nb_protection[s] = sum(data_predict$LRE[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE)
      thresholds$nb_passes[s] = sum(data_predict$LRE[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
      thresholds$nb_rotation_loss[s] = sum(data_predict$nb_rotations[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE) 
      thresholds$nb_rotation[s] = sum(data_predict$nb_rotations[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
      thresholds$nb_energy_loss[s] = sum(data_predict$production[which(data_predict$LRE_pred_inclusion >= thresholds$proba[s] & data_predict$rpm_10min_next > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))+1], na.rm = TRUE) 
      thresholds$nb_energy[s] = sum(data_predict$production[which(data_predict$rpm_10min > 1 & data_predict$rpm_10min_next != "NA" & (data_predict$jday > 265 & data_predict$jday<325))], na.rm = TRUE)
    }

  }
  
  # computing the final proportion of bat detections protected, rotations and energy lost for each probability threshold
  thresholds$protection = as.numeric(thresholds$nb_protection) / as.numeric(thresholds$nb_passes)
  thresholds$rotation_loss = as.numeric(thresholds$nb_rotation_loss) / as.numeric(thresholds$nb_rotation)
  thresholds$energy_loss = as.numeric(thresholds$nb_energy_loss) / as.numeric(thresholds$nb_energy)
  
  # saving
  write.csv(thresholds, paste0("./simul_",i,".csv"), row.names=F)
  
  print(paste0("Curtailment simulation done - ", format(round(difftime(Sys.time(),start.time),2))))
  
}

difftime(Sys.time(),start.time.all)


