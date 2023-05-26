#-------------------------------------------------------------------------------
# Title: ACTG 5248 analysis of CD8 T cell Memory Restoration using MSI data
# MSI = mean signal intensity
#
# Purpose: 
# Interested in looking at CD8 T cell memory restoration in people 
# living with HIV once they receive ART.
# 
# Question:
# How does MSI change as a function of time?
# 
# Two cohorts:
# People with HIV durably suppressed on ART (PWH-ART), referred to as DS
#    - 10 participants
#    - 3 DIFFERENT timepoints for each.
# ACTG5248
#   - 10 participants
#   - 12 SAME timepoints for each
#   - Days: 0, 2, 7, 10, 14, 21, 28, 56, 
#           84, 140, 252, 504
# Missing data:
# 291374 day 28
# 214273 day 2, 7, 14, 28, 84, 140
# 611175 day 2, 10, 56
# 611183 day 14, 21
# 611051 day 2, 10
# 21929 day 2, 7, 10
# 112185 day 2, 10
# 611172 day 2, 10
#
# Analysis: Gamma fixed-effects model with log link with indicators for PID.
# Statistician: Ann Marie Weideman
# Date Created: 4 May 2023
# Last Modified: 26 May 2023
# Notes:
# We chose this analysis over a mixed-effects model because fixed-effects models
# tend to be less susceptible to bias when we have only a few clusters 
# (in this case, patients). We have evidence from literature to suggest that 
# this is an appropriate approach. One citation can be found here:
# https://biol609.github.io/Readings/McNeish_Kelley_PsychMethods_2019.pdf. 
#-------------------------------------------------------------------------------

#-------------------------------------
# Import and preprocess data
#-------------------------------------
curr_dir<-'G:/projects/cfar/Nilu/ACTG_5248/'

#Modify here as data is updated
ACTG <- read.csv(file = paste0(curr_dir,'Data/msi_data/CD8+ MSI A5248_V2.csv'))
#Modify here as data is updated
DS <- read.csv(file = paste0(curr_dir,'Data/msi_data/CD8+ MSI LT-ART_V3.csv'))

# Generate indicator variables
#Modify here as data is updated
ACTG$PID291374 <- ifelse(ACTG$PID == 291374, 1, 0)
ACTG$PID363043 <- ifelse(ACTG$PID == 363043, 1, 0)
ACTG$PID363044 <- ifelse(ACTG$PID == 363044, 1, 0)
ACTG$PID214273 <- ifelse(ACTG$PID == 214273, 1, 0)
ACTG$PID611175 <- ifelse(ACTG$PID == 611175, 1, 0)
ACTG$PID611183 <- ifelse(ACTG$PID == 611183, 1, 0)
ACTG$PID611051 <- ifelse(ACTG$PID == 611051, 1, 0)
ACTG$PID21929 <- ifelse(ACTG$PID == 21929, 1, 0)
ACTG$PID112185 <- ifelse(ACTG$PID == 112185, 1, 0)
ACTG$PID611172 <- ifelse(ACTG$PID == 611172, 1, 0)

#Modify here as data is updated
DS$PID546 <- ifelse(DS$PID == 546, 1, 0)
DS$PID861 <- ifelse(DS$PID  == 861, 1, 0)
DS$PID728 <-  ifelse(DS$PID  == 728, 1, 0)
DS$PID674 <- ifelse(DS$PID  == 674, 1, 0)
DS$PID1095 <- ifelse(DS$PID  == 1095, 1, 0)
DS$PID929 <- ifelse(DS$PID  == 929, 1, 0)
DS$PID749 <- ifelse(DS$PID == 749, 1, 0)
DS$PID673 <- ifelse(DS$PID == 673, 1, 0)
DS$PID870 <- ifelse(DS$PID  == 870, 1, 0)
DS$PID1036 <- ifelse(DS$PID  == 1036, 1, 0)

# Create months variable
ACTG$ThirtyDays <- ACTG$Day/30

#-----------------------------------------------------------
# ACTG cohort: Fixed effects models with indicators for PID
#-----------------------------------------------------------
require(broom)
require(dplyr)

ACTG$ThirtyDays <- ACTG$Day/30
by_population_ACTG <- group_by(ACTG, Marker)

#Modify here as data is updated
require(margins)
fems_ACTG<-do(by_population_ACTG, 
              tidy(glm(MSI~ -1 + ThirtyDays + PID291374 + PID363043
                       + PID363044 + PID214273 + PID611175 + PID611183 
                       + PID611051  + PID21929 + PID112185 + PID611172, 
                       family=Gamma(link = "log"), data = .))[1,])

# Estimate on the decimal and percent scale
fems_ACTG$estimate_decimal<-(exp(fems_ACTG$estimate)-1)
fems_ACTG$estimate_percent<-(exp(fems_ACTG$estimate)-1)*100

# CIs on the decimal and percent scale
fems_ACTG$CI_lower_decimal<-(exp(fems_ACTG$estimate-1.96*fems_ACTG$std.error)-1)
fems_ACTG$CI_upper_decimal<-(exp(fems_ACTG$estimate+1.96*fems_ACTG$std.error)-1)
fems_ACTG$CI_lower_percent<-fems_ACTG$CI_lower_decimal*100
fems_ACTG$CI_upper_percent<-fems_ACTG$CI_upper_decimal*100

# Add descriptions of results
fems_ACTG$description<- ifelse(fems_ACTG$p.value<0.05,
                               paste0("For every 30 days post ART initiation, there is an estimated ",
                                      ifelse(sign(fems_ACTG$estimate_percent)==-1,"decrease","increase"),
                                      " of ",round(abs(fems_ACTG$estimate_percent),2),"% in MSI (",
                                      ifelse(fems_ACTG$p.value<0.0001,"p<0.0001",
                                             ifelse(fems_ACTG$p.value<0.001,"p<0.001",
                                                    ifelse(fems_ACTG$p.value<0.01,"p<0.01",
                                                           ifelse(fems_ACTG$p.value<0.05,"p<0.05",
                                                                  paste0("p=",round(fems_ACTG$p.value,2)))))),
                                      
                                      "), on average, for ",
                                      fems_ACTG$Marker,"."),
                               paste0("There is no evidence of a change in ", fems_ACTG$Marker,
                                      " across time.")
)

fems_ACTG<-dplyr::rename(fems_ACTG, p_value = p.value)

# Format estimates and p-values
fems_ACTG$slope<-format(fems_ACTG$estimate_percent, scientific=FALSE)
fems_ACTG$p_value<-format(fems_ACTG$p_value, scientific=FALSE)

fems_ACTG$term<-"Every 30 days"

# Drop unexponentiated/untransformed estimates
fems_ACTG<-subset(fems_ACTG,select=-c(std.error))

#---------------------------------------------------------------
# Graph deviance residuals to determine if normal
#---------------------------------------------------------------

#Deviance residuals for ACTG
require(dplyr)
ACTG_glms<-ACTG %>% nest_by(Marker) %>%
  mutate(mod = list(glm(MSI~ -1 + ThirtyDays + PID291374 + PID363043
                        + PID363044 + PID214273 + PID611175 + PID611183 
                        + PID611051  + PID21929 + PID112185 + PID611172, 
                        data=data,family=Gamma(link = "log")))) 

pdf(file = paste0(curr_dir,"data/output/CD8_ACTG_MSI_residuals.pdf"))
for(i in 1:length(ACTG_glms$mod)){
  ypred = predict(ACTG_glms$mod[[i]])
  res = residuals(ACTG_glms$mod[[i]], type = 'deviance')
  par(mfrow = c(2, 1)) 
  plot(ypred,res, ylab="residuals",xlab="predicted")
  hist(res, xlab="Residuals",main=NULL)
  mtext(ACTG_glms$Marker[[i]],                   
        side = 3,
        line = - 2,
        cex=1.5,
        outer = TRUE)
}
dev.off()

#-----------------------------------------------------------
# DS cohort: Fixed effects models with indicators for PID
#-----------------------------------------------------------

by_population_DS <- group_by(DS, Marker)

#Modify here as data is updated
require(margins)
fems_DS<-do(by_population_DS, 
            tidy((glm(MSI ~ -1 + Months.post.timepoint.1 +
                        + PID546 + PID861 + PID728 + PID674 + PID1095 
                      + PID929 + PID749 + PID673 + PID870 + PID1036, 
                      family=Gamma(link = "log"), data = .)))[1,])

# Estimate on the decimal and percent scale
fems_DS$estimate_decimal<-(exp(fems_DS$estimate)-1)
fems_DS$estimate_percent<-(exp(fems_DS$estimate)-1)*100

# CIs on the decimal and percent scale
fems_DS$CI_lower_decimal<-(exp(fems_DS$estimate-1.96*fems_DS$std.error)-1)
fems_DS$CI_upper_decimal<-(exp(fems_DS$estimate+1.96*fems_DS$std.error)-1)
fems_DS$CI_lower_percent<-fems_DS$CI_lower_decimal*100
fems_DS$CI_upper_percent<-fems_DS$CI_upper_decimal*100

# Add descriptions of results
fems_DS$description<- ifelse(fems_DS$p.value<0.05,
                             paste0("For every 30 days while durably suppressed, there is an estimated ",
                                    ifelse(sign(fems_DS$estimate_percent)==-1,"decrease","increase"), 
                                    " of ",round(abs(fems_DS$estimate_percent),2), "% in MSI (",
                                    ifelse(fems_DS$p.value<0.0001,"p<0.0001",
                                           ifelse(fems_DS$p.value<0.001,"p<0.001",
                                                  ifelse(fems_DS$p.value<0.01,"p<0.01",
                                                         ifelse(fems_DS$p.value<0.05,"p<0.05",
                                                                paste0("p=",round(fems_DS$p.value,2)))))),
                                    
                                    "), on average, for ",
                                    fems_DS$Marker,"."),
                             paste0("There is no evidence of a change in ",fems_DS$Marker,
                                    " across time.")
)

fems_DS<-dplyr::rename(fems_DS, p_value = p.value)

# Format estimates and p-values
fems_DS$slope<-format(fems_DS$estimate_percent, scientific=FALSE)
fems_DS$p_value<-format(fems_DS$p_value, scientific=FALSE)

fems_DS$term<-"Every 30 days"

# Drop unexponentiated/untransformed estimates
fems_DS<-subset(fems_DS,select=-c(std.error))

# Reorder columns prior to export
fems_ACTG<-fems_ACTG[,c("Marker", "term", "estimate_decimal", 
                        "estimate_percent", "statistic", "p_value",
                        "CI_lower_decimal", "CI_upper_decimal",
                        "CI_lower_percent", "CI_upper_percent", "description")]
fems_DS<-fems_DS[,c("Marker", "term", "estimate_decimal", 
                    "estimate_percent", "statistic", "p_value",
                    "CI_lower_decimal", "CI_upper_decimal",
                    "CI_lower_percent", "CI_upper_percent", "description")]

#---------------------------------------------------------------
# Graph deviance residuals to determine if normal
#---------------------------------------------------------------

#Deviance residuals for DS
require(dplyr)
DS_glms<-DS %>% nest_by(Marker) %>%
  mutate(mod = list(glm(MSI ~ -1 + Months.post.timepoint.1 +
                          + PID546 + PID861 + PID728 + PID674 + PID1095 
                        + PID929 + PID749 + PID673 + PID870 + PID1036, 
                        data=data, family=Gamma(link = "log")))) 

pdf(file = paste0(curr_dir,"data/output/CD8_DS_MSI_residuals.pdf"))
for(i in 1:length(DS_glms$mod)){
  ypred = predict(DS_glms$mod[[i]])
  res = residuals(DS_glms$mod[[i]], type = 'deviance')
  par(mfrow = c(2, 1)) 
  plot(ypred,res, ylab="residuals",xlab="predicted")
  hist(res, xlab="Residuals",main=NULL)
  mtext(DS_glms$Marker[[i]],                   
        side = 3,
        line = - 2,
        cex=1.5,
        outer = TRUE)
}
dev.off()

#-------------------------------------
# Write to csv
#-------------------------------------
require(writexl)
# #Modify here as data is updated
# write_xlsx(fems_ACTG, paste0(curr_dir,"data/output/MSI_CD8_ACTG_analysis_AMW_10participants_03May2023.xlsx"))
# #Modify here as data is updated
# write_xlsx(fems_DS, paste0(curr_dir,"data/output/MSI_CD8_DS_analysis_AMW_10participants_03May2023.xlsx"))