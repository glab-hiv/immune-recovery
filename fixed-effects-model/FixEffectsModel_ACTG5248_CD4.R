#-------------------------------------------------------------------------------
# Title: ACTG 5248 analysis of CD4 T cell Memory Restoration
#
# Purpose: 
# Interested in looking at CD4 T cell memory restoration in people 
# living with HIV once they receive ART.
# 
# Question:
# How do percentages of these markers change as a function of time?
# 
# Two cohorts:
# Chronic HIV (CHI) Durably-Suppressed (DS)
#    - 10 participants
#    - 3 DIFFERENT timepoints for each.
# ACTG 5248
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
# Analysis: Quasi-binomial fixed-effects model with indicators for PID.
# Statistician: Ann Marie Weideman
# Date Created: 17 March 2022
# Last Modified: 13 Dec 2022
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
ACTG <- read.csv(file = paste0(curr_dir,'Data/updated_v2/ACTG CD4 subpopulations Event Counts_10 participants_v10_with FC.csv'))
DS <- read.csv(file = paste0(curr_dir,'Data/updated_v2/CD4 CHI Durably Suppressed Event Counts_v9.csv'))

# Generate indicator variables
#ACTG
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

#DS
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

#-------------------------------------
# Graph data to inspect trends
#-------------------------------------

library(ggplot2)
pdf(file = paste0(curr_dir,"data/output/CD4_ACTG_PT_plots.pdf"))
for(i in unique(ACTG$Population)) {
  print(ggplot(ACTG[which(ACTG$Population==i),], aes(x = Day/30, 
          y = qlogis(Percentages), colour = as.factor(PID))) +
          geom_point() + facet_wrap( ~ PID) + 
          labs(title="ACTG group", 
               x ="Time (months)", y = paste0("logit(Proportion of ",i,")")))
  #Sys.sleep(2)
}
dev.off()

for(i in unique(DS$Population)) {
  print(ggplot(DS[which(DS$Population==i),], aes(x = Months.post.timepoint.1, 
          y = qlogis(Percentages), colour = as.factor(PID))) +
          geom_point() + facet_wrap( ~ PID) + 
          labs(title="DS group", 
               x ="Time (months)", y = paste0("logit(Proportion of ",i,")")))
  #Sys.sleep(2)
}

ggplot(ACTG, aes(x = Day/30, y = qlogis(Percentages), 
                 colour = as.factor(PID))) + geom_point() + 
  facet_wrap( ~ Population) + 
  labs(title="ACTG Group",
       x ="Time (months)", y = "logit(Proportion of Marker)")

ggplot(DS, aes(x = Months.post.timepoint.1, y = qlogis(Percentages), 
               colour = as.factor(PID))) + geom_point() + 
  facet_wrap( ~ Population) + 
  labs(title="Durably Suppressed (DS) Group",
       x ="Time (months)", y = "logit(Proportion of Marker)")

ACTG$ThirtyDays <- ACTG$Day/30

#-----------------------------------------------------------
# ACTG cohort: Fixed-effects models with indicators for PID
#-----------------------------------------------------------
require(broom)
require(dplyr)

ACTG$ThirtyDays <- ACTG$Day/30
by_population_ACTG <- group_by(ACTG, Population)

require(margins)
fems_ACTG<-do(by_population_ACTG, 
              tidy(margins(glm(cbind(Event.count, (Denominator - Event.count))~
                      -1 + ThirtyDays + PID291374 + PID363043
                      + PID363044 + PID214273 + PID611175 + PID611183 
                      + PID611051  + PID21929 + PID112185 + PID611172, 
                      family=quasibinomial, data = .),
                      variables="ThirtyDays")))
# Compute CIs
fems_ACTG$CI_lower_decimal<-fems_ACTG$estimate-1.96*fems_ACTG$std.error
fems_ACTG$CI_upper_decimal<-fems_ACTG$estimate+1.96*fems_ACTG$std.error

# Add descriptions of results
fems_ACTG$description<- ifelse(fems_ACTG$p.value<0.05,
      paste0("For every 30 days post ART initiation, there is an estimated ",
      ifelse(sign(fems_ACTG$estimate)==-1,"decrease","increase"),
      " of ",round(abs(fems_ACTG$estimate)*100,2)," percentage points (",
      ifelse(fems_ACTG$p.value<0.0001,"p<0.0001",
             ifelse(fems_ACTG$p.value<0.001,"p<0.001",
                    ifelse(fems_ACTG$p.value<0.01,"p<0.01",
                      ifelse(fems_ACTG$p.value<0.05,"p<0.05",
                             paste0("p=",round(fems_ACTG$p.value,2)))))),
      
      "), on average, in ",
      fems_ACTG$Population,"."),
      paste0("There is no evidence of a change in ", fems_ACTG$Population,
             " across time.")
)

# Rename "estimate" to "slope_decimal" so PI understands
fems_ACTG<-dplyr::rename(fems_ACTG, slope_decimal = estimate)
fems_ACTG<-dplyr::rename(fems_ACTG, standard_error_decimal = std.error)
fems_ACTG<-dplyr::rename(fems_ACTG, p_value = p.value)
# Create slope percent, se percent, and CI percent
fems_ACTG$slope_percentage_points<-fems_ACTG$slope_decimal*100
fems_ACTG$standard_error_percentage_points<-fems_ACTG$standard_error_decimal*100
fems_ACTG$CI_lower_percentage_points<-fems_ACTG$CI_lower_decimal*100
fems_ACTG$CI_upper_percentage_points<-fems_ACTG$CI_upper_decimal*100

# Format estimates and p-values
fems_ACTG$slope_decimal<-format(fems_ACTG$slope_decimal, scientific=FALSE)
fems_ACTG$slope_percentage_points<-format(fems_ACTG$slope_percentage_points, scientific=FALSE)
fems_ACTG$p_value<-format(fems_ACTG$p_value, scientific=FALSE)

fems_ACTG$term<-"Every 30 days"
                          
#-----------------------------------------------------------
# DS cohort: Fixed-effects models with indicators for PID
#-----------------------------------------------------------

by_population_DS <- group_by(DS, Population)

require(margins)
fems_DS<-do(by_population_DS, 
            tidy(margins(glm(cbind(Event.count, (Denominator - Event.count)) ~ 
                             -1 + Months.post.timepoint.1 +
                              + PID546 + PID861 + PID728 + PID674 + PID1095 
                              + PID929 + PID749 + PID673 + PID870 + PID1036, 
                              family=quasibinomial, data = .),
                              variables="Months.post.timepoint.1")))

# Compute CIs
fems_DS$CI_lower_decimal<-fems_DS$estimate-1.96*fems_DS$std.error
fems_DS$CI_upper_decimal<-fems_DS$estimate+1.96*fems_DS$std.error

# Add descriptions of results
fems_DS$description<- ifelse(fems_DS$p.value<0.05,
                               paste0("For every 30 days while durably suppressed, there is an estimated ",
                                      ifelse(sign(fems_DS$estimate)==-1,"decrease","increase"),
                                      " of ",round(abs(fems_DS$estimate)*100,2)," percentage points (",
                                      ifelse(fems_DS$p.value<0.0001,"p<0.0001",
                                             ifelse(fems_DS$p.value<0.001,"p<0.001",
                                                    ifelse(fems_DS$p.value<0.01,"p<0.01",
                                                           ifelse(fems_DS$p.value<0.05,"p<0.05",
                                                                  paste0("p=",round(fems_DS$p.value,2)))))),
                                      
                                      "), on average, in ",
                                      fems_DS$Population,"."),
              paste0("There is no evidence of a change in ",fems_DS$Population,
                     " across time.")
)

# Rename "estimate" to "slope_decimal" so PI understands
fems_DS<-dplyr::rename(fems_DS, slope_decimal = estimate)
fems_DS<-dplyr::rename(fems_DS, standard_error_decimal = std.error)
fems_DS<-dplyr::rename(fems_DS, p_value = p.value)

# Create slope percent, se percent, and CI percent
fems_DS$slope_percentage_points<-fems_DS$slope_decimal*100
fems_DS$standard_error_percentage_points<-fems_DS$standard_error_decimal*100
fems_DS$CI_lower_percentage_points<-fems_DS$CI_lower_decimal*100
fems_DS$CI_upper_percentage_points<-fems_DS$CI_upper_decimal*100

# Format estimates and p-values
fems_DS$slope_decimal<-format(fems_DS$slope_decimal, scientific=FALSE)
fems_DS$slope_percentage_points<-format(fems_DS$slope_percentage_points, scientific=FALSE)
fems_DS$p_value<-format(fems_DS$p_value, scientific=FALSE)

fems_DS$term<-"Every 30 days"

#-------------------------------------
# Write to csv
#-------------------------------------
require(writexl)
write_xlsx(fems_ACTG, paste0(curr_dir,"data/output/CD4_ACTG_analysis_AMW_10participants_13Dec2022.xlsx"))
write_xlsx(fems_DS, paste0(curr_dir,"data/output/CD4_DS_analysis_AMW_10participants_13Dec2022.xlsx"))
