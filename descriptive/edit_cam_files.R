##########################################################################################
# Program : edit_cam_files.R
#
# Project: Use the most recent v704_survival_wk80_tau_sieve.csv file to update
#          v704_survival_wk80_tau_sieve_v5_cam.csv with new survival data.
#    
# Location: /trials/vaccine/p704/analysis/sieve/code
#
# Input:
#   ../adata/v704_survival_wk80_tau_sieve.csv
#   ../adata/v704_survival_wk80_tau_sieve_v5_cam.csv
#
# Output:
#   ../adata/v704_survival_wk80_tau_sieve_v5_cam_er.csv
#
# Specifications: ../specs/readme_v70x_survival_wk80_tau_sieve.txt
#
# Code History
# ---------------------------------------------------------------------------------------
# Date        Programmer        Details 
# ---------------------------------------------------------------------------------------
# 2022Jun16   Erika Rudnicki    Version 1.0
###########################################################################################

# setwd("/trials/vaccine/p704/analysis/sieve/code")

library(dplyr)
library(lubridate)
options(dplyr.print_max = 1e9, stringsAsFactors=FALSE)

dat1 <- read.csv("../adata/v704_survival_wk80_tau_sieve.csv")
dat2 <- read.csv("../adata/v704_survival_wk80_tau_sieve_v5_cam.csv")

colnames(dat1)
colnames(dat2)

dat_out <-
  dat2 %>%
  select(-hiv1survday, -hiv1fpday) %>%
  left_join(dat1[,c("pub_id", "hiv1survday", "hiv1fpday")], by="pub_id")

write.csv(dat_out, "../adata/v704_survival_wk80_tau_sieve_v5_cam_er.csv", row.names=FALSE)

q(save="no")

