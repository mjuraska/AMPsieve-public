##########################################################################################
# Program : d_survival_ccwk80wk104_sieve.R
#
# Project: 
# - Consort diagrame has:
#   - 140 infections (137 /efficacy/adata/v704_survival.csv)
#   - 130 MITT infections (127 /efficacy/adata/v704_survival.csv)
#     - 98 primary
#     - 32 non-primar (27 /efficacy/adata/v704_survival.csv)
# - Why is this dicrepancy observed? 
#   - Pull FSR data and infection tallies and inspect agains the efficacy data
#    
# Location: /trials/vaccine/p704/analysis/sieve/code
#
# Input:
#   ../adata/v704_survival_wk80_tau_sieve_r4p.csv 
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_survival.csv
#   /trials/vaccine/p704/analysis/fsr/closed/adata/survival.csv
#   /trials/vaccine/p704/analysis/fsr/closed/adata/hiv_infections.sas7bdat
#
# Output:
#
# Specifications: ../specs/readme_v704_survival_wk80_tau_sieve_r4p.txt
#
# Code History
# ---------------------------------------------------------------------------------------
# Date        Programmer        Details 
# ---------------------------------------------------------------------------------------
# 2021May26   Kevin Gillespie    Version 1.0
###########################################################################################

library(dplyr)
library(haven)


fsrdir703 <- "/trials/vaccine/p703/analysis/fsr/closed/adata/"
fsrdir704 <- "/trials/vaccine/p704/analysis/fsr/closed/adata/"

effdir703 <- "/trials/vaccine/p703/analysis/efficacy/adata/"
effdir704 <- "/trials/vaccine/p704/analysis/efficacy/adata/"


surv_new_703 <- read.csv(paste0(fsrdir703, "survival.csv"))
surv_new_704 <- read.csv(paste0(fsrdir704, "survival.csv"))
                         
                         
surv_new <- bind_rows(surv_new_703, surv_new_704) %>% 
  mutate(
    ptid  = gsub("-", "", ptid),
    enrdt = as.Date(enrdt, "%d%b%Y")
  )

surv_old_703 <- read.csv(paste0(effdir703, "v703_survival.csv"))
surv_old_704 <- read.csv(paste0(effdir704, "v704_survival.csv"))

surv_old <- bind_rows(surv_old_703, surv_old_704) %>% 
  mutate(
    ptid  = gsub("-", "", ptid),
    enrdt = as.Date(enrdt, "%d%b%Y")
  )


hiv703 <- read_sas(paste0(fsrdir703, "hiv_infections.sas7bdat"))
hiv704 <- read_sas(paste0(fsrdir704, "hiv_infections.sas7bdat"))

hiv <- bind_rows(hiv703, hiv704) %>% 
  mutate(
    ptid  = gsub("-", "", ptid),
    enrdt = as.Date(enrdt, "%d%b%Y")
  )


# How many hiv records as of 2020-04-03 (when consort constructed)
hiv %>% filter(
  dxdt < "2020-04-03" 
)


# How many records after then and before the efficacy report went live
hiv %>% filter(
  dxdt >= "2020-04-03" & dxdt <= "2020-08-08"
)

# How many were after 04/03/2020
hiv %>% filter(
  dxdt >= "2020-04-03"
)

# Understanding discrepancy -----------------------------------------------

# Records in FSR dataset that were primary endpoints and had a record in hiv infections
test_new <- surv_new %>% 
  filter(ptid %in% hiv$ptid )
# 127 infections

test_old <- surv_old %>% 
  filter(ptid %in% hiv$ptid)


dim(test_new)


table(test_new$hiv1_event)

# Look at the rows that were eventually determined to be infected,
# that were not at the time the two survival datasets were constructed
test_new %>% filter(hiv1_event != 1)
test_old %>% filter(hiv1_event != 1)

# The collection of ppts who did become infected, but were not when 'survival'
# was created
# - The are HIV infected based on HIV infections data from FSR
# - Were not HIV infected based on efficacy survival
uninf <- test_new %>% 
  filter(hiv1_event == 1)

head(hiv)

# Now pull out the 142-127 = 15 who were not infected in survival
xtra_ids <- hiv[which(!(hiv$ptid %in% uninf$ptid)),]

# List important info about the 15
xtra_ids %>% select(pub_id, ptid, rx_code, dxday, efficacy_flag, safety_flag, MITT_flag)

# Focus on the ones that were also MITT
xtra_ids %>% filter(MITT_flag == 1) %>% select(pub_id, ptid, rx_code, dxday, efficacy_flag, safety_flag)

# Look at those 5 from survival perspective
surv_old %>% filter(ptid %in% xtra_ids$ptid & MITT_flag == 1)
surv_new %>% filter(ptid %in% xtra_ids$ptid & is.na(MITT_flag))
# Important notes:
# - the 5 that are not infections in survival are all dxdt > 04-03-2020
# - 1 is control and the other 4 are treated
#   - the control is at end of week 80 window, thus a primary endpoint
# - 2 were lost to follow up
#   - so 127 + (5 - 2) = 130
#   - also 27 + 5 added hiv infections == 32