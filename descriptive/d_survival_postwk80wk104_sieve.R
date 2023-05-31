##########################################################################################
# Program : d_survival_ccwk80wk104_sieve.R
#
# Project: 
# - The below-mentioned dataset was created using T:\vaccine\p704\analysis\sieve\code\d_survival_wk80_tau_sieve_r4p.R 
#   and its dataset specifications are saved as T:\vaccine\p704\analysis\sieve\specs\readme_v704_survival_wk80_tau_sieve_r4p.txt.
# - Per Michal’s suggestion, let’s create a new program in the same directory 
#   (name it ‘d_survival_ccwk80wk104_sieve.R’) that takes the existing v704_survival_wk80_tau_sieve_r4p.csv 
#   dataset, adds on the new variables, and outputs as ../adata/v704_survival_ ccwk80wk104_sieve.csv
#   - Please also create a new dataset specifications file (starting from readme_v704_survival_wk80_tau_sieve_r4p.txt)
#     for v704_survival_ ccwk80wk104_sieve.csv. 
# - These are the new variables w/ suggested naming conventions (Michal feel free to tweak):
#   - ccprefl = MITT complete case week 80-104 preliminary cohort flag, where we call it 
#     ‘preliminary’ because it considers only the first two bullet points mentioned below 
#     (1=meets those two criteria and 0=does not meet those two criteria)
#     - Use T:\vaccine\p704\analysis\efficacy\adata\v704_survival.csv as the source data for this derivation. 
#       - It’s a little confusing, but the specifications for this dataset are saved in the 
#         ‘cuminc_hiv’ tab of T:\vaccine\p704\analysis\dsmb\2020_08\docs\AMP_specs_704_703.Xlsx. 
#         More to the point:
#         - Efficacy_flag = once-randomized MITT participants
#         - statuswk80 = primary endpoints observed by the week 80 visit
#         - statuswk104 = endpoints observed by the week 104 visit
#         - s1wk80dt = schedule 1 week 80 visit date, if completed
#         - s1wk80dttgt = target schedule 1 week 80 visit date, if week 80 visit was missed but 
#           there is evidence of HIV-1 negative status in schedule 1 post-week 80.
#         - s4wk80dt = hypothetical schedule 4 week 80 visit date (hypothetical because there 
#           is no week 80 visit in schedule 4)
#   - s4negpostwk80fl = flag indicating evidence of negative HIV-1 status at a post-hypothetical schedule 4 week 80 visit date
#   - ccpre_hiv1event = indicator of MITT non-primary HIV infected case (NA if ccprefl=0; otherwise 1=yes, 0=no)
#     - Can use same dataset/variables as mentioned for ccprefl (statuswk80 and statuswk104)
#   - ccpre_hiv1survday = follow-up time from week 80 visit to diagnosis date (NA if ccprefl=0; time unit is days)
#     - Can use same dataset/variables as mentioned for ccprefl (fudtwk104 gives end of follow-up date, 
#       and use the blue variables to calculate start of follow-up date)
#   - ccpre_hiv1fpday = follow-up time from week 80 visit to first RNA positive sample (NA if ccprefl=0; time unit is days)
#     - Can use same dataset/variables as mentioned for ccpre_hiv1survday; the first RNA positive sample 
#       date is calculated within ../code/d_survival_wk80_tau_sieve_r4p.R
#   - IC50/IC80 variables:
#     - I misspoke on our call, thinking that the IC50/IC80 data shouldn’t change. Rather, they will change
#       because now we are focused on endpoints who meet ccpre_hiv1event=1 criteria whereas before we only 
#       included IC50/IC80 data for primary endpoints. I think the new derivations should be straightforward 
#       though. Instead of referencing /trials/vaccine/p704/analysis/efficacy/adata/v704_survival_wk80_tau_neut.csv 
#       (as currently done in ../code/d_survival_wk80_tau_sieve_r4p.R), you can pull the information 
#       from /trials/vaccine/p704/analysis/efficacy/adata/v704_survival_wk104_tau_neut_gpdx.csv instead. 
#    
# Location: /trials/vaccine/p704/analysis/sieve/code
#
# Input:
#   ../adata/v704_survival_wk80_tau_sieve_r4p.csv 
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_survival.csv
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_survival_wk104_tau_neut_gpdx.csv 
#
# Output:
#   ../adata/v704_survival_ ccwk80wk104_sieve.csv
#
# Specifications: ../specs/readme_v704_survival_wk80_tau_sieve_r4p.txt
#
# Code History
# ---------------------------------------------------------------------------------------
# Date        Programmer        Details 
# ---------------------------------------------------------------------------------------
# 2021May26   Kevin Gillespie    Version 1.0
###########################################################################################

# setwd("/trials/vaccine/p704/analysis/sieve/code")
library(haven)
library(dplyr)

# Data import -------------------------------------------------------------

sm <- read_sas("/trials/vaccine/p704/analysis/efficacy/adata/subject_master.sas7bdat") %>% 
  mutate(
    ptid = gsub("-", "", ptid)
  )

df <- read.csv("../adata/v704_survival_wk80_tau_sieve_r4p_v2_cam.csv")

surv <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/v704_survival.csv") %>% 
  mutate(
    ptid  = gsub("-", "", ptid),
    enrdt = as.Date(enrdt, "%d%b%Y")
  )


ic <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/v704_survival_wk104_tau_neut_gpdx.csv") %>% 
  mutate(ptid = as.character(ptid))

vl <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/v704_viral_loads.csv")

low <- read.csv("/trials/vaccine/p704/s670/qdata/VTN704_iSCA01_20200901.txt", sep="\t")

mfppt <- read.csv("../adata/hvtn704_ren_pfitter_tf_lookup_v2.csv")

seq <- read.csv("../adata/amp_sequences_slapnap_vrc01_cruncher01_all_v1_28Dec2020_v704.csv")


# Data manipulation -------------------------------------------------------

low <-
  low %>%
  mutate(
    ptid = as.character(ptid),
    drawdt=as.Date(drawdt, "%d%b%Y"),
    lowresult=ifelse(is.na(copies_per_ml), "Negative", "Positive")
  ) %>%
  select(ptid, visitno, drawdt, copies_per_ml, lowresult)


vl <-
  vl %>%
  select(ptid, visitno, drawdt, result, resultc, vl)

# filter 'seq' to include only the MF founder per participant as identified in 'mfppt'
seqmf <-
  mfppt %>%
  rename(seqname=tf.sequence) %>%
  left_join(seq, by=c("pubid", "seqname")) %>%
  mutate(pub_id=gsub("_", "-", gsub("V", "", pubid))) %>%
  select(pub_id, seqname, pred.prob.sens, pred.ic80) 

# using 'seqmf' calculate par score variables based on MF founder
seqmfpar <-
  seqmf %>%
  rename(seqname.mf=seqname,
         pred.ic80.mf=pred.ic80, 
         pred.prob.sens.mf=pred.prob.sens) %>%
  mutate(pred.prob.res.mf=1-pred.prob.sens.mf,
         parscore1.mf=log(pred.prob.res.mf/(1-pred.prob.res.mf)),
         parscore2.mf=log10(pred.ic80.mf),
         parscore3.mf=case_when(pred.ic80.mf < 1 ~ "<1",
                                pred.ic80.mf >= 1 & pred.ic80.mf <= 3 ~ "[1,3]",
                                pred.ic80.mf > 3 ~ ">3"))

# switch over to calculating follow-up time based on first positive visit
vlfp <-
  vl %>%
  mutate(
    drawdt = as.Date(drawdt),
    ptid = as.character(ptid)
  ) %>% 
  left_join(low, by=c("ptid", "drawdt", "visitno")) %>%
  filter(result %in% c(1,10,11) | !is.na(copies_per_ml))  %>%
  arrange(ptid, drawdt) %>%
  group_by(ptid) %>%
  filter(drawdt==min(drawdt)) %>%
  ungroup() %>%
  select(ptid, drawdt) %>%
  rename(drawdtfp=drawdt)

vlfpsurv <-
  vlfp %>%
  left_join(surv[,c("ptid", "enrdt")], by="ptid") %>%
  mutate(hiv1fpday=as.numeric(drawdtfp-enrdt))


# Construct the additional variables Michal has requested
blu <- surv %>% 
  left_join(vlfpsurv, by = c("ptid" = "ptid", "enrdt")) %>%
  mutate(
    ltfu = ifelse(s1wk80dt == "" & s1wk80dttgt == "" & s4negpostwk80fl == 0, 1, 0),
    postwk80fl = ifelse(MITT_flag == 1 & statuswk80 == 0 & ltfu == 0, 1, 0),
    postwk80_hiv1event = case_when(
      postwk80fl == 0                    ~ "0",
      postwk80fl == 1 & statuswk104 == 1 ~ "1",
      postwk80fl == 1 & statuswk104 == 0 ~ "0"
    ),
    origin_time = case_when(
      s1wk80dt != "" & s1wk80dttgt == "" & ltfu == 0 ~ as.character(.$s1wk80dt),
      s1wk80dt == "" & s1wk80dttgt != "" & ltfu == 0 ~ as.character(.$s1wk80dttgt),
      s4negpostwk80fl == 1 & ltfu == 0 ~ as.character(.$s4wk80dt)
    ),
    origin_time = as.Date(origin_time, "%d%b%Y"),
    fudtwk80   = as.Date(as.character(fudtwk80), "%d%b%Y"),
    fudtwk104   = as.Date(as.character(fudtwk104), "%d%b%Y"),
    dxdt        = as.Date(as.character(dxdt), "%d%b%Y"),
    postwk80_hiv1survday = dxdt - origin_time,
    postwk80_hiv1fpday   = drawdtfp - origin_time
  )  %>% 
  left_join(select(sm, ptid, pub_id)) %>% 
  select(pub_id, postwk80fl, postwk80_hiv1event, postwk80_hiv1survday, postwk80_hiv1fpday)



# Combine the data with what already exists -------------------------------

untuch <- df %>% 
  select(
    protocol, southAmerica, pub_id, tx, tx_pool, rx, hiv1event, 
    hiv1survday, hiv1fpday, seqname.mf, pred.ic80.mf, pred.prob.sens.mf,
    pred.prob.res.mf, parscore1.mf, parscore2.mf, parscore3.mf,
    hiv1eventparscore3.mf
  )

cam_df <- df %>% 
  select(pub_id, 31:ncol(df))

ic_merge <- ic %>% 
  select( -hiv1event, -hiv1survday, -ptid, -gmt50mscat, -gmt80mscat,
          -gmt50lscat, -gmt80lscat, -gmt50mfcat, -gmt80mfcat, 
          -p_1st4_overall, -hiv1eventType.8, -hiv1eventType.9)


df_new <- untuch %>% 
  left_join(ic_merge, by = c("protocol", "southAmerica", "pub_id", "tx", "tx_pool", "rx")) %>% 
  left_join(cam_df) %>% 
  left_join(blu)



df_new %>% 
  filter(postwk80_hiv1event == 1)

id_list <- df_new %>% 
  filter(hiv1event == 1) %>% 
  head()

old_list <- df %>% 
  filter(pub_id %in% id_list$pub_id)


dim(df)
dim(df_new)

all_equal(df, df_new[,1:(ncol(df))])

# write out the results ---------------------------------------------------

write.csv(df_new, "../adata/v704_survival_postwk80wk104_sieve.csv")
