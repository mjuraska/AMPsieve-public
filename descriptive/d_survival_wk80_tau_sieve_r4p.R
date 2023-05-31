##########################################################################################
# Program : d_survival_wk80_tau_sieve_r4p.R
#
# Project: Create survival dataset with sieve PAR score and neut sensitivity data 
#          based on the following parameters:
#           - HVTN 704.
#           - Enrollment through the week 80 visit, with follow-up time for endpoints
#             defined as (1) diagnosis date, and (2) first positive result.
#           - Censored at Tau.
#           - This dataset to be used for January 2021 HIV R4P presentations, hence 'r4p'
#             is appended to its name; there will likely be future versions too.
#    
# Location: /trials/vaccine/p704/analysis/sieve/code
#
# Input:
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_survival_wk80_tau.csv
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_survival.csv
#   ../adata/hvtn704_ren_pfitter_tf_lookup_v2.csv
#   ../adata/amp_sequences_slapnap_vrc01_cruncher01_all_v1_28Dec2020_v704.csv
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_viral_loads.csv
#   /trials/vaccine/p704/s670/qdata/VTN704_iSCA01_20200901.txt
#
# Output:
#   ../adata/v704_survival_wk80_tau_sieve_r4p.csv
#
# Specifications: ../specs/readme_v704_survival_wk80_tau_sieve_r4p.txt
#
# Code History
# ---------------------------------------------------------------------------------------
# Date        Programmer        Details 
# ---------------------------------------------------------------------------------------
# 2021Jan19   Erika Rudnicki    Version 1.0
###########################################################################################

# setwd("/trials/vaccine/p704/analysis/sieve/code")

library(dplyr)
library(lubridate)
options(dplyr.print_max = 1e9, stringsAsFactors=FALSE)

# load input data
# -- 'survneut' is the primary neutralization sensitivity analysis dataset
# -- 'survprim' is the primary analysis dataset with enrollment date
# -- 'mfppt' is sent by C. Magaret and identifies the most frequent founder per participant
# -- 'seq' is sent by C. Magaret and includes IC80 predictions for all founders
# -- 'vl' includes our standard RNA PCR viral load data
# -- 'low' includes low copy RNA viral load data  
# -- 'master' includes enrollment date

survneut <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/v704_survival_wk80_tau_neut.csv")
survprim <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/v704_survival.csv")
mfppt <- read.csv("../adata/hvtn704_ren_pfitter_tf_lookup_v2.csv")
seq <- read.csv("../adata/amp_sequences_slapnap_vrc01_cruncher01_all_v1_28Dec2020_v704.csv")
vl <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/v704_viral_loads.csv")
low <- read.csv("/trials/vaccine/p704/s670/qdata/VTN704_iSCA01_20200901.txt", sep="\t")

# process data as needed - will need updated for 703
survprim <-
  survprim %>%
  mutate(ptid=as.numeric(gsub("-", "", ptid)),
         enrdt=dmy(enrdt))

vl <-
  vl %>%
  mutate(drawdt=ymd(drawdt)) %>%
  select(ptid, visitno, drawdt, result, resultc, vl)

low <-
  low %>%
  mutate(drawdt=dmy(drawdt),
         lowresult=ifelse(is.na(copies_per_ml), "Negative", "Positive")) %>%
  select(ptid, visitno, drawdt, copies_per_ml, lowresult)

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
  left_join(survprim[,c("ptid", "enrdt")], by="ptid") %>%
  mutate(hiv1fpday=as.numeric(drawdtfp-enrdt))

# finally combine all these new data with 'survneut'
# note that we limit our new data (which includes all MITT cases) 
# to only the primary endpoint cases as indicated by 'survneut' 
dat_inf <-
  survneut %>%
  filter(hiv1event==1) %>%
  left_join(seqmfpar, by="pub_id") %>%
  left_join(vlfpsurv, by="ptid")
  
dat_uninf <-
  survneut %>%
  filter(hiv1event==0)

dat_final <-
  bind_rows(dat_inf, dat_uninf) %>%
  mutate(hiv1eventparscore3.mf=case_when(hiv1event==0 ~ 0,
                                         hiv1event==1 & parscore3.mf=="<1" ~ 1,
                                         hiv1event==1 & parscore3.mf=="[1,3]" ~ 2,
                                         hiv1event==1 & parscore3.mf==">3" ~ 3),
         hiv1fpday=ifelse(hiv1event==1, hiv1fpday, hiv1survday)) %>%
  select(protocol, southAmerica, pub_id, tx, tx_pool, rx, hiv1event, hiv1survday, hiv1fpday,
         seqname.mf, pred.ic80.mf, pred.prob.sens.mf, pred.prob.res.mf, parscore1.mf, parscore2.mf,
         parscore3.mf,  hiv1eventparscore3.mf, nisolates, gmt50ms, gmt80ms, gmt50ls, gmt80ls, gmt50mf, gmt80mf,
         hiv1event50ms, hiv1event80ms, hiv1event50ls, hiv1event80ls, hiv1event50mf, hiv1event80mf) %>%
  arrange(pub_id)

# Output final dataset
write.csv(dat_final, "../adata/v704_survival_wk80_tau_sieve_r4p.csv", row.names=FALSE)

q(save="no")
