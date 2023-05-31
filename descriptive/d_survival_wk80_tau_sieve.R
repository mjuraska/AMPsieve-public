##########################################################################################
# Program : d_survival_wk80_tau_sieve.R
#
# Project: Create survival dataset with sieve PAR score and neut sensitivity data 
#          based on the following parameters:
#           - HVTN 704.
#           - Enrollment through the week 80 visit, with follow-up time for endpoints
#             defined as (1) diagnosis date, and (2) first positive result.
#           - Censored at Tau.
#           - This dataset is an update to the one used for January 2021 HIV R4P presentations, 
#             (those programs/dataset are named identically with 'r4p' appended).
#    
# Location: /trials/vaccine/p704/analysis/sieve/code
#
# Input:
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_survival_wk80_tau.csv
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_survival.csv
#   ../adata/hvtn704_ren_tf_lookup_v2.csv
#   ../adata/par_score_v704_22Nov2021.csv
#   /trials/vaccine/p704/analysis/efficacy/adata/v704_viral_loads.csv
#   /trials/vaccine/p704/s670/qdata/VTN704_iSCA01_20200901.txt
#
# Output:
#   ../adata/v704_survival_wk80_tau_sieve.csv
#
# Specifications: ../specs/readme_v704_survival_wk80_tau_sieve_r4p.txt
#
# Code History
# ---------------------------------------------------------------------------------------
# Date        Programmer        Details 
# ---------------------------------------------------------------------------------------
# 2021Jan19   Erika Rudnicki    Version 1.0
# 2021Dec09   Erika Rudnicki    Updated 'seq' and 'mfppt' input data from Craig (see below)
# 2022Jun16   Erika Rudnicki    After discussion with Peter and Michal, new right-censoring
#                               methods are to be applied to this dataset. Namely,
#                               we no longer consider specific visits for primary follow-up,
#                               so uninfected ppts are censored at 
#                               min(last HIV negative visit - enrollment dt, tau=595 days/85 weeks).
###########################################################################################

# setwd("/trials/vaccine/p704/analysis/sieve/code")

library(dplyr)
library(lubridate)
options(dplyr.print_max = 1e9, stringsAsFactors=FALSE)

# load input data
# -- 'survneut' is the primary neutralization sensitivity analysis dataset
# -- 'survprim' is the primary analysis dataset with enrollment date
# -- 'mfppt' is sent by C. Magaret and identifies the most frequent founder per participant
# -- 'seq' is sent by C. Magaret and includes IC80 predictions for all founders (PAR Scores)
# -- 'vl' includes our standard RNA PCR viral load data
# -- 'low' includes low copy RNA viral load data  
# -- 'master' includes enrollment date

survneut <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/v704_survival_wk80_tau_neut.csv")
survprim <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/v704_survival.csv")
mfppt <- read.csv("../adata/hvtn704_ren_tf_lookup_v2.csv")
seq <- read.csv("../adata/par_score_v704_22Nov2021.csv")
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
# per email from Craig, replace 'REN' portion of mfppt.tf.sequence with 'Env' in order  
# to match with seq.seqname values. 
seqmf <-
  mfppt %>%
  mutate(seqname=gsub("REN", "Env", tf.sequence)) %>%
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

# here is where we can make the 16JUN2022 right-censoring adjustment
dat_final2 <-
  survprim %>%
  select(ptid, enrdt, lastnegdt, dxdt, statuswk80, fudayswk80) %>%
  right_join(survneut[,c("ptid", "pub_id")], by="ptid") %>%
  select(-ptid) %>%
  mutate(lastnegdt=dmy(lastnegdt),
         dxdt=dmy(dxdt),
         hiv1survday_v2=pmin(as.numeric(if_else(statuswk80==1, dxdt-enrdt, lastnegdt-enrdt)), 595)) %>% 
  right_join(dat_final, by="pub_id") %>%
  mutate(diff=hiv1survday_v2-hiv1survday,
         hiv1fpday_v2=ifelse(hiv1event==1, hiv1fpday, hiv1survday_v2))

dat_final2 %>%
  group_by(statuswk80, hiv1event) %>%
  summarise(min_diff=min(diff),
            max_diff=max(diff),
            mean_diff=mean(diff),
            median_diff=median(diff),
            min_time=min(hiv1survday),
            max_time=max(hiv1survday),
            min_time_v2=min(hiv1survday_v2),
            max_time_v2=max(hiv1survday_v2),
            n_at_tau=sum(hiv1survday==595),
            n_at_tau_v2=sum(hiv1survday_v2==595))

dat_final2 %>%
  group_by(statuswk80, hiv1event) %>%
  summarise(min_diff=min(hiv1fpday_v2-hiv1fpday),
            max_diff=max(hiv1fpday_v2-hiv1fpday),
            mean_diff=mean(hiv1fpday_v2-hiv1fpday),
            median_diff=median(hiv1fpday_v2-hiv1fpday),
            min_time=min(hiv1fpday),
            max_time=max(hiv1fpday),
            min_time_v2=min(hiv1fpday_v2),
            max_time_v2=max(hiv1fpday_v2),
            n_at_tau=sum(hiv1fpday==595),
            n_at_tau_v2=sum(hiv1fpday_v2==595))

nrow(filter(dat_final2, diff < 0)) # ppts previously censored at tau=609 who are now censored earlier
nrow(filter(dat_final2, diff == 0)) # ppts who dropped out early and data doesn't change, also includes all primary endpoints
nrow(filter(dat_final2, diff > 0)) # ppts who remained uninfected after week 80 and we're extending their follow-up time now

dat_final2 <-
  dat_final2 %>%
  mutate(hiv1survday=hiv1survday_v2,
         hiv1fpday=hiv1fpday_v2) %>%
  select(protocol, southAmerica, pub_id, tx, tx_pool, rx, hiv1event, hiv1survday, hiv1fpday,
         seqname.mf, pred.ic80.mf, pred.prob.sens.mf, pred.prob.res.mf, parscore1.mf, parscore2.mf,
         parscore3.mf,  hiv1eventparscore3.mf, nisolates, gmt50ms, gmt80ms, gmt50ls, gmt80ls, gmt50mf, gmt80mf,
         hiv1event50ms, hiv1event80ms, hiv1event50ls, hiv1event80ls, hiv1event50mf, hiv1event80mf) %>%
  arrange(pub_id)

# Output final (v2) dataset
write.csv(dat_final2, "../adata/v704_survival_wk80_tau_sieve.csv", row.names=FALSE)

source("../macro/cuminc_functions.R")  

# groupings to be compared
grpVar_pool  <- "tx_pool"
grpVar_ind <- "tx"

# follow-up time information
timeVar <- "hiv1survday"

# event indicator
eventIndVar <- "hiv1event"

# unique identifier
idVar <- "pub_id"

# the strata variable
strataVar <- "tx"

# reference level of your group variable
refLvl <- "C3"

# comparison level of your group variable
cmpLvl_pool <- c("T1+T2")
cmpLvl_ind <- c("T1", "T2")

strataWts <- c(C3 = 1, T2=0.5, T1=0.5)

mitt.cuminc_pool <- 
  naCumInc( 
    data = dat_final2,
    futimeVar = timeVar, 
    eventVar = eventIndVar, 
    groupVar = grpVar_pool, 
    strataVar = strataVar,
    idVar = idVar,
    stratifiedEstimator = TRUE,
    strataWeights = strataWts)

mitt.CIR_pool <- EffCIR( mitt.cuminc_pool, refLvl = refLvl, cmpLvl=cmpLvl_pool, nullHypEff=0)

mitt.CIR_pool$eff

q(save="no")

