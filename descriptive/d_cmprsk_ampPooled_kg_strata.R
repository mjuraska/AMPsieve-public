#----------------------------------------------------------------------------------------
# PROGRAM: d_cmprsk_ampPooled.R
#
# DESCRIPTION: calculate CIR-based efficacy estimates using sieve categorical PAR score 3
#              (1) AMP pooled
#              (2) scale of time since enrollment through week 80 visit
#              (3) pooled VRC01 vs. control
#              (4) censoring at tau
#              (5) competing risks analysis by PAR score 3
#              (6) for cases, endpoint is first positive date rather than diagnosis date
#
# CODED BY: Erika Rudnicki
#
# INPUT:  
#      ../adata/v704and703_survival_wk80_tau_sieve_v5_cam_er.csv
#
# OUTPUT: 
#      ../adata/v704and703_cmprisk_cuminc.csv
#      ../adata/v704and703_cmprisk_cir.csv
#
# MAINTENANCE HISTORY:
# Date           Programmer         Description
# ***************************************************************************************
# 21Jan2021      Erika Rudnicki     Version 1.0
# 24Mar2022      Erika Rudnicki     Drop 'r4p' references, new source data
# 24Mar2022      Erika Rudnicki     Update from v704 to AMP pooled
# 21Jun2022      Erika Rudnicki     Re-run on new input data
# 26Aug2022      Kevin Gillespie    Re-run on new input data & modified plots
# 20Sep2022      Kevin Gillespie    Re-run with updated sap filtering
#----------------------------------------------------------------------------------------

# setwd("/trials/vaccine/p704/analysis/sieve/code")

# reset library paths to use validated R in batch mode
.libPaths("")

# load macros
source("../macro/cuminc_functions.R")  

# data directory 
dataDir <- "../adata"

# input file name
dataFile <- "amp_sieve_pooled_marks_final_v5.csv"

# output file names
cuminc.csvFile <- c("v703and704_cmprisk_cuminc_mf_strata_25102022.csv", "v703and704_cmprisk_cuminc_ms_strata_25102022.csv", "v703and704_cmprisk_cuminc_ls_strata_25102022.csv")
cir.csvFile <- c("v703and704_cmprisk_cir_mf_strata_25102022.csv", "v703and704_cmprisk_cir_ms_strata_25102022.csv", "v703and704_cmprisk_cir_ls_strata_25102022.csv")

# specify variables names from dataset

  # name of variable containing groupings to be compared
  grpVar <- "tx_pool" 
  
  # name of variable containing follow-up time information
  timeVar <- "hiv1fpday"

  # name of variable containing info on which type of event has occurred
  # It's best to make this into a factor after it's read in, being careful to set the
  # eventType that corresponds to 'censoring' to be the first level of the factor.
  eventTypeVar <- c("hiv1eventparscore3.mf", "hiv1eventparscore3.ms", "hiv1eventparscore3.ls")   

  # idVar - name of variable containing a unique identifier
  idVar <- "pub_id"

  # reference level of your group variable (e.g reference trt group) 
  refLvl <- "C3"

  # comparison level of your group variable
  cmpLvl <- c("T1+T2")
  
  # Stratification variable?
  strataVar <- "strata"

# read in time-to-event dataset
dat <- read.csv( file.path( dataDir, dataFile), stringsAsFactors = FALSE )

# delete infected ppts without neut data
N_b4 <- nrow(dat)
nrow(subset(dat, is.na(hiv1fpday) & (is.na(parscore3.mf)))) # Michal edit
dat <- subset(dat, !is.na(hiv1fpday) & !(hiv1event==1 & is.na(parscore3.mf))) 
N_aft <- nrow(dat)

N_b4 - N_aft

dat$parscore3.mf <- ifelse(dat$hiv1event==0, NA, dat$parscore3.mf)


dat$strata <- sapply(1:nrow(dat), function(x) {
  strata1 <- gsub('HVTN ', '', dat[x,'protocol'])
  strata2 <- ifelse(strata1 == '703', dat[x,'southAfrica'], dat[x,'southAmerica'])
  
  strata  <- paste0(strata1, strata2)
  factor(
    strata, 
    levels = c('7030', '7031', '7040', '7041'), 
    labels = c('703/notSAfrica', '703/SAfrica', '704/notSAmerica', '704/SAmerica')
  )
  
})

for(i in 1:length(eventTypeVar)){
  
  typeVar <- eventTypeVar[i]
  
  ## Create eventTypeVar
  table(dat[[typeVar]])
  dat[[typeVar]] <- factor(dat[[typeVar]], 
                                levels= c(0, 1, 2, 3), 
                                labels = c("Censored", "< 1", "1 to 3", "> 3") )
  table(dat[[typeVar]])
  
  ## Calculate cumulative incidence via competing risks.
  ## Note our input dataset is already censored how we need it to be.
  
  cuminc_list <- cmpRisks(data=dat, 
                          futimeVar=timeVar, 
                          groupVar= grpVar,
                          # strataVar=strataVar,
                          eventVar= typeVar)
  
  ## get CIR/PE estimates for each cohort
  cir_list <- EffCIR(cuminc_list, refLvl= refLvl, cmpLvl= cmpLvl)
  
  ## output cuminc/cumhaz estimates
  write.csv(cuminc_list$cuminc, file= file.path(dataDir, cuminc.csvFile[i]), na="", row.names=FALSE, quote=FALSE)
  
  ## output CIR/PE estimates
  write.csv(cir_list$CIR, file= file.path(dataDir, cir.csvFile[i]), na="", row.names=FALSE, quote=FALSE)
}

q(save = "no")
  