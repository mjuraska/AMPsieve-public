#----------------------------------------------------------------------------------------
# PROGRAM: d_cmprsk.R
#
# DESCRIPTION: calculate CIR-based efficacy estimates using sieve categorical PAR score 3
#              (1) HVTN 704 
#              (2) scale of time since enrollment through week 80 visit
#              (3) pooled VRC01 vs. control
#              (4) censoring at tau
#              (5) competing risks analysis by PAR score 3
#              (6) for cases, endpoint is first positive date rather than diagnosis date
#
# CODED BY: Erika Rudnicki
#
# INPUT:  
#      ../adata/v704_survival_wk80_tau_sieve_v5_cam.csv
#
# OUTPUT: 
#      ../adata/v704_cmprisk_cuminc.csv
#      ../adata/v704_cmprisk_cir.csv
#
# MAINTENANCE HISTORY:
# Date           Programmer         Description
# ***************************************************************************************
# 21Jan2021      Erika Rudnicki     Version 1.0
# 24Mar2022      Erika Rudnicki     Drop 'r4p' references, new source data
#----------------------------------------------------------------------------------------

# setwd("/trials/vaccine/p704/analysis/sieve/code")

# reset library paths to use validated R in batch mode
.libPaths("")

# load macros
source("../macro/cuminc_functions.R")  

# data directory 
dataDir <- "../adata"

# input file name
dataFile <- "v704_survival_wk80_tau_sieve_v5_cam.csv"

# output file names
cuminc.csvFile <- "v704_cmprisk_cuminc.csv"
cir.csvFile <- "v704_cmprisk_cir.csv"

# specify variables names from dataset

  # name of variable containing groupings to be compared
  grpVar <- "tx_pool" 
  
  # name of variable containing follow-up time information
  timeVar <- "hiv1fpday"

  # name of variable containing info on which type of event has occurred
  # It's best to make this into a factor after it's read in, being careful to set the
  # eventType that corresponds to 'censoring' to be the first level of the factor.
  eventTypeVar <- "hiv1eventparscore3.mf"  

  # idVar - name of variable containing a unique identifier
  idVar <- "pub_id"

  # reference level of your group variable (e.g reference trt group) 
  refLvl <- "C3"

  # comparison level of your group variable
  cmpLvl <- c("T1+T2")

# read in time-to-event dataset
dat <- read.csv( file.path( dataDir, dataFile), stringsAsFactors = FALSE )

## Create eventTypeVar
table(dat[[eventTypeVar]])
dat[[eventTypeVar]] <- factor(dat[[eventTypeVar]], 
                              levels= c(0, 1, 2, 3), 
                              labels = c("Censored", "< 1", "1 to 3", "> 3") )
table(dat[[eventTypeVar]])

## Calculate cumulative incidence via competing risks.
## Note our input dataset is already censored how we need it to be.

cuminc_list <- cmpRisks(data=dat, 
                   futimeVar=timeVar, 
                   groupVar= grpVar,
                   eventVar= eventTypeVar)

print( str(cuminc_list ) )
print( levels(  cuminc_list$data$hiv1eventparscore3.mf ) )

print (table(cuminc_list$cuminc$eventType))


q(save = "no")

## get CIR/PE estimates for each cohort
cir_list <- EffCIR(cuminc_list, refLvl= refLvl, cmpLvl= cmpLvl)

## output cuminc/cumhaz estimates
write.csv(cuminc_list$cuminc, file= file.path(dataDir, cuminc.csvFile), na="", row.names=FALSE, quote=FALSE)

## output CIR/PE estimates
write.csv(cir_list$CIR, file= file.path(dataDir, cir.csvFile), na="", row.names=FALSE, quote=FALSE)

q(save = "no")
