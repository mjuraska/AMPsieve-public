# Purpose: Primary endpoint analysis of PE by predicted IC80 category (PAR score 3)
#          Point estimates and 95% CI for the mark-specific hazard-based PE(j), j=1,2,3, 
#          for the trichotomous mark "< 1 ug/ml" (j=1) vs. "[1, 3] ug/ml" (j=2) vs. "> 3 ug/ml" (j=3)
#          Two-sided p-value from the Wald test of H0: PE(j)=0, separately for j=1,2,3
#          Modified Lunn and McNeil test (Biometrics, 1995) p-value from the test of H0: PE(1)=PE(2)=PE(3)
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Competing risks Cox model
#          Modified Lunn and McNeil (1995, Biometrics)
# Input:   AMP adjudicated primary endpoints
#          Failure types calculated using the continuous predicted IC80 (PAR score 2)
# Output:  A single PDF file for 704, 703, and the trial-pooled analysis
# Author:  Michal Juraska
# Date:    Jan 21, 2021

rm(list=ls(all=TRUE))

# Global variables and sourced code ---------------------------------------

datDir <- paste0("t:/vaccine/p704/analysis/sieve/adata")
outDir <- "t:/vaccine/p704/analysis/sieve/code/markHazPE/Routput"
figDir <- paste0("t:/vaccine/p704/analysis/sieve/figures/markHazPE")

source("t:/vaccine/p704/analysis/sieve/code/myFunctions.R")

source("t:/vaccine/p704/analysis/sieve/code/common.R")

variant <- c("mf", "ms", "ls")

# tags that are part of the output .RData and .pdf file names
trialFileString <- c("704", "703", "704and703")

# Run markHazTE() ---------------------------------------------------------

for (trial in 1:length(datFile)){
  data <- read.csv(file.path(datDir, datFile[trial])) %>%
    # two 703 ppts with missing sequences have also a missing time-to-event
    filter(!is.na(hiv1fpday)) %>%
    # stratification variable per SAP Section 5
    mutate(strataVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                               protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                               protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                               protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))
  
  # this exclusion criterion is outdated
  # discard cases with missing neutralization sensitivity data per Section 1.1 in the AMP Genotypic Sieve SAP (analysis cohort definitions)
  # data <- subset(data, !(hiv1event==1 & is.na(gmt80mf)))
  
  for (vt in variant){
    data1 <- subset(data, select=c("tx_pool", "hiv1fpday", "hiv1event", paste0("hiv1eventparscore3.", vt), "strataVar"))
    
    colnames(data1) <- c("tx", "ftime", "fstatus", "ftype", "strataVar")
    data1$tx <- ifelse(as.character(data1$tx)=="C3", 0, 1)
    
    data1$ftype <- ifelse(data1$fstatus==0, 0, data1$ftype)
    
    # discard all rows with a missing failure type
    data1 <- subset(data1, !is.na(ftype))
    
    saveFile <- paste0(trialFileString[trial], "_PEbyPredIC80categ_", vt, "_MITT_VRC01poooled_placebo.RData")
    markHazTE(data1, ftypeCodeVector=1:3, stratified=TRUE, saveDir=outDir, saveFile=saveFile)
  }
}

# Make a single plot of all estimates and p-values ------------------------

for (trial in 1:length(datFile)){
  if (!dir.exists(figDir)){ dir.create(figDir, recursive=TRUE) }
  
  for (vt in variant){
    pdf(file.path(figDir, paste0(trialFileString[trial], "_PEbyPredIC80categ_", vt, "_MITT_VRC01poooled_placebo.pdf")), width=7, height=6)
    plotMarkHazTE(loadFile=paste0(trialFileString[trial], "_PEbyPredIC80categ_", vt, "_MITT_VRC01poooled_placebo.RData"), 
                  dir=outDir, 
                  xLab=expression("Predicted" ~ IC[80] ~ "Category"),
                  xLabLine=4,
                  yLab="Prevention Efficacy (%)",
                  xTickLab=c(expression("<" ~ 1 ~ mu * "g/ml"), expression("[1, 3]" ~ mu * "g/ml"), expression(">" ~ 3 ~ mu * "g/ml")),
                  panelLab=NULL,
                  yLim=c(-1, 1.13),
                  yPvalAt=1.03,
                  parMar=c(5.8, 6.2, 1, 1.5))
    dev.off()
  }
}


# Auxiliary analyses ------------------------------------------------------

# data <- read.csv(file.path(datDir, datFile[1]))
# data$tx <- ifelse(as.character(data$tx_pool)=="C3", 0, 1)
# fit <- coxph(Surv(hiv1survday, hiv1event) ~ tx, data=data)
# 1 - summary(fit)$coef[1, 2]

# load(file.path(outDir, paste0("markHazTE_trial=", trialFileString[1], ".RData")), verbose=TRUE)
# PE1st4 <- 1 - exp(sapply(fitMI, function(fit1imp){ fit1imp[[1]]$logHazRatio }))
# PEnot1st4 <- 1 - exp(sapply(fitMI, function(fit1imp){ fit1imp[[2]]$logHazRatio }))
# d <- data.frame("PE_1st4"=paste0(round(PE1st4 * 100, digits=0), "%"), "PE_not1st4"=paste0(round(PEnot1st4 * 100, digits=0), "%"))
# rownames(d) <- paste0("Imputation ", 1:10)
# d