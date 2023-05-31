# Purpose: predicting treatment assignment using pairwise variables among epitope distance, hamming distances, gmt80, parscores
# Input:  Master file (.csv) with time to event data, treatment group assignment, region (south america vs. other) and binary marks to be analyzed
# Output: .csv and forest-plots containing results table with estimates of odds ratio, and associated confidence intervals   
# Author:  Li Li
# Date:    Jan 20, 2023


rm(list=ls(all=TRUE))
repoDir <- "/Users/lili/AMPsieve"
figuresOutDir <- file.path(repoDir, "figures/DPE/tier2")

library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)

source(file.path(repoDir, "code/common.R")) #read in a datFile
source(file.path(repoDir,"code/DPE/DPEutils.R"))
source(file.path(repoDir, "code/DPE/forest.R"))

master <-  read.csv( file.path(repoDir,"adata", datFile[3])) 
master$gmt80ls[master$gmt80ls == ">100"] <- 100
master$gmt80ls <- log10(as.numeric(master$gmt80ls))

format.or <- function(x)format(x,nsmall=2,digits=2)
paste.ci <- function(x){return(paste(format.or(x[1])," (" , format.or(x[2]),", ",  format.or(x[3]), ")", sep = ""))}
predictors.all <- c("parscore1.ls", "epitope.dist.b.ls","epitope.dist.c.ls" ,
                    "hdist.zspace.sites.binding.all.ls","gmt80ls","")
predictors.all.lab <- c("Pred prob IC80 > 1 µg/ml", "Epitope dist to subtype B ref",
                           "Epitope dist to subtype C ref", "PC-weighted Hamming dist",
                           "Observed IC80", "")
names(predictors.all.lab) <- predictors.all
colors.all <- c("royalblue","darkred","purple","turquoise","gold","black")
names(colors.all) <- predictors.all.lab
for(trial in c("704", "703")){
  master.trial <- trial_dose_data (master, trial, "T1+T2")
  master.trial$treatment <- 1*(master.trial$tx_pool=="T1+T2")
  if(trial == "704"){
    predictors <- c("parscore1.ls", "epitope.dist.b.ls","hdist.zspace.sites.binding.all.ls","gmt80ls")
  }else{
    predictors <- c("parscore1.ls", "epitope.dist.c.ls","hdist.zspace.sites.binding.all.ls","gmt80ls")
  }
  #standardize the covariates
  for(i in 1:length(predictors)){
    master[,predictors[i]] <- scale(master[,predictors[i]])
  }
  results <- tibble("model" = character(), "predictor"= character(), "odds.estimate" = character(),
                    "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),"p" = character())
  k = 0
  for(i in 1:3){
    for(j in (i+1):4){
      k = k+1
      covx <- c(predictors[i],predictors[j], "stratVar")
      formula = as.formula(paste0("treatment~",paste0(covx,collapse = "+")))
      fit.summary <- data.frame(summary(glm(formula = formula, family = "binomial", data = master.trial))$coefficients)[,c(1,2,4)]
      colnames(fit.summary) <- c("Estimate", "Std.error", "P")
      fit.summary$odds.estimate <- round(exp(fit.summary$Estimate),2)
      fit.summary$odds.CIL <- round(exp(fit.summary$Estimate - 1.96*sqrt(fit.summary$Std.error)),2)
      fit.summary$odds.CIU <- round(exp(fit.summary$Estimate + 1.96*sqrt(fit.summary$Std.error)),2)
      results <- add_row(.data = results, "model" = as.character(k), "predictor"= "", 
                         "odds.estimate" = "", "mean" = NA, "lower" = NA, "upper" = NA, "p" = "")
      results <- add_row(.data = results, "model" = "", "predictor"= predictors.all.lab[row.names(fit.summary)[2]], 
                         "odds.estimate" = paste.ci(fit.summary[2,4:6]),
                         "mean" = fit.summary[2,4], "lower" = fit.summary[2,5], 
                         "upper" = fit.summary[2,6], "p" = format.p(fit.summary$P[2]))
      results <- add_row(.data = results, "model" = "", "predictor"= predictors.all.lab[row.names(fit.summary)[3]], 
                         "odds.estimate" = paste.ci(fit.summary[3,4:6]),
                         "mean" = fit.summary[3,4], "lower" = fit.summary[3,5], 
                         "upper" = fit.summary[3,6], "p" = format.p(fit.summary$P[3]))
    }
  }
  
  #model with all four predictors
  k = k+1
  covx <- c(predictors,"stratVar")
  formula = as.formula(paste0("treatment~",paste0(covx,collapse = "+")))
  fit.summary <- data.frame(summary(glm(formula = formula, family = "binomial", data = master.trial))$coefficients)[,c(1,2,4)]
  colnames(fit.summary) <- c("Estimate", "Std.error", "P")
  fit.summary$odds.estimate <- round(exp(fit.summary$Estimate),2)
  fit.summary$odds.CIL <- round(exp(fit.summary$Estimate - 1.96*sqrt(fit.summary$Std.error)),2)
  fit.summary$odds.CIU <- round(exp(fit.summary$Estimate + 1.96*sqrt(fit.summary$Std.error)),2)
  results <- add_row(.data = results, "model" = as.character(k), "predictor"= "", 
                     "odds.estimate" = "", "mean" = NA, "lower" = NA, "upper" = NA, "p" = "")
  results <- add_row(.data = results, "model" = "", "predictor"= predictors.all.lab[row.names(fit.summary)[2]], 
                     "odds.estimate" = paste.ci(fit.summary[2,4:6]),
                     "mean" = fit.summary[2,4], "lower" = fit.summary[2,5], 
                     "upper" = fit.summary[2,6], "p" = format.p(fit.summary$P[2]))
  results <- add_row(.data = results, "model" = "", "predictor"= predictors.all.lab[row.names(fit.summary)[3]], 
                     "odds.estimate" = paste.ci(fit.summary[3,4:6]),
                     "mean" = fit.summary[3,4], "lower" = fit.summary[3,5], 
                     "upper" = fit.summary[3,6], "p" = format.p(fit.summary$P[3]))
  results <- add_row(.data = results, "model" = "", "predictor"= predictors.all.lab[row.names(fit.summary)[4]], 
                     "odds.estimate" = paste.ci(fit.summary[4,4:6]),
                     "mean" = fit.summary[4,4], "lower" = fit.summary[4,5], 
                     "upper" = fit.summary[4,6], "p" = format.p(fit.summary$P[4]))
  
  results <- add_row(.data = results, "model" = "", "predictor"= predictors.all.lab[row.names(fit.summary)[5]], 
                     "odds.estimate" = paste.ci(fit.summary[5,4:6]),
                     "mean" = fit.summary[5,4], "lower" = fit.summary[5,5], 
                     "upper" = fit.summary[5,6], "p" = format.p(fit.summary$P[5]))
  
 

  forestplot.logisticReg( results,
                          xlog = TRUE, ref_line = 1,ticks_digits = 1,   xlim.x = c(0.3, 10),  ticks_at.x = c(0.5, 1, 2, 4, 8), 
                          ci_cols = colors.all[results$predictor],
                          figuresOutDir, paste0(trial,"logisticReg","_ls",".pdf"),
                          header = c("Model","Seq\nFeature","Odds Ratio (95% CI)",
                                             "mean", "lower", "upper","P-value")
                                 )  
}