# Purpose: Primary endpoint analysis; estimation of hazard ratio-based PE by logit of predicted probability that IC80 > 1 ug/ml and by log10 of predicted IC80; 
#          hypothesis testing; plotting
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Juraska and Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on CRAN
# Input:   AMP adjudicated primary endpoints, superlearner-based neutralization predictions
# Output:  PDF files, each containing a single plot
# Author:  Michal Juraska
# Date:    Jan 19, 2021

rm(list=ls(all=TRUE))

datDir <- paste0("t:/vaccine/p704/analysis/sieve/adata")
figDir <- paste0("t:/vaccine/p704/analysis/sieve/figures/sievePH_JuraskaGilbert2013")
tabDir <- paste0("t:/vaccine/p704/analysis/sieve/tables/sievePH_JuraskaGilbert2013")
outDir <- "t:/vaccine/p704/analysis/sieve/code/sievePH_JuraskaGilbert2013/Routput"

# sievePH version 1.0.3 required for using a stratified Cox model
library(sievePH)
library(tidyverse)
source("t:/vaccine/p704/analysis/sieve/code/sievePH_JuraskaGilbert2013/plot.summary.sievePH.R")

logit <- function(p){
  return(log(p/(1-p)))
}

# trial-specific data files created by Craig Magaret, with clinical and neutralization data from Erika Rudnicki
# earlier versions of the trial-pooled data file created by t:\vaccine\p704\analysis\sieve\code\createTrialPooledData.R
source("t:/vaccine/p704/analysis/sieve/code/common.R")

# parscore1.xx: logit predicted probability IC80 >= 1 ug/ml
# parscore2.xx: log base 10 predicted IC80 (ug/ml)
variant <- c("mf", "ms", "ls")
dist <- c(paste0("parscore1.", variant),
          paste0("parscore2.", variant))

# tags that are part of the output PDF file names
trialFileString <- c("704", "703", "704and703")
doseFileString <- paste0("VRC01", c("pooled", "lowdose", "highdose"))
markFileString <- c(paste0("logitPredProbResIC80_", variant), 
                    paste0("log10predIC80_", variant))

# tags that are part of the plot titles
doseTitleString <- c("PE", "PE(10)", "PE(30)")
markType <- rep(c("ProbResIC80", "IC80"), each=3)
# markTitleString <- c("Logit Pred Prob IC80 > 1 ug/ml of Most Freq Variant",
#                      "Log10 Pred IC80 of Most Freq Variant")

# plot labels
VRC01lab <- c("VRC01", "10 mg/kg", "30 mg/kg")
xLabel <- rep(c(expression("Predicted Probability of" ~ IC[80] > 1 ~ mu * "g/ml"), expression("Predicted" ~ IC[80] ~ "(" * mu * "g/ml)")), each=3)

# the output data frame of p-values
pVal <- data.frame(trial=rep(trialFileString, each=3), variant=rep(variant, 3),
                   HRbyProbResIC80unity=NA, HRbyIC80unity=NA, 
                   HRbyProbResIC80constant=NA, HRbyIC80constant=NA)

# used for getting a common xLim for a given mark type across all variants
dataPooled <- read.csv(file.path(datDir, datFile[3]))

# cycle over 704, 703, trial-pooled
for (trial in 1:length(datFile)){
  data <- read.csv(file.path(datDir, datFile[trial])) %>%
    # two 703 ppts with missing sequences have also a missing time-to-event
    filter(!is.na(hiv1fpday)) %>%
    # discard cases with missing IC80 because IC80 is part of the bivariate mark
    filter(!(hiv1event==1 & is.na(gmt80mf))) %>%
    # stratification variable per SAP Section 5
    mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                              protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                              protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                              protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"),
           ic80g1mf=as.numeric(as.numeric(gsub(">|<", "", gmt80mf))>1),
           ic80g1ms=as.numeric(as.numeric(gsub(">|<", "", gmt80ms))>1),
           ic80g1ls=as.numeric(as.numeric(gsub(">|<", "", gmt80ls))>1))
  
  # cycle over VRC01 pooled-dose (T1+T2), low-dose (T1), high-dose (T2)
  for (dose in 1:length(doseFileString)){
    
    removeDose <- c("T2", "T1")
    if (dose==1){
      data1 <- data
    } else {
      # remove the arm not analyzed in this comparison
      data1 <- filter(data, tx!=removeDose[dose - 1])
    }
    
    # for each predicted neutralization sensitivity mark
    for (mark in 1:length(dist)){
      # identify variant of the selected mark
      markVariant <- substring(dist[mark], first=nchar(dist[mark]) - 1)
      
      # get common xLim for a given mark type
      if (mark %in% c(1, 4)){
        if (mark<4){
          xLim <- range(select(dataPooled, all_of(dist[1:3])), na.rm=TRUE)
        } else {
          xLim <- range(select(dataPooled, all_of(dist[4:6])), na.rm=TRUE)
        }
      }
      
      data2 <- subset(data1, select=c("tx_pool", "hiv1fpday", "hiv1event", dist[mark], paste0("ic80g1", markVariant), "stratVar"))
      
      colnames(data2) <- c("tx", "eventTime", "eventInd", "mark", "ic80g1", "stratVar")
      data2$tx <- ifelse(as.character(data2$tx)=="C3", 0, 1)
      
      # convert mark values for non-primary endpoints through tau to NA
      data2$mark <- ifelse(data2$eventInd==0, NA, data2$mark)
      data2$ic80g1 <- ifelse(data2$eventInd==0, NA, data2$ic80g1)
      
      # complete-case analysis, i.e., discard cases with a missing mark
      data2 <- subset(data2, !(eventInd==1 & is.na(mark)))
      
      # fit the mark-specific HR model
      markRng <- range(data2$mark, na.rm=TRUE)
      markGrid <- seq(markRng[1], markRng[2], length.out=100)
      
      # stratified Cox model
      fit <- sievePH(eventTime=data2$eventTime, eventInd=data2$eventInd, mark=select(data2, mark, ic80g1), tx=data2$tx, strata=data2$stratVar)
      
      sfit <- summary(fit, markGrid=matrix(c(rep(markGrid, 2), rep(0:1, each=100)), nrow=200), sieveAlternative="twoSided")
      
      sfitFile <- paste0(trialFileString[trial], "_sievePH_PEby", markFileString[mark], "_and_IC80_", markVariant, "_MITT_", doseFileString[dose], "_placebo.RData")
      save(sfit, file=file.path(outDir, sfitFile))
      
      # TO DO
      # compute p-values only for the comparison of pooled VRC01 vs. placebo
      # if (dose==1){
      #   colName <- paste0("HRby", markType[mark], "unity")
      #   # 1-sided weighted Wald test of {H0: PE(v)=0 for all v} sensitive to PE>0 and PE(v) decreasing in v
      #   pVal[pVal$trial==trialFileString[trial] & pVal$variant==markVariant, colName] <- sfit$pWtWald.HRunity.1sided
      #   
      #   colName <- paste0("HRby", markType[mark], "constant")
      #   # 1-sided Wald test of {H0: PE(v)=PE for all v} against {H1: PE(v) decreasing in v}, i.e., we reject H0 when the Wald test statistic is large
      #   pVal[pVal$trial==trialFileString[trial] & pVal$variant==markVariant, colName] <- sfit$pWald.HRconstant.1sided  
      # }
      
      # if (mark %in% 1:3){ # for predicted probability of IC80 > 1
      #   xTickLab <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
      #   xTickAt <- logit(xTickLab)
      #   sec.xTick <- c(0.03, 0.3, 1, 3, 5, 10)
      # } else if (mark %in% 4:6){  # for predicted IC80
      #   xTickLab <- c(0.62, 1, 2, 3, 5, 10, 20, 30, 50, 80)
      #   xTickAt <- log10(xTickLab)
      # }
      # 
      # #xLab <- bquote(paste("IC", .(substr(markFileString[mark], start=3, stop=4)), " (", mu, "g/mL)"))
      # xLab <- xLabel[mark]
      # 
      
      # 
      # if (dose==1){
      #   p <- sfit$pWald.HRconstant.1sided
      #   fmt.p <- ifelse(p<0.001, "< 0.001", paste0("= ", format(p, digits=2, nsmall=2)))
      #   subtitle <- paste0("1-sided Sieve P ", fmt.p)  
      # } else {
      #   subtitle <- NULL
      # }
      # 
      # 
      
      # pdf(file.path(figDir, paste0(trialFileString[trial], "_sievePH_PEby", markFileString[mark], "_and_IC80_", markVariant, "_MITT_", doseFileString[dose], "_placebo.pdf")), width=0.88*7, height=0.88*6.3)
      data2$tx <- ifelse(data2$tx==0, "Placebo", "VRC01")
      fit.te <- split(sfit$te, sfit$te$ic80g1)
      title <- paste0(switch(markVariant, "mf"="Most Frequent", "ms"="Predicted Most Sensitive", "ls"="Predicted Most Resistant"), " Founder")
      p <- ggplotSievePH(fit.te, data=drop_na(select(data2, markValue=mark, ic80g1, treatment=tx)), top.panel.type="box", 
                         xlab=xLabel[mark], lineColor=c("#9ACD32", "#CD3333"),
                         breaks.x=logit(seq(0.3, 0.9, 0.1)), labels.x=seq(0.3, 0.9, 0.1), title=title)
      ggsave(file.path(figDir, paste0(trialFileString[trial], "_sievePH_PEby", markFileString[mark], "_and_IC80_", markVariant, "_MITT_", doseFileString[dose], "_placebo.pdf")),
             plot=p, width=0.9 * 7, height=0.9 * 8)
      
      # plot(sfit, 
      #      mark=data2$mark, 
      #      tx=data2$tx, 
      #      xlim=xLim,
      #      ylim=c(-0.4, 1),
      #      xtickAt=xTickAt,
      #      xtickLab=xTickLab,
      #      ytickAt=seq(-0.4, 1, by=0.2),
      #      ytickLab=seq(-40, 100, by=20),
      #      xlab=xLab,
      #      ylab="Prevention Efficacy (%)             ",
      #      txLab=c("Placebo", VRC01lab[dose]),
      #      title=title,
      #      subtitle=subtitle)
      # # mtext(paste0(doseTitleString[dose], " by ", markTitleString[mark], ": ", trialTitleString[trial]), side=3, font=2, line=0.5, cex=1.4)
      # dev.off()
    }
  }
}

# save(pVal, file=file.path(tabDir, "704_703_704and703_pValues_PARscores.RData"))



# # test validity of the method's assumptions ------------------------------------------------------------------
# # tests for the mark logit of the predicted probability of the most frequent founder's IC80 > 1 ug/ml
# 
# # PE(v) of pooled VRC01 arms vs. placebo arm
# data1 <- subset(data, select=c("tx_pool","hiv1fpday","hiv1event","parscore1.mf"))
# 
# colnames(data1) <- c("tx", "eventTime", "eventInd", "mark")
# data1$tx <- ifelse(as.character(data1$tx)=="C3", 0, 1)
# 
# # complete-case analysis, i.e., discard cases with a missing mark
# data1 <- subset(data1, !(eventInd==1 & is.na(mark)))
# 
# # test validity of the assumption that T and V are conditionally independent given Z
# testIndepTimeMark(subset(data1, tx==0, select=c("eventTime", "eventInd", "mark")), iter=1000)
# 0.295
# testIndepTimeMark(subset(data1, tx==1, select=c("eventTime", "eventInd", "mark")), iter=1000)
# 0.419
# # conclusion: we do not reject validity of the conditional independence assumption
# 
# # test validity of the specified mark density ratio model
# testDensRatioGOF(data1$eventInd, data1$mark, data1$tx)
# 0.608
# # conclusion: we do not reject validity of the mark density ratio model
