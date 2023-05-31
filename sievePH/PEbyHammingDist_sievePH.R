# Purpose: Primary endpoint analysis; estimation of hazard ratio-based PE by Hamming distance; 
#          hypothesis testing; plotting
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Juraska and Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on CRAN
# Input:   AMP adjudicated primary endpoints, Hamming distances calculated by Craig Magaret
# Output:  PDF files, each containing a single plot
# Author:  Michal Juraska
# Date:    August 29, 2022

rm(list=ls(all=TRUE))

datDir <- paste0("t:/vaccine/p704/analysis/sieve/adata")
figDir <- paste0("t:/vaccine/p704/analysis/sieve/figures/sievePH_JuraskaGilbert2013")
tabDir <- paste0("t:/vaccine/p704/analysis/sieve/tables/sievePH_JuraskaGilbert2013")
tabDir2 <- paste0("t:/vaccine/p704/analysis/sieve/tables/DVE/tier2")
outDir <- "t:/vaccine/p704/analysis/sieve/code/sievePH_JuraskaGilbert2013/Routput"

# sievePH version 1.0.3 required for using a stratified Cox model
library(sievePH)
library(tidyverse)
source("t:/vaccine/p704/analysis/sieve/code/sievePH_JuraskaGilbert2013/plot.summary.sievePH.R")
source("h:/SCHARP/sieveMethods/sievePH/R/ggplot.summary.sievePH.R")

# trial-specific data files created by Craig Magaret, with clinical and neutralization data from Erika Rudnicki
# earlier versions of the trial-pooled data file created by t:\vaccine\p704\analysis\sieve\code\createTrialPooledData.R
source("t:/vaccine/p704/analysis/sieve/code/common.R")

variant <- c("mf", "ms", "ls")
dist <- c(outer(c("hdist.zspace.sites.preselect.all.", "hdist.zspace.sites.binding.all."), variant, FUN=paste0))

# tags that are part of the output PDF file names
trialFileString <- c("704", "703", "704and703")
doseFileString <- paste0("VRC01", c("pooled", "lowdose", "highdose"))
markType <- c("HammingDistPreselectAll", "HammingDistBindingAll")
markFileString <- c(outer(paste0(markType, "_"), variant, FUN=paste0))

# tags that are part of the plot titles
trialTitleString <- c("704/085", "703/081", "703+704")
doseTitleString <- c("PE", "PE(10)", "PE(30)")
markTitleString <- dist

# plot labels
VRC01lab <- c("VRC01", "10 mg/kg", "30 mg/kg")

# the output data frame of p-values
pVal <- data.frame(trial=rep(trialFileString, each=3), variant=rep(variant, 3))
pVal <- cbind(pVal, as.data.frame(matrix(NA, nrow=9, ncol=2 * length(markType))))
colnames(pVal)[-(1:2)] <- c(paste0("HRby", markType, "Unity"), paste0("HRby", markType, "Constant"))

dataPooled <- read.csv(file.path(datDir, datFile[3]))  

# cycle over 704, 703, trial-pooled
for (trial in 1:length(datFile)){
  data <- read.csv(file.path(datDir, datFile[trial])) %>%
    # two 703 ppts with missing sequences have also a missing time-to-event
    filter(!is.na(hiv1fpday)) %>%
    # stratification variable per SAP Section 5
    mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                              protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                              protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                              protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))
  
  # this cohort exclusion criterion is outdated
  # discard cases with missing neutralization sensitivity data per Section 1.1 in the AMP Genotypic Sieve SAP (analysis cohort definitions)
  # filter(!(hiv1event==1 & is.na(gmt80mf)))
  
  # cycle over VRC01 pooled-dose (T1+T2), low-dose (T1), high-dose (T2)
  for (dose in 1:length(doseFileString)){
    
    removeDose <- c("T2", "T1")
    if (dose==1){
      data1 <- data
    } else {
      # remove the arm not analyzed in this comparison
      data1 <- filter(data, tx!=removeDose[dose - 1])
    }
    
    # for each sequence feature
    for (mark in 1:length(dist)){
      # identify variant of the selected mark
      markVariant <- substring(dist[mark], first=nchar(dist[mark]) - 1)
      
      # get common xLim for a given mark type
      distType <- substring(dist[mark], first=1, last=nchar(dist[mark]) - 3)
      xLim <- range(select(dataPooled, starts_with(distType)), na.rm=TRUE)
      
      data2 <- subset(data1, select=c("tx_pool", "hiv1fpday", "hiv1event", dist[mark], "stratVar"))
      
      colnames(data2) <- c("tx", "eventTime", "eventInd", "mark", "stratVar")
      data2$tx <- ifelse(as.character(data2$tx)=="C3", 0, 1)
      
      # convert mark values for non-primary endpoints through tau to NA
      data2$mark <- ifelse(data2$eventInd==0, NA, data2$mark)
      
      # complete-case analysis, i.e., discard cases with a missing mark
      data2 <- subset(data2, !(eventInd==1 & is.na(mark)))
      
      # fit the mark-specific HR model
      markRng <- range(data2$mark, na.rm=TRUE)
      markGrid <- seq(markRng[1], markRng[2], length.out=200)
      
      # stratified Cox model
      fit <- sievePH(eventTime=data2$eventTime, eventInd=data2$eventInd, mark=data2$mark, tx=data2$tx, strata=data2$stratVar)
      
      sfit <- summary(fit, markGrid=markGrid, sieveAlternative="oneSided")
      
      sfitFile <- paste0(trialFileString[trial], "_sievePH_PEby", markFileString[mark], "_MITT_", doseFileString[dose], "_placebo.RData")
      save(sfit, file=file.path(outDir, sfitFile))
      
      # compute p-values only for the comparison of pooled VRC01 vs. placebo
      if (dose==1){
        markType1 <- substring(markFileString[mark], first=1, last=nchar(markFileString[mark]) - 2)
        colName <- paste0("HRby", markType1, "Unity")
        # 2-sided Wald test of {H0: PE(v)=0 for all v}
        pVal[pVal$trial==trialFileString[trial] & pVal$variant==markVariant, colName] <- sfit$pWtWald.HRunity.1sided
        
        colName <- paste0("HRby", markType1, "Constant")
        # 2-sided Wald test of {H0: PE(v) constant for all v}
        pVal[pVal$trial==trialFileString[trial] & pVal$variant==markVariant, colName] <- sfit$pWald.HRconstant.1sided  
      }
      
      #xLab <- bquote(paste("IC", .(substr(markFileString[mark], start=3, stop=4)), " (", mu, "g/mL)"))
      xLab <- switch(distType, 
                     "hdist.zspace.sites.preselect.all"="PC-Weighted Hamming Distance in\n27 Positions Predictive of Neutralization",
                     "hdist.zspace.sites.binding.all"="PC-Weighted Hamming Distance in\n50 VRC01 or CD4 Binding Positions")
      
      title <- paste0(switch(markVariant, "mf"="Most Frequent", "ms"="Predicted Most Sensitive", "ls"="Predicted Most Resistant"), " Founder")
      
      # annotate the figure with sieve test unadjusted and adjusted p-values
      if (dose==1){
        pval.filename <- paste0(trialFileString[trial], "_WestfallYoungAdjPvalues_tier2_", markVariant, ".csv")
        pval.filepath <- file.path(tabDir2, pval.filename)
        if (file.exists(pval.filepath)){
          p.df <- read_csv(pval.filepath, show_col_types=FALSE)
          p <- as.numeric(p.df[p.df$mark==dist[mark], c("p.unadj", "p.FWER", "p.FDR")])
          # p <- sfit$pWald.HRconstant.1sided  # this line used when adjusted p-values were unavailable
          fmt.p <-  sapply(p, function(x){ ifelse(x<0.001, "< 0.001", paste0("= ", format(x, digits=2, nsmall=2))) })
          subtitle <- paste0("One-Sided Unadjusted Sieve P ", fmt.p[1], "\nFWER P ", fmt.p[2], ", Q ", fmt.p[3])
        }
        
        plotHeights <- c(0.38, 0.62)
        ggsave.width <- 0.7 * 7
        ggsave.height <- 0.7 * 6.3
        
      } else {
        subtitle <- NULL
        
        plotHeights <- c(0.32, 0.68)
        ggsave.width <- 0.7 * 7.3
        ggsave.height <- 0.7 * 6
      }
      
      # pdf(file.path(figDir, paste0(trialFileString[trial], "_sievePH_PEby", markFileString[mark], "_MITT_", doseFileString[dose], "_placebo.pdf")), 
      #     width=0.88*7, height=0.88*6.3)
      p <- ggplot(sfit, 
                  mark=data2$mark, 
                  tx=data2$tx, 
                  xlim=xLim,
                  ylim=c(-0.4, 1),
                  xtickAt=NULL,
                  xtickLab=NULL,
                  ytickAt=seq(-0.4, 1, by=0.2),
                  ytickLab=seq(-40, 100, by=20),
                  xlab=xLab,
                  ylab="Prevention Efficacy (%)",
                  axisLabSize=15.5,
                  legendLabSize=12,
                  txLab=c("Placebo", VRC01lab[dose]),
                  jitterFactor=0.1,
                  title=title,
                  subtitle=subtitle,
                  subtitleSize=13,
                  estLineSize=1.8,
                  ciLineSize=1.4,
                  pointSize=2.1,
                  plotHeights=plotHeights)
      ggsave(file.path(figDir, paste0(trialFileString[trial], "_sievePH_PEby", markFileString[mark], "_MITT_", doseFileString[dose], "_placebo.pdf")),
             plot=p, width=ggsave.width, height=ggsave.height)
      
      # mtext(paste0(doseTitleString[dose], " by ", markTitleString[mark], ": ", trialTitleString[trial]), side=3, font=2, line=0.5, cex=1.4)
      # dev.off()
    }
  }
}

save(pVal, file=file.path(tabDir, "704_703_704and703_pValues_hammingDist.RData"))


# test validity of the method's assumptions ------------------------------------------------------------------
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
