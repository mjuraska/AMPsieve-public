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
source("h:/SCHARP/sieveMethods/sievePH/R/ggplot.summary.sievePH.R")

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
    # stratification variable per SAP Section 5
    mutate(stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                              protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                              protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                              protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))
  
  # this exclusion criterion is outdated
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
        colName <- paste0("HRby", markType[mark], "unity")
        # 1-sided weighted Wald test of {H0: PE(v)=0 for all v} sensitive to PE>0 and PE(v) decreasing in v
        pVal[pVal$trial==trialFileString[trial] & pVal$variant==markVariant, colName] <- sfit$pWtWald.HRunity.1sided
        
        colName <- paste0("HRby", markType[mark], "constant")
        # 1-sided Wald test of {H0: PE(v)=PE for all v} against {H1: PE(v) decreasing in v}, i.e., we reject H0 when the Wald test statistic is large
        pVal[pVal$trial==trialFileString[trial] & pVal$variant==markVariant, colName] <- sfit$pWald.HRconstant.1sided  
      }
      
      if (mark %in% 1:3){ # for predicted probability of IC80 > 1
        xTickLab <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
        xTickAt <- logit(xTickLab)
        sec.xTick <- c(0.03, 0.3, 1, 3, 5, 10)
      } else if (mark %in% 4:6){  # for predicted IC80
        xTickLab <- c(0.62, 1, 2, 3, 5, 10, 20, 30, 50, 80)
        xTickAt <- log10(xTickLab)
      }

      #xLab <- bquote(paste("IC", .(substr(markFileString[mark], start=3, stop=4)), " (", mu, "g/mL)"))
      xLab <- xLabel[mark]
      
      title <- paste0(switch(markVariant, "mf"="Most Frequent", "ms"="Predicted Most Sensitive", "ls"="Predicted Most Resistant"), " Founder")
      
      if (dose==1){
        p <- sfit$pWald.HRconstant.1sided
        fmt.p <- ifelse(p<0.001, "< 0.001", paste0("= ", format(p, digits=2, nsmall=2)))
        subtitle <- paste0("One-Sided Sieve P ", fmt.p)  
        
        plotHeights <- c(0.38, 0.62)
        ggsave.width <- 0.7 * 7
        ggsave.height <- 0.7 * 6.3
        
      } else {
        subtitle <- NULL
        
        plotHeights <- c(0.33, 0.67)
        ggsave.width <- 0.7 * 7.2
        ggsave.height <- 0.7 * 6.1
      }
      

      # pdf(file.path(figDir, paste0(trialFileString[trial], "_sievePH_PEby", markFileString[mark], "_MITT_", doseFileString[dose], "_placebo.pdf")), width=0.88*7, height=0.88*6.3)
      p <- ggplot(sfit, 
                  mark=data2$mark, 
                  tx=data2$tx, 
                  xlim=xLim,
                  ylim=c(-0.4, 1),
                  xtickAt=xTickAt,
                  xtickLab=xTickLab,
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
                  subtitleSize=16,
                  estLineSize=1.8,
                  ciLineSize=1.4,
                  pointSize=2.1,
                  plotHeights=plotHeights)
      ggsave(file.path(figDir, paste0(trialFileString[trial], "_sievePH_PEby", markFileString[mark], "_MITT_", doseFileString[dose], "_placebo.pdf")),
             plot=p, width=ggsave.width, height=ggsave.height)
      # dev.off()
      
      # plotting with overlaid PE estimates from neutralization sieve analysis
      if (dose==1 && mark<4){
        
        # a data frame named 'out' with PE estimates on the grid for both marks
        load(file.path(outDir, paste0(trialFileString[trial], "_sievePH_compPE_logitPredProbResIC80_IC80_", markVariant, "_MITT_VRC01pooled_placebo.RData")), verbose=TRUE)
        
        # this transformation is needed in 'data1' to get tickmark labels of the secondary x-axis
        # transform to the scale on which the mark was analyzed
        gmt80name <- paste0("gmt80", markVariant)
        data1$gmt80 <- as.character(data1[, gmt80name])
        data1$gmt80 <- ifelse(data1$gmt80==">100", "100", data1$gmt80)
        data1$gmt80 <- as.numeric(data1$gmt80)
        # right-censor at 10
        data1$gmt80 <- pmin(10, data1$gmt80)
        # analyze the gmt80mf on the log10 scale
        # data1$gmt80mf <- log10(data1$gmt80mf)
        
        data1 <- filter(data1, !(hiv1event==1 & (is.na(gmt80) | is.na(parscore1.mf))))
        
        # the next two lines delete all but one value censored at 10 to stretch PE by measured IC80 across the x-axis
        # idx <- which(data1$gmt80mf==10)
        # data1$gmt80mf[idx[-1]] <- NA
        
        
        # tickmark labels for the secondary x-axis ('sec.xTickLab')
        # sec.xTickLab <- numeric(length(xTickLab) - ifelse(trial==1, 1, 0))
        # for (i in 1:length(sec.xTickLab)){
        #   prob <- mean(data1$parscore1.mf <= xTickAt[i], na.rm=TRUE)
        #   sec.xTickLab[i] <- quantile(data1$gmt80mf, probs=prob, na.rm=TRUE)
        # }
        # sec.xTickLab <- round(sec.xTickLab, digits=2)
        xTickGrid <- seq(min(data1[, dist[mark]], na.rm=TRUE), max(data1[, dist[mark]], na.rm=TRUE), length.out=500)
        sec.xTickGrid <- numeric(length(xTickGrid))
        Fn <- ecdf(data1[, dist[mark]])
        for (i in 1:length(xTickGrid)){
          prob <- Fn(xTickGrid[i])
          # quantile(..., type=1) is the inverse of ecdf()
          sec.xTickGrid[i] <- quantile(data1$gmt80, probs=prob, na.rm=TRUE, type=1)
        }
        idx <- sapply(sec.xTick, function(x){ which.min(abs(sec.xTickGrid - x)) })  
        
        
        pdf(file.path(figDir, paste0(trialFileString[trial], "_sievePH_compPE_logitPredProbResIC80_IC80_",  substring(dist[mark], first=nchar(dist[mark]) - 1),
                                            "_MITT_VRC01pooled_placebo.pdf")), width=0.86*7, height=0.86*6.7)
        plot(sfit, 
             mark=data2$mark, 
             tx=data2$tx, 
             xlim=xLim,
             ylim=c(-0.7, 1),
             xtickAt=xTickAt,
             xtickLab=xTickLab,
             ytickAt=seq(-0.4, 1, by=0.2),
             ytickLab=seq(-40, 100, by=20),
             xlab=xLab,
             ylab="Prevention Efficacy (%)             ",
             txLab=c("Placebo", VRC01lab[dose]),
             compTE=out,
             sec.xtickAt=xTickGrid[idx],
             sec.xtickLab=sec.xTick,
             sec.xlab=expression("Measured" ~ IC[80] ~ "(" * mu * "g/ml)            "),
             title=title,
             parMar=c(8, 6, 2, 1))
        # mtext(paste0(doseTitleString[dose], " by ", markTitleString[mark], ": ", trialTitleString[trial]), side=3, font=2, line=0.5, cex=1.4)
        dev.off()
      }
    }
  }
}

save(pVal, file=file.path(tabDir, "704_703_704and703_pValues_PARscores.RData"))



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
