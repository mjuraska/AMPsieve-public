# Purpose: Westfall and Young permutation-based multiplicity adjustment for sieve test p-values
#          The SAP states that p-values are calculated for the VRC01 pooled vs. placebo comparison only.
# Method:  Westfall and Young (1993)
#          Juraska and Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on CRAN
#          Lunn and McNeil (1995)
# Input:   AMP adjudicated primary endpoints, neutralization-associated sequence features specified by Craig Magaret in the sieve SAP
# Output:  A list with each component being a vector of sieve test p-values for the analyzed marks from a single permutation of the mark variable.
#          Only p-values for the VRC01 pooled vs. placebo comparison are calculated, separately for each trial.
# Author:  Michal Juraska

rm(list=ls(all=TRUE))

datDir <- paste0("t:/vaccine/p704/analysis/sieve/adata")
tabDir <- paste0("t:/vaccine/p704/analysis/sieve/tables/sievePH_JuraskaGilbert2013")
outDir <- "t:/vaccine/p704/analysis/sieve/code/sievePH_JuraskaGilbert2013/Routput"

# sievePH version 1.0.3 required for using a stratified Cox model
library(sievePH)
library(tidyverse)
source("t:/vaccine/p704/analysis/sieve/code/lunnMcneil.R")

# data files created by Craig Magaret, with clinical and neutralization data from Erika Rudnicki
# earlier versions of the trial-pooled data file created by t:\vaccine\p704\analysis\sieve\code\createTrialPooledData.R
source("t:/vaccine/p704/analysis/sieve/code/common.R")

# mark variable names in 'datFile'
variant <- c("mf", "ms", "ls")
dist <- c(outer(c("hxb2.60.A.", "hxb2.170.Q.", "hxb2.230.D.", "hxb2.279.N.", "hxb2.280.N.", "hxb2.317.F.",
                  "hxb2.365.S.", "hxb2.429.E.", "hxb2.456.R.", "hxb2.458.G.", "hxb2.459.G.", "hxb2.471.G.", 
                  "hxb2.156.pngs.", "hxb2.229.pngs.", "hxb2.234.pngs.", "hxb2.616.pngs.", "hxb2.824.pngs.",  
                  "length.gp120.", "length.v1v2.", "length.v5.", "num.pngs.gp120.", "num.pngs.v1v2.", 
                  "num.pngs.v5.", "num.cysteine.gp120."),
                variant, FUN=paste0))

# number of permutations of the mark variable
nPerm <- 1000

# tags that are part of the output file names
trialFileString <- c("704", "703", "704and703")


# Compute p-values from data sets with resampled marks --------------------

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

  # turn 'tx_pool' into VRC01 indicator (format needed by 'sievePH')
  data$tx_pool <- ifelse(as.character(data$tx_pool)=="C3", 0, 1)
  
  # dichotomize two integer-valued mark types as follows:
  # num.pngs.v5.xx: (0, 1) vs. (2, 3)
  vars <- grep("num.pngs.v5", dist, value=TRUE)
  for (v in vars){ data[, v] <- as.numeric(data[, v] >= 2) }
  
  # num.cysteine.gp120.xx: 18 vs. (19, 20, 21)
  vars <- grep("num.cysteine.gp120", dist, value=TRUE)
  for (v in vars){ data[, v] <- as.numeric(data[, v] >= 19) }
  
  # apply the variability filter
  # 'distVarFilter' is a subset of 'dist' that includes only marks that pass the variability filter
  distVarFilter <- sapply(dist, function(dist1){
    dist1Type <- substring(dist1, first=1, last=nchar(dist1) - 3)
    
    # all quantitative marks have sufficient variability and this filter doesn't apply
    if (dist1Type %in% c("length.gp120", "length.v1v2", "length.v5", "num.pngs.gp120", "num.pngs.v1v2")){
      return(dist1)
    } else {
      counts <- as.vector(table(data[data$hiv1event==1, dist1], useNA="no"))
      if (length(counts)==2 && all(counts>=6)){
        return(dist1)
      } else {
        return(NA)
      }
    }
  })
  distVarFilter <- na.omit(distVarFilter)
  
  # get the p-values for individual permutations of the observed marks
  # the first vector in 'pvals' stores unadjusted p-values based on original data
  pvals <- lapply(1:(nPerm + 1), function(seed){
    set.seed(seed)
    
    # permute the marks observed or missing in cases
    idx <- sample(1:sum(data$hiv1event))
    
    # 'pvals1' is a vector of p-values (one for each mark variable) for the single permutation
    pvals1 <- sapply(1:length(distVarFilter), function(i){
      data1 <- subset(data, select=c("tx_pool", "hiv1fpday", "hiv1event", distVarFilter[i], "stratVar"))
      colnames(data1) <- c("tx", "eventTime", "eventInd", "mark", "stratVar")
      
      # convert mark values for non-primary endpoints through tau to NA
      data1$mark <- ifelse(data1$eventInd==0, NA, data1$mark)
      
      # apply the permutation
      if (seed>1){
        data1$mark[data1$eventInd==1] <- data1$mark[data1$eventInd==1][idx]  
      }
      
      # complete-case analysis, i.e., discard cases with a missing mark
      data1 <- subset(data1, !(eventInd==1 & is.na(mark)))
      
      # if the mark is quantitative
      distType <- substring(distVarFilter[i], first=1, last=nchar(distVarFilter[i]) - 3)
      if (distType %in% c("length.gp120", "length.v1v2", "length.v5", "num.pngs.gp120", "num.pngs.v1v2")){
        # fit the mark-specific HR model
        markRng <- range(data1$mark, na.rm=TRUE)
        markGrid <- seq(markRng[1], markRng[2], length.out=200)
        
        # stratified Cox model
        fit <- with(data1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx, strata=stratVar))

        sfit <- summary(fit, markGrid=markGrid, sieveAlternative="twoSided")
        
        # 2-sided Wald test of {H0: PE(v) constant for all v}
        return(sfit$pWald.HRconstant.2sided)
      } else {
        sfit <- with(data1, lunnMcneilTestS(eventTime, eventInd, mark + 1, tx, stratVar=stratVar))
        return(sfit$coef[3, 5])
      }
    })
    
    names(pvals1) <- distVarFilter
    return(pvals1)
  })
  
  pvals <- do.call(rbind, pvals)
  save(pvals, file=file.path(tabDir, paste0(trialFileString[trial], "_WestfallYoungPermPvalues_neutAssocSeqFeatures.RData")))
}


# Apply Westfall and Young (1993) to obtain adjusted p-values -------------

# a modification of kyotil::p.adj.perm()
source("t:/vaccine/p704/analysis/sieve/code/p.adj.perm2.R")

for (trial in 1:length(datFile)){
  if (file.exists(file.path(tabDir, paste0(trialFileString[trial], "_WestfallYoungPermPvalues_neutAssocSeqFeatures.RData")))){
    # a matrix named 'pvals'
    # the first row contains unadjusted p-values based on original data
    load(file=file.path(tabDir, paste0(trialFileString[trial], "_WestfallYoungPermPvalues_neutAssocSeqFeatures.RData")))
    
    # perform multiplicity adjustment separately for each founder variant
    for (vt in variant){
      vtCols <- grep(vt, colnames(pvals))
      pvals.adj <- p.adj.perm2(p.unadj=pvals[1, vtCols], p.perms=pvals[-1, vtCols], alpha=1)
      write.csv(pvals.adj, file=file.path(tabDir, paste0(trialFileString[trial], "_WestfallYoungAdjPvalues_neutAssocSeqFeatures_", vt, ".csv")), row.names=FALSE)  
    }
  }
}
