# Purpose: Westfall and Young permutation-based multiplicity adjustment for sieve test p-values
#          The SAP states that p-values are calculated for the VRC01 pooled vs. placebo comparison only.
# Method:  Westfall and Young (1993)
#          Juraska and Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on CRAN
#          Lunn and McNeil (1995)
# Input:   AMP adjudicated primary endpoints, AA site scanning sequence features
# Output:  A list with each component being a vector of sieve test p-values for the analyzed marks from a single permutation of the mark variable.
#          Only p-values for the VRC01 pooled vs. placebo comparison are calculated, separately for each trial.
# Author:  Li Li

rm(list=ls(all=TRUE))


repoDir <- "/Users/lili/AMPsieve"
adjPvaluesDir <- file.path(repoDir, "tables", "sievePH_JuraskaGilbert2013")
tablesOutDir <- file.path(repoDir, "tables/DVE")
figuresOutDir <- file.path(repoDir, "figures/DVE/updated")
outputDir <- file.path(repoDir, "code","DVE/Routput")
#sievePH version 1.0.3 required for using a stratified Cox model
library(sievePH)
library(tidyverse)

source(file.path(repoDir, "code/common.R")) #read in a datFile
source(file.path(repoDir, "code/lunnMcneil.R"))


### Read master file with covariates and time to event data for both trials (note the column southAmerica now has the combined trial/location info as follows:  value 0 = 704, not south america, 1 = 704, south america, 2 = 703)
master <-  read.csv( file.path(repoDir,"adata", datFile[3])) 

#obtain the site scanning residues
tier2_AAsiteScan <- master %>% colnames(.) %>% .[grepl("\\.is\\.", .)]%>%.[grepl(".ls",.)]%>%
  plyr::laply(., function(x)paste(strsplit(x, split = "\\.")[[1]][1:4],collapse="."))

#obtain the Hamming distance variables: preselect.all or binding.all
tier2_HammingDist <- c("hdist.zspace.sites.preselect.all", "hdist.zspace.sites.binding.all")

#obtain the epitope.distance
tier2_epitopeDist <- c("epitope.dist.any", "epitope.dist.b", "epitope.dist.c")

tier2_numFounders <- c("num.founders.all", "num.founders.tfl")
# mark variable names in 'datFile'
variant <- c(".mf", ".ms", ".ls")
dist <- c(outer(tier2_AAsiteScan,variant, FUN=paste0))

# number of permutations of the mark variable
nPerm <- 1000



# Compute p-values from data sets with resampled marks --------------------

# cycle over 704, 703, trial-pooled
for (trial in c("703", "704", "704and703")){
  if(trial == "704and703"){
    data = master
  }else{
    data <- master %>% filter(., protocol == paste0("HVTN ", trial) )
  }
  data <- filter(data, !is.na(hiv1fpday))  # two 703 ppts with missing sequences have also a missing time-to-event
  
  # turn 'tx_pool' into VRC01 indicator (format needed by 'sievePH')
  data$tx_pool <- ifelse(as.character(data$tx_pool)=="C3", 0, 1)
  
  # stratification variable per SAP Section 5
  data <- mutate(data, stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                                                          protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                                                          protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                                                          protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))
  
  # apply the variability filter
  # 'distVarFilter' is a subset of 'dist' that includes only marks that pass the variability filter
  distVarFilter <- sapply(dist, function(dist1){
    counts <- as.vector(table(data[data$hiv1event==1, dist1], useNA="no"))
    if (length(counts)==2 && all(counts>=6)){
      return(dist1)
    } else {
      return(NA)
    }
    
  })
  distVarFilter <- na.omit(distVarFilter)
  hammingDistVar <-outer(tier2_HammingDist,variant, FUN=paste0)
  epitopeDistVar <- paste0(ifelse(trial == "704", "epitope.dist.b", 
                                  ifelse(trial == "703", "epitope.dist.c","epitope.dist.any")), variant)
  allVar <- c(distVarFilter, hammingDistVar, epitopeDistVar, tier2_numFounders)
  # get the p-values for individual permutations of the observed binary marks
  # the first vector in 'pvals' stores unadjusted p-values based on original data
  pvals <- lapply(1:(nPerm + 1), function(seed){
    set.seed(seed)
    
    # permute the marks observed or missing in cases
    idx <- sample(1:sum(data$hiv1event))
    
    # 'pvals1' is a vector of p-values (one for each mark variable) for the single permutation
    pvals1 <- sapply(1:length(allVar), function(i){
      data1 <- subset(data, select=c("tx_pool", "hiv1fpday", "hiv1event", allVar[i], "stratVar"))
      colnames(data1) <- c("tx", "eventTime", "eventInd", "mark", "stratVar")
      
      # convert mark values for non-primary endpoints through tau to NA
      data1$mark <- ifelse(data1$eventInd==0, NA, data1$mark)
      
      # apply the permutation
      if (seed>1){
        data1$mark[data1$eventInd==1] <- data1$mark[data1$eventInd==1][idx]  
      }
      
      # complete-case analysis, i.e., discard cases with a missing mark
      data1 <- subset(data1, !(eventInd==1 & is.na(mark)))
      if(allVar[i] %in% c(distVarFilter, tier2_numFounders )){
        if(allVar[i] %in% tier2_numFounders){
          #Dichotomize  number of founders
          data1$mark <- 1*(data1$mark >1)
        }
        sfit <- with(data1, lunnMcneilTestS(eventTime, eventInd, mark + 1, tx, stratVar=stratVar))
        return(sfit$coef[3, 5])
      }else{
        markRng <- range(data1$mark, na.rm=TRUE)
        markGrid <- seq(markRng[1], markRng[2], length.out=200)
        
        if (trial =="704" | trial =="703"){  # unstratified Cox model
          fit <- with(data1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx)) 
        } else {         # stratified Cox model
          fit <- with(data1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx, strata=stratVar))
        }
        
        sfit <- summary(fit, markGrid=markGrid, sieveAlternative="oneSided")
        
        # 1-sided Wald test of {H0: PE(v) constant for all v}
        return(sfit$pWald.HRconstant.1sided)
      }
      
    })
    
    names(pvals1) <- allVar 
    return(pvals1)
  })
  
  
  
  pvals<- do.call(rbind, pvals)
  save(pvals, file=file.path(outputDir, paste0(trial, "_WestfallYoungPermPvalues_tier2.RData")))
  
}


# Apply Westfall and Young (1993) to obtain adjusted p-values -------------

# a modification of kyotil::p.adj.perm()
source(file.path(repoDir, "code/p.adj.perm2.R")) #read in a datFile
# tags that are part of the output file names
trialFileString <- c("704", "703", "704and703")

for (trial in 1:length(datFile)){
  if (file.exists(file.path(outputDir, paste0(trialFileString[trial], "_WestfallYoungPermPvalues_tier2.RData")))){
    # a matrix named 'pvals'
    # the first row contains unadjusted p-values based on original data
    load(file=file.path(outputDir, paste0(trialFileString[trial], "_WestfallYoungPermPvalues_tier2.RData")))
    
    # perform multiplicity adjustment separately for each founder variant
    for (vt in c("mf","ms","ls")){
      pvals_quantitative <- pvals[,colnames(pvals)%in% c(paste0(tier2_HammingDist,".",vt), paste0(tier2_epitopeDist,".",vt) )]
      pvals_binary <- pvals[,colnames(pvals)%in% paste0(tier2_AAsiteScan,".",vt) ]
     
      pvals_quantitative  <- p.adj.perm2(p.unadj= pvals_quantitative[1, ], p.perms= pvals_quantitative[-1, ], alpha=1)
      pvals_quantitative $featureType = "tier2_quantitative"
      
      pvals.adj_binary <- p.adj.perm2(p.unadj= pvals_binary[1, ], p.perms= pvals_binary[-1, ], alpha=1)
      pvals.adj_binary$featureType = "tier2_3"
      
      pvals.adj = rbind(pvals.adj_binary, pvals_quantitative)
      write.csv(pvals.adj, file=file.path(tablesOutDir,"tier2", paste0(trialFileString[trial], "_WestfallYoungAdjPvalues_tier2_", vt, ".csv")), row.names=FALSE)  
    }
    
    pvals_numFounders <- pvals[,colnames(pvals) %in% tier2_numFounders]
    pvals.adj_numFounders <- p.adj.perm2(p.unadj= pvals_numFounders[1, ], p.perms= pvals_numFounders[-1, ], alpha=1)
    write.csv(pvals.adj_numFounders, file=file.path(tablesOutDir,"tier2", paste0(trialFileString[trial], "_WestfallYoungAdjPvalues_tier2_numFounders", ".csv")), row.names=FALSE)  
    
  }
}
