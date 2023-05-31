# Purpose: Estimation of hazard ratio-based PE by tier2 AA position-specific features and Hamming distances; 
# hypothesis testing Method:  Lunn & McNeil; Custom script allowing for stratification of baseline hazards, see dve3() below
# Input:  Master file (.csv) with time to event data, treatment group assignment, region (south america vs. other) and binary marks to be analyzed
# Output: .csv and forest-plots containing results table with estimates of mark-specific PE, hazard rate ratio, and associated confidence intervals   
# Author:  Li Li
# Date:    Sep 19, 2022


rm(list=ls(all=TRUE))
repoDir <- "/Users/lili/AMPsieve"
adjPvaluesDir <- file.path(repoDir, "tables", "DPE/tier2")
tablesOutDir <- file.path(repoDir, "tables/DPE/tier2")
figuresOutDir <- file.path(repoDir, "figures/DPE/tier2")
outputDir <- file.path(repoDir, "code","DPE/Routput")

library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)

source(file.path(repoDir, "code/common.R")) #read in a datFile
source(file.path(repoDir, "code/lunnMcneil.R"))
source(file.path(repoDir,"code/DPE/DPEutils.R"))
source(file.path(repoDir, "code/DPE/forest.R"))
source(file.path(repoDir, "code/p.adj.perm2.R")) #read in a datFile

### Read master file with covariates and time to event data for both trials (note the column southAmerica now has the combined trial/location info as follows:  value 0 = 704, not south america, 1 = 704, south america, 2 = 703)
master <-  read.csv( file.path(repoDir,"adata", datFile[3])) 

#Sieve analyses for site scanning of residues
tier2_AAsiteScan <- master %>% colnames(.) %>% .[grepl("\\.is\\.", .)]%>%.[grepl(".ls",.)]%>%
                    plyr::laply(., function(x)paste(strsplit(x, split = "\\.")[[1]][1:4],collapse="."))


#Obtain results separately for each trial and pooled, and separately for mf, ms, ls
for(sequence in c("mf", "ms", "ls")){
  for(trial in c("703", "704", "704and703")){
    for(dose in c("T1", "T2", "T1+T2")){
      
      master.trial <- trial_dose_data (master, trial, dose)
      features.to.analyze.all <- paste0(tier2_AAsiteScan,".",sequence)
      # apply the variability filter
      features.to.analyze.x <- sapply(features.to.analyze.all, function(dist1){
        counts <- as.numeric(table(master.trial[master.trial$hiv1event==1, dist1], useNA="no"))
        if (length(counts)==2 && all(counts>=6)){
          return(dist1)
        } else {
          return(NA)
        }
        
      })
      features.to.analyze.x <- na.omit(features.to.analyze.x)
     
      
      ### List of columns in output 
      results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                             "DPE.p.value", 
                             paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                             paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;
      
      # We'll store the results here:
      AAsites <- sapply( strsplit(features.to.analyze.x,"\\."), function(.feature){
        paste( .feature[2])
      }) 
      results <- data.frame( feature = features.to.analyze.x, feature.name = AAsites )
      results[ , 3:( length( results.colnames) ) ] <- NA;
      colnames( results ) <- results.colnames
      
      # Begin for loop 
      .results.row.i <- 0;
      for( .feature in features.to.analyze.x ){
        
        .results.row.i <- .results.row.i + 1;
        
        # Be verbose:
        cat( "Computing results for feature", results[ .results.row.i, 1 ], "\n" );
        
        # Grab the column containing feature values from the master.trial table
        .x <- master.trial[ , .feature]
        
        ### Run DVE
        .dPE.results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event, .x, master.trial$txLabel, master.trial$stratVar )
        

        # Insert into results table 
        results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
                                     paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                                     paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.results , dose, .x)
      }
      
     
      #add adjusted p-values
      if(dose == "T1+T2"){
        p.value.file <- read.csv(file.path(adjPvaluesDir, paste0(trial,"_WestfallYoungAdjPvalues_tier2_",sequence,".csv")))
        results <- full_join(results, p.value.file, by = "mark")
        PEtable <- tibble("feature" = character(),"mark" = character(), "Haplotype" = character(), "n" = character(),
                          "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                          "P" = character(), "DiffP" = character(), "FWER" = character(), "FDR" = character())
        for( .feature in features.to.analyze.x ){
          result.f <- filter(results, mark == .feature)
          annotatedMark <- result.f$mark.name
          
          #sig.unadj <- ifelse(result.f$p.unadj<=0.05, 1, 0) 
          #sig.fdr <- ifelse(result.f$p.FDR<=0.2 & result.f$p.unadj <=0.05, 1, 0) 
          #sig.holm <- ifelse(result.f$p.FWER<=0.05, 1, 0) 
          #annotatedMark = ifelse(sig.unadj, paste(annotatedMark,"$^{\\P}$",sep=""),annotatedMark)
          #annotatedMark = ifelse(sig.fdr, paste(annotatedMark,"$^{\\S}$",sep=""), annotatedMark)
          #annotatedMark = ifelse(sig.holm, paste(annotatedMark,"$^{\\dag}$",sep=""), annotatedMark)
          
          feature.x <- "tier2_3"
          
          PEtable <- add_row(.data = PEtable, "feature" = feature.x , "mark" = annotatedMark, "Haplotype" = "", "n" = "",
                               "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                               "DiffP" = format.p(result.f$p.unadj), "FWER" = format.p(result.f$p.FWER), "FDR" = format.p(result.f$p.FDR))
            
                     
          
          Haplotype <- table.seqFeatLabel.tier2(paste0(strsplit(.feature, split = "\\.")[[1]][1:4], collapse="."))
          
          
          PEtable <- add_row(.data = PEtable,"feature" = feature.x , "mark" = "", "Haplotype" = Haplotype$Haplotype1, 
                             "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo) , 
                             "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                             "mean" = as.numeric(result.f$PE.type.1.estimate),
                             "lower" = as.numeric(result.f$PE.type.1.CI.low), #to be able to plot -Inf
                             "upper" = as.numeric(result.f$PE.type.1.CI.high),
                             "P" = format.p(result.f$PE.type.1.p.value),
                             "DiffP" = "", "FWER" = "", "FDR" = "")
          PEtable <- add_row(.data = PEtable,"feature" = feature.x , "mark" = "", "Haplotype" = Haplotype$Haplotype0, 
                             "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                             "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                             "mean" = as.numeric(result.f$PE.type.0.estimate),
                             "lower" = as.numeric(result.f$PE.type.0.CI.low),
                             "upper" = as.numeric(result.f$PE.type.0.CI.high),
                             "P" = format.p(result.f$PE.type.0.p.value),
                             "DiffP" = "", "FWER" = "", "FDR" = "")
          
        }
      
        write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_AAsiteScan_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,".csv")))
        
        #forest plots
        PEtable$lower[PEtable$lower == -Inf] = -110
        PEtable$feature <- NULL
        maxMarkPerPage = 31
        nMarkTotal = length(PEtable$mark)/3
        nrowsPerMark = 3
        if(nMarkTotal <= maxMarkPerPage){
          page.list <- list(page1 = 1: (nMarkTotal*nrowsPerMark))
        }else if(nMarkTotal <= 2*maxMarkPerPage){
          
          page.list <-list( page1 = 1:(floor(nMarkTotal/2)*nrowsPerMark),  page2 = (floor(nMarkTotal/2)*nrowsPerMark+1):(nMarkTotal*nrowsPerMark))
        }else{
          page.list <- list( page1 = 1:(floor(nMarkTotal/3)*nrowsPerMark),  page2 = (floor(nMarkTotal/3)*nrowsPerMark+1):(floor(nMarkTotal/3)*2*nrowsPerMark), 
                             page3 = (floor(nMarkTotal/3)*2*nrowsPerMark+1):(nMarkTotal*nrowsPerMark))
        }
        
        for(page in 1:length(page.list)){
          PEtable.rows <- PEtable[min(page.list[[page]]): max(page.list[[page]]), ]
          PEtable.rows <- filter(PEtable.rows, rowSums(is.na(PEtable.rows)) != ncol(PEtable.rows))
          
          if(length(PEtable.rows$mark)>=1){
            PEtable.forestplot.withSieveT(PEtable.rows,
              xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
              figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_AAsiteScan_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,"_",page,".pdf"))
            
          }
          
        }
      }else{
        
        PEtable <- tibble("feature" = character(),"mark" = character(), "Haplotype" = character(), "n" = character(),
                          "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                          "P" = character())
        for( .feature in features.to.analyze.x ){
          result.f <- filter(results, mark == .feature)
          Haplotype <- table.seqFeatLabel.tier2(paste0(strsplit(.feature, split = "\\.")[[1]][1:4], collapse="."))
          feature.x <- ifelse(result.f$mark %in% features.to.analyze.all, "tier2_3", "glycosite230")
          
          PEtable <- add_row(.data = PEtable,"feature" = feature.x, "mark" = result.f$mark.name, "Haplotype" = Haplotype$Haplotype1, 
                             "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo),
                             "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                             "mean" = as.numeric(result.f$PE.type.1.estimate),
                             "lower" = as.numeric(result.f$PE.type.1.CI.low),
                             "upper" = as.numeric(result.f$PE.type.1.CI.high),
                             "P" = format.p(result.f$PE.type.1.p.value))
          PEtable <- add_row(.data = PEtable,"feature" = feature.x, "mark" = "", "Haplotype" = Haplotype$Haplotype0, 
                             "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                             "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                             "mean" = as.numeric(result.f$PE.type.0.estimate),
                             "lower" = as.numeric(result.f$PE.type.0.CI.low),
                             "upper" = as.numeric(result.f$PE.type.0.CI.high),
                             "P" = format.p(result.f$PE.type.0.p.value))
          
        }
       
        write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_AAsiteScan_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,".csv")))
        
        PEtable$lower[PEtable$lower == -Inf] = -110
        PEtable$feature <- NULL
        maxMarkPerPage = 34
        nMarkTotal = length(PEtable$mark)/2
        nrowsPerMark = 2
        if(nMarkTotal <= maxMarkPerPage){
          page.list <- list(page1 = 1: (nMarkTotal*nrowsPerMark))
        }else if(nMarkTotal <= 2*maxMarkPerPage){
          
          page.list <-list( page1 = 1:(floor(nMarkTotal/2)*nrowsPerMark),  page2 = (floor(nMarkTotal/2)*nrowsPerMark+1):(nMarkTotal*nrowsPerMark))
        }else{
          page.list <- list( page1 = 1:(floor(nMarkTotal/3)*nrowsPerMark),  page2 = (floor(nMarkTotal/3)*nrowsPerMark+1):(floor(nMarkTotal/3)*2*nrowsPerMark), 
                             page3 = (floor(nMarkTotal/3)*2*nrowsPerMark+1):(nMarkTotal*nrowsPerMark))
        }
        
        for(page in 1:length(page.list)){
          PEtable.rows <- PEtable[min(page.list[[page]]): max(page.list[[page]]), ]
          PEtable.rows <- filter(PEtable.rows, rowSums(is.na(PEtable.rows)) != ncol(PEtable.rows))
          if(length(PEtable.rows$mark)>=1){
          PEtable.forestplot( PEtable.rows,
                                        xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                                        figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_AAsiteScan_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,"_",page,".pdf"))
          }
        }
      }
      
      
    }
    
    
  }
}

#P-value tables
p.value.files <- c("703_WestfallYoungAdjPvalues_tier2_ls.csv", "703_WestfallYoungAdjPvalues_tier2_mf.csv", "703_WestfallYoungAdjPvalues_tier2_ms.csv",
                   "704_WestfallYoungAdjPvalues_tier2_ls.csv", "704_WestfallYoungAdjPvalues_tier2_mf.csv", "704_WestfallYoungAdjPvalues_tier2_ms.csv",
                   "704and703_WestfallYoungAdjPvalues_tier2_ls.csv", "704and703_WestfallYoungAdjPvalues_tier2_mf.csv", "704and703_WestfallYoungAdjPvalues_tier2_ms.csv")
for(i in 1:length(p.value.files)){
  p.value.table <- read.csv(file.path(adjPvaluesDir, p.value.files[i]))
  
  #For binary marks, unadjusted p value threshold is 0.05
  p.value.table.annotated <- data.frame()
  for(featureType.x in c("tier2_3", "tier2_quantitative")){
    p.value.table1 <- filter(p.value.table, featureType == featureType.x)
    annotatedMark <- p.value.table1$mark
    if(featureType.x == "tier2_3"){
      cutoff.unadj = 0.05; cutoff.FDR = 0.2; cutoff.FWER = 0.05
    }else{
      cutoff.unadj = 0.025; cutoff.FDR = 0.1; cutoff.FWER = 0.025
    }
    sig.unadj <- ifelse(p.value.table1$p.unadj<=cutoff.unadj, 1, 0) 
    sig.fdr <- ifelse(p.value.table1$p.FDR<=cutoff.FDR & p.value.table1$p.unadj <=cutoff.unadj, 1, 0) 
    sig.holm <- ifelse(p.value.table1$p.FWER<=cutoff.FWER, 1, 0) 
    annotatedMark = ifelse(sig.unadj, paste(annotatedMark,"$^{\\ast}$",sep=""),annotatedMark)
    annotatedMark = ifelse(sig.fdr, paste(annotatedMark,"$^{\\S}$",sep=""), annotatedMark)
    annotatedMark = ifelse(sig.holm, paste(annotatedMark,"$^{\\dag}$",sep=""), annotatedMark)
    p.value.table1$mark <- annotatedMark
    p.value.table1$p.unadj <- format.p(p.value.table1$p.unadj,3)
    p.value.table1$p.FDR <- format.p(p.value.table1$p.FDR,3)
    p.value.table1$p.FWER <- format.p(p.value.table1$p.FWER,3)
    p.value.table.annotated <- rbind(p.value.table.annotated,p.value.table1 )
  }
 write.csv(p.value.table.annotated, file.path(adjPvaluesDir, paste0("annotated_",p.value.files[i]))) 
}




#Founder type analyses
tier2_founderType <- c("all.founders.sensitive", "at.least.one.founder.resistant", "all.founders.resistant", "all.founders.all.sens.resist")
for(trial in c("703", "704", "704and703")){
  for(dose in c("T1", "T2", "T1+T2")){
    master.trial <- trial_dose_data (master, trial, dose)
    
    ### List of columns in output 
    results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                           "DPE.p.value", 
                           paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                           paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;
    
    results <- data.frame( feature = tier2_founderType, feature.name = tier2_founderType)
    results[ , 3:( length( results.colnames) ) ] <- NA;
    colnames( results ) <- results.colnames
    .results.row.i = 0
    for(.feature in tier2_founderType){
      .results.row.i <- .results.row.i + 1;
     
      # Be verbose:
      cat( "Computing results for feature", results[ .results.row.i, 1 ], "\n" );
      
      # Grab the column containing feature values from the master.trial table
      .x <- master.trial[ , .feature]
      
      #if mark is all.founders.all.sens.resist, censor all events with NA mark
      if(.feature == "all.founders.all.sens.resist"){
        master.trial$hiv1event[master.trial$hiv1event==1&is.na(.x)] <- 0
      }
      
      ### Run DVE
      .dPE.results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event, .x, master.trial$txLabel, master.trial$stratVar )
      
      
      # Insert into results table 
      results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
                                   paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                                   paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.results , dose, .x)
    }
    
    if(dose == "T1+T2"){
       PEtable <- tibble("mark" = character(), "n" = character(),
                        "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                        "P" = character(), "DiffP" = character())
       
       result1 <- filter(results, mark == "all.founders.sensitive")
       result2 <- filter(results, mark == "all.founders.resistant")
       result3 <- filter(results, mark == "all.founders.all.sens.resist")
       
       PEtable <- add_row(.data = PEtable ,"mark" = "", "n" = "",
                          "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                          "DiffP" = format.p(result1$DPE.p.value)) 
       PEtable <- add_row(.data = PEtable , "mark" = "All founders predicted sensitive", 
                          "n" = paste0(result1$.n.events.type.1.trt, " vs. ", result1$.n.events.type.1.placebo) , 
                          "PE" = paste0(format.PE(result1$PE.type.1.estimate)," (", format.PE(result1$PE.type.1.CI.low),", ",format.PE(result1$PE.type.1.CI.high),")" ),
                          "mean" = as.numeric(result1$PE.type.1.estimate),
                          "lower" = as.numeric(result1$PE.type.1.CI.low), #to be able to plot -Inf
                          "upper" = as.numeric(result1$PE.type.1.CI.high),
                          "P" = format.p(result1$PE.type.1.p.value),
                          "DiffP" = "")
       
       PEtable <- add_row(.data = PEtable , "mark" = "At least 1 founder predicted resistant",
                          "n" = paste0(result1$.n.events.type.0.trt, " vs. ", result1$.n.events.type.0.placebo) ,
                          "PE" = paste0(format.PE(result1$PE.type.0.estimate)," (", format.PE(result1$PE.type.0.CI.low),", ",format.PE(result1$PE.type.0.CI.high),")" ),
                          "mean" = as.numeric(result1$PE.type.0.estimate),
                          "lower" = as.numeric(result1$PE.type.0.CI.low),
                          "upper" = as.numeric(result1$PE.type.0.CI.high),
                          "P" = format.p(result1$PE.type.0.p.value),
                          "DiffP" = "")
       
       PEtable <- add_row(.data = PEtable ,"mark" = "", "n" = "",
                          "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                          "DiffP" = format.p(result3$DPE.p.value)) 
       PEtable <- add_row(.data = PEtable , "mark" = "All founders predicted sensitive", 
                          "n" = paste0(result1$.n.events.type.1.trt, " vs. ", result1$.n.events.type.1.placebo) , 
                          "PE" = paste0(format.PE(result1$PE.type.1.estimate)," (", format.PE(result1$PE.type.1.CI.low),", ",format.PE(result1$PE.type.1.CI.high),")" ),
                          "mean" = as.numeric(result1$PE.type.1.estimate),
                          "lower" = as.numeric(result1$PE.type.1.CI.low), #to be able to plot -Inf
                          "upper" = as.numeric(result1$PE.type.1.CI.high),
                          "P" = format.p(result1$PE.type.1.p.value),
                          "DiffP" = "")
       
       PEtable <- add_row(.data = PEtable , "mark" = "All founders predicted resistant", 
                          "n" = paste0(result2$.n.events.type.1.trt, " vs. ", result2$.n.events.type.1.placebo) , 
                          "PE" = paste0(format.PE(result2$PE.type.1.estimate)," (", format.PE(result2$PE.type.1.CI.low),", ",format.PE(result2$PE.type.1.CI.high),")" ),
                          "mean" = as.numeric(result2$PE.type.1.estimate),
                          "lower" = as.numeric(result2$PE.type.1.CI.low), #to be able to plot -Inf
                          "upper" = as.numeric(result2$PE.type.1.CI.high),
                          "P" = format.p(result2$PE.type.1.p.value),
                          "DiffP" = "")
       
      
      write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_founderType_",ifelse(dose=="T1+T2", "dosePooled", dose),".csv")))
      PEtable.forestplot.withSieveT2(PEtable,
                                    xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                                    figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_founderType_",ifelse(dose=="T1+T2", "dosePooled", dose),".pdf"))
      #forest plots
    }else{
      
      PEtable <- tibble("mark" = character(), "n" = character(),
                        "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                        "P" = character())
      result1 <- filter(results, mark == "all.founders.sensitive")
      result2 <- filter(results, mark == "all.founders.resistant")
     
     
      PEtable <- add_row(.data = PEtable , "mark" = "All founders predicted sensitive", 
                         "n" = paste0(result1$.n.events.type.1.trt, " vs. ", result1$.n.events.type.1.placebo) , 
                         "PE" = paste0(format.PE(result1$PE.type.1.estimate)," (", format.PE(result1$PE.type.1.CI.low),", ",format.PE(result1$PE.type.1.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.1.estimate),
                         "lower" = as.numeric(result1$PE.type.1.CI.low), #to be able to plot -Inf
                         "upper" = as.numeric(result1$PE.type.1.CI.high),
                         "P" = format.p(result1$PE.type.1.p.value))
      
      PEtable <- add_row(.data = PEtable , "mark" = "At least 1 founder predicted resistant",
                         "n" = paste0(result1$.n.events.type.0.trt, " vs. ", result1$.n.events.type.0.placebo) ,
                         "PE" = paste0(format.PE(result1$PE.type.0.estimate)," (", format.PE(result1$PE.type.0.CI.low),", ",format.PE(result1$PE.type.0.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.0.estimate),
                         "lower" = as.numeric(result1$PE.type.0.CI.low),
                         "upper" = as.numeric(result1$PE.type.0.CI.high),
                         "P" = format.p(result1$PE.type.0.p.value))
      
      
      PEtable <- add_row(.data = PEtable , "mark" = "All founders predicted resistant", 
                         "n" = paste0(result2$.n.events.type.1.trt, " vs. ", result2$.n.events.type.1.placebo) , 
                         "PE" = paste0(format.PE(result2$PE.type.1.estimate)," (", format.PE(result2$PE.type.1.CI.low),", ",format.PE(result2$PE.type.1.CI.high),")" ),
                         "mean" = as.numeric(result2$PE.type.1.estimate),
                         "lower" = as.numeric(result2$PE.type.1.CI.low), #to be able to plot -Inf
                         "upper" = as.numeric(result2$PE.type.1.CI.high),
                         "P" = format.p(result2$PE.type.1.p.value))
      
      write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_founderType_",ifelse(dose=="T1+T2", "dosePooled", dose),".csv")))
      PEtable.forestplot( PEtable,
                          xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                          figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_founderType_",ifelse(dose=="T1+T2", "dosePooled", dose),".pdf"),
                          header = c("Founder Multiplicity Type","No. Cases (VRC01 vs. P)\n(Incidence per PYR)",
                                     "PE (%) (95% CI)","mean", "lower", "upper","P-value"),
                          ncolors = 1)
      
      
    }
    
    
}}    




# Glycosite 230 related features
#glycosite230 <- c("hxb2.230.N", "hxb2.230.pngs")
#marks <- colnames(master)%>% .[grepl("230", .)]
marks <- c("hxb2.230.D", "hxb2.230.N", "hxb2.230.pngs")
for(sequence in c("mf", "ms", "ls")){
  features.to.analyze.x <- paste0(marks, ".", sequence)
  for(trial in c("703", "704", "704and703")){
    for(dose in c("T1", "T2", "T1+T2")){
      master.trial <- trial_dose_data (master, trial, dose)
      
      ### List of columns in output 
      results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                             "DPE.p.value", 
                             paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                             paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;
      
      results <- data.frame( feature = features.to.analyze.x, feature.name = features.to.analyze.x)
      results[ , 3:( length( results.colnames) ) ] <- NA;
      colnames( results ) <- results.colnames
      .results.row.i = 0
      for(.feature in features.to.analyze.x){
        .results.row.i <- .results.row.i + 1;
        
        # Be verbose:
        cat( "Computing results for feature", results[ .results.row.i, 1 ], "\n" );
        
        # Grab the column containing feature values from the master.trial table
        .x <- master.trial[ , .feature]
        
        #if mark is all.founders.all.sens.resist, censor all events with NA mark
        if(.feature == "all.founders.all.sens.resist"){
          master.trial$hiv1event[master.trial$hiv1event==1&is.na(.x)] <- 0
        }
        
        ### Run DVE
        .dPE.results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event, .x, master.trial$txLabel, master.trial$stratVar )
        
        
        # Insert into results table 
        results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
                                     paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                                     paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.results , dose, .x)
      }
      
      if(dose == "T1+T2"){
        PEtable <- tibble("mark" = character(), "Haplotype" = character(), "n" = character(),
                          "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                          "P" = character(), "DiffP" = character(), "FWER" = character(), "FDR" = character())
        adjPvaluesDirSievePH <- file.path(repoDir, "tables", "sievePH_JuraskaGilbert2013")
        p.value.file <- read.csv(file.path(adjPvaluesDirSievePH, paste0(trial,"_WestfallYoungAdjPvalues_neutAssocSeqFeatures_",sequence,".csv")))
        results <- left_join(results, p.value.file, by = "mark")
        for( .feature in features.to.analyze.x){
          result.f <- filter(results, mark == .feature)
          
          sig.unadj <- ifelse(result.f$DPE.p.value<=0.05, 1, 0) 
          
          #annotatedMark = ifelse(sig.unadj, paste(annotatedMark,"$^{\\P}$",sep=""),annotatedMark)
          pngs.ind <- strsplit( result.f$mark.name, split = "\\.")[[1]][3]
          PEtable <- add_row(.data = PEtable ,"mark" = ifelse(pngs.ind == "pngs", "230-232", "230"), "Haplotype" = "", "n" = "",
                             "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                             "DiffP" = format.p(result.f$DPE.p.value),"FWER" = format.p(result.f$p.FWER), "FDR" = format.p(result.f$p.FDR))          
          Haplotype.x <- table.seqFeatLabel.tier1(paste0(strsplit( result.f$mark.name, split = "\\.")[[1]][1:3], collapse="."))
          
          PEtable <- add_row(.data = PEtable , "mark" = "", "Haplotype" = Haplotype.x$Haplotype1, 
                             "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo) , 
                             "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                             "mean" = as.numeric(result.f$PE.type.1.estimate),
                             "lower" = as.numeric(result.f$PE.type.1.CI.low), #to be able to plot -Inf
                             "upper" = as.numeric(result.f$PE.type.1.CI.high),
                             "P" = format.p(result.f$PE.type.1.p.value),
                             "DiffP" = "","FWER" = "", "FDR" = "")
          PEtable <- add_row(.data = PEtable , "mark" = "", "Haplotype" = Haplotype.x$Haplotype0, 
                             "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                             "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                             "mean" = as.numeric(result.f$PE.type.0.estimate),
                             "lower" = as.numeric(result.f$PE.type.0.CI.low),
                             "upper" = as.numeric(result.f$PE.type.0.CI.high),
                             "P" = format.p(result.f$PE.type.0.p.value),
                             "DiffP" = "","FWER" = "", "FDR" = "")
          
        }
        
        write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_","glycosite230","_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,".csv")))
        PEtable.forestplot.withSieveT(PEtable,
                                       xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                                       figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_","glycosite230","_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,".pdf"))
        #forest plots
      }else{
        
        PEtable <- tibble("mark" = character(), "Haplotype" = character(), "n" = character(),
                          "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                          "P" = character())
        for( .feature in features.to.analyze.x ){
          result.f <- filter(results, mark == .feature)
          Haplotype.x <- table.seqFeatLabel.tier1(paste0(strsplit( result.f$mark.name, split = "\\.")[[1]][1:3], collapse="."))
          pngs.ind <- strsplit( result.f$mark.name, split = "\\.")[[1]][3]
          PEtable <- add_row(.data = PEtable,"mark" = ifelse(pngs.ind == "pngs", "230-232", "230"), "Haplotype" = Haplotype.x$Haplotype1, 
                             "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo),
                             "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                             "mean" = as.numeric(result.f$PE.type.1.estimate),
                             "lower" = as.numeric(result.f$PE.type.1.CI.low),
                             "upper" = as.numeric(result.f$PE.type.1.CI.high),
                             "P" = format.p(result.f$PE.type.1.p.value))
          PEtable <- add_row(.data = PEtable, "mark" = "", "Haplotype" = Haplotype.x$Haplotype0, 
                             "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                             "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                             "mean" = as.numeric(result.f$PE.type.0.estimate),
                             "lower" = as.numeric(result.f$PE.type.0.CI.low),
                             "upper" = as.numeric(result.f$PE.type.0.CI.high),
                             "P" = format.p(result.f$PE.type.0.p.value))
          
        }
        
        write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_","glycosite230","_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,".csv")))
        PEtable.forestplot( PEtable,
                            xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                            figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_","glycosite230","_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,".pdf"),
                            header = c("Seq Feature","Seq\nFeature","No. Cases (VRC01 vs. P)\n(Incidence per PYR)",
                                       "PE (%) (95% CI)","mean", "lower", "upper","P-value"))
        
        
      }
      
      
    }}    
}

####PNGS forest plots

for(sequence in c("mf", "ms", "ls")){
  features.to.analyze.x <- paste0("hxb2.230.pngs", ".", sequence)
  PEtable <- tibble("trial" = character(),"dose" = character(), "Haplotype" = character(), "n" = character(),
                    "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                    "P" = character(), "DiffP" = character())
  
  for(trial in c("704", "703")){
    for(dose in c("T1+T2","T2", "T1")){
      master.trial <- trial_dose_data (master, trial, dose)
      doseDesc <- ifelse(dose == "T1", "10mg/kg", ifelse(dose == "T2", "30 mg/kg", "Dose-pooled"))
      trialDesc <- ifelse(trial == "704", "HVTN 704/HPTN085", "HVTN703/HPTN081")
      ### List of columns in output 
      results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                             "DPE.p.value", 
                             paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                             paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;
      
      results <- data.frame( feature = features.to.analyze.x, feature.name = features.to.analyze.x)
      .feature <- features.to.analyze.x
      results[ , 3:( length( results.colnames) ) ] <- NA;
      colnames( results ) <- results.colnames
      .results.row.i = 1
      # Be verbose:
      cat( "Computing results for feature", results[ .results.row.i, 1 ], "\n" );
      
      # Grab the column containing feature values from the master.trial table
      .x <- master.trial[ , .feature]
      
      #if mark is all.founders.all.sens.resist, censor all events with NA mark
      if(.feature == "all.founders.all.sens.resist"){
        master.trial$hiv1event[master.trial$hiv1event==1&is.na(.x)] <- 0
      }
      
      ### Run DVE
      .dPE.results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event, .x, master.trial$txLabel, master.trial$stratVar )
      
      
      # Insert into results table 
      results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
                                   paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                                   paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.results , dose, .x)
      
      
      if(dose == "T1+T2"){
        
        result.f <- filter(results, mark == features.to.analyze.x)
        
        PEtable <- add_row(.data = PEtable ,"trial" = trialDesc,"dose" = doseDesc, "Haplotype" = "", "n" = "",
                           "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                           "DiffP" = format.p(result.f$DPE.p.value,2))          
        Haplotype.x <- table.seqFeatLabel.tier1(paste0(strsplit( result.f$mark.name, split = "\\.")[[1]][1:3], collapse="."))
        
        PEtable <- add_row(.data = PEtable,trial = "" , "dose" = "", "Haplotype" = Haplotype.x$Haplotype1, 
                           "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo) , 
                           "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                           "mean" = as.numeric(result.f$PE.type.1.estimate),
                           "lower" = as.numeric(result.f$PE.type.1.CI.low), #to be able to plot -Inf
                           "upper" = as.numeric(result.f$PE.type.1.CI.high),
                           "P" = format.p(result.f$PE.type.1.p.value,2),
                           "DiffP" = "")
        PEtable <- add_row(.data = PEtable,trial = "" , "dose" = "", "Haplotype" = Haplotype.x$Haplotype0, 
                           "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                           "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                           "mean" = as.numeric(result.f$PE.type.0.estimate),
                           "lower" = as.numeric(result.f$PE.type.0.CI.low),
                           "upper" = as.numeric(result.f$PE.type.0.CI.high),
                           "P" = format.p(result.f$PE.type.0.p.value,2),
                           "DiffP" = "")
        
           #forest plots
      }else{
        
        result.f <- filter(results, mark == features.to.analyze.x)
        PEtable <- add_row(.data = PEtable, trial = "" ,"dose" = doseDesc, "Haplotype" = "", "n" = "",
                           "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                           "DiffP" = "--")  
        Haplotype.x <- table.seqFeatLabel.tier1(paste0(strsplit( result.f$mark.name, split = "\\.")[[1]][1:3], collapse="."))
        
        PEtable <- add_row(.data = PEtable, trial = "" ,"dose" = "", "Haplotype" = Haplotype.x$Haplotype1, 
                           "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo),
                           "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                           "mean" = as.numeric(result.f$PE.type.1.estimate),
                           "lower" = as.numeric(result.f$PE.type.1.CI.low),
                           "upper" = as.numeric(result.f$PE.type.1.CI.high),
                           "P" = format.p(result.f$PE.type.1.p.value), "DiffP" = "")
        PEtable <- add_row(.data = PEtable, trial = "" , "dose" = "", "Haplotype" = Haplotype.x$Haplotype0, 
                           "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                           "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                           "mean" = as.numeric(result.f$PE.type.0.estimate),
                           "lower" = as.numeric(result.f$PE.type.0.CI.low),
                           "upper" = as.numeric(result.f$PE.type.0.CI.high),
                           "P" = format.p(result.f$PE.type.0.p.value), "DiffP" = "")
        
        
         
        
      }
      
      
    }
    
  } 
    write.csv(PEtable, file.path(tablesOutDir, paste0("DVE_tier2_","glycosite230PNGS","_",sequence,".csv")))
    PEtable.forestplot.withSieveT4( PEtable[,-1],
                        xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                        figuresOutDir, paste0("DVE_tier2_","glycosite230PNGS","_",sequence,".pdf"),
                        header = c("VRC01\nDose","Seq\nFeature","No. Cases (VRC01 vs. P)\n(Incidence per PYR)",
                                   "PE (%) (95% CI)","mean", "lower", "upper","P-value","Diff PE\nP-value"),
                        add.line = 9)
    
    
       
}




#number of founders
num.founders <- c("num.founders.all", "num.founders.tfl")
for(trial in c("703", "704", "704and703")){
  for(dose in c("T1", "T2", "T1+T2")){
    master.trial <- trial_dose_data (master, trial, dose)
    
    ### List of columns in output 
    results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                           "DPE.p.value", 
                           paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                           paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;
    
    results <- data.frame( feature = num.founders, feature.name = num.founders)
    results[ , 3:( length( results.colnames) ) ] <- NA;
    colnames( results ) <- results.colnames
    .results.row.i = 0
    for(.feature in num.founders){
      .results.row.i <- .results.row.i + 1;
      
      # Be verbose:
      cat( "Computing results for feature", results[ .results.row.i, 1 ], "\n" );
      
      # Grab the column containing feature values from the master.trial table
      .x <- master.trial[ , .feature]
      # Dichotomize the feature (=1 vs. >1)
      .x <- 1*(.x >1)
     
      
      ### Run DVE
      .dPE.results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event, .x, master.trial$txLabel, master.trial$stratVar )
      
      
      # Insert into results table 
      results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
                                   paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                                   paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.results , dose,.x)
    }
    
    if(dose == "T1+T2"){
      p.value.file <- read.csv(file.path(adjPvaluesDir, paste0(trial,"_WestfallYoungAdjPvalues_tier2_numFounders",".csv")))
      results <- full_join(results, p.value.file, by = "mark")
      
      PEtable <- tibble("mark" = character(), "n" = character(),
                        "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                        "P" = character(), "DiffP" = character(), "holm" = character(), "fdr" = character())
      
     
      result1 <- filter(results, mark == "num.founders.tfl")
      result2 <- filter(results, mark == "num.founders.all")
      
      PEtable <- add_row(.data = PEtable ,"mark" = "", "n" = "",
                         "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                         "DiffP" = format.p(result1$p.unadj), "holm" = format.p(result1$p.FWER), "fdr" = format.p(result1$p.FDR)) 
     
      
      PEtable <- add_row(.data = PEtable , "mark" = "Single early lineage",
                         "n" = paste0(result1$.n.events.type.0.trt, " vs. ", result1$.n.events.type.0.placebo) ,
                         "PE" = paste0(format.PE(result1$PE.type.0.estimate)," (", format.PE(result1$PE.type.0.CI.low),", ",format.PE(result1$PE.type.0.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.0.estimate),
                         "lower" = as.numeric(result1$PE.type.0.CI.low),
                         "upper" = as.numeric(result1$PE.type.0.CI.high),
                         "P" = format.p(result1$PE.type.0.p.value),
                         "DiffP" = "", "holm" = "", "fdr" = "")
      PEtable <- add_row(.data = PEtable , "mark" = "Multiple early lineages", 
                         "n" = paste0(result1$.n.events.type.1.trt, " vs. ", result1$.n.events.type.1.placebo) , 
                         "PE" = paste0(format.PE(result1$PE.type.1.estimate)," (", format.PE(result1$PE.type.1.CI.low),", ",format.PE(result1$PE.type.1.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.1.estimate),
                         "lower" = as.numeric(result1$PE.type.1.CI.low), #to be able to plot -Inf
                         "upper" = as.numeric(result1$PE.type.1.CI.high),
                         "P" = format.p(result1$PE.type.1.p.value),
                         "DiffP" = "", "holm" = "", "fdr" = "")
      
      
      PEtable <- add_row(.data = PEtable ,"mark" = "", "n" = "",
                         "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                         "DiffP" = format.p(result2$p.unadj), "holm" = format.p(result2$p.FWER), "fdr" = format.p(result2$p.FDR)) 
      
      
      PEtable <- add_row(.data = PEtable , "mark" = "No. unique lineages = 1",
                         "n" = paste0(result2$.n.events.type.0.trt, " vs. ", result2$.n.events.type.0.placebo) ,
                         "PE" = paste0(format.PE(result2$PE.type.0.estimate)," (", format.PE(result2$PE.type.0.CI.low),", ",format.PE(result2$PE.type.0.CI.high),")" ),
                         "mean" = as.numeric(result2$PE.type.0.estimate),
                         "lower" = as.numeric(result2$PE.type.0.CI.low),
                         "upper" = as.numeric(result2$PE.type.0.CI.high),
                         "P" = format.p(result2$PE.type.0.p.value),
                         "DiffP" = "", "holm" = "", "fdr" = "")
      PEtable <- add_row(.data = PEtable , "mark" = "No. unique lineages > 1", 
                         "n" = paste0(result2$.n.events.type.1.trt, " vs. ", result2$.n.events.type.1.placebo) , 
                         "PE" = paste0(format.PE(result2$PE.type.1.estimate)," (", format.PE(result2$PE.type.1.CI.low),", ",format.PE(result2$PE.type.1.CI.high),")" ),
                         "mean" = as.numeric(result2$PE.type.1.estimate),
                         "lower" = as.numeric(result2$PE.type.1.CI.low), #to be able to plot -Inf
                         "upper" = as.numeric(result2$PE.type.1.CI.high),
                         "P" = format.p(result2$PE.type.1.p.value),
                         "DiffP" = "", "holm" = "", "fdr" = "")
      
      
      write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_num.founders_",ifelse(dose=="T1+T2", "dosePooled", dose),".csv")))
      PEtable.forestplot.withSieveT3(PEtable,
                                     xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                                     figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_num_founders_",ifelse(dose=="T1+T2", "dosePooled", dose),".pdf"))
      #forest plots
    }else{
      
      PEtable <- tibble("mark" = character(), "n" = character(),
                        "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                        "P" = character())
      
      results$holm <- p.adjust(results$DPE.p.value, method = "holm")
      results$fdr<- p.adjust(results$DPE.p.value, method = "fdr")
      
      result1 <- filter(results, mark == "num.founders.tfl")
      result2 <- filter(results, mark == "num.founders.all")
      
   
      
      PEtable <- add_row(.data = PEtable , "mark" = "Single early lineage",
                         "n" = paste0(result1$.n.events.type.0.trt, " vs. ", result1$.n.events.type.0.placebo) ,
                         "PE" = paste0(format.PE(result1$PE.type.0.estimate)," (", format.PE(result1$PE.type.0.CI.low),", ",format.PE(result1$PE.type.0.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.0.estimate),
                         "lower" = as.numeric(result1$PE.type.0.CI.low),
                         "upper" = as.numeric(result1$PE.type.0.CI.high),
                         "P" = format.p(result1$PE.type.0.p.value))
      PEtable <- add_row(.data = PEtable , "mark" = "Multiple early lineages", 
                         "n" = paste0(result1$.n.events.type.1.trt, " vs. ", result1$.n.events.type.1.placebo) , 
                         "PE" = paste0(format.PE(result1$PE.type.1.estimate)," (", format.PE(result1$PE.type.1.CI.low),", ",format.PE(result1$PE.type.1.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.1.estimate),
                         "lower" = as.numeric(result1$PE.type.1.CI.low), #to be able to plot -Inf
                         "upper" = as.numeric(result1$PE.type.1.CI.high),
                         "P" = format.p(result1$PE.type.1.p.value))
      
      
     
      
      PEtable <- add_row(.data = PEtable , "mark" = "No. unique lineages = 1",
                         "n" = paste0(result2$.n.events.type.0.trt, " vs. ", result2$.n.events.type.0.placebo) ,
                         "PE" = paste0(format.PE(result2$PE.type.0.estimate)," (", format.PE(result2$PE.type.0.CI.low),", ",format.PE(result2$PE.type.0.CI.high),")" ),
                         "mean" = as.numeric(result2$PE.type.0.estimate),
                         "lower" = as.numeric(result2$PE.type.0.CI.low),
                         "upper" = as.numeric(result2$PE.type.0.CI.high),
                         "P" = format.p(result2$PE.type.0.p.value))
      PEtable <- add_row(.data = PEtable , "mark" = "No. unique lineages > 1", 
                         "n" = paste0(result2$.n.events.type.1.trt, " vs. ", result2$.n.events.type.1.placebo) , 
                         "PE" = paste0(format.PE(result2$PE.type.1.estimate)," (", format.PE(result2$PE.type.1.CI.low),", ",format.PE(result2$PE.type.1.CI.high),")" ),
                         "mean" = as.numeric(result2$PE.type.1.estimate),
                         "lower" = as.numeric(result2$PE.type.1.CI.low), #to be able to plot -Inf
                         "upper" = as.numeric(result2$PE.type.1.CI.high),
                         "P" = format.p(result2$PE.type.1.p.value))
      
      write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_num.founders_",ifelse(dose=="T1+T2", "dosePooled", dose),".csv")))
      PEtable.forestplot( PEtable,
                          xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                          figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_num_founders_",ifelse(dose=="T1+T2", "dosePooled", dose),".pdf"),
                          header = c("No. Lineages","No. Cases (VRC01 vs. P)\n(Incidence per PYR)",
                                     "PE (%) (95% CI)","mean", "lower", "upper","P-value"),
                          ncolors = 1, addLine = 2)
      
      
    }
    
    
  }}



#number of founders version 2
num.founders <- c("num.founders.tfl")
for(trial in c("703", "704", "704and703")){
  for(dose in c("T1", "T2", "T1+T2")){
    master.trial <- trial_dose_data (master, trial, dose)
    
    ### List of columns in output 
    results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                           "DPE.p.value", 
                           paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                           paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;
    
    results <- data.frame( feature = num.founders, feature.name = num.founders)
    results[ , 3:( length( results.colnames) ) ] <- NA;
    colnames( results ) <- results.colnames
    .results.row.i = 0
    for(.feature in num.founders){
      .results.row.i <- .results.row.i + 1;
      
      # Be verbose:
      cat( "Computing results for feature", results[ .results.row.i, 1 ], "\n" );
      
      # Grab the column containing feature values from the master.trial table
      .x <- master.trial[ , .feature]
      # Dichotomize the feature (=1 vs. >1)
      .x <- 1*(.x >1)
      
      
      ### Run DVE
      .dPE.results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event, .x, master.trial$txLabel, master.trial$stratVar )
      
      
      # Insert into results table 
      results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
                                   paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                                   paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.results , dose, .x)
    }
    
    if(dose == "T1+T2"){
      p.value.file <- read.csv(file.path(adjPvaluesDir, paste0(trial,"_WestfallYoungAdjPvalues_tier2_numFounders",".csv")))
      results <- full_join(results, p.value.file, by = "mark")
      
      PEtable <- tibble("mark" = character(), "n" = character(),
                        "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                        "P" = character(), "DiffP" = character())
      
      
      result1 <- filter(results, mark == "num.founders.tfl")
     
      PEtable <- add_row(.data = PEtable ,"mark" = "", "n" = "",
                         "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                         "DiffP" = format.p(result1$p.unadj)) 
      
      
      PEtable <- add_row(.data = PEtable , "mark" = "Single early lineage",
                         "n" = paste0(result1$.n.events.type.0.trt, " vs. ", result1$.n.events.type.0.placebo) ,
                         "PE" = paste0(format.PE(result1$PE.type.0.estimate)," (", format.PE(result1$PE.type.0.CI.low),", ",format.PE(result1$PE.type.0.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.0.estimate),
                         "lower" = as.numeric(result1$PE.type.0.CI.low),
                         "upper" = as.numeric(result1$PE.type.0.CI.high),
                         "P" = format.p(result1$PE.type.0.p.value),
                         "DiffP" = "")
      PEtable <- add_row(.data = PEtable , "mark" = "Multiple early lineages", 
                         "n" = paste0(result1$.n.events.type.1.trt, " vs. ", result1$.n.events.type.1.placebo) , 
                         "PE" = paste0(format.PE(result1$PE.type.1.estimate)," (", format.PE(result1$PE.type.1.CI.low),", ",format.PE(result1$PE.type.1.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.1.estimate),
                         "lower" = as.numeric(result1$PE.type.1.CI.low), #to be able to plot -Inf
                         "upper" = as.numeric(result1$PE.type.1.CI.high),
                         "P" = format.p(result1$PE.type.1.p.value),
                         "DiffP" = "")
      
     
     
      
      
      write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_num_early_founders_",ifelse(dose=="T1+T2", "dosePooled", dose),".csv")))
      PEtable.forestplot.withSieveT2(PEtable,
                                     xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                                     figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_num_early_founders_",ifelse(dose=="T1+T2", "dosePooled", dose),".pdf"))
      #forest plots
    }else{
      
      PEtable <- tibble("mark" = character(), "n" = character(),
                        "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                        "P" = character())
      
      results$holm <- p.adjust(results$DPE.p.value, method = "holm")
      results$fdr<- p.adjust(results$DPE.p.value, method = "fdr")
      
      result1 <- filter(results, mark == "num.founders.tfl")
      
      
      PEtable <- add_row(.data = PEtable , "mark" = "Single early lineage",
                         "n" = paste0(result1$.n.events.type.0.trt, " vs. ", result1$.n.events.type.0.placebo) ,
                         "PE" = paste0(format.PE(result1$PE.type.0.estimate)," (", format.PE(result1$PE.type.0.CI.low),", ",format.PE(result1$PE.type.0.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.0.estimate),
                         "lower" = as.numeric(result1$PE.type.0.CI.low),
                         "upper" = as.numeric(result1$PE.type.0.CI.high),
                         "P" = format.p(result1$PE.type.0.p.value))
      PEtable <- add_row(.data = PEtable , "mark" = "Multiple early lineages", 
                         "n" = paste0(result1$.n.events.type.1.trt, " vs. ", result1$.n.events.type.1.placebo) , 
                         "PE" = paste0(format.PE(result1$PE.type.1.estimate)," (", format.PE(result1$PE.type.1.CI.low),", ",format.PE(result1$PE.type.1.CI.high),")" ),
                         "mean" = as.numeric(result1$PE.type.1.estimate),
                         "lower" = as.numeric(result1$PE.type.1.CI.low), #to be able to plot -Inf
                         "upper" = as.numeric(result1$PE.type.1.CI.high),
                         "P" = format.p(result1$PE.type.1.p.value))
      
      
      
      
      
      
      write.csv(PEtable, file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_num_early_founders_",ifelse(dose=="T1+T2", "dosePooled", dose),".csv")))
      PEtable.forestplot( PEtable,
                          xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                          figuresOutDir, paste0("HVTN",trial,"_DVE_tier2_num_early_founders_",ifelse(dose=="T1+T2", "dosePooled", dose),".pdf"),
                          header = c("No. Lineages","No. Cases (VRC01 vs. P)\n(Incidence per PYR)",
                                     "PE (%) (95% CI)","mean", "lower", "upper","P-value"),
                          ncolors = 1, addLine = 2)
      
      
    }
    
    
  }}


#704 subtypeB PNGS230

for(sequence in c("ls")){
  features.to.analyze.x <- paste0("hxb2.230.pngs", ".", sequence)
  PEtable <- tibble("trial" = character(),"dose" = character(), "Haplotype" = character(), "n" = character(),
                    "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                    "P" = character(), "DiffP" = character())
  
  for(trial in c("704")){
    for(dose in c("T1+T2","T2", "T1")){
      master.trial <- trial_dose_data (master, trial, dose)
      #restrict cases to be endpoints being subtypeB and censor endpoints that are not type B
      master.trial$hiv1event[master.trial$subtype !="B"] = 0
     
      doseDesc <- ifelse(dose == "T1", "10mg/kg", ifelse(dose == "T2", "30 mg/kg", "Dose-pooled"))
      trialDesc <- ifelse(trial == "704", "HVTN 704/HPTN085", "HVTN703/HPTN081")
      ### List of columns in output 
      results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                             "DPE.p.value", 
                             paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                             paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;
      
      results <- data.frame( feature = features.to.analyze.x, feature.name = features.to.analyze.x)
      .feature <- features.to.analyze.x
      results[ , 3:( length( results.colnames) ) ] <- NA;
      colnames( results ) <- results.colnames
      .results.row.i = 1
      # Be verbose:
      cat( "Computing results for feature", results[ .results.row.i, 1 ], "\n" );
      
      # Grab the column containing feature values from the master.trial table
      .x <- master.trial[ , .feature]
      
      #if mark is all.founders.all.sens.resist, censor all events with NA mark
      if(.feature == "all.founders.all.sens.resist"){
        master.trial$hiv1event[master.trial$hiv1event==1&is.na(.x)] <- 0
      }
      
      ### Run DVE
      .dPE.results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event, .x, master.trial$txLabel, master.trial$stratVar )
      
      
      # Insert into results table 
      results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DPE.p.value", 
                                   paste( "PE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                                   paste( "PE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- result.summary (.dPE.results , dose, .x)
      
      
      if(dose == "T1+T2"){
        
        result.f <- filter(results, mark == features.to.analyze.x)
        
        PEtable <- add_row(.data = PEtable ,"trial" = trialDesc,"dose" = doseDesc, "Haplotype" = "", "n" = "",
                           "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                           "DiffP" = format.p(result.f$DPE.p.value,2))          
        Haplotype.x <- table.seqFeatLabel.tier1(paste0(strsplit( result.f$mark.name, split = "\\.")[[1]][1:3], collapse="."))
        
        PEtable <- add_row(.data = PEtable,trial = "" , "dose" = "", "Haplotype" = Haplotype.x$Haplotype1, 
                           "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo) , 
                           "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                           "mean" = as.numeric(result.f$PE.type.1.estimate),
                           "lower" = as.numeric(result.f$PE.type.1.CI.low), #to be able to plot -Inf
                           "upper" = as.numeric(result.f$PE.type.1.CI.high),
                           "P" = format.p(result.f$PE.type.1.p.value,2),
                           "DiffP" = "")
        PEtable <- add_row(.data = PEtable,trial = "" , "dose" = "", "Haplotype" = Haplotype.x$Haplotype0, 
                           "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                           "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                           "mean" = as.numeric(result.f$PE.type.0.estimate),
                           "lower" = as.numeric(result.f$PE.type.0.CI.low),
                           "upper" = as.numeric(result.f$PE.type.0.CI.high),
                           "P" = format.p(result.f$PE.type.0.p.value,2),
                           "DiffP" = "")
        
        #forest plots
      }else{
        
        result.f <- filter(results, mark == features.to.analyze.x)
        PEtable <- add_row(.data = PEtable, trial = "" ,"dose" = doseDesc, "Haplotype" = "", "n" = "",
                           "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                           "DiffP" = "--")  
        Haplotype.x <- table.seqFeatLabel.tier1(paste0(strsplit( result.f$mark.name, split = "\\.")[[1]][1:3], collapse="."))
        
        PEtable <- add_row(.data = PEtable, trial = "" ,"dose" = "", "Haplotype" = Haplotype.x$Haplotype1, 
                           "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo),
                           "PE" = paste0(format.PE(result.f$PE.type.1.estimate)," (", format.PE(result.f$PE.type.1.CI.low),", ",format.PE(result.f$PE.type.1.CI.high),")" ),
                           "mean" = as.numeric(result.f$PE.type.1.estimate),
                           "lower" = as.numeric(result.f$PE.type.1.CI.low),
                           "upper" = as.numeric(result.f$PE.type.1.CI.high),
                           "P" = format.p(result.f$PE.type.1.p.value), "DiffP" = "")
        PEtable <- add_row(.data = PEtable, trial = "" , "dose" = "", "Haplotype" = Haplotype.x$Haplotype0, 
                           "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                           "PE" = paste0(format.PE(result.f$PE.type.0.estimate)," (", format.PE(result.f$PE.type.0.CI.low),", ",format.PE(result.f$PE.type.0.CI.high),")" ),
                           "mean" = as.numeric(result.f$PE.type.0.estimate),
                           "lower" = as.numeric(result.f$PE.type.0.CI.low),
                           "upper" = as.numeric(result.f$PE.type.0.CI.high),
                           "P" = format.p(result.f$PE.type.0.p.value), "DiffP" = "")
        
        
        
        
      }
      
      
    }
    
  } 
  write.csv(PEtable, file.path(tablesOutDir, paste0("DVE_tier2_","glycosite230PNGS","_704SubtypeB",sequence,".csv")))
  PEtable.forestplot.withSieveT5( PEtable[,-1],
                                  xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                                  figuresOutDir, paste0("DVE_tier2_","glycosite230PNGS","_704SubtypeB",sequence,".pdf"),
                                  header = c("VRC01\nDose","Seq\nFeature","No. Cases (VRC01 vs. P)\n(Incidence per PYR)",
                                             "PE (%) (95% CI)","mean", "lower", "upper","P-value","Diff PE\nP-value"))
  
  
  
}









#viruses with both IC80 ≤ 1 µg/ml and a PNGS at 230–232
master.trial <- trial_dose_data (master, "704", "T1+T2")
master.trial$gmt80ls[master.trial$gmt80ls==">100"] = "100"
master.trial$mark <- 1*(as.numeric(master.trial$gmt80ls) <=1)*(master.trial$hxb2.230.pngs.ls==1)
.x <-  master.trial$mark
results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event,.x, master.trial$txLabel, master.trial$stratVar )

result.summary (results, "T1+T2", .x)

#test the risk ratio 
eventInd1 <- master.trial$hiv1event==1&.x==1
df <- table(master.trial$txLabel,eventInd1)

rownames(df) <- c("placebo","vaccine")
colnames(df) <- c("nonEvent", "Event")

df2 <- df[c(2,1),c(2,1)]
library(DescTools)
1-RelRisk(df2, method = "score", conf.level = 0.95)
1-RelRisk(df2, method = "wald", conf.level = 0.95)
1-RelRisk(df2, method = "wald", delta = 0.5, conf.level = 0.95)


pv = (0.5)/(1784+0.5)
pp = (3+0.5)/(3+892+0.5)
v = 1/3.5+1/0.5-1/895.5-1/1784.5
#confidence limits for VE
ci = c(1-pv/pp*exp(1.96*sqrt(v)), 1 - pv/pp*exp(-1.96*sqrt(v)))
est = 1-pv/pp



pv = -log(0.05)/1784
pp = 3/895
v = pv/pp




