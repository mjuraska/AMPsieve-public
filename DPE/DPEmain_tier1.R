# Purpose: Estimation of hazard ratio-based PE by tier1 AA position-specific features; 
# Hypothesis testing method:  Lunn & McNeil 
#          Custom script allowing for stratification of baseline hazards, see dve3() below
# Input:  Master file (.csv) with time to event data, treatment group assignment, region (south america vs. other) and binary marks to be analyzed
# Output: .Rdata containing results table with estimates of mark-specific PE, hazard rate ratio, and associated confidence intervals   
# Author:  Li Li
# Date:    Sep 6, 2022

rm(list=ls(all=TRUE))
repoDir <- "/Users/lili/AMPsieve"
adjPvaluesDir <- file.path(repoDir, "tables", "sievePH_JuraskaGilbert2013")
tablesOutDir <- file.path(repoDir, "tables/DPE")
figuresOutDir <- file.path(repoDir, "figures/DPE/tier1")
library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)

source(file.path(repoDir, "code/common.R")) #read in a datFile
source(file.path(repoDir, "code/lunnMcneil.R"))
source(file.path(repoDir,"code/DVE/DPEutils.R"))
source(file.path(repoDir, "code/DVE/forest.R"))

### Read master file with covariates and time to event data for both trials (note the column southAmerica now has the combined trial/location info as follows:  value 0 = 704, not south america, 1 = 704, south america, 2 = 703)
master <-  read.csv( file.path(repoDir,"adata", datFile[3])) 

master$num.pngs.v5.mf <- I( master$num.pngs.v5.mf > 1 )
master$num.cysteine.gp120.mf <- I( master$num.cysteine.gp120.mf > 18 )
master$num.pngs.v5.ms <- I( master$num.pngs.v5.ms > 1 )
master$num.cysteine.gp120.ms. <- I( master$num.cysteine.gp120.ms > 18 )
master$num.pngs.v5.ls <- I( master$num.pngs.v5.ls > 1 )
master$num.cysteine.gp120.ls <- I( master$num.cysteine.gp120.ls > 18 )

### List of features to analyze
features.to.analyze <- c( "hxb2.60.A", "hxb2.170.Q", "hxb2.230.D", "hxb2.279.N", "hxb2.280.N", "hxb2.317.F",
                          "hxb2.365.S", "hxb2.429.E", "hxb2.456.R", "hxb2.458.G", "hxb2.459.G", "hxb2.471.G",
                          "hxb2.156.pngs", "hxb2.229.pngs", "hxb2.234.pngs", "hxb2.616.pngs", "hxb2.824.pngs", 
                          "num.pngs.v5", "num.cysteine.gp120" )

### Make some more readable feature names
feature.names <- sapply( strsplit(features.to.analyze,"\\."), function(.feature){
  paste( .feature[2])
}) 
feature.names[ 18:19 ] <- c("# PNGS in v5", "# Cys in gp120")

dve_results <- list()
#Obtain results separately for each trial and pooled, and separately for mf, ms, ls
for(sequence in c("mf", "ms", "ls")){
  features.to.analyze.x <- paste0(features.to.analyze,".",sequence)
  for(trial in c("703", "704", "704and703")){
    for(dose in c("T1", "T2", "T1+T2")){
      
      if(trial == "704and703"){
        master.trial <- master
      }else{
        master.trial <- filter(master, protocol == paste0("HVTN ", trial) )
      }
      if(dose %in% c("T1", "T2")){
        master.trial <- filter(master.trial, tx %in% c("C3", dose))#dose-specific analysis exclude data from the other dosage
        master.trial$txLabel <- 1*(master.trial$tx == dose)
      }else{
        master.trial$txLabel <- 1*(master.trial$tx_pool == "T1+T2")
      } 
      master.trial <- filter(master.trial, !is.na(hiv1fpday))  # two 703 ppts with missing sequences have also a missing time-to-event
      
      
      # stratification variable per SAP Section 5
      master.trial <- mutate(master.trial, stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                                  protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                                  protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                                  protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))
      
     
      
      ### List of columns in output 
      results.colnames <- c( "mark", "mark.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), 
                             "DVE.p.value", 
                             paste( "VE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ),
                             paste( "VE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;
      
      # We'll store the results here:
      results <- data.frame( feature = features.to.analyze.x, feature.name = feature.names );
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
        .dve.results <- DVE.f( master.trial$hiv1fpday, master.trial$hiv1event, .x, master.trial$txLabel, master.trial$stratVar )
        if(dose == "T1+T2"){
          .DVE.p.value <- .dve.results$Diffp #Differential VE test is only done for dose pooled vs. Placebo
        }else{
          .DVE.p.value <- NA
        }
        
        
        
        # Extract VE values
        .type.1.VE <- ( 1 - summary( .dve.results$fit1 )$conf[ , -2 ] ) * 100;
        .type.1.VE[ c( 2, 3 ) ] <- .type.1.VE[ c( 3, 2 ) ];
        .type.1.VE <- c( .type.1.VE, summary( .dve.results$fit1 )$sc[ "pvalue" ] );
        .type.0.VE <- ( 1 - summary( .dve.results$fit0 )$conf[ , -2 ] ) * 100;
        .type.0.VE[ c( 2, 3 ) ] <- .type.0.VE[ c( 3, 2 ) ];
        .type.0.VE <- c( .type.0.VE, summary( .dve.results$fit0 )$sc[ "pvalue" ] );
        
        # Extract numbers of events
        .n.events.type.1.trt <- sum( .x == 1 & master.trial$txLabel == 1 & master.trial$hiv1event, na.rm = TRUE)
        .n.events.type.0.trt <- sum( .x == 0 & master.trial$txLabel == 1 & master.trial$hiv1event, na.rm = TRUE)
        .n.events.type.1.placebo <- sum( .x == 1 & master.trial$txLabel == 0 & master.trial$hiv1event, na.rm = TRUE)
        .n.events.type.0.placebo <- sum( .x == 0 & master.trial$txLabel == 0 & master.trial$hiv1event, na.rm = TRUE)
        # Obtain incidence rates (number of cases per person-years)
        .inc.1.trt <- .n.events.type.1.trt/(sum(master.trial$hiv1fpday[.x == 1 & master.trial$txLabel == 1], na.rm = TRUE)/365.25)
        .inc.0.trt <- .n.events.type.0.trt/(sum(master.trial$hiv1fpday[.x == 0 & master.trial$txLabel == 1], na.rm = TRUE)/365.25)
        .inc.1.placebo <- .n.events.type.1.placebo/(sum(master.trial$hiv1fpday[.x == 1 & master.trial$txLabel == 0], na.rm = TRUE)/365.25)
        .inc.0.placebo <- .n.events.type.0.placebo/(sum(master.trial$hiv1fpday[.x == 0 & master.trial$txLabel == 0], na.rm = TRUE)/365.25) 
        #paste the number of events and incidence rates together
        paste.n.inc <- function(n, inc){paste0(n, " (", ifelse(is.na(inc), "--", format(inc, nsmall=2, digits=2)), ")")}
        .n.events.type.1.trt <- paste.n.inc (.n.events.type.1.trt , .inc.1.trt)
        .n.events.type.0.trt <- paste.n.inc (.n.events.type.0.trt , .inc.0.trt)
        .n.events.type.1.placebo <- paste.n.inc (.n.events.type.1.placebo , .inc.1.placebo)
        .n.events.type.0.placebo <- paste.n.inc (.n.events.type.0.placebo , .inc.0.placebo)
        
        # Insert into results table 
        results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","placebo"),each=2), sep = "."), "DVE.p.value", 
                                     paste( "VE.type.1", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), 
                                     paste( "VE.type.0", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- c( .n.events.type.1.trt, 
                                                                                                                                    .n.events.type.0.trt, .n.events.type.1.placebo, .n.events.type.0.placebo, .DVE.p.value,.type.1.VE, .type.0.VE )
      }
      
      dve_results[[paste0(trial,"_",dose,".",sequence)]] <- results
      #add adjusted p-values
      if(dose == "T1+T2"){
        p.value.file <- read.csv(file.path(adjPvaluesDir, paste0(trial,"_WestfallYoungAdjPvalues_neutAssocSeqFeatures_",sequence,".csv")))
        results <- full_join(results, p.value.file, by = "mark")
        PEtable <- tibble("mark" = character(), "Haplotype" = character(), "n" = character(),
                          "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                          "P" = character(), "DiffP" = character(), "FWER" = character(), "FDR" = character())
        for( .feature in features.to.analyze.x ){
          result.f <- filter(results, mark == .feature)
          PEtable <- add_row(.data = PEtable, "mark" = result.f$mark.name, "Haplotype" = "", "n" = "",
                             "PE" = "","mean" = NA, "lower" = NA, "upper" = NA, "P" = "",
                             "DiffP" = format.p(result.f$p.unadj), "FWER" = format.p(result.f$p.FWER), "FDR" = format.p(result.f$p.FDR))
          Haplotype <- table.seqFeatLabel.tier1(paste0(strsplit(.feature, split = "\\.")[[1]][1:3], collapse="."))
          
          
          PEtable <- add_row(.data = PEtable, "mark" = "", "Haplotype" = Haplotype$Haplotype1, 
                             "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo) , 
                             "PE" = paste0(format.PE(result.f$VE.type.1.estimate)," (", format.PE(result.f$VE.type.1.CI.low),", ",format.PE(result.f$VE.type.1.CI.high),")" ),
                             "mean" = as.numeric(result.f$VE.type.1.estimate),
                             "lower" = as.numeric(result.f$VE.type.1.CI.low), #to be able to plot -Inf
                             "upper" = as.numeric(result.f$VE.type.1.CI.high),
                             "P" = format.p(result.f$VE.type.1.p.value),
                             "DiffP" = "", "FWER" = "", "FDR" = "")
          PEtable <- add_row(.data = PEtable, "mark" = "", "Haplotype" = Haplotype$Haplotype0, 
                             "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                             "PE" = paste0(format.PE(result.f$VE.type.0.estimate)," (", format.PE(result.f$VE.type.0.CI.low),", ",format.PE(result.f$VE.type.0.CI.high),")" ),
                             "mean" = as.numeric(result.f$VE.type.0.estimate),
                             "lower" = as.numeric(result.f$VE.type.0.CI.low),
                             "upper" = as.numeric(result.f$VE.type.0.CI.high),
                             "P" = format.p(result.f$VE.type.0.p.value),
                             "DiffP" = "", "FWER" = "", "FDR" = "")
          
        }
        #forest plots
        PEtable$lower[PEtable$lower == -Inf] = -110
        PEtable.forestplot.withSieveT(PEtable,xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                           figuresOutDir, paste0("HVTN",trial,"_DVE_tier1_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,".pdf"))
      
      }else{
          
        PEtable <- tibble("mark" = character(), "Haplotype" = character(), "n" = character(),
                          "PE" = character(), "mean" = numeric(), "lower" = numeric(), "upper" = numeric(),
                          "P" = character())
        for( .feature in features.to.analyze.x ){
          result.f <- filter(results, mark == .feature)
          Haplotype <- table.seqFeatLabel.tier1(paste0(strsplit(.feature, split = "\\.")[[1]][1:3], collapse="."))
          
          PEtable <- add_row(.data = PEtable, "mark" = result.f$mark.name, "Haplotype" = Haplotype$Haplotype1, 
                             "n" = paste0(result.f$.n.events.type.1.trt, " vs. ", result.f$.n.events.type.1.placebo),
                             "PE" = paste0(format.PE(result.f$VE.type.1.estimate)," (", format.PE(result.f$VE.type.1.CI.low),", ",format.PE(result.f$VE.type.1.CI.high),")" ),
                             "mean" = as.numeric(result.f$VE.type.1.estimate),
                             "lower" = as.numeric(result.f$VE.type.1.CI.low),
                             "upper" = as.numeric(result.f$VE.type.1.CI.high),
                             "P" = format.p(result.f$VE.type.1.p.value))
          PEtable <- add_row(.data = PEtable, "mark" = "", "Haplotype" = Haplotype$Haplotype0, 
                             "n" = paste0(result.f$.n.events.type.0.trt, " vs. ", result.f$.n.events.type.0.placebo) ,
                             "PE" = paste0(format.PE(result.f$VE.type.0.estimate)," (", format.PE(result.f$VE.type.0.CI.low),", ",format.PE(result.f$VE.type.0.CI.high),")" ),
                             "mean" = as.numeric(result.f$VE.type.0.estimate),
                             "lower" = as.numeric(result.f$VE.type.0.CI.low),
                             "upper" = as.numeric(result.f$VE.type.0.CI.high),
                             "P" = format.p(result.f$VE.type.0.p.value))
          
        }
        #forest plots
        PEtable$lower[PEtable$lower == -Inf] = -110
        PEtable.forestplot(PEtable,xlim.x = c(-100, 105), ticks_at.x = c(-100, -50, 0, 50, 100), 
                           figuresOutDir, paste0("HVTN",trial,"_DVE_tier1_",ifelse(dose=="T1+T2", "dosePooled", dose),"_",sequence,".pdf"))
        
        
        
        }
      
      
    }
   
   
  }
}

saveRDS( dve_results, file = file.path( tablesOutDir, "HVTN_DVE_tier1_results_by_feature.Rdata" )) 
