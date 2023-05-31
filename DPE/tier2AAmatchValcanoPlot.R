# Purpose: Volcano plot showing the difference in the log HR for the two values of the mark (as the effect size) 
# on the x-axis and –log10 differential PE unadjusted 2-sided p-value on the y-axis. 
# A dashed horizontal line could run through 0.05.
# Author:  Li Li
# Date:    Dec 23, 2022


rm(list=ls(all=TRUE))
repoDir <- "/Users/lili/AMPsieve"
adjPvaluesDir <- file.path(repoDir, "tables", "DPE/tier2")
tablesOutDir <- file.path(repoDir, "tables/DPE/tier2")
figuresOutDir <- file.path(repoDir, "figures/DPE/tier2")
outputDir <- file.path(repoDir, "code","DPE/Routput")
library(tidyverse)
source(file.path(repoDir, "code/common.R")) #read in a datFile
source(file.path(repoDir, "code/DPE/DPEutils.R")) #read in a datFile


for(sequence in c("mf", "ms", "ls")){
  for(trial in c("703", "704", "704and703")){
 
    PEtable <- read.csv(file.path(tablesOutDir, paste0("HVTN",trial,"_DVE_tier2_AAsiteScan_dosePooled","_",sequence,".csv")))
    p.value.file <- read.csv(file.path(adjPvaluesDir, paste0(trial,"_WestfallYoungAdjPvalues_tier2_",sequence,".csv")))
    plot.df <- tibble("mark" = character(),"haplotype"= character() ,"difflogHR" = numeric(), 
                      "neglog10unadjustedp" = numeric(), "unadjustedp" = character(), 
                      "fwer" = character(), "fdr" = character(),"labeltest" = character())
    for(i in 1: (length(PEtable$feature)/3)){
      resulti = PEtable[((i-1)*3+1): (i*3),]
      markPosi = as.character(resulti$mark[1])
      haplotypei = paste0(resulti$Haplotype[2], " vs. ", resulti$Haplotype[3])
      marki = paste0("hxb2.", markPosi, ".is.", resulti$Haplotype[2], ".", sequence)
      difflogHRi = log(1-resulti$mean[2]/100) -  log(1-resulti$mean[3]/100)
      pvaluesi = p.value.file[p.value.file$mark == marki, ]
      neglog10unadjustedpi = -log10(pvaluesi$p.unadj)
      plot.df <- add_row(.data = plot.df, "mark" = markPosi,"haplotype"= haplotypei ,"difflogHR" = difflogHRi, 
                         "neglog10unadjustedp" = neglog10unadjustedpi, 
                         "unadjustedp" = format.p(pvaluesi$p.unadj), "fwer" = format.p(pvaluesi$p.FWER), 
                         "fdr" = format.p(pvaluesi$p.FDR),
                         "labeltest" = paste0("Env AA Pos. ", markPosi,"; ", haplotypei,"\n",
                                              "FWER P", format.p2(pvaluesi$p.FWER),"; ",
                                              "Q", format.p2(pvaluesi$p.FDR)))
      
    }
    
    sequenceLabel <- case_when(sequence == "mf" ~ "Most Frequent Founder",
                               sequence == "ms" ~ "Predicted Most Sensitive Founder",
                               sequence == "ls" ~ "Predicted Most Resistant Founder")
    p <- ggplot()+
         geom_point(aes(x = difflogHR, y = neglog10unadjustedp), size = 4, alpha = 0.5, 
                    data = filter(plot.df, 10^(-neglog10unadjustedp) > 0.05))+
        geom_point(aes(x = difflogHR, y = neglog10unadjustedp), size = 4, alpha = 1, color = "red",
                 data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05))+
         #geom_text(aes(x = difflogHR, y = neglog10unadjustedp, label = labeltest), nudge_y = 0.25, 
          #        data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05))+
         ggrepel::geom_text_repel(aes(x = difflogHR, y = neglog10unadjustedp, label = labeltest),
                               data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05), size = ifelse(trial == "704and703", 5.1, 5.5), 
                               force = 10, seed = 0, lineheight = 0.85)+
        
         ylab("Differential PE Unadjusted 2-Sided P-value")+
         xlab("Difference in Log Hazard Ratio")+
         scale_x_continuous(limits = c(-3, 3), breaks = c(-2,-1,0,1, 2), minor_breaks = NULL)+
         scale_y_continuous(limits = c(0, ifelse(trial == "704and703", 3.05, 2)), breaks = -log10(c(1, 0.5, 0.2, 0.1,  0.05, 0.01, 0.02, 0.001)),
                            labels = c(1, 0.5, 0.2, 0.1,  0.05, 0.01, 0.02, 0.001), minor_breaks = NULL)+
         theme_bw()+
      ggtitle(sequenceLabel)+
      theme(plot.title = element_text(hjust = 0.5, vjust = 2, size = 20),
            axis.title.x = element_text(size = 19),
            axis.title.y = element_text(size = 19),
            axis.text.x = element_text(size = 18, colour = "black"),
            axis.text.y = element_text(size = 18, colour = "black"))
     
    ggplot2::ggsave(filename = paste0("tier2AAmatchValcano_", trial,"_", sequence,".pdf"), 
                    plot = p, 
                    path = figuresOutDir,
                    width=0.95*7.3, height=0.9*7.3)
     
    
}}
