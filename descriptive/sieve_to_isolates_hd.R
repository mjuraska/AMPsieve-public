# Purpose: Descriptive plots for the Hamming distance in primary endpoints between the most resistant synthesized sequence
#for measuring TZM-bl in vitro neutralization sensitivity and the predicted most resistant founder sequence 
#used for predicting TZM-bl in vitro neutralization sensitivity
# Author:  Li Li
# Date:    1/25, 2023

rm(list=ls(all=TRUE))
repoDir <- "/Users/lili/AMPsieve"
figuresOutDir <- file.path(repoDir, "figures/descriptive")
source(file.path(repoDir, "code/common.R")) #read in a datFile
hd_data <-  read.csv( file.path(repoDir,"adata", sieve_to_isolates_hd))
master <-  read.csv( file.path(repoDir,"adata", datFile[3])) 

library(tidyverse)

hd_data <- mutate(hd_data, HD.status = case_when (!HD %in% c("no matching ptid", "No seq to calc HD") ~ 1,
                                                           HD %in% c("no matching ptid", "No seq to calc HD") ~ 0))

hd_data <- mutate(hd_data, HD.marks.status = case_when (!HD.marks %in% c("no matching ptid", "No seq to calc HD") ~ 1,
                                                        HD.marks %in% c("no matching ptid", "No seq to calc HD") ~ 0))

#restrict to participants with sequence and IC80 data
parID <- filter(master, hiv1event ==1 & !is.na(hxb2.142.T.ls) & !is.na(gmt80ls))$pub_id
data <- filter(hd_data, hiv1event==1 & pub_id%in%parID)

table(data$protocol, data$HD.status ==0)
table(data$protocol, data$HD.marks.status==0)




for(sites in c("env", "neut")){
  if(sites == "env"){
    plot_data <- filter(data, HD.status == 1)
    markLab <- "Hamming Distance in Env"
    plot_data$mark <- as.numeric(plot_data$HD)
  }else{
    plot_data <- filter(data, HD.marks.status == 1)
    markLab <- "Hamming Distance in\nNeutralization-Associated Env Positions"
    plot_data$mark <- as.numeric(plot_data$HD.marks)
  }
  
  plot_data$protocol <- factor(plot_data$protocol, levels = c("HVTN 704", "HVTN 703"))
  set.seed(100)
  #scatter plot for Hamming distance in Env and the union of sites predicting neutralization
  
  p <- ggplot(plot_data) +
    geom_boxplot(aes(x = protocol, y = mark), width = 0.8, fill = "gray80", color = "black", lwd = 0.8, 
                 outlier.shape = NA, na.rm = TRUE) + 
    geom_jitter(aes(x = protocol, y = mark, color = as.factor(protocol)), alpha = 1, 
                size = 2,
                height = 0, fill = "white", stroke = 1, na.rm = TRUE)
  p <- p + scale_color_manual(name="", values = c("red3", "blue"), breaks = c("HVTN 703", "HVTN 704"))+
    guides(color = "none") 
  
  p <- p + xlab("") +
    ylab(markLab) +
    ylim(c(0,160))+
    scale_x_discrete(breaks = c("HVTN 703", "HVTN 704"),
                     labels = c("HVTN 703\nHPTN 081", "HVTN 704\nHPTN 085"))+
    theme_bw() +
    theme(axis.title.x = element_text(size = 19, margin = margin(t = 10)),
          axis.title.y = element_text(size = 19, margin = margin(r = 10)),
          axis.text.x = element_text(size = 19, colour = "black"),
          axis.text.y = element_text(size = 19, colour = "black"))
  ggplot2::ggsave(filename = paste0("hd_", sites, "_scatterplot.pdf"), 
                  plot = p, 
                  path = figuresOutDir,
                  width=6, height=6)  
  
  #RCDF plot 
  p <- ggplot(plot_data)+
    aes(x = mark, y = 1-stat(y), colour=protocol) + 
    stat_ecdf(alpha=0.8, geom = "step", size = 1.5, pad = "TRUE")+
    xlim(c(0,160))+
    scale_color_manual(name="", values = c("blue", "red3"), breaks = c("HVTN 704", "HVTN 703"),
                       labels = c("HVTN 704/HPTN 085", "HVTN 703/HPTN 081"))+
    ylab("RCDF") +
    xlab(markLab) +
    theme_bw() +
    theme(legend.position = c(0.5,0.8 ),
          legend.key.size = unit(0.45, "cm"),
          legend.title=element_blank(),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(size=19),
          legend.key = element_blank(),
          legend.key.width = unit(0.5,"cm"),
          legend.background = element_blank(),
          axis.title.x = element_text(size = 19, margin = margin(t = 10)),
          axis.title.y = element_text(size = 19, margin = margin(r = 10)),
          axis.text.x = element_text(size = 19, colour = "black"),
          axis.text.y = element_text(size = 19, colour = "black"))
  ggplot2::ggsave(filename = paste0("hd_", sites, "_rcdf.pdf"), 
                  plot = p, 
                  path = figuresOutDir,
                  width=6, height=6) 
}



   
