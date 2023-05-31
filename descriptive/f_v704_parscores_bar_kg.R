#----------------------------------------------------------------
# PROGRAM: f_parscores_violin.R
#
# DESCRIPTION: Generate ordered bar plots of HVTN 704 par score 3 
#              values by treatment group
# 
# CODED BY: Erika Rudnicki
#
# SAVED AS: /trials/vaccine/p704/analysis/sieve/code
#
# INPUT: 
#   ../adata/v704_survival_wk80_tau_sieve_v5_cam.csv
#
# OUTPUT: 
#   ../figures/descriptive/v704barparscore3.mf.pdf
#
# MAINTENANCE HISTORY:
# Date          Programmer       Description
# 2021Jan21     Erika Rudnicki   Version 1.0
# 24Mar2022     Erika Rudnicki   Drop 'r4p' references, new source data
# 20May2022     Erika Rudnicki   Output file names with '_' instead of '.'
# 21Jun2022     Erika Rudnicki   Re-run on new input data
# 26Aug2022     Kevin Gillespie  Re-run on new input data & modified plots 
# 20Sep2022     Kevin Gillespie  Re-run with updated sap filtering
#----------------------------------------------------------------

# setwd("/trials/vaccine/p704/analysis/sieve/code")

library(dplyr)
library(ggplot2)
library(tidyr)
library(Cairo)

# load data
dat <- read.csv("../adata/amp_sieve_v704_marks_final_v4.csv")

# delete infected ppts without neut data
N_b4 <- nrow(dat)
nrow(subset(dat, hiv1event==1 & !is.na(parscore1.mf) & is.na(gmt80mf))) # Michal edit
dat <- subset(dat, !(hiv1event==1 & is.na(parscore1.mf) & !is.na(gmt80mf)))
N_aft <- nrow(dat)

N_b4 - N_aft

datplot <- 
  dat %>%
  filter(hiv1event==1) %>%
  select(pub_id, tx, contains("parscore3"), -contains("hiv1event"))

datplot2 <-
  datplot %>%
  gather(foundertype, level, contains("parscore3")) %>%
  mutate(tx2 = factor(case_when(tx=="C3" ~ "Placebo",
                                tx=="T1" ~ "VRC01 \n10 mg/kg",
                                tx=="T2" ~ "VRC01 \n30 mg/kg"), levels=c("Placebo", "VRC01 \n10 mg/kg", "VRC01 \n30 mg/kg")),
         level=factor(level, levels=c("<1", "[1,3]", ">3"))) %>%
  group_by(tx2, foundertype, level) %>%
  arrange(tx2, foundertype) %>%
  tally()

# function to plot
myplot <- function(dat, x, y, fill,
                   ylim=c(0, 40), 
                   ybreaks=c(0,5,10,15,20,25,30,35,40), 
                   ylabels=c("0","5","10","15","20","25","30","35","40"), 
                   ytitle=NULL,
                   xtitle="Treatment Arm",
                   toptitle=NULL){
  ggplot(data=dat, aes_string(x=x, y=y, fill=fill)) +
    geom_hline(yintercept=c(0,5,10,15,20,25,30,35,40), color="gray") +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_y_continuous(limits=ylim, breaks=ybreaks, labels=ylabels) +
    labs(x=xtitle, y=ytitle, title=toptitle) +
    scale_fill_manual(values = c("<1" = "#0AB7C9", 
                                 "[1,3]" = "#FF6F1B",
                                 ">3" = "#810094"), 
                      name="",
                      labels = c("<1" = "< 1 \U003BCg/ml", 
                                 "[1,3]" = "[1,3] \U003BCg/ml",
                                 ">3" = ">3 \U003BCg/ml")) +    
    theme_bw() +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          legend.position=c(0.8, 0.925),
          legend.background = element_rect(fill="transparent"),
          legend.text=element_text(size=20),
          plot.title = element_text(hjust = 0.5),
          text=element_text(size=23),
          axis.text.x = element_text(size=24),
          axis.text.y = element_text(size=24),
          axis.title.x = element_text(margin = margin(t = 10, r=0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r=10, b = 0, l = 0)))
}

plots <- c("parscore3.mf", "parscore3.ms", "parscore3.ls")
plots2 <- c("parscore3_mf", "parscore3_ms", "parscore3_ls")
plots_ytitles <- rep("Number of Primary Endpoints", 3)
plots_titles <- c("Most Frequent Founder", "Predicted Most Sensitive Founder", "Predicted Most Resistant Founder")
for (i in 1:length(plots)){
  CairoPDF(paste0("../figures/descriptive/v704bar", plots2[i], "_20092022.pdf"), width=8, height=8)
    print(myplot(dat=filter(datplot2, foundertype==plots[i]), x="tx2", y="n", fill="level", ytitle=plots_ytitles[i],toptitle=plots_titles[i]))
  dev.off()
}

print( warnings() )

q(save = "no")

