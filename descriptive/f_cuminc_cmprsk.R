#----------------------------------------------------------------------------------------
# PROGRAM: f_cuminc_cmprsk.R
#
# DESCRIPTION: plot cumulative incidence based on the following parameters
#              (1) HVTN 704 
#              (2) scale of time since enrollment through week 80 visit
#              (3) pooled VRC01 vs. control
#              (4) censoring at tau
#              (5) competing risks analysis by PAR score 3
#              (6) for cases, endpoint is first positive date rather than diagnosis date
#
# CODED BY: Erika Rudnicki
#
# INPUT:  
#      ../adata/v704_survival_wk80_tau_sieve_v5_cam_er.csv
#      ../adata/v704_cmprisk_cir.csv
#
# OUTPUT: 
#     ../figures/cmprsk/v704_cuminc_cmprsk.pdf
#
# MAINTENANCE HISTORY:
# Date           Programmer         Description
# ***************************************************************************************
# 21Jan2021      Erika Rudnicki     Version 1.0
# 24Mar2022      Erika Rudnicki     Drop 'r4p' references, new source data
# 21Jun2022      Erika Rudnicki     Re-run on new input data
#----------------------------------------------------------------------------------------

# setwd("/trials/vaccine/p704/analysis/efficacy/code")

library(dplyr)
library(ggplot2)
library(grid)
library(Cairo)

# input directory and file names
adataDir <- "../adata"
adataFile1 <- "v704_survival_wk80_tau_sieve_v5_cam_er.csv"
adataFile2 <- "v704_cmprisk_cir.csv"

# output directory and file names
pdfDir <- "../figures/cmprsk"
pdfFile <- c("v704_cuminc_cmprsk.pdf")

pdfFileSave <- file.path(pdfDir, pdfFile)

# source input data
surv <- 
  read.csv(file.path(adataDir, adataFile1), stringsAsFactors = FALSE)

# delete infected ppts without neut data
nrow(surv)
subset(surv, hiv1event==1 & nisolates==0)
surv <- subset(surv, !(hiv1event==1 & nisolates==0))
nrow(surv)

cir <- 
  read.csv(file.path(adataDir, adataFile2), stringsAsFactors = FALSE) %>%
  mutate(eventType = factor(eventType, levels=c("< 1", "1 to 3", "> 3")))

# reformat 'CIR' data for plotting, taking steps to ensure plots
# end at tau (our censoring time point) rather than the last event time
(tau <- max(surv$hiv1fpday))
(tau_wks <- tau/7)

cuminc_v <-
  cir %>%
  group_by(eventType) %>%
  mutate(time = if_else(time==max(time), tau, time),
         level = cmpLvl,
         cuminc = cuminc.cmp) %>%
  select(eventType, time, level, cuminc)

cuminc_p <-
  cir %>%
  group_by(eventType) %>%
  mutate(time = if_else(time==max(time), tau, time),
         level = refLvl,
         cuminc = cuminc.ref) %>%
  select(eventType, time, level, cuminc)

cuminc_plot <- 
  bind_rows(cuminc_v, cuminc_p) %>%
  mutate(level = factor(level, levels=c("T1+T2", "C3")),
         time_wks = time/7)
    
# here is where we could create footnotes using 'surv' data
# but will omit for now because plots going in the report
# do not have room for the large number of rows we'd need to create
# ...
    
# define font sizes (larger if w/o footnotes)
font_size=32

# plot Pooled VRC01 vs. Control w/o footnotes
p <-
  
  ggplot(data=cuminc_plot, aes(x=time_wks, y=cuminc, color=eventType, linetype=level)) +
  scale_x_continuous(name="Weeks since Enrollment",
                     limits=c(0, 88),
                     breaks = seq(0, 88, by=16),
                     labels = seq(0, 88, by=16)) +
  scale_y_continuous(name="Cumulative Incidence (%)",
                     limits=c(0, 0.0425),
                     breaks=seq(0, 0.04, 0.01),
                     labels=seq(0, 4, 1)) +
  geom_step(lwd=2) +
  scale_color_manual(values = c("< 1" = "#0AB7C9", 
                                "1 to 3" = "#FF6F1B",
                                "> 3" = "#810094"), 
                     name="",
                     labels = c("< 1" = "< 1 \U003BCg/ml", 
                                "1 to 3" = "[1, 3] \U003BCg/ml",
                                "> 3" = "> 3 \U003BCg/ml")) +
  scale_linetype_manual(values = c("T1+T2" = "solid", "C3" = "dotted"),
                        name="",
                        labels = c("T1+T2" = "VRC01 Pooled", "C3" = "Control")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=font_size, margin=margin(t=20)),
        axis.title.y = element_text(size=font_size, margin=margin(r=20)),
        axis.text.x = element_text(size=font_size, margin=margin(t=10)),
        axis.text.y = element_text(size=font_size, margin=margin(r=10)),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title = element_text(hjust=0.5, size=font_size, face="bold"),
        plot.margin = unit(c(1,1,1,1), "cm"), #trbl
        legend.position = c(0.5, 0.925),
        legend.text=element_text(size=font_size*.75),
        legend.key.width = unit(2,"cm"),
        legend.box = "horizontal") +
  labs(title="Most Frequent Founders") +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))

ggsave(pdfFileSave[1], plot = p, width=8.5, height=8.5, device=cairo_pdf)

q(save = "no")
