#----------------------------------------------------------------------
# PROGRAM: f_violin_swg.R
#
# DESCRIPTION: Quickly extend work in f_parscores_violin_r4p to stratify
#              by D230 match vs. mismatch for SWG meeting tomorrow
# 
# CODED BY: Erika Rudnicki
#
# SAVED AS: /trials/vaccine/p704/analysis/sieve/code
#
# INPUT: 
#   ../adata/v704_survival_wk80_tau_sieve_r4p_v2_cam.csv
#   ../tables/DVE/HVTN704_site_scanning_results_DVE.csv
#
# OUTPUT: 
#   ../figures/descriptive/v704_violin_d230_gmt80mf.pdf
#
# MAINTENANCE HISTORY:
# Date          Programmer       Description
# 2021Jun01     Erika Rudnicki   Version 1.0, copied from 
#  /trials/vaccine/p704/analysis/sieve/code/f_parscores_violin_r4p.R
# 2022Sep09     Kevin Gillespie  Updated to use amp_sieve_v704_marks_final_v3b.csv
#                                  and changed the PE estimated per email from Michal
#                                  sent 2022Sep08. 
# 20Sep2022     Kevin Gillespie  Re-run with updated sap filtering
#----------------------------------------------------------------------

# setwd("/trials/vaccine/p704/analysis/sieve/code")

library(dplyr)
library(ggplot2)
options(dplyr.print_max = 100, stringsAsFactors=FALSE)


# load data
dat1 <- read.csv("../adata/amp_sieve_v704_marks_final_v4a.csv")
dat2 <- read.csv("../tables/DVE/ver 2021 May 24/HVTN704_site_scanning_results_DVE.csv")

colnames(dat1)
colnames(dat2)

dat1 %>% group_by(hiv1event, nisolates, hxb2.230.D.ls) %>% tally()

datplot <- subset(dat1, hiv1event==1 & nisolates>=1)

# prep facet labels
dat2sub <-
  dat2 %>%
  filter(X=="hxb2.230.D.ls") %>%
  mutate(pe_match = paste0("PE = ", round(VE.vs.Match.estimate), "% (", round(VE.vs.Match.CI.low), ", ", round(VE.vs.Match.CI.high), ")"),
         pe_mismatch = paste0("PE = ", round(VE.vs.Mismatch.estimate), "% (", round(VE.vs.Mismatch.CI.low), ", ", round(VE.vs.Mismatch.CI.high), ")"))
dat2sub$pe_match
dat2sub$pe_mismatch

# factor tx so boxplots ordered by control, low dose, high dose
datplot$tx2 <- with(datplot, case_when(tx=="C3" ~ "Placebo",
                                       tx=="T1" ~ "VRC01 \n10 mg/kg",
                                       tx=="T2" ~ "VRC01 \n30 mg/kg"))
datplot$tx2 <- factor(datplot$tx2, levels=c("Placebo", "VRC01 \n10 mg/kg", "VRC01 \n30 mg/kg"))

# datplot$facetvar <- with(datplot, case_when(hxb2.230.D.mf==1 ~ paste0("D230 Match \n", dat2sub$pe_match),
#                                             hxb2.230.D.mf==0 ~ paste0("D230 Mismatch \n", dat2sub$pe_mismatch)))
# (09SEP2022) KG: At the request of Michal these values are being recoded with udpated data from email 
# sent 08SEP2022. New values should read:
# - pe_match = '-23 (-116, 30)'
# - pe_mismatch = '64 (27, 82)'
datplot$facetvar <- with(datplot, case_when(hxb2.230.D.ls==1 ~ paste0("D230 Match \n", 'PE = -10 (-90, 36)'),
                                            hxb2.230.D.ls==0 ~ paste0("D230 Mismatch \n", 'PE = 53 (10, 75)')))

datplot$facetvar2 <- with(datplot, case_when(hxb2.230.D.ls==1 ~ paste0("D230 Match"),
                                             hxb2.230.D.ls==0 ~ paste0("D230 Mismatch")))

datplot %>% group_by(facetvar, facetvar2, tx2) %>% tally()

# summarize data
datplot$gmt80lsn = as.numeric(gsub(">", "", datplot$gmt80ls))
range(datplot$gmt80lsn)

# function to plot
myplot <- function(dat, x, y, facetvar,
                   ylim=NULL,
                   ybreaks=NULL,
                   ylabels=NULL,
                   ytitle=NULL,
                   ytrans=NULL,
                   xtitle="Treatment Group",
                   toptitle=NULL,
                   boxcol=c("blue", "red3", "red3", "blue", "red3", "red3")){
  set.seed(39573056)
  ggplot(data=dat, aes_string(x=x, y=y, color=x)) +
    facet_grid(cols = vars(facetvar)) +
    geom_violin() +
    geom_boxplot(fill=boxcol, width=0.15, lwd=1.25, alpha = 0.3, outlier.shape=NA) + #remove outlier points
    geom_point(size=2, position = position_jitter(w = 0.3, h = 0)) +
    scale_y_continuous(breaks=ybreaks, labels=ylabels, trans="log10") +
    labs(x=xtitle, y=ytitle, title=toptitle) +
    scale_color_manual(values=boxcol) +
    theme_bw() +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          legend.position="none",
          plot.title = element_text(hjust = 0.5),
          text=element_text(size=23),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(margin = margin(t = 10, r=0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r=10, b = 0, l = 0)))
}

plots <- c("gmt80lsn")
plots_ytitles <- expression("IC80 ("*mu*"g/mL)")
# plots_titles <- paste0("Differential p = ", round(dat2sub$DVE.p.value, 3))
plots_titles <- ''
plots_ybreaks=c(0.03, 0.1, 0.3, 1, 3, 10, 30, 100)
plots_ylabels=c("0.03", "0.1", "0.3", "1", "3", "10", "30", expression("">="100"))

for (i in 1:length(plots)){
  pdf(paste0("../figures/descriptive/v704_violin_d230_gmt80ls", "14OCT2022.pdf"), width=8, height=6)
    print(myplot(dat=datplot, 
                 x="tx2", 
                 y=plots[i], 
                 facetvar="facetvar",
                 ytitle=plots_ytitles[i],
                 toptitle=plots_titles[i], 
                 ylim=plots_ylim[i],
                 ybreaks=plots_ybreaks,
                 ylabels=plots_ylabels,
                 ytrans=plots_ytrans[i]))
  dev.off()
}

# function to plot
myplot2 <- function(dat, x, y, facetvar2,
                   ylim=NULL,
                   ybreaks=NULL,
                   ylabels=NULL,
                   ytitle=NULL,
                   ytrans=NULL,
                   xtitle="Treatment Group",
                   toptitle=NULL,
                   boxcol=c("blue", "red3", "red3", "blue", "red3", "red3")){
  set.seed(39573056)
  ggplot(data=dat, aes_string(x=x, y=y, color=x)) +
    facet_grid(cols = vars(facetvar2)) +
    geom_violin() +
    geom_boxplot(fill=boxcol, width=0.15, lwd=1.25, alpha = 0.3, outlier.shape=NA) + #remove outlier points
    geom_point(size=2, position = position_jitter(w = 0.3, h = 0)) +
    scale_y_continuous(breaks=ybreaks, labels=ylabels, trans="log10") +
    labs(x=xtitle, y=ytitle, title=toptitle) +
    scale_color_manual(values=boxcol) +
    theme_bw() +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          legend.position="none",
          plot.title = element_text(hjust = 0.5),
          text=element_text(size=25),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(margin = margin(t = 10, r=0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r=10, b = 0, l = 0)))
}

for (i in 1:length(plots)){
  pdf(paste0("../figures/descriptive/v704_violin_d230_gmt80ls_v2", "14OCT2022.pdf"), width=8, height=6)
  print(myplot2(dat=datplot, 
               x="tx2", 
               y=plots[i], 
               facetvar="facetvar2",
               ytitle=plots_ytitles[i],
               toptitle=" ", 
               ylim=plots_ylim[i],
               ybreaks=plots_ybreaks,
               ylabels=plots_ylabels,
               ytrans=plots_ytrans[i]))
  dev.off()
}

print( warnings() )

q(save = "no")

