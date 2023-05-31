#----------------------------------------------------------------
# PROGRAM: f_parscores_violin.R
#
# DESCRIPTION: Generate violin plot of HVTN 704 parscore1 and 
#              parscore2 values on log scale by treatment group
# 
# CODED BY: Erika Rudnicki
#
# SAVED AS: /trials/vaccine/p704/analysis/sieve/code
#
# INPUT: 
#   ../adata/v704_survival_wk80_tau_sieve_v5_cam.csv
#
# OUTPUT: 
#   ../figures/descriptive/v704_violin_pred.prob.res.mf.pdf
#   ../figures/descriptive/v704_violin_pred.ic80.mf.pdf
#
# MAINTENANCE HISTORY:
# Date          Programmer       Description
# 2021Jan19     Erika Rudnicki   Version 1.0
# 2021Jan22     Erika Rudnicki   Add boxplots within 
# 24Mar2022     Erika Rudnicki   Drop 'r4p' references, new source data
# 20May2022     Erika Rudnicki   Output file names with '_' instead of '.'
# 21Jun2022     Erika Rudnicki   Re-run on new input data
#----------------------------------------------------------------

# setwd("/trials/vaccine/p704/analysis/sieve/code")

library(dplyr)
library(ggplot2)
options(dplyr.print_max = 100, stringsAsFactors=FALSE)

# load data
dat <- read.csv("../adata/v704_survival_wk80_tau_sieve_v5_cam_er.csv")

# delete infected ppts without neut data
nrow(dat)
subset(dat, hiv1event==1 & nisolates==0)
dat <- subset(dat, !(hiv1event==1 & nisolates==0))
nrow(dat)

datplot <- subset(dat, hiv1event==1)

# factor tx so boxplots ordered by control, low dose, high dose
datplot$tx2 <- with(datplot, case_when(tx=="C3" ~ "Control",
                                       tx=="T1" ~ "VRC01 \n10 mg/kg",
                                       tx=="T2" ~ "VRC01 \n30 mg/kg"))
datplot$tx2 <- factor(datplot$tx2, levels=c("Control", "VRC01 \n10 mg/kg", "VRC01 \n30 mg/kg"))

# summarize data
colnames(datplot)
with(datplot, range(parscore1.mf))
with(datplot, range(parscore2.mf))
with(datplot, range(pred.prob.res.mf))
with(datplot, range(pred.ic80.mf))

# function to plot
myplot <- function(dat, x, y, 
                   ylim=NULL,
                   ybreaks=NULL,
                   ylabels=NULL,
                   ytitle=NULL,
                   ytrans=NULL,
                   xtitle="Treatment Group",
                   toptitle=NULL,
                   boxcol=c("blue", "red3", "red3")){
  set.seed(39573056)
  ggplot(data=dat, aes_string(x=x, y=y, color=x)) +
    geom_violin() +
    geom_boxplot(fill=boxcol, width=0.15, lwd=1.5, alpha = 0.3, outlier.shape=NA) + #remove outlier points
    geom_point(size=3, position = position_jitter(w = 0.3, h = 0)) +
    scale_y_continuous(trans=ytrans) + #limits=ylim, breaks=ybreaks, labels=ylabels, 
    labs(x=xtitle, y=ytitle, title=toptitle) +
    scale_color_manual(values=boxcol) +
    theme_bw() +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          legend.position="none",
          plot.title = element_text(hjust = 0.5),
          text=element_text(size=25),
          axis.text.x = element_text(size=24),
          axis.text.y = element_text(size=24),
          axis.title.x = element_text(margin = margin(t = 10, r=0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r=10, b = 0, l = 0)))
}

plots <- c("pred.prob.res.mf", "pred.ic80.mf")
plots2 <- c("pred_prob_res_mf", "pred_ic80_mf")
plots_ytitles <- c(expression("Predicted Probability of IC80 > 1 "*mu*"g/mL"), expression("Predicted IC80 ("*mu*"g/mL)"))
plots_titles <- rep(c("Most Frequent Founder"), 2)
plots_ylim=c(c(0.1, 1), c(0.5, 15.2))
plots_ybreaks=c(c(0.3, 0.5, 0.7, 0.9), c(0.63, 1, 2, 3, 5, 10, 15)) 
plots_ylabels=c(c("0.3", "0.5", "0.7", "0.86"), c("0.63", "1", "2", "3", "5", "10", "15"))
plots_ytrans <- c("logit", "log10")

for (i in 1:length(plots)){
  pdf(paste0("../figures/descriptive/v704violin", plots2[i], ".pdf"), width=8, height=8)
    print(myplot(dat=datplot, 
                 x="tx2", 
                 y=plots[i], 
                 ytitle=plots_ytitles[i],
                 toptitle=plots_titles[i], 
                 ylim=plots_ylim[i],
                 ybreaks=plots_ybreaks[i],
                 ylabels=plots_ylabels[i],
                 ytrans=plots_ytrans[i]))
  dev.off()
}

print( warnings() )

q(save = "no")

