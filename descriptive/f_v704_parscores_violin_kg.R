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
# 26Aug2022     Kevin Gillespie  Re-run on new input data & modified plots 
# 20Sep2022     Kevin Gillespie  Re-run with updated sap filtering
#----------------------------------------------------------------

# setwd("/trials/vaccine/p704/analysis/sieve/code")

library(dplyr)
library(ggplot2)
library(reshape2)
options(dplyr.print_max = 100, stringsAsFactors=FALSE)

# load data
dat <- read.csv("../adata/amp_sieve_v704_marks_final_v4a.csv")

# delete infected ppts without neut data
N_b4 <- nrow(dat)
nrow(subset(dat, hiv1event==1 & !is.na(parscore1.mf) & is.na(gmt80mf))) # Michal edit
dat <- subset(dat, !(hiv1event==1 & is.na(parscore1.mf) & !is.na(gmt80mf)))
N_aft <- nrow(dat)

N_b4 - N_aft

datplot <- subset(dat, hiv1event==1)

# factor tx so boxplots ordered by control, low dose, high dose
datplot$tx2 <- with(datplot, case_when(tx=="C3" ~ "Placebo",
                                       tx=="T1" ~ "VRC01 \n10 mg/kg",
                                       tx=="T2" ~ "VRC01 \n30 mg/kg"))
datplot$tx2 <- factor(datplot$tx2, levels=c("Placebo", "VRC01 \n10 mg/kg", "VRC01 \n30 mg/kg"))

datplot <- datplot %>% 
  select(tx2, contains('pred.prob.res'), contains('pred.ic80'), contains('parscore'))


# summarize data
colnames(datplot)

# Most freq
with(datplot, range(parscore1.mf))
with(datplot, range(parscore2.mf))
with(datplot, range(pred.prob.res.mf))
with(datplot, range(pred.ic80.mf))

# Most sens
with(datplot, range(parscore1.ms))
with(datplot, range(parscore2.ms))
with(datplot, range(pred.prob.res.ms))
with(datplot, range(pred.ic80.ms))

# Least sens
with(datplot, range(parscore1.ls))
with(datplot, range(parscore2.ls))
with(datplot, range(pred.prob.res.ls))
with(datplot, range(pred.ic80.ls))

# inv. logit
logit <- function(f) {
  log(f/(1-f))
}


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
    # ylim(ylim) + 
    scale_y_continuous(trans=ytrans, limits=ylim, breaks=ybreaks, labels=ylabels) +
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


plots <- c("pred.prob.res.mf", "pred.ic80.mf", "pred.prob.res.ms", "pred.ic80.ms", "pred.prob.res.ls", "pred.ic80.ls")
plots2 <- c("pred_prob_res_mf", "pred_ic80_mf","pred_prob_res_ms", "pred_ic80_ms" ,"pred_prob_res_ls", "pred_ic80_ls")
plots_ytitles <- rep(c(expression("Predicted Probability of IC80 > 1 "*mu*"g/mL"), expression("Predicted IC80 ("*mu*"g/mL)")),3)
plots_titles <- c(rep(c("Most Frequent Founder"), 2), rep(c("Predicted Most Sensitive Founder"), 2), rep(c("Predicted Most Resistant Founder"), 2))
plots_ylim=list('prob' = c(0.2, 0.95), 'ic80' = c(0.5, 100))
plots_ybreaks=list('prob' = c(0.3, 0.5, 0.7, 0.9), 'ic80' = c(0.1, 1, 10, 100)) #c(0.63, 1, 2, 3, 5, 10, 15))
plots_ylabels=list('prob' = c("0.3", "0.5", "0.7", "0.9"), 'ic80' = c('0.1', '1', '10', '100')) #c("0.63", "1", "2", "3", "5", "10", "15"))
plots_ytrans <- list('prob' = "logit", 'ic80' = "log10")

y_list <- list(plots_ylim, plots_ybreaks, plots_ylabels, plots_ytrans)

for (i in 1:length(plots)){
  
  val <- ifelse(i %in% seq(1, length(plots), 2), 'prob', 'ic80')
  
  if (val == "prob") {
    ylim <- y_list[[1]][[val]]
    ybreaks <- y_list[[2]][[val]]
    # ylim <- sort(1/(1+exp(c(-3.66, 3.66))))
  } else {
    ylim <- y_list[[1]][[val]]
    ybreaks <- y_list[[2]][[val]]
  }
  
  pdf(paste0("../figures/descriptive/v704violin", plots2[i], "_20092022.pdf"), width=8, height=8)
    print(myplot(dat=datplot, 
                 x="tx2", 
                 y=plots[i], 
                 ytitle=plots_ytitles[i],
                 toptitle=plots_titles[i], 
                 ylim=ylim,
                 ybreaks=ybreaks,
                 ylabels=y_list[[3]][[val]],
                 ytrans=y_list[[4]][[val]]))
  dev.off()
}

print( warnings() )

q(save = "no")

