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
# 14Oct2022     Michal Juraska   Revised for manuscript
#----------------------------------------------------------------

rm(list=ls(all=TRUE))

datDir <- paste0("t:/vaccine/p704/analysis/sieve/adata")
figDir <- paste0("t:/vaccine/p704/analysis/sieve/figures/manuscript/descriptive")

library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)
options(dplyr.print_max = 100, stringsAsFactors=FALSE)

source("t:/vaccine/p704/analysis/sieve/code/common.R")

logit <- function(f) {
  log(f/(1-f))
}


# Plotting functions ------------------------------------------------------

myplot <- function(dat, x, y,
                   ylim=NULL,
                   ybreaks=NULL,
                   ylabels=NULL,
                   ytitle=NULL,
                   ytrans=NULL,
                   xtitle=NULL,
                   toptitle=NULL,
                   boxcol=c("blue", "red3", "red3")){
  set.seed(39573056)
  ggplot(data=dat, aes_string(x=x, y=y, color=x)) +
    geom_violin() +
    geom_boxplot(fill=boxcol, width=0.15, lwd=1.5, alpha=0.3, outlier.shape=NA) + #remove outlier points
    geom_point(size=3, position=position_jitter(w=0.3, h=0)) +
    scale_y_continuous(trans=ytrans, breaks=ybreaks, labels=ylabels) +
    coord_cartesian(ylim=ylim) +
    labs(x=xtitle, y=ytitle, title=toptitle) +
    scale_color_manual(values=boxcol) +
    scale_shape_identity() +
    theme_bw() +
    theme(plot.margin = unit(c(0.1,0.2,0.1,0.1), "in"),
          legend.position="none",
          plot.title = element_text(hjust = 0.5),
          text=element_text(size=25),
          axis.text.x = element_text(size=24, lineheight=0.7),
          axis.text.y = element_text(size=24)
          #axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
          #axis.title.y = element_text(margin=margin(t=0, r=10, b=0, l=0))
          )
}

myplot2 <- function(dat, x, y,
                    ylim=NULL,
                    ybreaks=NULL,
                    ylabels=NULL,
                    ytitle=NULL,
                    ytrans=NULL,
                    xtitle=expression("Measured" ~ IC[80] ~ "(" * mu * "g/ml)"),
                    toptitle=NULL,
                    boxcol=c("blue", "red3", "red3")){
  ggplot(data=dat, aes_string(x=x, y=y)) +
    #stat_smooth(method="loess", color="black", size=2, se=FALSE) +
    stat_smooth(method=lm, color="black") +
    geom_point(size=3) +
    scale_x_continuous(trans="log10", breaks=c(0.1, 0.3, 1, 3, 10, 30, 100), 
                       labels=c("0.1", "0.3", "1", "3", "10", "30", expression(" " >= "100"))) +
    scale_y_continuous(trans=ytrans, breaks=ybreaks, labels=ylabels) +
    coord_cartesian(xlim=c(0.05, 110), ylim=ylim) +
    labs(x=xtitle, y=ytitle, title=toptitle) +
    #scale_color_manual(values=boxcol) +
    annotate("text", x=0, y=1, hjust=-0.3, vjust=2, size=7, 
             label=as.character(as.expression(substitute("corr"==r, 
                                                         list(r=format(cor(dat[, x], dat[, y], method="spearman"), nsmall=2, digits=2))))), 
             parse=TRUE) +
    theme_bw() +
    theme(plot.margin = unit(c(0.1,0.2,0.1,0.1), "in"),
          legend.position="none",
          plot.title = element_blank(),
          text=element_text(size=25),
          axis.text.x = element_text(size=24),
          axis.text.y = element_text(size=24),
          axis.title.x = element_text(margin=margin(t=4, r=0, b=0, l=0))
          #axis.title.y = element_text(margin=margin(t=0, r=5, b=0, l=0))
          )
}


# Set plotting parameters -------------------------------------------------

trial <- c("704", "703")

plots <- c("pred.prob.res.mf", "pred.ic80.mf", "pred.prob.res.ms", "pred.ic80.ms", "pred.prob.res.ls", "pred.ic80.ls")
plots2 <- c("pred_prob_res_mf", "pred_ic80_mf","pred_prob_res_ms", "pred_ic80_ms" ,"pred_prob_res_ls", "pred_ic80_ls")
plots_ytitles <- rep(c(expression("Predicted Probability of" ~ IC[80] ~ "> 1" ~ mu * "g/ml"), expression("Predicted" ~ IC[80] ~ "(" * mu * "g/ml)")),3)
plots_titles <- c("HVTN 704/HPTN 085", "HVTN 703/HPTN 081")
plots_ylim=list('prob' = c(0.25, 0.95), 'ic80' = c(0.5, 100))
plots_ybreaks=list('prob' = c(0.3, 0.5, 0.7, 0.9), 'ic80' = c(0.1, 1, 10, 100)) #c(0.63, 1, 2, 3, 5, 10, 15))
plots_ylabels=list('prob' = c("0.3", "0.5", "0.7", "0.9"), 'ic80' = c('0.1', '1', '10', '100')) #c("0.63", "1", "2", "3", "5", "10", "15"))
plots_ytrans <- list('prob' = "logit", 'ic80' = "log10")

y_list <- list(plots_ylim, plots_ybreaks, plots_ylabels, plots_ytrans)



# Plot single-panel figures -----------------------------------------------

ggplotList1 <- vector("list", length=2 * length(plots))
names(ggplotList1) <- c(t(outer(trial, plots, FUN=paste, sep=".")))


plots.ic80 <- rep(c("gmt80mf", "gmt80ms", "gmt80ls"), each=2)
plots2.ic80 <- rep(c("ic80_mf", "ic80_ms", "ic80_ls"), each=2)

ggplotList2 <- vector("list", length=2 * length(plots))
names(ggplotList2) <- c(t(outer(trial, plots, FUN=paste, sep=".")))

for (j in 1:2){
  # load data
  datplot <- read.csv(file.path(datDir, datFile[j])) %>%
    filter(hiv1event==1 & !is.na(pred.prob.res.mf)) %>%
    mutate(tx2=case_when(tx=="C3" ~ "Placebo",
                         tx=="T1" ~ "VRC01 \n10 mg/kg",
                         tx=="T2" ~ "VRC01 \n30 mg/kg"),
           gmt80mf=as.numeric(ifelse(gmt80mf==">100", "100", gmt80mf)),
           gmt80ms=as.numeric(ifelse(gmt80ms==">100", "100", gmt80ms)),
           gmt80ls=as.numeric(ifelse(gmt80ls==">100", "100", gmt80ls))) %>%
    select(tx2, contains("pred.prob.res"), contains("pred.ic80"), contains("gmt80"))
  
  # factor tx so boxplots ordered by control, low dose, high dose
  datplot$tx2 <- factor(datplot$tx2, levels=c("Placebo", "VRC01 \n10 mg/kg", "VRC01 \n30 mg/kg"))
  
  
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
    
    # scatter/box/violin plot of PAR score
    pdf(file=file.path(figDir, paste0(trial[j], "_violin_", plots2[i], ".pdf")), width=8, height=8)
    print(ggplotList1[[length(plots) * (j - 1) + i]] <- myplot(dat=datplot, 
                                                              x="tx2", 
                                                              y=plots[i], 
                                                              ytitle=if (j==1) plots_ytitles[i] else NULL,
                                                              toptitle=plots_titles[j], 
                                                              ylim=ylim,
                                                              ybreaks=ybreaks,
                                                              ylabels=y_list[[3]][[val]],
                                                              ytrans=y_list[[4]][[val]]))
    dev.off()
    
    # bivariate scatter plot of PAR score vs. IC80
    datplot2 <- filter(datplot, !is.na(gmt80mf))
    
    pdf(file=file.path(figDir, paste0(trial[j], "_scatter_", plots2[i], "_", plots2.ic80[i], ".pdf")), width=8, height=8)
    print(ggplotList2[[length(plots) * (j - 1) + i]] <- myplot2(dat=datplot2,
                                                               x=plots.ic80[i],
                                                               y=plots[i],
                                                               ytitle=if (j==1) plots_ytitles[i] else NULL,
                                                               toptitle=NULL,
                                                               ylim=ylim,
                                                               ybreaks=ybreaks,
                                                               ylabels=y_list[[3]][[val]],
                                                               ytrans=y_list[[4]][[val]]))
    dev.off()
  }
}


# Combine panels into multi-panel manuscript figures ----------------------

pdf(file=file.path(figDir, "fig2_predProbResPARscore_ls.pdf"), width=13.4, height=13)
ggplotList1[["704.pred.prob.res.ls"]] + 
  ggplotList1[["703.pred.prob.res.ls"]] + 
  ggplotList2[["704.pred.prob.res.ls"]] +
  ggplotList2[["703.pred.prob.res.ls"]] +
  plot_layout(ncol=2, byrow=TRUE)
  #plot_annotation(tag_levels="A")
dev.off()


# multiPanelPlot <- ggpubr::ggarrange(plotlist=c(ggplotList1[c("704.pred.prob.res.mf", "703.pred.prob.res.mf")],
#                                                ggplotList2[c("704.pred.prob.res.mf", "703.pred.prob.res.mf")]), nrow=2, ncol=2, align="hv")
# ggpubr::ggexport(multiPanelPlot, filename=file.path(figDir, "fig2_predProbResPARscore_mf.pdf"), width=16, height=16)
