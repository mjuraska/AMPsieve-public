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


# Plotting function ------------------------------------------------------

myplot <- function(dat, x, y,
                   ylim=NULL,
                   ybreaks=NULL,
                   ylabels=NULL,
                   ytitle=NULL,
                   xtitle=NULL,
                   toptitle=NULL,
                   boxcol=c("blue", "red3", "red3"),
                   plotMargin=c(0.1,0.2,0.1,0.1),
                   axis.text.x=TRUE){
  set.seed(39573056)
  ggplot(data=dat, aes_string(x=x, y=y, color=x)) +
    geom_violin() +
    geom_boxplot(fill=boxcol, width=0.15, lwd=1.5, alpha=0.3, outlier.shape=NA) + #remove outlier points
    geom_point(size=3, position=position_jitter(w=0.3, h=0)) +
    scale_y_continuous(breaks=ybreaks, labels=ylabels) +
    coord_cartesian(ylim=ylim) +
    labs(x=xtitle, y=ytitle, title=toptitle) +
    scale_color_manual(values=boxcol) +
    scale_shape_identity() +
    theme_bw() +
    theme(plot.margin = unit(plotMargin, "in"),
          legend.position="none",
          plot.title = element_text(hjust = 0.5),
          text=element_text(size=25),
          axis.text.x = if (axis.text.x) element_text(size=24, lineheight=0.7) else element_blank(),
          axis.text.y = element_text(size=24)
          #axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
          #axis.title.y = element_text(margin=margin(t=0, r=10, b=0, l=0))
    )
}


# Plot single-panel figures -----------------------------------------------

datPooled <- read.csv(file.path(datDir, datFile[3]))

trial <- c("704", "703")

ggplotList <- vector("list", length=8)
# names(ggplotList) <- c(t(outer(trial, plots, FUN=paste, sep=".")))


for (j in 1:2){
  plots <- c("parscore1.ls", "parscore2.ls", paste0("epitope.dist.", ifelse(j==1, "b", "c"), ".ls"), "hdist.zspace.sites.binding.all.ls")
  plots2 <- gsub("\\.", "\\_", plots)
  plots_ytitles <- c(bquote(atop("Predicted Probability of", IC[80] ~ "> 1" ~ mu * "g/ml")), 
                     expression("Predicted" ~ IC[80] ~ "(" * mu * "g/ml)"),
                     "VRC01 Epitope Distance",
                     "PC-Weighted Hamming Distance in\n50 VRC01 or CD4 Binding Positions")
  plots_titles <- c("HVTN 704/HPTN 085", "HVTN 703/HPTN 081")
  plots_ylim=lapply(plots, function(x, datPooled){ 
    if (grepl("epitope", x)){
      range(datPooled[, c("epitope.dist.b.ls", "epitope.dist.c.ls")], na.rm=TRUE)
    } else {
      range(datPooled[, x], na.rm=TRUE)   
    }
    }, datPooled=datPooled)
    
  plots_ybreaks=list("prob"=logit(seq(0.3, 0.9, by=0.1)), 
                     "predIC80"=log10(c(0.62, 1, 2, 3, 5, 10, 20, 30, 50, 80)),
                     "epiDist"=seq(0, 3, by=0.5),
                     "hammingDist"=seq(2, 12, by=2))
  plots_ylabels=list("prob"=seq(0.3, 0.9, by=0.1), 
                     "predIC80"=c(0.62, 1, 2, 3, 5, 10, 20, 30, 50, 80),
                     "epiDist"=seq(0, 3, by=0.5),
                     "hammingDist"=seq(2, 12, by=2))

  # y_list <- list(plots_ylim, plots_ybreaks, plots_ylabels, plots_ytrans)
  
  # load data
  datplot <- read.csv(file.path(datDir, datFile[j])) %>%
    filter(hiv1event==1 & !is.na(parscore1.ls)) %>%
    mutate(tx2=case_when(tx=="C3" ~ "Placebo",
                         tx=="T1" ~ "VRC01 \n10 mg/kg",
                         tx=="T2" ~ "VRC01 \n30 mg/kg"),
           tx2=factor(tx2, levels=c("Placebo", "VRC01 \n10 mg/kg", "VRC01 \n30 mg/kg"))) %>%
    select(tx2, all_of(plots))
  
  
  for (i in 1:length(plots)){
    
    # val <- ifelse(i %in% seq(1, length(plots), 2), 'prob', 'ic80')
    # 
    # if (val == "prob") {
    #   ylim <- y_list[[1]][[val]]
    #   ybreaks <- y_list[[2]][[val]]
    #   # ylim <- sort(1/(1+exp(c(-3.66, 3.66))))
    # } else {
    #   ylim <- y_list[[1]][[val]]
    #   ybreaks <- y_list[[2]][[val]]
    # }
    
    # scatter/box/violin plot
    # pdf(file=file.path(figDir, paste0(trial[j], "_violinBoxScatter_", plots2[i], ".pdf")), width=8, height=8)
    ggplotList[[length(plots) * (j - 1) + i]] <- myplot(dat=datplot, 
                                                        x="tx2", 
                                                        y=plots[i], 
                                                        ytitle=if (j==1) plots_ytitles[i] else NULL,
                                                        toptitle=if (i==1) plots_titles[j] else NULL, 
                                                        ylim=plots_ylim[[i]],
                                                        plotMargin=if (j==1) c(0.1,0.1,0.1,0.5) else c(0.1,0.1,0.1,0.1),
                                                        axis.text.x=TRUE,
                                                        ybreaks=plots_ybreaks[[i]],
                                                        ylabels=plots_ylabels[[i]])
    # dev.off()
  }
}
ggplotList <- ggplotList[c(1,5,2,6,3,7,4,8)]

# Combine panels into multi-panel manuscript figures ----------------------

pdf(file=file.path(figDir, "signaturesViolinBoxScatterByTxAndTrial_ls.pdf"), width=0.9 * 20, height=0.9 * 26)
# p <- ggpubr::ggarrange(plotlist=ggplotList, widths=c(1, 1), heights = rep(1, 4), ncol=2, nrow=4, align="v")
# p <- annotate_figure(p, top=text_grob("HVTN 704/HPTN 085", size=40, hjust=0))
# p <- annotate_figure(p, top=text_grob("HVTN 703/HPTN 081", size=40, hjust=1))
p <- ggplotList[[1]] + ggplotList[[2]] + ggplotList[[3]] + ggplotList[[4]] + 
  ggplotList[[5]] + ggplotList[[6]] + ggplotList[[7]] + ggplotList[[8]] +
  patchwork::plot_layout(ncol=2, height=rep(1, 4))
print(p)
dev.off()


# multiPanelPlot <- ggpubr::ggarrange(plotlist=c(ggplotList1[c("704.pred.prob.res.mf", "703.pred.prob.res.mf")],
#                                                ggplotList2[c("704.pred.prob.res.mf", "703.pred.prob.res.mf")]), nrow=2, ncol=2, align="hv")
# ggpubr::ggexport(multiPanelPlot, filename=file.path(figDir, "fig2_predProbResPARscore_mf.pdf"), width=16, height=16)
