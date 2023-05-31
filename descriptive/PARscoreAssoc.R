# Purpose: Association of PAR scores 1 and 2
# Method:  Scatter plot and LOESS
# Input:   PAR scores 1 and 2 for AMP adjudicated primary endpoints
# Output:  PDF files, each containing a single plot
# Author:  Michal Juraska
# Date:    Mar 31, 2022

rm(list=ls(all=TRUE))

datDir <- paste0("t:/vaccine/p704/analysis/sieve/adata")
figDir <- paste0("t:/vaccine/p704/analysis/sieve/figures/descriptive")

library(ggplot2)

# trial-specific data files created by Craig Magaret, with clinical and neutralization data from Erika Rudnicki
# earlier versions of the trial-pooled data file created by t:\vaccine\p704\analysis\sieve\code\createTrialPooledData.R
source("t:/vaccine/p704/analysis/sieve/code/common.R")

marks <- data.frame(x=paste0("parscore2.", c("mf", "ms", "ls")), 
                    y=paste0("pred.prob.res.", c("mf", "ms", "ls")))

# tags that are part of the output PDF file names
trialFileString <- c("704", "703", "704and703")

# plot labels
xLab <- expression("Predicted" ~ IC[80] ~ "(" * mu * "g/ml)")
yLab <- expression("Predicted Probability of" ~ IC[80] > 1 ~ mu * "g/ml")

xTickLabs <- c(0.62, 1, 2, 3, 5, 10, 20, 30, 50, 80)

# for computing common 'xLim' across 704, 703, and pooled
datPooled <- read.csv(file.path(datDir, datFile[3]))

# compute 'xLim' based on trial-pooled data
parscore2.pooled <- datPooled[datPooled$hiv1event==1 & !is.na(datPooled$gmt80mf), marks[ , "x"]]
xLim <- range(parscore2.pooled, na.rm=TRUE)

# cycle over 704, 703, trial-pooled
for (trial in 1:length(datFile)){
  dat <- read.csv(file.path(datDir, datFile[trial]))
  
  # for each row in 'marks'
  for (m in 1:NROW(marks)){
    # keep only cases the two PAR scores to be plotted
    dat1 <- dat[dat$hiv1event==1 & !is.na(dat[, marks[m, "x"]]) & !is.na(dat[, marks[m, "y"]]), ]
    
    title <- paste0(switch(substring(marks[m, "x"], first=nchar(marks[m, "x"]) - 1), "mf"="Most Frequent", "ms"="Pred. Most Sensitive", "ls"="Pred. Most Resistant"), " Founder")
    
    p <- ggplot(dat1, aes_string(x=marks[m, "x"], y=marks[m, "y"])) +
      geom_point(shape=21, size=1.3) +
      labs(x=xLab, y=yLab, title=title) +
      coord_cartesian(xlim=xLim, ylim=0:1) +   # x-axis range based on the support of PAR score 2 in trial-pooled data
      scale_x_continuous(breaks=log10(xTickLabs), labels=xTickLabs) +
      scale_y_continuous(breaks=seq(0, 1, by=0.2)) +
      stat_smooth(method="loess", color="blue", se=FALSE) +
      annotate("text", x=Inf, y=-Inf, hjust=1.1, vjust=-0.3, 
               label=as.character(as.expression(substitute("Spearman's corr"==r, 
                                                           list(r=format(cor(dat1[, marks[m, "x"]], dat1[, marks[m, "y"]], method="spearman"), nsmall=2, digits=2))))), 
               parse=TRUE) +
      theme_bw()
    
    if (!dir.exists(figDir)){ dir.create(figDir) }
    ggsave(file.path(figDir, paste0(trialFileString[trial], "_predIC80vsPredProbResIC80_", substring(marks[m, "x"], first=nchar(marks[m, "x"]) - 1), ".pdf")),
           plot=p, width=3, height=3)
  }
}
