#----------------------------------------------------------------------
# PROGRAM: DPEplots.R
#
# DESCRIPTION: Violin plots stratify by PNGS at 230 vs not PNGS at 230
# 
# Originally CODED BY: Erika Rudnicki
# Revised BY: Li Li

# INPUT: 
#   ../adata/v704_survival_wk80_tau_siePE_r4p_v2_cam.csv
#   ../tables/DPE/DPE_tier2_glycosite230PNGS_ls.csv
#   ../tables/DPE/DPE_tier2_glycosite230PNGS_mf.csv
#   ../tables/DPE/DPE_tier2_glycosite230PNGS_ms.csv
# OUTPUT: 
#   ../figures/descriptiPE/v704_violin_dPNGS230_gmt80mf.pdf
#   ../figures/descriptiPE/v704_violin_dPNGS230_gmt80ms.pdf
#   ../figures/descriptiPE/v704_violin_dPNGS230_gmt80ls.pdf
#   ../figures/descriptiPE/v703_violin_dPNGS230_gmt80mf.pdf
#   ../figures/descriptiPE/v703_violin_dPNGS230_gmt80ms.pdf
#   ../figures/descriptiPE/v703_violin_dPNGS230_gmt80ls.pdf

rm(list=ls(all=TRUE))
repoDir <- "/Users/lili/AMPsieve"
tablesOutDir <- file.path(repoDir, "tables/DPE/tier2")
figuresOutDir <- file.path(repoDir, "figures/DPE/tier2")
source(file.path(repoDir, "code/common.R")) #read in a datFile
source(file.path(repoDir,"code/DPE/DPEutils.R"))
### Read master file with covariates and time to ePEnt data for both trials (note the column southAmerica now has the combined trial/location info as follows:  value 0 = 704, not south america, 1 = 704, south america, 2 = 703)
master <-  read.csv( file.path(repoDir,"adata", datFile[3])) 

library(dplyr)
library(ggplot2)
options(dplyr.print_max = 100, stringsAsFactors=FALSE)

# function to plot
myplot <- function(dat, x, y, facetvar,
                   ylim=NULL,
                   ybreaks=NULL,
                   ylabels=NULL,
                   ytitle=NULL,
                   xtitle="Treatment Group",
                   toptitle=NULL,
                   boxcol=c("blue", "red3", "red3", "blue", "red3", "red3")){
  set.seed(39573056)
  ggplot(data=dat, aes_string(x=x, y=y, color=x)) +
    facet_grid(cols = vars(facetvar)) +
    geom_violin() +
    geom_boxplot(fill=boxcol, width=0.15, lwd=1.25, alpha = 0.3, outlier.shape=NA) + #remoPE outlier points
    geom_point(size=2, position = position_jitter(w = 0.3, h = 0)) +
    scale_y_continuous(breaks=ybreaks, labels=ylabels, trans="log10") +
    labs(x=xtitle, y=ytitle, title=toptitle) +
    scale_color_manual(values=boxcol) +
    theme_bw() +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          legend.position="none",
          plot.title = element_text(hjust = 0.5),
          text=element_text(size=23),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          axis.title.y = element_text(margin = margin(t = 0, r=0, b = 0, l = -2)))
}

for(sequence in c("mf", "ms", "ls")){
  for(trial in c("704", "703")){
    feature <- paste0("hxb2.230.pngs.",sequence)
    master.trial <- trial_dose_data (master, trial, "T1+T2")
    dat2 <- read.csv(file.path(tablesOutDir, paste0("DVE_tier2_glycosite230PNGS_",sequence,".csv")))
    if(trial == "704"){
      PE <- round(dat2[2:3, c("mean", "lower", "upper")], 0)
    }else{
      PE <- round(dat2[11:12, c("mean", "lower", "upper")], 0)
    }
    PE_PNGS <- paste0(PE$mean[1], "% (", PE$lower[1], ", ", PE$upper[1], ")")
    PE_noPNGS <- paste0(PE$mean[2], "% (", PE$lower[2], ", ", PE$upper[2], ")")
    
    datplot <- subset(master.trial, hiv1event==1 & nisolates>=1)
    datplot$tx2 <- with(datplot, case_when(tx=="C3" ~ "Placebo",
                                           tx=="T1" ~ "VRC01 \n10 mg/kg",
                                           tx=="T2" ~ "VRC01 \n30 mg/kg"))
    datplot$tx2 <- factor(datplot$tx2, levels=c("Placebo", "VRC01 \n10 mg/kg", "VRC01 \n30 mg/kg"))
    datplot$feature <- datplot[,feature]
    datplot$gmt80 = as.numeric(gsub(">", "", datplot[,paste0("gmt80",sequence)]))
    
    datplot$facetvar <- with(datplot, case_when(feature==1 ~ paste0("PNGS at 230-232 \n", "PE = ", PE_PNGS),
                                                feature==0 ~ paste0("No PNGS at 230-232 \n", "PE = ", PE_noPNGS)))
    
    datplot$facetvar <- factor(datplot$facetvar, levels = c(paste0("PNGS at 230-232 \n", "PE = ", PE_PNGS), paste0("No PNGS at 230-232 \n", "PE = ", PE_noPNGS)))
    plots_ytitles <- expression("Measured IC80 ("*mu*"g/ml)")
    plots_ybreaks=c(0.03, 0.1, 0.3, 1, 3, 10, 30, 100)
    plots_ylabels=c("0.03", "0.1", "0.3", "1", "3", "10", "30", expression("">="100"))
    pdf(file.path(repoDir,"figures/descriptive/",paste0("v",trial,"_violin_pngs230_gmt80_", sequence, ".pdf")), width=8, height=6)
    print(myplot(dat=datplot, 
                 x="tx2", 
                 y="gmt80", 
                 facetvar="facetvar",
                 ytitle=plots_ytitles,
                 xtitle = "Treatment Group",
                 toptitle="", 
                 ylim=plots_ylim,
                 ybreaks=plots_ybreaks,
                 ylabels=plots_ylabels))
    dev.off()
  }
}



