# Purpose: Plotting of AA residue bar plots using analyzed viral variants only (each of mf, ms, ls separately)
#          Bar plots are generated for Tier 1 AA positions pre-identified to be associated with VRC01 neutralization (AMP SAP feature set 2)
# Method:  Brian Simpkins' plotting code
#          Original source code at https://github.com/brianSimpkins/AMPsite_scanning/blob/master/AMP_aa_dist_plotter.Rmd
#          Modified code used herein available in \trials\vaccine\p704\analysis\sieve\code\aaResidueBarPlots\aa_dist_plotter_functions.R
# Input:   Measured Env AA sequences in the AMP trials; 1 fasta file per primary endpoint
#          HXB2 lookup table pertaining to the used AA alignment in /data/VTN/VTN703_704/sieve/dat/seq/work/aa/ref/hxb2.map
# Output:  A PDF figure stacking bar plots of AA residue distributions for Tier 1 pre-identified neutralization-associated positions
# Author:  Michal Juraska

rm(list=ls(all=TRUE))

figDir <- paste0("t:/vaccine/p704/analysis/sieve/figures/aaResidueBarPlots")

library(ggpubr)
library(RColorBrewer)
library(gridExtra)
source("t:/vaccine/p704/analysis/sieve/code/aaResidueBarPlots/aa_dist_plotter_functions.R")

# read in sequence data
seq_path <- "o:/VTN/VTN703_704/sieve/analysis/cmagaret/13_mkSieveData/results/seqs_analysis"
seq_filename <- paste0(c("v704", "v703", "pooled"), "_sieve_env_v1.fasta")

# read in treatment mapping
source("t:/vaccine/p704/analysis/sieve/code/common.R")
seqMaster <- read.csv(file.path("t:/vaccine/p704/analysis/sieve/adata", datFile[3]))

# read in hxb2 mapping
map <- read.table("o:/VTN/VTN703_704/sieve/dat/seq/work/aa/ref/hxb2.map", sep="|", header=TRUE)

protocol <- paste0("HVTN ", c(704, 703))
variant <- c("mf", "ms", "ls")

# read in VRC01 sites
#VRC01sites <- read.delim("u:/bsimpkin/Site-Scanning/bsites_vrc01.dat", header=FALSE)$V1
# read in CD4 sites
#CD4sites <- read.delim("u:/bsimpkin/Site-Scanning/bsites_cd4.dat", header=FALSE)$V1

# neutralization-associated AA positions
sites <- c(60, 170, 230, 279, 280, 317, 365, 429, 456, 458, 459, 471)

for (trial in 1:length(seq_filename)){
  for (v in variant){
    # load in the data
    seqData <- gatherSequences2(seq_path, seq_filename[trial], v)
    
    # clean up seqMaster
    seqMaster1 <- seqMaster[seqMaster$protocol==protocol[trial], c("protocol", "pub_id", "tx", "hiv1event", "parscore1.mf")]
    seqMaster1$pub_id <- gsub("-", "_", seqMaster1$pub_id)
    seqMaster1 <- seqMaster1[order(seqMaster1$pub_id), ]
    seqMaster1$tx <- factor(seqMaster1$tx, levels=c("C3", "T1", "T2"), labels=c("Placebo", "VRC01 10mg/kg", "VRC01 30mg/kg"))
    
    # this cohort exclusion criterion is outdated
    # restrict to primary endpoints with in vitro neutralization sensitivity data
    # seqMaster1 <- subset(seqMaster1, hiv1event==1 & !is.na(gmt80ls))
    seqMaster1 <- subset(seqMaster1, hiv1event==1 & !is.na(parscore1.mf))
    
    # subset seqData
    seqData <- subset(seqData, pubid %in% seqMaster1$pub_id)
    
    # Site distribution plot
    # get list of plots
    plots <- dist.plot.builder(sequenceData=seqData, sites=sites, map=map, seqMaster=seqMaster1)
    
    # alternate backgrounds of the plots
    for(i in 1:length(plots)){
      if(i%%2==1){
        plots[[i]] <- plots[[i]] + theme(plot.background=element_rect(fill="grey85", color="grey85"))
      }
    }
    
    # Give the first page a facet label by creating an empty plot and appending it to the beginning
    firstP <- ggplot(seqMaster1, aes(x=pub_id, y=pub_id, fill=tx)) +
      geom_col(aes(pub_id, NA), alpha=0) +
      facet_grid(. ~ tx, scales="free", space="fixed", switch="x") +
      theme_void() +
      theme(legend.title=element_blank(), legend.text=element_blank(), strip.text=element_text(size=12, vjust=1),
            plot.margin=margin(b=.1, t=-.5, unit='cm'), plot.background = element_rect(fill="grey99", color="grey99"),
            axis.title.x = element_blank(), axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    plots <- c(list(firstP), plots)
    
    # save them to an output file
    multiPanelPlot <- ggpubr::ggarrange(plotlist = plots, nrow = ceiling(length(plots)/1), ncol = 1, align='v', heights=c(.3, rep(1, length(plots)-1)))
    multiPanelPlot <- annotate_figure(multiPanelPlot, bottom=text_grob("Primary Endpoint Cases", size=14), left=text_grob("Env AA Position", rot=90, size=14))
    
    if (!dir.exists(figDir)){ dir.create(figDir, recursive=TRUE) }
    # pdf(file.path(figDir, paste0(substring(protocol[trial], first=6), "_AAresidueBarPlot_neutAssocPos.pdf")), width=8.5, height=11)
    # ggpubr::ggarrange(plotlist = plots, nrow = ceiling(length(plots)/1), ncol = 1, align='v', heights=c(.3, rep(1, length(plots)-1)))
    # dev.off()
    
    ggpubr::ggexport(multiPanelPlot, filename=file.path(figDir, paste0(substring(protocol[trial], first=6), "_AAresidueBarPlot_", v, "_neutAssocPos.pdf")), width=8.5*0.9, height=11*0.9)
  }
}
