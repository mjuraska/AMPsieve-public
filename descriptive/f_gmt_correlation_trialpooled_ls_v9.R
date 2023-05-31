##########################################################################################
# Program : f_gmt_correlation.R
#
# Project: Create pairs plots of neutralization sensitivity results with panels defined
#          by IC50, IC80, and IIP cross classified with most sensitive, least sensitive,
#          and majority founder variant. Include covariate-adjusted spearman rhos.
#
#          This program plots combined AMP trials data, adjusting for treatment arm and trial.
#    
# Location: /trials/vaccine/p704/analysis/efficacy/code
#
# Input:
#   ../adata/amp_survival_wk80_tau_neut.csv
#
# Output:
#   ../figures/adhoc/amp_gmt_correlation.pdf
#
# Code History
# ---------------------------------------------------------------------------------------
# Date        Programmer        Details 
# ---------------------------------------------------------------------------------------
# 30Jun2020   Erika Rudnicki    Version 1.0, updated from 704
# 08Jul2020   Erika Rudnicki    Add IIP
# 30Jul2020   Erika Rudnicki    Remove title
# 31Oct2022   Kevin Gillespie   Updated to work with seive data
# 09Jan2023   Kevin Gillespie   Updated source data and targets for Michal
###########################################################################################

# setwd("/trials/vaccine/p704/analysis/efficacy/code")

library(PResiduals)
library(plyr)
source("../macro/pairs2b.R")

options(stringsAsFactors=F)

if(T){#Extra Functions
  
  panel.cor.adj <- function(x, y, z=datplot$tx, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- partial_Spearman(x|y ~ z)$TS$TB$ts
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    # test <- partial_Spearman(x|y ~ z)$TS$TB$pval
    # Signif <- symnum(test, corr = FALSE, na = FALSE, 
    #                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
    #                                                                           "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    # text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  
  hist.panel = function(x, i, ...) {
    par(new = TRUE)
    if (i == 4) {
      hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
           main = "", breaks = "FD", xlim = c(0,13))
    } else if (i ==3) {
      hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
           main = "", breaks = "FD", xlim = c(0,3))
    } else if (i == 2) {
      hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
           main = "", breaks = "FD", xlim = c(-1,2))
    } else if (i == 1) {
      hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
           main = "", breaks = "FD", xlim = c(-1,3))
    } else {
      hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
           main = "", breaks = "FD", xlim = c(-2,3))
    }
    
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  
  panel.smooth = function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                           cex = 1, col.smooth = "red", lwd.smooth=2.25, span = 7/8, iter = 3, ...) {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, lwd=lwd.smooth, ...)
    
  }
  
  panel.diag = function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                         cex = 1, col.diag = "red", lwd.diag=1, span = 2/3, iter = 3, ...) {
    # require(reshape2)
    # z <- merge(data.frame(x,y), melt(table(x ,y)),sort =F)$value
    # print(data.frame(x,y))
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    abline(a=0, b=1, col=col.diag, lwd=lwd.diag)
  }
  
  
  draw.plot=function(dat,vars,covars,labels,lgnd,title){  
    
    z=dat[,c(covars)]
    df=dat[,c(vars)]
    
    type.matrix <- matrix(2, nrow=ncol(df), ncol=ncol(df))
    # type.matrix[4:9,1:3] <- 2
    
    pairs2b(df, gap = 0, lower.panel1 = panel.diag, lower.panel2 = panel.smooth, type.matrix = type.matrix,
           upper.panel = panel.cor.adj, pch=19,
           diag.panel = hist.panel, text.panel=NULL, oma=c(10,11,5,4)) #bltr
    
    # title(paste0(title, " (n=", nrow(dat), ")"))
    vars.n=sapply(vars,function(x)length(grep(x,names(df))))
    # cols=c(rep("blue", length(lgnd)), rep("forestgreen", length(lgnd)), rep("purple", length(lgnd)))
    cols = c('blue', 'blue', 'forestgreen', 'forestgreen', 'purple')
    
    start.x=par('usr')[1]+(par('usr')[2]-par('usr')[1])/7.5
    end.x=par('usr')[2]-(par('usr')[2]-par('usr')[1])/20
    start.x1=start.x+(end.x-start.x)/2/ncol(df)
    end.x1=end.x-(end.x-start.x)/2/ncol(df)
    fix.y=par('usr')[3]+(par('usr')[4]-par('usr')[3])/10
    
    start.y=par('usr')[4]-(par('usr')[4]-par('usr')[3])/8.5*0.65
    end.y=par('usr')[3]+(par('usr')[4]-par('usr')[3])/6.5
    start.y1=start.y+(end.y-start.y)/2/ncol(df)
    end.y1=end.y-(end.y-start.y)/2/ncol(df)
    fix.x=par('usr')[1]+(par('usr')[2]-par('usr')[1])/17
    
    par(font=2)
    text(x=seq(start.x1,end.x1,length.out=ncol(df)),y=fix.y,
         col=cols,xpd=T,cex=0.9, #adj=0,
         labels=labels)
    
    text(y=seq(start.y1,end.y1,length.out=ncol(df)),x=fix.x,
         col=cols,xpd=T,cex=0.9, #adj=1,
         labels=labels)
    
    par(xpd=T)
    
    # legend(x=par('usr')[1]+(par('usr')[2]-par('usr')[1])/8.5*1.5,
    #        y=par('usr')[3]+(par('usr')[4]-par('usr')[3])/14,
    #        ncol=length(lgnd),cex=1,bty="n",text.width=rep.int(max(strwidth(lgnd)),(1+length(lgnd))),
    #        c(lgnd),
    #        fill=c(c("blue","forestgreen","purple","yellow","orange")[1:length(vars)]),
    #        border=NA)
    
    par(font=1)
  }
  
}# End of Extra Functions


# load data
dat <- read.csv("../adata/amp_sieve_pooled_marks_final_v9.csv")

colnames_pairs <- c(
  "parscore1.ls", "parscore2.ls",
  "epitope.dist.any.ls", 
  "hdist.zspace.sites.binding.all.ls",
  "gmt80ls"
)

datplot <- dat[,c('protocol', 'tx', colnames_pairs)]
datplot[,colnames_pairs[which(grepl("gmt", colnames_pairs))]] <- sapply(colnames_pairs[which(grepl("gmt", colnames_pairs))], function(y) log10(as.numeric(gsub(">|<", "", datplot[,y]))))

datplot$covar = paste0(datplot$protocol, datplot$tx)

cairo_pdf(paste0('../figures/descriptive/v703and704_gmt_correlation_ls_v9a.pdf'),width=11,height=8.5)

draw.plot(dat=datplot,
          covars=c("covar"),
          vars=c("parscore1.ls", "parscore2.ls", "epitope.dist.any.ls", 
                 "hdist.zspace.sites.binding.all.ls", "gmt80ls"),
          labels=c(
            "Pred Prob\nIC80 > 1 \u03bcg/ml",
            # expression(atop('Logit Pred Prob','IC80 > 1'*mu*"g/ml")),
            "Pred IC80",
            "Epitope Dist to\nSubtype Any Ref",
            "PC-Weighted\nHamming Dist",
            "IC80"),
          lgnd=c("Predicted Prob\nIC80 > 1 mcg/ml", "Predicted IC80 (mcg/ml)", 
                 "Epitope Dist", "PC-Weighted", 
                 "Observed IC80 (mcg/ml)"),
          title="Pooled AMP Trials Primary Cases with Breakthrough Neutralization Data")

dev.off()

q(save = "no")


