##########################################################################################
# Program : f_gmt_correlation.R
#
# Project: Create pairs plots of neutralization sensitivity results with panels defined
#          by IC50, IC80, and IIP cross classified with most sensitive, least sensitive,
#          and Most Frequent variant. Include covariate-adjusted spearman rhos.
#
#          This program plots HVTN 704 data, adjusting for treatment arm.
#    
# Location: /trials/vaccine/p704/analysis/efficacy/code
#
# Input:
#   ../adata/v704_survival_wk80_tau_neut.csv
#
# Output:
#   ../figures/adhoc/v704_gmt_correlation.pdf
#   ../tables/v704_primary_cases_nisolates.csv
#
# Code History
# ---------------------------------------------------------------------------------------
# Date        Programmer        Details 
# ---------------------------------------------------------------------------------------
# 30Jun2020   Erika Rudnicki    Version 1.0, heavily influenced by 
#             T:\vaccine\p100\analysis\adhoc\2016_12_CoR\correlation_analysis.R
# 08Jul2020   Erika Rudnicki    Add IIP
# 17Jul2020   Erika Rudnicki    Use y=x instead of linear smoother
# 22Jul2020   Erika Rudnicki    Use pairs2 macro to allow different line types by panel
# 30Jul2020   Erika Rudnicki    Remove title
# 31Oct2022   Kevin Gillespie   Updated to work with seive data
###########################################################################################

# setwd("/trials/vaccine/p704/analysis/efficacy/code")

library(PResiduals)
library(plyr)
source("../macro/pairs2.R")

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
  
  hist.panel = function(x, ...) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
         main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  
  panel.smooth = function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                           cex = 1, col.smooth = "red", lwd.smooth=1.5, span = 2/3, iter = 3, ...) {
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
    type.matrix[7:9,1:6] <- 1
    
    pairs2(df, gap = 0, lower.panel1 = panel.smooth, lower.panel2 = panel.diag, type.matrix = type.matrix,
           upper.panel = panel.cor.adj, pch=19,
          diag.panel = hist.panel, text.panel=NULL, xlim=c(-2,5), ylim=c(-2,5),oma=c(10,11,5,4)) #bltr
    
    # title(paste0(title, " (n=", nrow(dat), ")"))
    vars.n=sapply(vars,function(x)length(grep(x,names(df))))
    cols=c(rep("blue", length(lgnd)), rep("forestgreen", length(lgnd)), rep("purple", length(lgnd)))

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
         col=cols,xpd=T,cex=.7, #adj=0,
         labels=labels)
    
    text(y=seq(start.y1,end.y1,length.out=ncol(df)),x=fix.x,
         col=cols,xpd=T,cex=.7, #adj=1,
         labels=labels)
    
    par(xpd=T)
    
    legend(x=par('usr')[1]+(par('usr')[2]-par('usr')[1])/8.5*3.5,
           y=par('usr')[3]+(par('usr')[4]-par('usr')[3])/14,
           ncol=length(lgnd),cex=1,bty="n",text.width=rep.int(max(strwidth(lgnd)),(1+length(lgnd))),
           c(lgnd),
           fill=c(c("blue","forestgreen","purple","yellow","orange")[1:length(vars)]),
           border=NA)
    
    par(font=1)
  }
  
}# End of Extra Functions

# load data
dat <- read.csv("../adata/v704_survival_wk80_tau_neut.csv")
datplot <- subset(dat, nisolates > 0)
datplot[5,11:16]
datplot[,11:16] <- sapply(datplot[,11:16], function(y) log10(as.numeric(gsub(">|<", "", y))))

# plot
pdf(paste0('../figures/adhoc/v704_gmt_correlation.pdf'),width=11,height=8.5)

draw.plot(dat=datplot,
          covars=c("tx"),
          vars=c("gmt50ms","gmt50mf", "gmt50ls", "gmt80ms","gmt80mf", "gmt80ls", "IIPmsPool", "IIPmfPool", "IIPlsPool"),
          labels=c("Log10 IC50 \nMost Sensitive", "Log10 IC50 \nMost Frequent", "Log10 IC50 \nLeast Sensitive", 
                   "Log10 IC80 \nMost Sensitive", "Log10 IC80 \nMost Frequent", "Log10 IC80 \nLeast Sensitive", 
                   "IIP \nMost Sensitive", "IIP \nMost Frequent", "IIP \nLeast Sensitive"),
          lgnd=c("IC50", "IC80", "IIP"),
          title="HVTN 704/HPTN 085 Primary Endpoint Cases with Breakthrough Neutralization Data")

dev.off()

q(save = "no")


