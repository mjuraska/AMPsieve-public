library(ggpubr) # for combining two panels into one
library(ggpmisc) # for adding a table to the top panel
library(scales) # for pretty x-axis breaks
library(patchwork)

#' Plotting Mark-Specific Proportional Hazards Model Fits
#'
#' \code{plot} method for class \code{summary.sievePH}. For univariate marks, it plots point and interval estimates of the mark-specific treatment effect parameter specified by \code{contrast} in \code{\link{summary.sievePH}}, and,
#' optionally, scatter/box plots of the observed mark values by treatment. For bivariate marks, plotting is restricted to the point estimate, which is displayed as a surface. No plotting is provided for marks of higher dimensions.
#'
#' @param x an object returned by \code{\link{summary.sievePH}}
#' @param mark either a numeric vector specifying a univariate continuous mark or a data frame specifying a multivariate continuous mark.
#' For subjects with a right-censored time-to-event, the value(s) in \code{mark} should be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo)
#' @param xlim a numeric vector of length 2 specifying the x-axis range (\code{NULL} by default)
#' @param ylim a numeric vector of length 2 specifying the y-axis range (\code{NULL} by default)
#' @param zlim a numeric vector of length 2 specifying the z-axis range in a 3-dimensional plot (\code{NULL} by default)
#' @param xtickAt a numeric vector specifing the position of x-axis tickmarks (\code{NULL} by default)
#' @param xtickLab a numeric vector specifying labels for tickmarks listed in \code{xtickAt}. If \code{NULL} (default), the labels are determined by \code{xtickAt}.
#' @param ytickAt a numeric vector specifing the position of y-axis tickmarks (\code{NULL} by default)
#' @param ytickLab a numeric vector specifying labels for tickmarks listed in \code{ytickAt}. If \code{NULL} (default), the labels are determined by \code{ytickAt}.
#' @param xlab a character string specifying the x-axis label (\code{NULL} by default)
#' @param ylab a character string specifying the y-axis label (\code{NULL} by default)
#' @param zlab a character string specifying the z-axis label in a 3-dimensional plot (\code{NULL} by default)
#' @param txLab a character vector of length 2 specifying the placebo and treatment labels (in this order). The default labels are \code{placebo} and \code{treatment}.
#' @param title a character string specifying the plot title (\code{NULL} by default)
#' @param ... other arguments to be passed to plotting functions
#'
#' @details
#' For bivariate marks, \code{markGrid} in \code{\link{summary.sievePH}} must have equally spaced values for each component.
#'
#' @return None. The function is called solely for plot generation.
#'
#' @examples
#' n <- 500
#' tx <- rep(0:1, each=n/2)
#' tm <- c(rexp(n/2, 0.2), rexp(n/2, 0.2 * exp(-0.4)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' mark <- ifelse(eventInd==1, c(rbeta(n/2, 2, 5), rbeta(n/2, 2, 2)), NA)
#' markRng <- range(mark, na.rm=TRUE)
#'
#' # fit a model with a univariate mark
#' fit <- sievePH(eventTime, eventInd, mark, tx)
#' sfit <- summary(fit, markGrid=seq(markRng[1], markRng[2], length.out=10))
#' plot(sfit, mark, tx)
#'
#' @seealso \code{\link{sievePH}}, \code{\link{sievePHipw}}, \code{\link{sievePHaipw}} and \code{\link{summary.sievePH}}
#'
#' @export
plot.summary.sievePH <- function(x, mark=NULL, tx=NULL, xlim=NULL, ylim=NULL, zlim=NULL, xtickAt=NULL, xtickLab=NULL, ytickAt=NULL, ytickLab=NULL, 
                                 xlab=NULL, xLabLine=3, ylab=NULL, zlab=NULL, txLab=c("Placebo", "Treatment"), title=NULL, subtitle=NULL,
                                 compTE=NULL, sec.xtickAt=NULL, sec.xtickLab=NULL, sec.xlab=NULL,
                                 parMar=c(4.5, 6, ifelse(is.null(subtitle), 2, 3), 1), ...){
  contrast <- names(x)[length(names(x))]
  
  cexAxis <- 1.5
  cexLab <- 1.7
  cexTitle <- 2
  cexSubtitle <- 1.5
  cexText <- 1.3
  cexLegend <- 1.3
  
  par(mar=parMar, oma=rep(0, 4), cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexTitle)
  
  # a 2-dimensional plot only when the mark is univariate
  if (NCOL(x[[contrast]])==4){
    if (is.null(xlim)){
      xlim <- range(x[[contrast]][, 1])
    }
    
    if (is.null(ylim)){
      ylim <- range(x[[contrast]][, -1], na.rm=TRUE)
    }
    ySplit <- ylim[2]
    
    # need extra room for box plots on the top
    if (!any(c(is.null(mark), is.null(tx)))){
      ylim <- c(ylim[1], ylim[2] + ifelse(is.null(sec.xtickLab), 0.35, 0.35) * (ylim[2] - ylim[1]))
    }
    
    if (is.null(xlab)){ xlab <- colnames(x[[contrast]])[1] }
    if (is.null(ylab)){ ylab <- switch(colnames(x[[contrast]])[2], TE="Treatment Efficacy", HR="Hazard Ratio", LogHR="Log Hazard Ratio") }
    
    plot(x[[contrast]][, 1], x[[contrast]][, 2], xlim=xlim, ylim=ylim, type="n", xlab="", ylab="", xaxt=ifelse(is.null(xtickAt), "s", "n"), yaxt="n", bty="l", ...)
    
    if (!is.null(xtickAt)){
      if (is.null(xtickLab)){ xtickLab <- xtickAt }
      axis(side=1, at=xtickAt[seq(1, length(xtickAt), by=2)], labels=xtickLab[seq(1, length(xtickLab), by=2)], cex.axis=cexAxis)
      axis(side=1, at=xtickAt[seq(2, length(xtickAt), by=2)], labels=xtickLab[seq(2, length(xtickLab), by=2)], cex.axis=cexAxis)
    }
    
    if (!is.null(sec.xtickAt) && !is.null(sec.xtickLab)){
      axis(side=1, line=4.5, at=sec.xtickAt[seq(1, length(sec.xtickAt), by=2)], labels=sec.xtickLab[seq(1, length(sec.xtickLab), by=2)], cex.axis=cexAxis, 
           col="magenta3", col.ticks="magenta3", col.axis="magenta3")
      axis(side=1, line=4.5, at=sec.xtickAt[seq(2, length(sec.xtickAt), by=2)], labels=sec.xtickLab[seq(2, length(sec.xtickLab), by=2)], cex.axis=cexAxis, 
           col="magenta3", col.ticks="magenta3", col.axis="magenta3")
    }
    
    if (!is.null(ytickAt)){
      if (is.null(ytickLab)){ ytickLab <- ytickAt }
      axis(side=2, at=ytickAt, labels=ytickLab, las=1, cex.axis=cexAxis)
    } else {
      # to avoid overlapping tickmarks
      if (!any(c(is.null(mark), is.null(tx)))){
        axis(side=2, at=axTicks(2)[axTicks(2) <= 1.1 * max(x[[contrast]][, 4], na.rm=TRUE)])
      }
    }
    
    if (is.null(sec.xlab)){
      mtext(xlab, side=1, line=xLabLine, cex=cexLab)
    } else { 
      mtext(xlab, side=1, line=2.7, cex=cexLab)
      mtext(sec.xlab, side=1, line=7, cex=cexLab, col="magenta3") 
    }
    
    mtext(ylab, side=2, line=3.5, las=3, cex=cexLab)
    
    if (!is.null(title)){ mtext(title, side=3, font=1, line=ifelse(is.null(subtitle), 0, 1.5), cex=cexTitle) }
    if (!is.null(subtitle)){ mtext(subtitle, side=3, font=1, line=0, cex=cexSubtitle) }
    
    abline(h=ifelse(colnames(x[[contrast]])[2]=="HR", 1, 0), col="gray70", lwd=2)
    
    ylimGap <- 0
    if (!is.null(compTE)){
      ylimGap <- 0.3
      # colCI <- "darkgoldenrod2"
      # colRGB <- c(col2rgb(colCI))
      # colRGB <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha=255*0.55, maxColorValue=255)
      # polygon(c(x[[contrast]][, 1], rev(x[[contrast]][, 1])), 
      #         c(ifelse(compTE$LB2 >= ylim[1] & compTE$LB2 <= ySplit, compTE$LB2, ylim[1]), rev(ifelse(compTE$UB2 <= 1, compTE$UB2, 1))), col=colRGB, border=NA)
      lines(x[[contrast]][, 1], ifelse(compTE$TE2 >= ylim[1] + ylimGap & compTE$TE2 <= ySplit, compTE$TE2, NA), lwd=6, col="magenta3")
    }
    
    lines(x[[contrast]][, 1], ifelse(x[[contrast]][, 2] >= ylim[1] + ylimGap & x[[contrast]][, 2] <= ySplit, x[[contrast]][, 2], NA), lwd=6)
    lines(x[[contrast]][, 1], ifelse(x[[contrast]][, 3] >= ylim[1] + ylimGap & x[[contrast]][, 3] <= ySplit, x[[contrast]][, 3], NA), lwd=5.5, lty="dashed")
    lines(x[[contrast]][, 1], ifelse(x[[contrast]][, 4] >= ylim[1] + ylimGap & x[[contrast]][, 4] <= ySplit, x[[contrast]][, 4], NA), lwd=5.5, lty="dashed")
    
    # text(min(out$v), -0.6, paste0("Marginal Sieve Test P ",ifelse(pMarginalSieve<0.001,"< 0.001",paste0("= ",format(pMarginalSieve,digits=2))),marginalSignifMark), pos=4, cex=cexText)
    
    #legend("bottomleft", fill=colCI, border=colCI, legend="95% Pointwise CI", cex=cexLegend, bty="n")
    #legend("bottomleft", lwd=c(3.5,2), lty=c("dashed","longdash"), legend=c("95% Pointwise CI","Overall Hazard-Ratio PE"), col=c("black","darkorange"), cex=cexLegend, bty="n")
    if (is.null(compTE)){
      legend(x=xlim[1], y=ifelse(colnames(x[[contrast]])[2]=="TE", ylim[1] + ifelse(parMar[1]<5, 0.05, 0.07) * (ylim[2] - ylim[1]), ySplit), lwd=5.5, lty="dashed", 
             legend="95% Pointwise CI", col="black", cex=cexLegend, bty="n")  
    } else {
      legend(x=xlim[1], y=ifelse(colnames(x[[contrast]])[2]=="TE", ylim[1] + 0.16 * (ylim[2] - ylim[1]), ySplit), lwd=c(5.5, 6), lty=c("dashed", "solid"), 
             legend=c("95% Pointwise CI", expression("Est. PE by Measured" ~ IC[80])), col=c("black", "magenta3"), cex=cexLegend, bty="n",
             x.intersp=0.5)
    }
    
    # add scatter/box plots of the observed mark values by treatment
    if (!any(c(is.null(mark), is.null(tx)))){
      if (is.null(sec.xtickLab)){
        par(fig=c(0,1,0.7,1), new=TRUE)  
      } else {
        par(fig=c(0,1,0.75,1), new=TRUE)  
      }
      
      data <- na.omit(cbind(mark, tx))
      plotMarkHoriz(data[, 1], data[, 2], parMar=c(0.5, parMar[-1]), yLim=xlim, txLab=txLab)
      
      par(fig=c(0,1,0,1), new=TRUE)
    }
    
    # a 3-dimensional plot (a surface) when the mark is bivariate
  } else if (NCOL(x[[contrast]])==5){
    if (is.null(xlim)){
      xlim <- range(x[[contrast]][, 1])
    }
    
    if (is.null(ylim)){
      ylim <- range(x[[contrast]][, 2])
    }
    
    if (is.null(zlim)){
      zlim <- range(x[[contrast]][, 3])
    }
    
    if (is.null(xlab)){ xlab <- colnames(x[[contrast]])[1] }
    if (is.null(ylab)){ ylab <- paste0("\n", colnames(x[[contrast]])[2]) }
    if (is.null(zlab)){ zlab <- switch(colnames(x[[contrast]])[3], TE="\n\nTreatment Efficacy", HR="\n\nHazard Ratio", LogHR="\n\nLog Hazard Ratio") }
    
    # the first two arguments must be vectors with equally spaced values in ascending order
    persp(sort(unique(x[[contrast]][, 1])), sort(unique(x[[contrast]][, 2])), getOuterProduct(x[[contrast]], zlim),
          xlab=xlab, ylab=ylab, zlab=zlab, col="lightgreen", theta=150, phi=20, ticktype="detailed",
          nticks=5, xlim=xlim, ylim=ylim, zlim=zlim, r=3, expand=0.8, main=title)
  } else {
    stop("Plotting of results is available for univariate and bivariate marks only.")
  }
}

plotMarkHoriz <- function(mark, tx, parMar, yLim, txLab=c("Placebo", "Treatment")){
  cexAxis <- 1.5
  cexLab <- 1.7
  
  par(mar=parMar, oma=c(0,0,0,0), cex.axis=cexAxis, cex.lab=cexLab)
  # vioplot(mark[tx==0], mark[tx==1], at=c(0.5, 1.5), names=NA, horizontal=TRUE, drawRect=TRUE, col="white", border=c("blue", "red3"), rectCol=c("blue", "red3"), 
  #         lineCol=c("blue", "red3"), frame.plot=FALSE, yaxt="n")
  boxplot(mark ~ as.factor(tx), at=c(0.5,1.5), xlim=c(0,2), ylim=c(yLim[1], yLim[2]), frame.plot=FALSE, xaxt="n", yaxt="n",
          xlab="", ylab="", boxwex=0.75, outline=FALSE, border="black", lwd=3, horizontal=TRUE)
  axis(side=2, at=c(0.5,1.5), labels=txLab, cex.axis=cexAxis, las=1)
  points(mark, jitter(tx + 0.5, factor=0.7), col=ifelse(tx==1, "red3", "blue"), pch=ifelse(tx==1, 24, 21), lwd=1.9, cex=ifelse(tx==1, 1.4, 1.4), bg="white")
  # points(mark, jitter(tx + 0.5, factor=0.6), col=ifelse(tx==1, "red3", "blue"), pch=20, lwd=2, cex=0.6)
}

getOuterProduct <- function(df, zlim){
  mark1 <- sort(unique(df[, 1]))
  mark2 <- sort(unique(df[, 2]))
  
  out <- matrix(NA, nrow=length(mark1), ncol=length(mark2))
  for (i in 1:length(mark1)){
    for (j in 1:length(mark2)){
      idx <- which(df[, 1]==mark1[i] & df[, 2]==mark2[j])
      if (length(idx) > 1){
        stop("There are replicates on the marker grid.")
      } else if (length(idx)==1){
        out[i, j] <- df[idx, 3]
      }
      
    }
  }
  out <- ifelse(out >= zlim[1] & out <= zlim[2], out, NA)
  return(out)
}


#' Sieve PH plot with two panels: a descriptive panel and a vaccine efficacy panel
#' @param fit.te a data frame returned by summary.sievePH or a list of data frames, each formatted as the output from summary.sievePH and used for
#' plotting a single VE curve
#' @param data a data frame specifying a continuous mark under column "markValue" and a categorical variable describing treatment (Vaccine/Placebo).
#' For subjects with a right-censored time-to-event, the value(s) in mark should be set to NA.
#' @param top.panel.type options are a scatter/box plot ("box"), violin plot ("violin") and rcdf plot ("rcdf)
#' @param limits.x the range of the mark to be plotted
#' @param breaks.x the values where ticks will be placed for the mark
#' @param labels.x x-axis tick mark labels. If left unspecified, then values in \code{breaks.x} will be used.
#' @param xlab the x-axis label
#' @param sieveTestLabel a character string describing sieve test p-values
#' @param title title for the plot
#' @param lineColor a character vector of the same length as \code{fit.te} specifying colors of the VE curves
ggplotSievePH <- function(fit.te, data, top.panel.type = c("box", "violin", "rcdf"), limits.x=NULL, breaks.x=NULL, labels.x=breaks.x, 
                          xlab, limits.y=NULL, sieveTestLabel=NULL, title=NULL, lineColor=NULL){
  if (is.null(limits.x)){
    limits.x <- range(data$markValue, na.rm=TRUE)
  }
  
  if (is.null(limits.y)){
    limits.y <- c(-0.4, 1)
  }
  
  legend.x.just.p1 <- 0:1
  
  
  if (is.data.frame(fit.te)){
    p1 <- ggplot(fit.te) +
      geom_line(aes(x=mark, y=TE), size=1.7, colour=ifelse(is.null(lineColor), "black", lineColor)) +
      geom_line(aes(x=mark, y=LB, linetype='95% Pointwise CI', size='95% Pointwise CI'), colour=ifelse(is.null(lineColor), "black", lineColor)) +
      geom_line(aes(x=mark, y=UB, linetype='95% Pointwise CI', size='95% Pointwise CI'), colour=ifelse(is.null(lineColor), "black", lineColor)) +
      scale_linetype_manual(name="", labels=c('95% Pointwise CI'), values=c('95% Pointwise CI'= "dashed")) +
      scale_size_manual(name="", labels=c('95% Pointwise CI'), values=c('95% Pointwise CI'= 1.5)) +
      scale_y_continuous(name="Prevention Efficacy (%)", breaks=seq(-0.4, 1, 0.2), labels=seq(-0.4, 1, 0.2) * 100) +
      coord_cartesian(ylim=limits.y) +
      xlab(xlab) +
      theme_bw() +
      #geom_hline(yintercept=0, colour="grey", size = 1) +
      # annotate(geom = "text", x = min(limits.x), y = text.y, 
      #          label = sieveTestLabel, 
      #          size = 4, hjust = 0)+
      theme(legend.key.size = unit(0.65, "cm"),
            legend.title=element_text(size = 15),
            legend.margin=margin(grid::unit(0,"cm")),
            legend.text=element_text(size=13),
            legend.position = c(0.5, 1),
            legend.justification = c(0.5, 1),
            legend.key = element_blank(),
            legend.key.width = unit(1.4,"cm"),
            legend.background = element_blank(),
            plot.title = element_text(hjust = 0.5, vjust = 2, size = 20),
            axis.title.x = element_text(size = 15, vjust = -0.5),
            axis.title.y = element_text(size = 15, hjust = 0.2),
            axis.text.x = element_text(size = 14, hjust = 0.5,vjust = 0.5, colour = "black"),
            axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),
            plot.margin=unit(c(-0.62,0.2,0.2,0.2), "cm"))
    
  } else if (is.list(fit.te)){
    p1 <- ggplot() +
      geom_hline(yintercept=0, colour="gray50")
    
    for (i in 1:length(fit.te)){
      fit.te[[i]]$ic80g1 <- factor(fit.te[[i]]$ic80g1, levels=0:1)
      
      p1 <- p1 + 
        geom_ribbon(aes(x=mark, ymin=LB, ymax=UB, colour=ic80g1, fill=ic80g1), data=fit.te[[i]], alpha=0.2) +
        geom_line(aes(x=mark, y=TE, colour=ic80g1), data=fit.te[[i]], na.rm=TRUE, size=1.9)
    }
    
    legendLabels <- c(expression("" <= 1 ~ mu * "g/ml"), expression("" > 1 ~ mu * "g/ml"))
    
    p1 <- p1 +
      scale_colour_manual(values=c("#159724", "#CD3333"), name=expression("Measured" ~ IC[80]), labels=legendLabels, guide=guide_legend(ncol=2, title.position="left")) +
      scale_fill_manual(values=c("#159724", "#CD3333"), name=expression("Measured" ~ IC[80]), labels=legendLabels, guide=guide_legend(ncol=2, title.position="left")) +
      scale_y_continuous(name="Prevention Efficacy (%)", breaks=seq(-0.4, 1, 0.2), labels=seq(-0.4, 1, 0.2) * 100) +
      coord_cartesian(ylim=limits.y) +
      xlab(xlab) +
      theme_bw() +
      #geom_hline(yintercept=0, colour="grey", size = 1) +
      # annotate(geom = "text", x = min(limits.x), y = text.y, 
      #          label = sieveTestLabel, 
      #          size = 4, hjust = 0)+
      theme(legend.key.size = unit(0.65, "cm"),
            legend.title=element_text(size = 15),
            legend.margin=margin(grid::unit(0,"cm")),
            legend.text=element_text(size=14),
            legend.position = "bottom",
            #legend.justification = c(0.5, 1),
            legend.key = element_blank(),
            legend.key.width = unit(1.4,"cm"),
            legend.background = element_blank(),
            plot.title = element_text(hjust = 0.5, vjust = 2, size = 20),
            axis.title.x = element_text(size = 19, vjust = -0.5),
            axis.title.y = element_text(size = 19, hjust = 0.5),
            axis.text.x = element_text(size = 17, hjust = 0.5,vjust = 0.5, colour = "black"),
            axis.text.y = element_text(size = 17, colour = "black", hjust = 0.5),
            plot.margin=unit(c(0,0.2,0.2,0.2), "cm"))
    
  } else {
    stop("'fit.te' must be a data frame or a list of data frames.")
  }

  if (is.null(breaks.x)){
    p1 <- p1 + scale_x_continuous(limits = limits.x, breaks = pretty_breaks(n = 4))
  } else {
    p1 <- p1 + scale_x_continuous(limits = limits.x, breaks = breaks.x, labels=labels.x)
  }
  
  data_c <- drop_na(data) #Remove observations with missing mark
  nCasesVaccine <- length(data_c$treatment[data_c$treatment!= "Placebo"])
  nCasesPlacebo <- length(data_c$treatment[data_c$treatment == "Placebo"])
  quantilesVaccine <- quantile(data_c$markValue[data_c$treatment != "Placebo"], probs = c(0, 0.25, 0.5, 0.75, 1))
  quantilesPlacebo <- quantile(data_c$markValue[data_c$treatment== "Placebo"], probs = c(0, 0.25, 0.5, 0.75, 1))
  quantilesVaccine <- sapply(quantilesVaccine, function(x)format(as.numeric(x), digits = 2, nsmall = 2))
  quantilesPlacebo <- sapply(quantilesPlacebo, function(x)format(as.numeric(x), digits = 2, nsmall = 2))
  table = data.frame(rbind(c("VRC01",nCasesVaccine,quantilesVaccine ), c("Placebo",nCasesPlacebo, quantilesPlacebo)))
  colnames(table) <- c("Treatment", "No. of Cases", "Min", "Q1", "Median", "Q3", "Max")
  
  if (top.panel.type=="box"){
    data <- data %>% 
      mutate(ic80g1=factor(ic80g1, levels=0:1),
             txIC80=factor(interaction(treatment, ic80g1), levels=c("Placebo.0", "VRC01.0", "Placebo.1", "VRC01.1")))
      #add_row(markValue = NA, ic80g1=NA, treatment=NA, txIC80="Numerical Summary") %>%
  
    p2 <- ggplot(data) +
      ggtitle(title) +
      geom_boxplot(aes(x=txIC80, y=markValue, colour=ic80g1), outlier.shape=NA, na.rm=TRUE) +
      geom_jitter(aes(x=txIC80, y=markValue, colour=ic80g1), na.rm=TRUE, width=0.2, height=0, shape=21, fill="white", size=2, stroke=1.3) +
      scale_colour_manual(values=c("#159724", "#CD3333"), guide="none") +
      coord_flip() +
      #annotate(geom="table", x="Numerical Summary", y=min(data$markValue, na.rm=TRUE), label=list(table), vjust=0, fill="gray95") +
      #scale_x_discrete(expand= expansion(add = c(0.6, 2)))+
      scale_x_discrete(name=NULL, labels=c("Placebo.0"="Placebo", "VRC01.0"="VRC01", "Placebo.1"="Placebo", "VRC01.1"="VRC01")) +
      
      theme_bw() +
      theme(plot.title=element_text(hjust=0.5, size=22),
            axis.text.y=element_text(size=12, colour="black"),
            axis.ticks.x=element_blank(),
            plot.margin=unit(c(0.2,0.2,0,0.2), "cm")) 
    
    if (is.null(breaks.x)){
      p2 <- p2 + scale_y_continuous(name=NULL, limits = limits.x, breaks = pretty_breaks(n = 4), labels=NULL)
    } else {
      p2 <- p2 + scale_y_continuous(name=NULL, limits = limits.x, breaks = breaks.x, labels=NULL)
    }
    
  }

  # p <- ggarrange(p2, p1, heights=c(1, 2), ncol=1, nrow=2, align="v")
  # 
  # if (!is.null(title)){
  #   p <- annotate_figure(p, top=text_grob(title, color="black", face="bold", size=14, vjust=2.3))  
  # }
  
  p <- p2 + p1 + plot_layout(ncol=1, heights=c(1, 2))
  
  return(p)
}

# this function is currently broken
# there are two calls of scale_colour_manual(), one for the curves and one for the box plot, which doesn't work
ggplotSievePH2 <- function(fit.te, data, top.panel.type = c("box", "violin", "rcdf"), limits.x=NULL, breaks.x=NULL, labels.x=breaks.x, 
                          xlab, limits.y=NULL, sieveTestLabel=NULL, title=NULL, lineColor=NULL){
  if (is.null(limits.x)){
    limits.x <- range(data$markValue, na.rm=TRUE)
  }
  
  if (is.null(limits.y)){
    limits.y <- c(-0.4, 2)
  }
  
  legend.x.just.p1 <- 0:1
  
  
  if (is.data.frame(fit.te)){
    p1 <- ggplot(fit.te) +
      geom_line(aes(x=mark, y=TE), size=1.7, colour=ifelse(is.null(lineColor), "black", lineColor)) +
      geom_line(aes(x=mark, y=LB, linetype='95% Pointwise CI', size='95% Pointwise CI'), colour=ifelse(is.null(lineColor), "black", lineColor)) +
      geom_line(aes(x=mark, y=UB, linetype='95% Pointwise CI', size='95% Pointwise CI'), colour=ifelse(is.null(lineColor), "black", lineColor)) +
      scale_linetype_manual(name="", labels=c('95% Pointwise CI'), values=c('95% Pointwise CI'= "dashed")) +
      scale_size_manual(name="", labels=c('95% Pointwise CI'), values=c('95% Pointwise CI'= 1.5)) +
      scale_y_continuous(name="Prevention Efficacy (%)", breaks=seq(-0.4, 1, 0.2), labels=seq(-0.4, 1, 0.2) * 100) +
      coord_cartesian(ylim=limits.y) +
      xlab(xlab) +
      theme_bw() +
      #geom_hline(yintercept=0, colour="grey", size = 1) +
      # annotate(geom = "text", x = min(limits.x), y = text.y, 
      #          label = sieveTestLabel, 
      #          size = 4, hjust = 0)+
      theme(legend.key.size = unit(0.65, "cm"),
            legend.title=element_text(size = 15),
            legend.margin=margin(grid::unit(0,"cm")),
            legend.text=element_text(size=13),
            legend.position = c(0.5, 1),
            legend.justification = c(0.5, 1),
            legend.key = element_blank(),
            legend.key.width = unit(1.4,"cm"),
            legend.background = element_blank(),
            plot.title = element_text(hjust = 0.5, vjust = 2, size = 20),
            axis.title.x = element_text(size = 15, vjust = -0.5),
            axis.title.y = element_text(size = 15, hjust = 0.2),
            axis.text.x = element_text(size = 14, hjust = 0.5,vjust = 0.5, colour = "black"),
            axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),
            plot.margin=unit(c(-0.62,0.2,0.2,0.2), "cm"))
    
  } else if (is.list(fit.te)){
    p1 <- ggplot() +
      geom_hline(yintercept=0, colour="gray50")
    
    for (i in 1:length(fit.te)){
      fit.te[[i]]$ic80g1 <- factor(fit.te[[i]]$ic80g1, levels=0:1)
      
      p1 <- p1 + 
        geom_ribbon(aes(x=mark, ymin=LB, ymax=UB, colour=ic80g1, fill=ic80g1), data=fit.te[[i]], alpha=0.2) +
        geom_line(aes(x=mark, y=TE, colour=ic80g1), data=fit.te[[i]], na.rm=TRUE, size=1.9)
    }
    
    legendLabels <- c(expression("" <= 1 ~ mu * "g/ml"), expression("" > 1 ~ mu * "g/ml"))
    
    p1 <- p1 +
      scale_colour_manual(values=c("#159724", "#CD3333"), name=expression("Measured" ~ IC[80]), labels=legendLabels, guide=guide_legend(ncol=2, title.position="left")) +
      scale_fill_manual(values=c("#159724", "#CD3333"), name=expression("Measured" ~ IC[80]), labels=legendLabels, guide=guide_legend(ncol=2, title.position="left")) +
      scale_y_continuous(name="Prevention Efficacy (%)", breaks=c(seq(-0.4, 1, 0.2), 1.2, 1.4, 1.6, 1.8), 
                         labels=c(seq(-0.4, 1, 0.2) * 100, "Placebo", "VRC01", "Placebo", "VRC01")) +
      coord_cartesian(ylim=limits.y) +
      xlab(xlab) +
      theme_bw() +
      #geom_hline(yintercept=0, colour="grey", size = 1) +
      # annotate(geom = "text", x = min(limits.x), y = text.y, 
      #          label = sieveTestLabel, 
      #          size = 4, hjust = 0)+
      theme(legend.key.size = unit(0.65, "cm"),
            legend.title=element_text(size = 15),
            legend.margin=margin(grid::unit(0,"cm")),
            legend.text=element_text(size=14),
            legend.position = "bottom",
            #legend.justification = c(0.5, 1),
            legend.key = element_blank(),
            legend.key.width = unit(1.4,"cm"),
            legend.background = element_blank(),
            plot.title = element_text(hjust = 0.5, vjust = 2, size = 22),
            axis.title.x = element_text(size = 19, vjust = -0.5),
            axis.title.y = element_text(size = 19, hjust = 0.5),
            axis.text.x = element_text(size = 17, hjust = 0.5,vjust = 0.5, colour = "black"),
            axis.text.y = element_text(size = 17, colour = "black", hjust = 0.5),
            plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"))
    
  } else {
    stop("'fit.te' must be a data frame or a list of data frames.")
  }
  
  if (is.null(breaks.x)){
    p1 <- p1 + scale_x_continuous(limits = limits.x, breaks = pretty_breaks(n = 4))
  } else {
    p1 <- p1 + scale_x_continuous(limits = limits.x, breaks = breaks.x, labels=labels.x)
  }
  
  data_c <- drop_na(data) #Remove observations with missing mark
  nCasesVaccine <- length(data_c$treatment[data_c$treatment!= "Placebo"])
  nCasesPlacebo <- length(data_c$treatment[data_c$treatment == "Placebo"])
  quantilesVaccine <- quantile(data_c$markValue[data_c$treatment != "Placebo"], probs = c(0, 0.25, 0.5, 0.75, 1))
  quantilesPlacebo <- quantile(data_c$markValue[data_c$treatment== "Placebo"], probs = c(0, 0.25, 0.5, 0.75, 1))
  quantilesVaccine <- sapply(quantilesVaccine, function(x)format(as.numeric(x), digits = 2, nsmall = 2))
  quantilesPlacebo <- sapply(quantilesPlacebo, function(x)format(as.numeric(x), digits = 2, nsmall = 2))
  table = data.frame(rbind(c("VRC01",nCasesVaccine,quantilesVaccine ), c("Placebo",nCasesPlacebo, quantilesPlacebo)))
  colnames(table) <- c("Treatment", "No. of Cases", "Min", "Q1", "Median", "Q3", "Max")
  
  if (top.panel.type=="box"){
    data <- data %>% 
      mutate(ic80g1=factor(ic80g1, levels=0:1),
             txIC80=factor(interaction(treatment, ic80g1), levels=c("Placebo.0", "VRC01.0", "Placebo.1", "VRC01.1")))
    #add_row(markValue = NA, ic80g1=NA, treatment=NA, txIC80="Numerical Summary") %>%
    
    p1 <- p1 + 
      ggtitle(title) +
      geom_boxplot(aes(x=markValue, y=1 + 0.2 * as.numeric(txIC80), colour=txIC80), data=data, outlier.shape=NA, na.rm=TRUE) +
      geom_jitter(aes(x=markValue, y=1 + 0.2 * as.numeric(txIC80), colour=txIC80), data=data, na.rm=TRUE, width=0, height=0.08, shape=21, fill="white", size=2, stroke=1.3) 
      # scale_colour_manual(values=c("#159724", "#CD3333"), guide="none")
      #coord_flip() +
      #annotate(geom="table", x="Numerical Summary", y=min(data$markValue, na.rm=TRUE), label=list(table), vjust=0, fill="gray95") +
      #scale_x_discrete(expand= expansion(add = c(0.6, 2)))+
      #scale_y_discrete(name=NULL, labels=c("Placebo.0"="Placebo", "VRC01.0"="VRC01", "Placebo.1"="Placebo", "VRC01.1"="VRC01")) +
      # theme(plot.title=element_text(hjust=0.5, size=22),
      #       axis.text.y=element_text(size=12, colour="black"),
      #       axis.ticks.x=element_blank(),
      #       plot.margin=unit(c(0.2,0.2,0,0.2), "cm")) 
    

        # if (is.null(breaks.x)){
    #   p2 <- p2 + scale_y_continuous(name=NULL, limits = limits.x, breaks = pretty_breaks(n = 4), labels=NULL)
    # } else {
    #   p2 <- p2 + scale_y_continuous(name=NULL, limits = limits.x, breaks = breaks.x, labels=NULL)
    # }
  }
  
  # p <- ggarrange(p2, p1, heights=c(1, 2), ncol=1, nrow=2, align="v")
  # 
  # if (!is.null(title)){
  #   p <- annotate_figure(p, top=text_grob(title, color="black", face="bold", size=14, vjust=2.3))  
  # }
  
  # p <- p2 + p1 + plot_layout(ncol=1, heights=c(1, 2))
  
  return(p1)
}



#' Sieve PH plot with variants and variant VE annotations
#' @param fit.te an data frame returned by summary.sievePH 
#' @param data a data frame specifying a continuous mark under column "markValue", 
#'  a categorical variable describing treatment (Vaccine/Placebo), and a column under "variant" specifying variant
#' For subjects with a right-censored time-to-event, the value(s) in mark should be set to NA.
#' @param top.panel.type options are violin plot ("violin") and rcdf plot ("rcdf)
#' @param limits.x the range of the mark to be plotted
#' @param breaks.x the values where ticks will be placed for the mark
#' @param n.pretty_breaks if breaks.x is NULL, the number of breaks to be placed using function pretty_breaks() from R package "scales"
#' @param xlab the x-axis label
#' @param sieveTestLabel a character string describing sieve test p-values
#' @param title title for the plot
#' @variantFit a list vith variant VEs; if null, no variant VE will be overlaid with the VE curves
#' @pad.x.left a numerical value to padd space to the left of the x-axis so that the confidence interval bar of VE can be plotted at the smallest mark value
#' @pad.x.right a numerical value to padd space to the right of the x-axis so that the confidence interval bar of VE can be plotted at the smallest mark value

plot.sievePH.w.variant <- function(fit.te, data,limits.x, breaks.x, n.pretty_breaks = 4,lab.x, sieveTestLabel,title, variantFit = NULL, 
                                   pad.x.left, pad.x.right = 0, jitter.width = 0.3){
  cols <-  c("gray50", "gold", "purple", "turquoise", "green", "blue", "salmon", "darkred", "red", "black")
  #cols <- c("gray20", "orange", "darkgreen", "skyblue4", "yellow", "midnightblue", "darkred", "olivedrab4", "purple4","green3")
  variants.list <- c("Ancestral\nLineage", "Alpha", "Beta", "Gamma", "Delta","Zeta", "Iota", "Epsilon", "Lambda", "Mu")
  names(cols) <- variants.list 
  
  lengend.y = max(fit.te$UB) +0.01
  text.y = max(fit.te$UB)*100 +30
  
  limits.x2 = limits.x
  limits.x2[1] = limits.x2[1]-pad.x.left
  limits.x2[2] = limits.x2[2]+pad.x.right
  
  
  legend.y.p1 = 0.92
  legend.x.p1 = 0.05
  legend.x.just.p1=c(0, 1)
  text.y = 127
  
  p1 <- ggplot(fit.te)+
    geom_line(aes(x = markValues, y = TE*100), linetype = "solid", size = 1.7)+
    geom_line(aes(x = markValues, y = LB*100, linetype = '95% Pointwise CI', size = '95% Pointwise CI'))+
    geom_line(aes(x = markValues, y = UB*100, linetype = '95% Pointwise CI', size = '95% Pointwise CI'))+
    theme_bw()+
    scale_linetype_manual(name="",
                          labels = c('95% Pointwise CI'),
                          values = c('95% Pointwise CI'= "dashed"))+
    scale_size_manual(name="",
                      labels = c('95% Pointwise CI'),
                      values = c('95% Pointwise CI'= 1.5))+
    scale_y_continuous(limits = c(-40, text.y+10), breaks = c(-20, 0, 20, 40, 60, 80, 100))+
    xlab (lab.x)+
    ylab ("Vaccine Efficacy (%)")+
    #geom_hline(yintercept=0, colour="grey", size = 1) +
    annotate(geom = "text", x = min(limits.x), y = text.y, 
             label = sieveTestLabel, 
             size = 4, hjust = 0)+
    theme(legend.key.size = unit(0.65, "cm"),
          legend.title=element_text(size = 15),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(size=11),
          legend.position = c(legend.x.p1,legend.y.p1),
          legend.justification = c(legend.x.just.p1, 0),
          legend.key = element_blank(),
          legend.key.width = unit(1.4,"cm"),
          legend.background = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 2, size = 20),
          axis.title.x = element_text(size = 15, vjust = -0.5),
          axis.title.y = element_text(size = 15, hjust = 0.35),
          axis.text.x = element_text(size = 14, hjust = 0.5,vjust = 0.5, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),
          plot.margin=unit(c(-0.5,1,1,1), "cm"))
  if(is.null(breaks.x)){
    
    p1 <- p1 +scale_x_continuous(limits = limits.x2, breaks = pretty_breaks(n = n.pretty_breaks))
  }else{
    p1 <- p1 +scale_x_continuous(limits = limits.x2, breaks = breaks.x)
  }
  if(!is.null(variantFit)){
    summary = variantFit$summary
    VEtable <- data.frame(variant = summary$variants.fitted$variant, label = summary$variants.fitted$label, 
                          VE = summary$VE, CIL = summary$CIL, CIU = summary$CIU)
    VEtable <- filter(VEtable, variant != "other")
    #median of the mark in the placebo group restricting to each variant
    medianMark <- drop_na(ddply(filter(data, treatment == "Placebo"), .(variant), function(df) median(df$markValue,na.rm = TRUE)))
    colnames(medianMark) <- c("variant", "medianMark")
    VEtable <- left_join(VEtable, medianMark, by = "variant")
    VEtable$variant[VEtable$variant == "Ancestral.Lineage"] = "Ancestral\nLineage"
    variants.list.plot <- variants.list [variants.list %in% VEtable$variant]
    cols.plot <- cols[variants.list.plot]
    p1 <- p1 + geom_segment(aes(x = medianMark, y = CIL*100, xend = medianMark, yend = CIU*100, colour = variant), size = 1, alpha = 1, data = VEtable)+
      scale_colour_manual(name = "", labels = variants.list.plot , values = cols.plot, guide="none")
    barwidth = 0.22/30*max(limits.x)
    p1 <- p1 + geom_segment(aes(x = medianMark-barwidth, y = CIL*100, xend = medianMark+barwidth, yend = CIL*100, colour = variant), size = 1, alpha = 1, data = VEtable)
    p1 <- p1 + geom_segment(aes(x = medianMark-barwidth, y = CIU*100, xend = medianMark+barwidth, yend = CIU*100, colour = variant), size = 1, alpha = 1, data = VEtable)
    p1 <- p1 + geom_point(aes(x = medianMark, y = VE*100, colour = variant), size = 1.5, alpha = 1, data = VEtable)
  }
  
  
  
  data_c <- drop_na(data) #Remove observations with missing mark
  nCasesVaccine <- length(data_c$treatment[data_c$treatment!= "Placebo"])
  nCasesPlacebo <- length(data_c$treatment[data_c$treatment == "Placebo"])
  table = data.frame(rbind(c("Vaccine",nCasesVaccine), c("Placebo",nCasesPlacebo)))
  colnames(table) <- c("Treatment", "No. of Cases")
  data_c$variant[data_c$variant == "Ancestral.Lineage"] <- "Ancestral\nLineage"
  
  
  data_c <- add_row(. = data_c, markValue = NA, treatment = "Numerical Summary")
  data_c$treatment <- factor(data_c$treatment, levels = c("Placebo", "Vaccine", "Numerical Summary"))
  data_c2 <- drop_na(data_c)
  
  #browser()
  variants.list.plot <- variants.list[variants.list %in%unique(data_c2$variant)]
  cols.plot <- cols[variants.list.plot]
  p2 <- ggplot(data_c)+
    geom_boxplot(aes(x = treatment, y = markValue), outlier.shape = NA, data = data_c)+
    geom_jitter(aes(x = treatment, y = markValue, col = variant), shape=1, width = jitter.width, height = 0, fill = "white",  stroke = 0.8,data = data_c2)+
    coord_flip()+
    annotate(geom="table", x = "Numerical Summary", y = min(limits.x), label = list(table), vjust = 0.5, fill = "gray95")+
    xlab("")+
    ylab("")+
    #scale_x_discrete(expand= expansion(add = c(0.6, 2)))+
    #scale_y_continuous(limits = limits.x2, breaks = breaks.x, labels = NULL)+
    scale_x_discrete(labels = c("Numerical Summary" = "", "Vaccine" = "Vaccine", "Placebo" = "Placebo"))+
    theme_bw()+
    scale_colour_manual(name = "", labels = variants.list.plot , values = cols.plot)+
    guides(color=guide_legend(ncol=3, bycol = TRUE,override.aes = list(size=1.5, shape = 1, stroke = 2)))+
    
    theme(axis.text.y = element_text(size=15, colour = "black", angle = 90, hjust = 0.5),
          #axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin=unit(c(1,1,-0.12,1), "cm"),
          legend.position = c(0.7,0.85),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.key.size = unit(0.5, "cm"),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(size=10),
          legend.background = element_blank()
    )
  if(is.null(breaks.x)){
    p2 <- p2 +scale_y_continuous(limits = limits.x2, breaks = pretty_breaks(n = n.pretty_breaks),labels = NULL)
  }else{
    p2 <- p2 +scale_y_continuous(limits = limits.x2, breaks = breaks.x, labels = NULL)
  }
  
  p <- ggarrange(p2,p1,heights = c(1.5,2),ncol=1,nrow=2,align = "v")
  p <- annotate_figure(p, top = text_grob(title, 
                                          color = "black", face = "bold", size = 14, vjust = 2.3))
  
  return(p)
  
}
