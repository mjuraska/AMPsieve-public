# Purpose: Collection of R functions sourced by multiple other R files in subfolders of t:/vaccine/p704/analysis/efficacy/code
# Author:  Michal Juraska

library(survival)
source("t:/vaccine/p704/analysis/sieve/code/lunnMcneil.R")

# 'markHazVE' estimates either the unstratified or stratified mark-specific hazard-based TE, TE(j), j=1,2, from the competing risks Cox model; 
# it returns a point estimate, 95% CI, and Wald test p-value for {H0: lambda(t,j|Z=1) = lambda(t,j|Z=0), t>0}
# INPUT ARGUMENTS:
# data            - a data frame with columns ftime, fstatus, ftype (positive integers; cencode=0), tx (tx=1, P=0), and, optionally, strataVar
# ftypeCodeVector - a numeric vector of labels of marks, e.g., 1:2 for a binary mark
# stratified      - shall the stratified Cox model be used?
markHazTE <- function(data, ftypeCodeVector=1:2, stratified=FALSE, saveDir=NULL, saveFile=NULL){
  out <- list()
  runLunnMcneilTest <- TRUE
  
  for (i in 1:length(ftypeCodeVector)){
    ftypeCode <- ftypeCodeVector[i]
    data$fstatusGt <- as.numeric(data$ftype==ftypeCode)
    
    # numbers of tx and placebo cases for plotting purposes
    nCases <- as.numeric(tapply(data$fstatusGt, data$tx, sum, na.rm=TRUE))
    out[[i]] <- list(ftypeCode=ftypeCode)
    out[[i]]$nCasesPla <- nCases[1]
    out[[i]]$nCasesTx <- nCases[2]
    
    if (sum(nCases)==0){
      runLunnMcneilTest <- FALSE
      out[[i]]$logHazRatio <- NA
      out[[i]]$varLogHazRatio <- NA
      out[[i]]$hazTE <- NA
      out[[i]]$LBhazTE <- NA
      out[[i]]$UBhazTE <- NA
      out[[i]]$pWald <- NA
    } else {
      if (stratified){
        cox <- coxph(Surv(ftime, fstatusGt) ~ tx + strata(strataVar), data=data)
      } else {
        cox <- coxph(Surv(ftime, fstatusGt) ~ tx, data=data)
      }      
      scox <- summary(cox)$coef
      out[[i]]$logHazRatio <- scox[1]
      out[[i]]$varLogHazRatio <- scox[3]^2
      out[[i]]$hazTE <- 1 - scox[2]
      out[[i]]$LBhazTE <- 1 - exp(scox[1] + qnorm(0.975)*scox[3])
      out[[i]]$UBhazTE <- 1 - exp(scox[1] - qnorm(0.975)*scox[3])
      out[[i]]$pWald <- scox[5]
    }          
  }
  
  if (length(ftypeCodeVector)>1){
    if (runLunnMcneilTest){
      sfit <- lunnMcneilTest3gt(data$ftime, data$fstatus, ifelse(data$ftype>0, data$ftype - 1, data$ftype), data$tx, stratVar=data$strataVar, stratified=stratified)
      
      out$interactionCoefLunnMcNeil <- sfit$coef[3, 1]
      out$varInteractionCoefLunnMcNeil <- sfit$coef[3, 3]^2
      out$pLunnMcNeil <- sfit$coef[3, 5]
    } else {
      out$interactionCoefLunnMcNeil <- NA
      out$varInteractionCoefLunnMcNeil <- NA
      out$pLunnMcNeil <- NA
    }
  }
  
  if (!is.null(saveDir) && !is.null(saveFile)){ 
    save(out, file=file.path(saveDir, saveFile))    
  }
  
  return(invisible(out))
}

# 'markHazTE.MI' returns a list with MI inference about mark-specific hazard-based TE, TE(j), including the MI version of the Wald test p-value 
# evaluating {H0: lambda(time,j|Z=0)=lambda(time,j|Z=1), time>0} for each j=1,...,J, where J is an ordered categorical mark; 
# the output includes the MI version of the Lunn and McNeil test p-value;
# MI for the Lunn and McNeil test of {H0: hazTE(j) constant in j} is applied at the Cox model parameter estimator stage (using the data duplication method);
# INPUT ARGUMENTS
# ftypeCodeVector           - a numeric vector of labels of marks, e.g., 1:2 for a binary mark
# loadFile                  - an .RData file containing a list, each component of which is the output list from markHazTE()
# saveFile                  - if NULL and saveDir is provided, saveFile will be generated
# dir                       - the input and output directory
markHazTE.MI <- function(ftypeCodeVector=1:2, loadFile, saveFile=NULL, dir){
  # a list named 'fitMI', each component of which is the output from markHazTE()
  load(file.path(dir, loadFile))
  
  # number of imputations
  M <- length(fitMI)
  
  out <- as.list(NULL)
  for (i in 1:length(ftypeCodeVector)){    
    nCasesPlaVector <- sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$nCasesPla }, i=i)
    nCasesTxVector <- sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$nCasesTx }, i=i)
    
    # MI on the log hazard ratio scale
    logHazRatio.MI <- mean(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$logHazRatio }, i=i)) 
    
    # within-imputation variance
    meanVarLogHazRatio <- mean(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$varLogHazRatio }, i=i))
    
    # between-imputation variance
    sampleVarLogHazRatio <- var(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$logHazRatio }, i=i))
    
    # Rubin's rule for the estimated total variance (e.g., Schafer, SMMR 1999)
    varLogHazRatio.MI <- meanVarLogHazRatio + (1 + 1/M) * sampleVarLogHazRatio
    
    # degrees of freedom of the t-distribution of the overall estimate (e.g., Schafer, SMMR 1999)
    df <- (M - 1) * (1 + meanVarLogHazRatio / ((1 + 1/M) * sampleVarLogHazRatio))^2
    
    # 95% confidence bounds for the log hazard ratio
    LBlogHazRatio.MI <- logHazRatio.MI - qt(0.975, df=df) * sqrt(varLogHazRatio.MI)
    UBlogHazRatio.MI <- logHazRatio.MI + qt(0.975, df=df) * sqrt(varLogHazRatio.MI)
    
    # 2-sided Wald test (using t-distribution of the overall estimate)
    waldStat.MI <- logHazRatio.MI / sqrt(varLogHazRatio.MI)
    pWald.MI <- 2 * pt(-abs(waldStat.MI), df=df)

    out[[i]] <- list(ftypeCode=ftypeCodeVector[i])
    out[[i]]$nCasesPla <- round(mean(nCasesPlaVector), digits=0)
    out[[i]]$nCasesTx <- round(mean(nCasesTxVector), digits=0)    
    out[[i]]$logHazRatio <- logHazRatio.MI
    out[[i]]$varLogHazRatio <- varLogHazRatio.MI
    out[[i]]$hazTE <- 1-exp(logHazRatio.MI)
    out[[i]]$LBhazTE <- 1 - exp(UBlogHazRatio.MI)
    out[[i]]$UBhazTE <- 1 - exp(LBlogHazRatio.MI)    
    out[[i]]$pWald <- pWald.MI
  }
  
  # Lunn and McNeil test (1995, Biometrics)
  if (length(ftypeCodeVector)>1){
    # MI on the log hazard ratio scale for the interaction term
    logHazRatio.MI <- mean(sapply(fitMI, "[[", "interactionCoefLunnMcNeil")) 
    
    # within-imputation variance
    meanVarLogHazRatio <- mean(sapply(fitMI, "[[", "varInteractionCoefLunnMcNeil"))
    
    # between-imputation variance
    sampleVarLogHazRatio <- var(sapply(fitMI, "[[", "interactionCoefLunnMcNeil"))
    
    # Rubin's rule for the estimated total variance (e.g., Schafer, SMMR 1999)
    varLogHazRatio.MI <- meanVarLogHazRatio + (1 + 1/M) * sampleVarLogHazRatio
    
    # degrees of freedom of the t-distribution of the overall estimate (e.g., Schafer, SMMR 1999)
    df <- (M - 1) * (1 + meanVarLogHazRatio / ((1 + 1/M) * sampleVarLogHazRatio))^2
    
    lunnMcNeilStat.MI <- logHazRatio.MI / sqrt(varLogHazRatio.MI)
    out$pLunnMcNeil <- 2 * pt(-abs(lunnMcNeilStat.MI), df=df)
  } else {
    out$pLunnMcNeil <- NA
  }
  
  if (!is.null(saveFile)){ 
    save(out, file=file.path(dir, saveFile))
    cat("Output saved in:\n", file.path(dir, saveFile), "\n\n")
  }
  
  return(invisible(out))
}

# plot1est() is an internal function for plotMarkHazTE()
# INPUT ARGUMENTS
# fitTE      - a list outputted by markHazTE.MI()
# ftypeCode  - specifies the component (failure type) of fitTE
plot1est <- function(fitTE, ftypeCode, xBase, xSpread, xSide="left", yLim1=-1, yPvalAt=0.9, color="black"){
  if (fitTE[[ftypeCode]]$LBhazTE >= yLim1){
    segments(x0=xBase, y0=fitTE[[ftypeCode]]$LBhazTE, y1=fitTE[[ftypeCode]]$UBhazTE, lwd=2.5, col=color)
    segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$LBhazTE, lwd=2.5, col=color)
    segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$UBhazTE, lwd=2.5, col=color)
    
    text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$UBhazTE, 
         format(round(fitTE[[ftypeCode]]$UBhazTE * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
    text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$LBhazTE, 
         format(round(fitTE[[ftypeCode]]$LBhazTE * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
    text(xBase, yPvalAt, ifelse(fitTE[[ftypeCode]]$pWald<0.001, "p<0.001", paste0("p=", format(round(fitTE[[ftypeCode]]$pWald, 2), nsmall=2))), pos=1, cex=1.1)
    
  } else {
    segments(x0=xBase, y0=yLim1, y1=fitTE[[ftypeCode]]$UBhazTE, lwd=2.5, col=color)
    segments(x0=xBase - xSpread/3, x1=xBase + xSpread/3, y0=yLim1 - 0.05 * abs(yLim1), y1=yLim1 + 0.05 * abs(yLim1), lwd=2.5, col=color)
    segments(x0=xBase - xSpread/3, x1=xBase + xSpread/3, y0=yLim1 - 0.1 * abs(yLim1), y1=yLim1, lwd=2.5, col=color)
    segments(x0=xBase, y0=yLim1 - 0.05 * abs(yLim1), y1=yLim1 - 0.15 * abs(yLim1), lwd=2.5, col=color)
    segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=yLim1 - 0.15 * abs(yLim1), lwd=2.5, col=color)
    segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$UBhazTE, lwd=2.5, col=color)
    
    text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$UBhazTE, 
         format(round(fitTE[[ftypeCode]]$UBhazTE * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
    text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, yLim1 - 0.15 * abs(yLim1), 
         format(round(fitTE[[ftypeCode]]$LBhazTE * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
    text(xBase, yPvalAt, ifelse(fitTE[[ftypeCode]]$pWald<0.001, "p<0.001", paste0("p=", format(round(fitTE[[ftypeCode]]$pWald, 2), nsmall=2))), pos=1, cex=1.1)
  }
  
  if (fitTE[[ftypeCode]]$hazTE >= yLim1){
    points(xBase, fitTE[[ftypeCode]]$hazTE, pch=19, cex=1.3, col=color)  
    text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$hazTE, 
         format(round(fitTE[[ftypeCode]]$hazTE*100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
  }
}

# loadFile   - a character vector
# xTickLab   - a character vector of failure type labels along the x-axis
# panelLab   - a character vector of the same length as 'loadFile'
plotMarkHazTE <- function(loadFile, dir, xLab, yLab, xTickLab, panelLab=NULL, panelLabLine=-1, yLim=c(-0.4, 1), xLabLine=4.7, parMar=c(6.2, 6.2, 1.5, 1.5), 
                          xTxLabAt=0.55, yPvalAt=0.8){
  nPanels <- length(loadFile)
  
  xGap <- 0.25
  xSpread <- 0.05
  
  par(mar=parMar, las=1, cex.axis=1.3, cex.lab=1.5)
  plot(-1, 0, xlim=c(0.5, nPanels + 0.5), xaxt="n", ylim=yLim, yaxt="n", type="n", xlab="", ylab="", bty="l")
  
  axis(side=1, at=c(1:nPanels - xGap, 1:nPanels, 1:nPanels + xGap), labels=FALSE, tick=TRUE)
  axis(side=1, at=1:nPanels - xGap, labels=rep(xTickLab[1], nPanels), tick=FALSE, line=0.7, cex.axis=1.1)
  axis(side=1, at=1:nPanels, labels=rep(xTickLab[2], nPanels), tick=FALSE, line=0.7, cex.axis=1.1)
  axis(side=1, at=1:nPanels + xGap, labels=rep(xTickLab[3], nPanels), tick=FALSE, line=0.7, cex.axis=1.1)
  mtext(xLab, side=1, line=xLabLine, cex=1.5)
  axis(side=2, at=seq(-0.4, 1, by=0.2), labels=seq(-40, 100, by=20), cex.axis=1.3)
  mtext(yLab, side=2, line=3.3, cex=1.5, las=0)
  # axis(side=3, at=1:nPanels, labels=panelLab, tick=FALSE, line=panelLabLine, cex.axis=1.4)

  axis(side=2, at=yLim[1] + 0.07 * abs(yLim[1]), labels="No. of", tick=FALSE, cex.axis=1.1)
  axis(side=2, at=yLim[1] - 0.03 * abs(yLim[1]), labels="Cases", tick=FALSE, cex.axis=1.1)
  text(xTxLabAt, yLim[1] + 0.07 * abs(yLim[1]), labels="VRC01", cex=1.1)
  text(xTxLabAt, yLim[1] - 0.03 * abs(yLim[1]), labels="Placebo", cex=1.1)

  abline(h=seq(-0.4, 1, by=0.2), col="gray80", lwd=1)

  for (i in 1:nPanels){
    # a list named 'out', outputted by markHazTE()
    load(file.path(dir, loadFile[i]))

    if (out[[1]]$nCasesPla > 0 && out[[1]]$nCasesTx > 0){
      plot1est(fitTE=out, ftypeCode=1, xBase=i - xGap, xSpread=xSpread, xSide="left", yLim1=-0.6, yPvalAt=yPvalAt)
    }
    if (out[[2]]$nCasesPla > 0 && out[[2]]$nCasesTx > 0){
      plot1est(fitTE=out, ftypeCode=2, xBase=i, xSpread=xSpread, xSide="left", yLim1=-0.6, yPvalAt=yPvalAt)
    }
    if (out[[3]]$nCasesPla > 0 && out[[3]]$nCasesTx > 0){
      plot1est(fitTE=out, ftypeCode=3, xBase=i + xGap, xSpread=xSpread, xSide="right", yLim1=-0.6, yPvalAt=yPvalAt)
    }

    magF1 <- 1.01
    magF2 <- 1.06
    segments(x0=i - xGap, x1=i + xGap, y0=magF2 * yPvalAt)
    segments(x0=i - xGap, y0=magF1 * yPvalAt, y1=magF2 * yPvalAt)
    segments(x0=i, y0=magF1 * yPvalAt, y1=magF2 * yPvalAt)
    segments(x0=i + xGap, y0=magF1 * yPvalAt, y1=magF2 * yPvalAt)
    text(i, 1.1 * yPvalAt, ifelse(out$pLunnMcNeil<0.001, "2-sided p<0.001", paste0("2-sided p=", round(out$pLunnMcNeil, 2))), pos=3, offset=0.01, cex=1.1)

    text(i - xGap, yLim[1] + 0.07 * abs(yLim[1]), labels=out[[1]]$nCasesTx, cex=1.1)
    text(i - xGap, yLim[1] - 0.03 * abs(yLim[1]), labels=out[[1]]$nCasesPla, cex=1.1)
    text(i, yLim[1] + 0.07 * abs(yLim[1]), labels=out[[2]]$nCasesTx, cex=1.1)
    text(i, yLim[1] - 0.03 * abs(yLim[1]), labels=out[[2]]$nCasesPla, cex=1.1)
    text(i + xGap, yLim[1] + 0.07 * abs(yLim[1]), labels=out[[3]]$nCasesTx, cex=1.1)
    text(i + xGap, yLim[1] - 0.03 * abs(yLim[1]), labels=out[[3]]$nCasesPla, cex=1.1)
  }
}

# 'markHazTE.byDose' estimates either the unstratified or stratified mark-specific hazard-based TE30(j) and TE10(j), j=1,2, from the competing risks Cox model; 
# it returns a point estimate, 95% CI, and Wald test p-value for {H0: lambda(t,j|Z=1) = lambda(t,j|Z=0), t>0}
# INPUT ARGUMENTS:
# data            - a data frame with columns ftime, fstatus, ftype (positive integers; cencode=0), tx (tx=1, P=0), and, optionally, strataVar
# ftypeCodeVector - a numeric vector of labels of marks, e.g., 1:2 for a binary mark
# stratified      - shall the stratified Cox model be used?
markHazTE.byDose <- function(data, ftypeCodeVector=1:2, stratified=FALSE, saveDir=NULL, saveFile=NULL){
  out <- list()

  for (i in 1:length(ftypeCodeVector)){
    ftypeCode <- ftypeCodeVector[i]
    # numbers of tx and placebo cases for plotting purposes
    data$fstatusGt <- as.numeric(data$ftype==ftypeCode)
    
    nCases <- as.numeric(tapply(data$fstatusGt, data$tx, sum, na.rm=TRUE))
    out[[i]] <- list(ftypeCode=ftypeCode)
    out[[i]]$nCasesPla <- nCases[1]
    out[[i]]$nCasesLowDose <- nCases[2]
    out[[i]]$nCasesHighDose <- nCases[3]
    
    if (sum(nCases)==0){
      out[[i]]$logHazRatioLowDose <- NA
      out[[i]]$logHazRatioHighDose <- NA
      out[[i]]$varLogHazRatioLowDose <- NA
      out[[i]]$varLogHazRatioHighDose <- NA
      out[[i]]$hazTElowDose <- NA
      out[[i]]$hazTEhighDose <- NA
      out[[i]]$LBhazTElowDose <- NA
      out[[i]]$LBhazTEhighDose <- NA
      out[[i]]$UBhazTElowDose <- NA
      out[[i]]$UBhazTEhighDose <- NA
      out[[i]]$pWaldLowDose <- NA
      out[[i]]$pWaldHighDose <- NA
    } else {
      if (stratified){
        cox <- coxph(Surv(ftime, fstatusGt) ~ txLowDose + txHighDose + strata(strataVar), data=data)
      } else {
        cox <- coxph(Surv(ftime, fstatusGt) ~ txLowDose + txHighDose, data=data)
      }      
      scox <- summary(cox)$coef
      out[[i]]$logHazRatioLowDose <- scox[1, 1]
      out[[i]]$logHazRatioHighDose <- scox[2, 1]
      out[[i]]$varLogHazRatioLowDose <- scox[1, 3]^2
      out[[i]]$varLogHazRatioHighDose <- scox[2, 3]^2
      out[[i]]$hazTElowDose <- 1 - scox[1, 2]
      out[[i]]$hazTEhighDose <- 1 - scox[2, 2]
      out[[i]]$LBhazTElowDose <- 1 - exp(scox[1, 1] + qnorm(0.975) * scox[1, 3])
      out[[i]]$LBhazTEhighDose <- 1 - exp(scox[2, 1] + qnorm(0.975) * scox[2, 3])
      out[[i]]$UBhazTElowDose <- 1 - exp(scox[1, 1] - qnorm(0.975) * scox[1, 3])
      out[[i]]$UBhazTEhighDose <- 1 - exp(scox[2, 1] - qnorm(0.975) * scox[2, 3])
      out[[i]]$pWaldLowDose <- scox[1, 5]
      out[[i]]$pWaldHighDose <- scox[2, 5]
    }          
  }
  
  if (!is.null(saveDir) & !is.null(saveFile)){ 
    save(out, file=file.path(saveDir, saveFile))    
  }
  
  return(invisible(out))
}

# 'markHazTE.byDose.MI' returns a list with MI inference about mark-specific hazard-based TE30(j) and TE10(j), including the MI version of the Wald test p-value 
# evaluating {H0: lambda(time,j|Z=0)=lambda(time,j|Z=1), time>0} for each j=1,...,J, where J is an ordered categorical mark; 
# INPUT ARGUMENTS
# ftypeCodeVector           - a numeric vector of labels of marks, e.g., 1:2 for a binary mark
# loadFile                  - an .RData file containing a list, each component of which is the output list from markHazTE()
# saveFile                  - if NULL and saveDir is provided, saveFile will be generated
# dir                       - the input and output directory
markHazTE.byDose.MI <- function(ftypeCodeVector=1:2, loadFile, saveFile=NULL, dir){
  # a list named 'fitMI', each component of which is the output from markHazTE.byDose()
  load(file.path(dir, loadFile))
  
  # number of imputations
  M <- length(fitMI)
  
  out <- vector("list", length(ftypeCodeVector))
  for (i in 1:length(ftypeCodeVector)){    
    nCasesPlaVector <- sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$nCasesPla }, i=i)
    nCasesLowDoseVector <- sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$nCasesLowDose }, i=i)
    nCasesHighDoseVector <- sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$nCasesHighDose }, i=i)
    
    # MI on the log hazard ratio scale
    logHazRatioLowDose.MI <- mean(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$logHazRatioLowDose }, i=i))
    logHazRatioHighDose.MI <- mean(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$logHazRatioHighDose }, i=i)) 
    
    # within-imputation variance
    meanVarLogHazRatioLowDose <- mean(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$varLogHazRatioLowDose }, i=i))
    meanVarLogHazRatioHighDose <- mean(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$varLogHazRatioHighDose }, i=i))
    
    # between-imputation variance
    sampleVarLogHazRatioLowDose <- var(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$logHazRatioLowDose }, i=i))
    sampleVarLogHazRatioHighDose <- var(sapply(fitMI, function(fit1imp, i){ fit1imp[[i]]$logHazRatioHighDose }, i=i))
    
    # Rubin's rule for the estimated total variance (e.g., Schafer, SMMR 1999)
    varLogHazRatioLowDose.MI <- meanVarLogHazRatioLowDose + (1 + 1/M) * sampleVarLogHazRatioLowDose
    varLogHazRatioHighDose.MI <- meanVarLogHazRatioHighDose + (1 + 1/M) * sampleVarLogHazRatioHighDose
    
    # degrees of freedom of the t-distribution of the overall estimate (e.g., Schafer, SMMR 1999)
    dfLowDose <- (M - 1) * (1 + meanVarLogHazRatioLowDose / ((1 + 1/M) * sampleVarLogHazRatioLowDose))^2
    dfHighDose <- (M - 1) * (1 + meanVarLogHazRatioHighDose / ((1 + 1/M) * sampleVarLogHazRatioHighDose))^2
    
    # 95% confidence bounds for the log hazard ratio
    LBlogHazRatioLowDose.MI <- logHazRatioLowDose.MI - qt(0.975, df=dfLowDose) * sqrt(varLogHazRatioLowDose.MI)
    LBlogHazRatioHighDose.MI <- logHazRatioHighDose.MI - qt(0.975, df=dfHighDose) * sqrt(varLogHazRatioHighDose.MI)
    UBlogHazRatioLowDose.MI <- logHazRatioLowDose.MI + qt(0.975, df=dfLowDose) * sqrt(varLogHazRatioLowDose.MI)
    UBlogHazRatioHighDose.MI <- logHazRatioHighDose.MI + qt(0.975, df=dfHighDose) * sqrt(varLogHazRatioHighDose.MI)
    
    # 2-sided Wald test (using t-distribution of the overall estimate)
    waldStatLowDose.MI <- logHazRatioLowDose.MI / sqrt(varLogHazRatioLowDose.MI)
    waldStatHighDose.MI <- logHazRatioHighDose.MI / sqrt(varLogHazRatioHighDose.MI)
    pWaldLowDose.MI <- 2 * pt(-abs(waldStatLowDose.MI), df=dfLowDose)
    pWaldHighDose.MI <- 2 * pt(-abs(waldStatHighDose.MI), df=dfHighDose)
    
    out[[i]] <- list(ftypeCode=ftypeCodeVector[i])
    out[[i]]$nCasesPla <- round(mean(nCasesPlaVector), digits=0)
    out[[i]]$nCasesLowDose <- round(mean(nCasesLowDoseVector), digits=0)
    out[[i]]$nCasesHighDose <- round(mean(nCasesHighDoseVector), digits=0)
    out[[i]]$logHazRatioLowDose <- logHazRatioLowDose.MI
    out[[i]]$logHazRatioHighDose <- logHazRatioHighDose.MI
    out[[i]]$varLogHazRatioLowDose <- varLogHazRatioLowDose.MI
    out[[i]]$varLogHazRatioHighDose <- varLogHazRatioHighDose.MI
    out[[i]]$hazTElowDose <- 1 - exp(logHazRatioLowDose.MI)
    out[[i]]$hazTEhighDose <- 1 - exp(logHazRatioHighDose.MI)
    out[[i]]$LBhazTElowDose <- 1 - exp(UBlogHazRatioLowDose.MI)
    out[[i]]$LBhazTEhighDose <- 1 - exp(UBlogHazRatioHighDose.MI)
    out[[i]]$UBhazTElowDose <- 1 - exp(LBlogHazRatioLowDose.MI)
    out[[i]]$UBhazTEhighDose <- 1 - exp(LBlogHazRatioHighDose.MI)
    out[[i]]$pWaldLowDose <- pWaldLowDose.MI
    out[[i]]$pWaldHighDose <- pWaldHighDose.MI
  }
  
  if (!is.null(saveFile)){ 
    save(out, file=file.path(dir, saveFile))
    cat("Output saved in:\n", file.path(dir, saveFile), "\n\n")
  }
  
  return(invisible(out))
}

# plot1est.byDose() is an internal function for plotMarkHazTE.byDose()
# INPUT ARGUMENTS
# fitTE      - a list outputted by markHazTE.byDose.MI()
# ftypeCode  - specifies the component (failure type) of fitTE
plot1est.byDose <- function(fitTE, tx=c("lowdose", "highdose"), ftypeCode, xBase, xSpread, xSide="left", yLim1=-1, yPvalAt=0.9, color="black"){
  tx <- match.arg(tx)
  
  if (tx=="lowdose"){
    
    if (fitTE[[ftypeCode]]$LBhazTElowDose >= yLim1){
      segments(x0=xBase, y0=fitTE[[ftypeCode]]$LBhazTElowDose, y1=fitTE[[ftypeCode]]$UBhazTElowDose, lwd=2.5, col=color)
      segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$LBhazTElowDose, lwd=2.5, col=color)
      segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$UBhazTElowDose, lwd=2.5, col=color)
      
      text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$UBhazTElowDose, 
           format(round(fitTE[[ftypeCode]]$UBhazTElowDose * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
      text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$LBhazTElowDose, 
           format(round(fitTE[[ftypeCode]]$LBhazTElowDose * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
      text(xBase, yPvalAt, ifelse(fitTE[[ftypeCode]]$pWaldLowDose<0.001, "p<0.001", paste0("p=", format(round(fitTE[[ftypeCode]]$pWaldLowDose, 2), nsmall=2))), pos=1, cex=1.1)
      
    } else {
      segments(x0=xBase, y0=yLim1, y1=fitTE[[ftypeCode]]$UBhazTElowDose, lwd=2.5, col=color)
      segments(x0=xBase - xSpread/3, x1=xBase + xSpread/3, y0=yLim1 - 0.05 * abs(yLim1), y1=yLim1 + 0.05 * abs(yLim1), lwd=2.5, col=color)
      segments(x0=xBase - xSpread/3, x1=xBase + xSpread/3, y0=yLim1 - 0.1 * abs(yLim1), y1=yLim1, lwd=2.5, col=color)
      segments(x0=xBase, y0=yLim1 - 0.05 * abs(yLim1), y1=yLim1 - 0.15 * abs(yLim1), lwd=2.5, col=color)
      segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=yLim1 - 0.15 * abs(yLim1), lwd=2.5, col=color)
      segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$UBhazTElowDose, lwd=2.5, col=color)
      
      text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$UBhazTElowDose, 
           format(round(fitTE[[ftypeCode]]$UBhazTElowDose * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
      text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, yLim1 - 0.15 * abs(yLim1), 
           format(round(fitTE[[ftypeCode]]$LBhazTElowDose * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
      text(xBase, yPvalAt, ifelse(fitTE[[ftypeCode]]$pWaldLowDose<0.001, "p<0.001", paste0("p=", format(round(fitTE[[ftypeCode]]$pWaldLowDose, 2), nsmall=2))), pos=1, cex=1.1)
    }
    
    points(xBase, fitTE[[ftypeCode]]$hazTElowDose, pch=19, cex=1.3, col=color)
    text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$hazTElowDose, 
         format(round(fitTE[[ftypeCode]]$hazTElowDose*100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
    
  } else {
    
    if (fitTE[[ftypeCode]]$LBhazTEhighDose >= yLim1){
      segments(x0=xBase, y0=fitTE[[ftypeCode]]$LBhazTEhighDose, y1=fitTE[[ftypeCode]]$UBhazTEhighDose, lwd=2.5, col=color)
      segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$LBhazTEhighDose, lwd=2.5, col=color)
      segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$UBhazTEhighDose, lwd=2.5, col=color)
      
      text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$UBhazTEhighDose, 
           format(round(fitTE[[ftypeCode]]$UBhazTEhighDose * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
      text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$LBhazTEhighDose, 
           format(round(fitTE[[ftypeCode]]$LBhazTEhighDose * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
      text(xBase, yPvalAt, ifelse(fitTE[[ftypeCode]]$pWaldHighDose<0.001, "p<0.001", paste0("p=", format(round(fitTE[[ftypeCode]]$pWaldHighDose, 2), nsmall=2))), pos=1, cex=1.1)
      
    } else {
      segments(x0=xBase, y0=yLim1, y1=fitTE[[ftypeCode]]$UBhazTEhighDose, lwd=2.5, col=color)
      segments(x0=xBase - xSpread/3, x1=xBase + xSpread/3, y0=yLim1 - 0.05 * abs(yLim1), y1=yLim1 + 0.05 * abs(yLim1), lwd=2.5, col=color)
      segments(x0=xBase - xSpread/3, x1=xBase + xSpread/3, y0=yLim1 - 0.1 * abs(yLim1), y1=yLim1, lwd=2.5, col=color)
      segments(x0=xBase, y0=yLim1 - 0.05 * abs(yLim1), y1=yLim1 - 0.15 * abs(yLim1), lwd=2.5, col=color)
      segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=yLim1 - 0.15 * abs(yLim1), lwd=2.5, col=color)
      segments(x0=xBase - xSpread, x1=xBase + xSpread, y0=fitTE[[ftypeCode]]$UBhazTEhighDose, lwd=2.5, col=color)
      
      text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$UBhazTEhighDose, 
           format(round(fitTE[[ftypeCode]]$UBhazTEhighDose * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
      text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, yLim1 - 0.15 * abs(yLim1), 
           format(round(fitTE[[ftypeCode]]$LBhazTEhighDose * 100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
      text(xBase, yPvalAt, ifelse(fitTE[[ftypeCode]]$pWaldHighDose<0.001, "p<0.001", paste0("p=", format(round(fitTE[[ftypeCode]]$pWaldHighDose, 2), nsmall=2))), pos=1, cex=1.1)
    }
    
    points(xBase, fitTE[[ftypeCode]]$hazTEhighDose, pch=19, cex=1.3, col=color)
    text(xBase + ifelse(xSide=="left", -1, 1) * xSpread, fitTE[[ftypeCode]]$hazTEhighDose, 
         format(round(fitTE[[ftypeCode]]$hazTEhighDose*100, 1), nsmall=1), pos=ifelse(xSide=="left", 2, 4), cex=1.1)
  }
}


# loadFile   - a character string
# xTickLab   - a character vector of failure type labels along the x-axis
# panelLab   - a character vector of the same length as 'loadFile'
plotMarkHazTE.byDose <- function(loadFile, dir, xLab, yLab, xTickLab, xTickLabLine=0.7, panelLab=NULL, panelLabLine=-1, yLim=c(-0.4, 1), xLabLine=4.7, yLabLine=3.3,
                                 parMar=c(6.2, 6.2, 1.5, 1.5), xTxLabAt=0.55, yPvalAt=0.8, yTxLabLine=1){
  # a list named 'out', outputted by markHazTE.byDose.MI()
  load(file.path(dir, loadFile))
  
  # the number of panels is the number of treatment arms times the number of failure types
  nPanels <- 2 * length(out)
  
  xGap <- 0.18
  xSpread <- 0.05
  
  par(mar=parMar, las=1, cex.axis=1.3, cex.lab=1.5)
  plot(-1, 0, xlim=c(0.5, nPanels + 0.5), xaxt="n", ylim=yLim, yaxt="n", type="n", xlab="", ylab="", bty="l")
  
  axis(side=1, at=1:nPanels, labels=FALSE, tick=TRUE)
  axis(side=1, at=1:nPanels, labels=xTickLab, tick=FALSE, line=xTickLabLine, cex.axis=1.1)
  mtext(xLab, side=1, line=xLabLine, cex=1.5)
  axis(side=2, at=seq(-0.4, 1, by=0.2), labels=seq(-40, 100, by=20), cex.axis=1.3)
  mtext(yLab, side=2, line=yLabLine, cex=1.5, las=0)
  mtext(panelLab, side=3, line=panelLabLine, cex=1.5)
  
  axis(side=2, at=yLim[1] + yTxLabLine * 0.07 * abs(yLim[1]), labels="No. of", tick=FALSE, cex.axis=1.1)
  axis(side=2, at=yLim[1] - 0.03 * abs(yLim[1]), labels="Cases", tick=FALSE, cex.axis=1.1)
  text(xTxLabAt, yLim[1] + yTxLabLine * 0.07 * abs(yLim[1]), labels="VRC01", cex=1.1)
  text(xTxLabAt, yLim[1] - 0.03 * abs(yLim[1]), labels="Control", cex=1.1)
  
  abline(h=seq(-0.4, 1, by=0.2), col="gray80", lwd=1)
  
  # 'txLab', 'ftypeCode' and 'nCasesLab' determine the order of plotted results from left to right
  txLab <- rep(c("highdose", "lowdose"), 2)
  ftypeCode <- c(1, 1, 2, 2)
  nCasesTxLab <- rep(c("nCasesHighDose", "nCasesLowDose"), 2)
  for (i in 1:nPanels){
    plot1est.byDose(fitTE=out, tx=txLab[i], ftypeCode=ftypeCode[i], xBase=i, xSpread=xSpread, xSide="left", yLim1=-0.5, yPvalAt=yPvalAt)

    text(i, yLim[1] + yTxLabLine * 0.07 * abs(yLim[1]), labels=out[[ftypeCode[i]]][[nCasesTxLab[i]]], cex=1.1)
    text(i, yLim[1] - 0.03 * abs(yLim[1]), labels=out[[ftypeCode[i]]]$nCasesPla, cex=1.1)
  }
}


censorData <- function(data, cutoff){
  data$infustatus <- ifelse(data$infufudays > cutoff & data$infustatus==1, 0, data$infustatus)
  data$infufudays <- pmin(data$infufudays, cutoff)

  return(data)
}

censorData2 <- function(data, cutoff){
  data$hiv1event <- ifelse(data$hiv1survday > cutoff & data$hiv1event==1, 0, data$hiv1event)
  data$hiv1survday <- pmin(data$hiv1survday, cutoff)
  
  return(data)
}

censorData3 <- function(data, cutoff){
  data$statuswk104 <- ifelse(data$fudayswk104 > cutoff & data$statuswk104==1, 0, data$statuswk104)
  data$fudayswk104 <- pmin(data$fudayswk104, cutoff)
  
  return(data)
}


max2 <- function(x){
  x <- sort(x)
  return(x[length(x) - 1])
}

plotSmoothHaz <- function(data, txVar, txName, txLabel, xMax=NULL, yMax=NULL, showLegend=TRUE){
  source("t:/vaccine/p704/analysis/efficacy/code/confbandhazplusindiv.r")
  
  # standardize variable names
  data$tx <- data[, txVar]
  data$fTime <- data$infufudays
  data$fInd <- data$infustatus
  
  ngrid <- 100
  
  # set u1=t1 to the maximum of the smallest treatment-specific failure times
  t1 <- max(tapply(data$fTime, data$tx, min))
  
  # set u2=t2 to the minimum of the 2nd largest treatment-specific failure times
  t2 <- min(tapply(data$fTime, data$tx, max2))
  
  N <- 10000
  band2 <- band1 <- 0
  band11 <- (t2 - t1) / 2
  biasadjust <- tailsl <- tailsu <- TRUE
  Nv <- 0
  
  # revise format of 'tx' as required by 'Confbandhazplusindiv'
  # idx <- which(levels(data$tx)==txName)
  # if (length(idx)==0){ stop("There is no treatment label ", txName, " in the treatment variable ", txVar, ".") }
  # # reorder factor levels so that the factor level of interest is first
  # data$tx <- factor(data$tx, levels=c(levels(data$tx)[idx], levels(data$tx)[-idx]))
  # levels(data$tx) <- 1:length(levels(data$tx))
  # data$tx <- as.numeric(data$tx)
  
  data$tx <- 2 - as.numeric(data$tx==txName)
  
  # if there are 2 treatment arms
  # if (max(data$tx)==2){
  hazFit <- Confbandhazplusindiv(data$fTime, data$fInd, data$tx, N, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  # extract results for data$tx==1
  hazFit <- hazFit[[2]]
  # if there are 3 treatment arms
  # } else if (max(data$tx)==3){
  #   hazFit1 <- Confbandhazplusindiv(data$fTime[data$tx<3], data$fInd[data$tx<3], data$tx[data$tx<3], nBoot, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  #   hazFit2 <- Confbandhazplusindiv(data$fTime[data$tx>1], data$fInd[data$tx>1], data$tx[data$tx>1] - 1, nBoot, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  #   hazFit <- hazFit1[2:3]
  #   hazFit[[3]] <- hazFit2[[3]]
  # }

  #xMax <- max(hazFit[[1]][, 1], hazFit[[2]][, 1])
  if (is.null(xMax)){ xMax <- max(hazFit[, 1]) }
  if (is.null(yMax)){ yMax <- max(hazFit[, 6]) }
  # if (max(data$tx)==2){
  #   txLabels <- c("Placebo", "VRC01 10 + 30 mg/kg")
  # } else if (max(data$tx==3)){
  #   txLabels <- c("Placebo", "VRC01 10 mg/kg", "VRC01 30 mg/kg")
  # }
  
  par(mar=c(5, 6.5, 3, 1), las=1, cex.axis=1.2)
  cexPlotTitle <- 1.4
  cexAxisLabel <- 1.3
  cexLegend <- 1.2
  
  # for each treatment arm (in the order in 'txLabels')
  #for (i in ifelse(includePlacebo, 1, 2):max(data$tx)){
  timeGrid <- hazFit[, 1]
  haz <- pmax(hazFit[timeGrid<=xMax, 2], 0)
  ptLB <- pmax(hazFit[timeGrid<=xMax, 3], 0)
  ptUB <- hazFit[timeGrid<=xMax, 4]
  smLB <- pmax(hazFit[timeGrid<=xMax, 5], 0)
  smUB <- hazFit[timeGrid<=xMax, 6]
  timeGrid <- timeGrid[timeGrid<=xMax]
  
  plot(timeGrid, haz, xlab="", ylab="", type="n", xlim=c(0, xMax), ylim=c(0, yMax), xaxt="n", yaxt="n", bty="l")
  
  axis(1, at=seq(0, xMax, by=14), labels=seq(0, xMax, by=14) / 7)
  axis(1, at=seq(7, xMax, by=14), labels=seq(7, xMax, by=14) / 7)
  axis(2, at=seq(0, 6e-4, by=5e-5), labels=format(seq(0, 6e-4, by=5e-5), scientific=FALSE))
  
  mtext("Weeks since Last Infusion", side=1, line=2.5, cex=cexAxisLabel)
  mtext("Hazard of Primary Endpoint", side=2, line=5, las=0, cex=cexAxisLabel)
  mtext(txLabel, side=3, line=0.5, cex=cexPlotTitle)
  
  lines(timeGrid, haz, lwd=3.5)
  lines(timeGrid, ptLB, lty="dashed", lwd=3)
  lines(timeGrid, ptUB, lty="dashed", lwd=3)
  lines(timeGrid, smLB, lty="dotted", lwd=3)
  lines(timeGrid, smUB, lty="dotted", lwd=3)
  
  if (showLegend){
    legend("topleft", legend=c("95% pointwise CI","95% simultaneous CI"), cex=cexLegend, lwd=3, lty=c("dashed", "dotted"), bty="n", x.intersp=0.7)
  }
}

# same analysis as in plotSmoothHaz() except for time since enrollment as the failure time
# the input data set v704_survival_wk80_tau_neut.csv has different variable names than the input data set for plotSmoothHaz()
# the x-axis has a different range
plotSmoothHaz2 <- function(data, txVar, txName, txLabel, xMax=NULL, yMax=NULL, showLegend=TRUE){
  source("t:/vaccine/p704/analysis/efficacy/code/confbandhazplusindiv.r")
  
  # standardize variable names
  data$tx <- data[, txVar]
  data$fTime <- data$hiv1survday
  data$fInd <- data$hiv1event
  
  ngrid <- 100
  
  # set u1=t1 to the maximum of the smallest treatment-specific failure times
  t1 <- max(tapply(data$fTime, data$tx, min))
  
  # set u2=t2 to the minimum of the 2nd largest treatment-specific failure times
  t2 <- min(tapply(data$fTime, data$tx, max2))
  
  N <- 10000
  band2 <- band1 <- 0
  band11 <- (t2 - t1) / 2
  biasadjust <- tailsl <- tailsu <- TRUE
  Nv <- 0
  
  # revise format of 'tx' as required by 'Confbandhazplusindiv'
  # idx <- which(levels(data$tx)==txName)
  # if (length(idx)==0){ stop("There is no treatment label ", txName, " in the treatment variable ", txVar, ".") }
  # # reorder factor levels so that the factor level of interest is first
  # data$tx <- factor(data$tx, levels=c(levels(data$tx)[idx], levels(data$tx)[-idx]))
  # levels(data$tx) <- 1:length(levels(data$tx))
  # data$tx <- as.numeric(data$tx)
  
  data$tx <- 2 - as.numeric(data$tx==txName)
  
  # if there are 2 treatment arms
  # if (max(data$tx)==2){
  hazFit <- Confbandhazplusindiv(data$fTime, data$fInd, data$tx, N, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  # extract results for data$tx==1
  hazFit <- hazFit[[2]]
  # if there are 3 treatment arms
  # } else if (max(data$tx)==3){
  #   hazFit1 <- Confbandhazplusindiv(data$fTime[data$tx<3], data$fInd[data$tx<3], data$tx[data$tx<3], nBoot, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  #   hazFit2 <- Confbandhazplusindiv(data$fTime[data$tx>1], data$fInd[data$tx>1], data$tx[data$tx>1] - 1, nBoot, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  #   hazFit <- hazFit1[2:3]
  #   hazFit[[3]] <- hazFit2[[3]]
  # }
  
  #xMax <- max(hazFit[[1]][, 1], hazFit[[2]][, 1])
  if (is.null(xMax)){ xMax <- max(hazFit[, 1]) }
  if (is.null(yMax)){ yMax <- max(hazFit[, 6]) }
  # if (max(data$tx)==2){
  #   txLabels <- c("Placebo", "VRC01 10 + 30 mg/kg")
  # } else if (max(data$tx==3)){
  #   txLabels <- c("Placebo", "VRC01 10 mg/kg", "VRC01 30 mg/kg")
  # }
  
  par(mar=c(5, 6.5, 3, 1), las=1, cex.axis=1.2)
  cexPlotTitle <- 1.4
  cexAxisLabel <- 1.3
  cexLegend <- 1.2
  
  # for each treatment arm (in the order in 'txLabels')
  #for (i in ifelse(includePlacebo, 1, 2):max(data$tx)){
  timeGrid <- hazFit[, 1]
  haz <- hazFit[timeGrid<=xMax, 2]
  ptLB <- pmax(hazFit[timeGrid<=xMax, 3], 0)
  ptUB <- hazFit[timeGrid<=xMax, 4]
  smLB <- pmax(hazFit[timeGrid<=xMax, 5], 0)
  smUB <- hazFit[timeGrid<=xMax, 6]
  timeGrid <- timeGrid[timeGrid<=xMax]
  
  plot(timeGrid, haz, xlab="", ylab="", type="n", xlim=c(0, xMax), ylim=c(0, yMax), xaxt="n", yaxt="n", bty="l")
  
  axis(1, at=seq(0, xMax, by=16 * 7), labels=seq(0, xMax, by=16 * 7) / 7)
  #axis(1, at=seq(8 * 7, xMax, by=16 * 7), labels=seq(8 * 7, xMax, by=16 * 7) / 7)
  axis(2, at=seq(0, 3e-4, by=1e-4), labels=format(seq(0, 3e-4, by=1e-4), scientific=FALSE))
  
  mtext("Weeks since Enrollment", side=1, line=2.5, cex=cexAxisLabel)
  mtext("Hazard of Primary Endpoint", side=2, line=5, las=0, cex=cexAxisLabel)
  mtext(txLabel, side=3, line=0.5, cex=cexPlotTitle)
  
  lines(timeGrid, haz, lwd=3.5)
  lines(timeGrid, ptLB, lty="dashed", lwd=3)
  lines(timeGrid, ptUB, lty="dashed", lwd=3)
  lines(timeGrid, smLB, lty="dotted", lwd=3)
  lines(timeGrid, smUB, lty="dotted", lwd=3)
  
  if (showLegend){
    legend("topleft", legend=c("95% pointwise CI","95% simultaneous CI"), cex=cexLegend, lwd=3, lty=c("dashed", "dotted"), bty="n", x.intersp=0.7)
  }
}

# same analysis as in plotSmoothHaz2() except for use with the data set amp_survival.csv (Week 104 analyses), which has different variable names
# the x-axis has a different range
plotSmoothHaz3 <- function(data, txVar, txName, txLabel, xMax=NULL, yMax=NULL, showLegend=TRUE){
  source("t:/vaccine/p704/analysis/efficacy/code/confbandhazplusindiv.r")
  
  # standardize variable names
  data$tx <- data[, txVar]
  data$fTime <- data$fudayswk104
  data$fInd <- data$statuswk104
  
  ngrid <- 100
  
  # set u1=t1 to the maximum of the smallest treatment-specific failure times
  t1 <- max(tapply(data$fTime, data$tx, min))
  
  # set u2=t2 to the minimum of the 2nd largest treatment-specific failure times
  t2 <- min(tapply(data$fTime, data$tx, max2))
  
  N <- 10000
  band2 <- band1 <- 0
  band11 <- (t2 - t1) / 2
  biasadjust <- tailsl <- tailsu <- TRUE
  Nv <- 0
  
  # revise format of 'tx' as required by 'Confbandhazplusindiv'
  # idx <- which(levels(data$tx)==txName)
  # if (length(idx)==0){ stop("There is no treatment label ", txName, " in the treatment variable ", txVar, ".") }
  # # reorder factor levels so that the factor level of interest is first
  # data$tx <- factor(data$tx, levels=c(levels(data$tx)[idx], levels(data$tx)[-idx]))
  # levels(data$tx) <- 1:length(levels(data$tx))
  # data$tx <- as.numeric(data$tx)
  
  data$tx <- 2 - as.numeric(data$tx==txName)
  
  # if there are 2 treatment arms
  # if (max(data$tx)==2){
  hazFit <- Confbandhazplusindiv(data$fTime, data$fInd, data$tx, N, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  # extract results for data$tx==1
  hazFit <- hazFit[[2]]
  # if there are 3 treatment arms
  # } else if (max(data$tx)==3){
  #   hazFit1 <- Confbandhazplusindiv(data$fTime[data$tx<3], data$fInd[data$tx<3], data$tx[data$tx<3], nBoot, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  #   hazFit2 <- Confbandhazplusindiv(data$fTime[data$tx>1], data$fInd[data$tx>1], data$tx[data$tx>1] - 1, nBoot, 0, t1, t2, t1, t2, band11, ngrid, band1, band2, 0, 1, biasadjust, tailsl, tailsu, Nv)
  #   hazFit <- hazFit1[2:3]
  #   hazFit[[3]] <- hazFit2[[3]]
  # }
  
  #xMax <- max(hazFit[[1]][, 1], hazFit[[2]][, 1])
  if (is.null(xMax)){ xMax <- max(hazFit[, 1]) }
  if (is.null(yMax)){ yMax <- max(hazFit[, 6]) }
  # if (max(data$tx)==2){
  #   txLabels <- c("Placebo", "VRC01 10 + 30 mg/kg")
  # } else if (max(data$tx==3)){
  #   txLabels <- c("Placebo", "VRC01 10 mg/kg", "VRC01 30 mg/kg")
  # }
  
  par(mar=c(5, 6.5, 3, 1), las=1, cex.axis=1.2)
  cexPlotTitle <- 1.4
  cexAxisLabel <- 1.3
  cexLegend <- 1.2
  
  # for each treatment arm (in the order in 'txLabels')
  #for (i in ifelse(includePlacebo, 1, 2):max(data$tx)){
  timeGrid <- hazFit[, 1]
  haz <- hazFit[timeGrid<=xMax, 2]
  ptLB <- pmax(hazFit[timeGrid<=xMax, 3], 0)
  ptUB <- hazFit[timeGrid<=xMax, 4]
  smLB <- pmax(hazFit[timeGrid<=xMax, 5], 0)
  smUB <- hazFit[timeGrid<=xMax, 6]
  timeGrid <- timeGrid[timeGrid<=xMax]
  
  plot(timeGrid, haz, xlab="", ylab="", type="n", xlim=c(0, xMax), ylim=c(0, yMax), xaxt="n", yaxt="n", bty="l")
  
  axis(1, at=seq(0, xMax, by=16 * 7), labels=seq(0, xMax, by=16 * 7) / 7)
  axis(1, at=108 * 7, labels=108)
  axis(2, at=seq(0, 6e-4, by=1e-4), labels=format(seq(0, 6e-4, by=1e-4), scientific=FALSE))
  
  mtext("Weeks since Enrollment", side=1, line=2.5, cex=cexAxisLabel)
  mtext("Hazard of Primary Endpoint", side=2, line=5, las=0, cex=cexAxisLabel)
  mtext(txLabel, side=3, line=0.5, cex=cexPlotTitle)
  
  lines(timeGrid, haz, lwd=3.5)
  lines(timeGrid, ptLB, lty="dashed", lwd=3)
  lines(timeGrid, ptUB, lty="dashed", lwd=3)
  lines(timeGrid, smLB, lty="dotted", lwd=3)
  lines(timeGrid, smUB, lty="dotted", lwd=3)
  
  if (showLegend){
    legend("topleft", legend=c("95% pointwise CI","95% simultaneous CI"), cex=cexLegend, lwd=3, lty=c("dashed", "dotted"), bty="n", x.intersp=0.7)
  }
}


# txNames = character vector; if specified, the first component is assumed to be the placebo group label
plotSmoothHazPE <- function(data, txVar, txNames, title, xMax=NULL, showLegend=TRUE, nBoot=1000){
  source("t:/vaccine/p704/analysis/efficacy/code/confidenceband.loghazardratio.r")
  
  # standardize variable names
  data$tx <- data[, txVar]
  data$fTime <- data$infufudays
  data$fInd <- data$infustatus
  
  # keep only data for the specified comparison
  data <- subset(data, tx %in% txNames)
  data$tx <- factor(data$tx, levels=txNames)
  
  # set input parameters
  N <- 10000
  
  # set u1=t1 to the maximum of the smallest treatment-specific failure times
  t1 <- max(tapply(data$fTime, data$tx, min))
  
  # set u2=t2 to the minimum of the 2nd largest treatment-specific failure times
  t2 <- min(tapply(data$fTime, data$tx, max2))
  
  band12 <- band11 <- (t2 - t1) / 2
  
  # default settings for the rest of the parameters
  band1 <- 0
  band2 <- 0
  ngrid <- 50
  biasadjust <- TRUE
  tailsl <- TRUE
  tailsu <- TRUE
  Nv <- 0
  
  # revise format of 'tx' as required by 'confidenceband.loghazardratio'
  levels(data$tx) <- 1:length(levels(data$tx))
  data$tx <- as.numeric(data$tx)
  
  # if there are 2 treatment arms
  # if (max(data$tx)==2){
  logHRfit <- confidenceband.loghazardratio(data$fTime, data$fInd, data$tx, t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nBoot, Nv)
    # if there are 3 treatment arms
  # } else if (max(data$tx)==3){
  #   logHRfit <- list(confidenceband.loghazardratio(data$fTime[data$tx<3], data$fInd[data$tx<3], data$tx[data$tx<3], t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nboot, Nv))
  #   logHRfit[[2]] <- confidenceband.loghazardratio(data$fTime[data$tx!=2], data$fInd[data$tx!=2], as.numeric(data$tx[data$tx!=2]==3) + 1, t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nboot, Nv)
  # }
  
  timeGrid <- logHRfit[, 1]
  if (is.null(xMax)){ xMax <- max(timeGrid) }
  
  logHR <- logHRfit[timeGrid<=xMax, 2]
  ptLB <- logHRfit[timeGrid<=xMax, 3]
  ptUB <- logHRfit[timeGrid<=xMax, 4]
  smLB <- logHRfit[timeGrid<=xMax, 5]
  smUB <- logHRfit[timeGrid<=xMax, 6]
  timeGrid <- timeGrid[timeGrid<=xMax]
  
  yLim <- c(-1, 1.2)
  
  par(mar=c(5, 6.5, 3, 1), las=1, cex.axis=1.2)
  cexPlotTitle <- 1.4
  cexAxisLabel <- 1.3
  cexLegend <- 1.2
  
  plot(timeGrid, 1 - exp(logHR), xlab="", ylab="", type="n", xaxt="n", yaxt="n", xlim=c(0, xMax), ylim=yLim, bty="l")
  
  axis(1, at=seq(0, xMax, by=14), labels=seq(0, xMax, by=14) / 7)
  axis(1, at=seq(7, xMax, by=14), labels=seq(7, xMax, by=14) / 7)
  axis(2, at=seq(-1, 1, by=0.25), labels=seq(-100, 100, by=25))
  
  mtext("Weeks since Last Infusion", side=1, line=2.8, cex=cexAxisLabel)
  mtext("Instantaneous Hazard-Ratio PE (%)", side=2, line=3.3, cex=cexAxisLabel, las=0)
  mtext(title, side=3, line=0.5, cex=cexPlotTitle, las=0)
  
  abline(h=0, col="gray50", lwd=2)
  lines(timeGrid, 1 - exp(logHR), lwd=3.5)
  lines(timeGrid, 1 - exp(ptLB), lty="dashed", lwd=3)
  lines(timeGrid, 1 - exp(ptUB), lty="dashed", lwd=3)
  lines(timeGrid, 1 - exp(smLB), lty="dotted", lwd=3)
  lines(timeGrid, 1 - exp(smUB), lty="dotted", lwd=3)
  
  if (showLegend){
    legend("topleft", legend=c("95% pointwise CI","95% simultaneous CI"), cex=cexLegend, lwd=3, lty=c("dashed", "dotted"), bty="n", x.intersp=0.7)
  }
}

# txNames = character vector; if specified, the first component is assumed to be the placebo group label
# loadFile = the output .RData from confidenceband.loghazardratio() to be loaded
# saveFile = the output .RData from confidenceband.loghazardratio() to be saved
plotSmoothHazPE2 <- function(data, txVar, txNames, title, xMax=NULL, showLegend=TRUE, nBoot=NULL, loadFile=NULL, saveFile=NULL, saveDir=NULL){
  source("t:/vaccine/p704/analysis/efficacy/code/confidenceband.loghazardratio.r")
  
  # standardize variable names
  data$tx <- data[, txVar]
  data$fTime <- data$hiv1survday
  data$fInd <- data$hiv1event
  
  # keep only data for the specified comparison
  data <- subset(data, tx %in% txNames)
  data$tx <- factor(data$tx, levels=txNames)
  
  # set input parameters
  N <- 10000
  
  # set u1=t1 to the maximum of the smallest treatment-specific failure times
  t1 <- max(tapply(data$fTime, data$tx, min))
  
  # set u2=t2 to the minimum of the 2nd largest treatment-specific failure times
  t2 <- min(tapply(data$fTime, data$tx, max2))
  
  if (is.null(nBoot)){
    band12 <- band11 <- (t2 - t1) / 2  
  } else {
    band12 <- band11 <- 0  
  }
  
  # default settings for the rest of the parameters
  band1 <- 0
  band2 <- 0
  ngrid <- 100
  biasadjust <- TRUE
  tailsl <- TRUE
  tailsu <- TRUE
  Nv <- 0
  
  # revise format of 'tx' as required by 'confidenceband.loghazardratio'
  levels(data$tx) <- 1:length(levels(data$tx))
  data$tx <- as.numeric(data$tx)
  
  if (is.null(loadFile)){
    # if there are 2 treatment arms
    # if (max(data$tx)==2){
    logHRfit <- confidenceband.loghazardratio(data$fTime, data$fInd, data$tx, t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nBoot, Nv)  
    
    # save the output from confidenceband.loghazardratio() which is the bottleneck of this code
    if (!(is.null(saveFile) && is.null(saveDir))){
      save(logHRfit, file=file.path(saveDir, saveFile))  
    }
  } else if (!is.null(saveDir)){
    load(file.path(saveDir, loadFile))
  } else {
    stop("'saveDir' is missing and required to locate 'loadFile'.")
  }
  
  # if there are 3 treatment arms
  # } else if (max(data$tx)==3){
  #   logHRfit <- list(confidenceband.loghazardratio(data$fTime[data$tx<3], data$fInd[data$tx<3], data$tx[data$tx<3], t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nboot, Nv))
  #   logHRfit[[2]] <- confidenceband.loghazardratio(data$fTime[data$tx!=2], data$fInd[data$tx!=2], as.numeric(data$tx[data$tx!=2]==3) + 1, t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nboot, Nv)
  # }
  
  timeGrid <- logHRfit[, 1]
  if (is.null(xMax)){ xMax <- max(timeGrid) }
  
  logHR <- logHRfit[timeGrid<=xMax, 2]
  ptLB <- logHRfit[timeGrid<=xMax, 3]
  ptUB <- logHRfit[timeGrid<=xMax, 4]
  smLB <- logHRfit[timeGrid<=xMax, 5]
  smUB <- logHRfit[timeGrid<=xMax, 6]
  timeGrid <- timeGrid[timeGrid<=xMax]
  
  yLim <- c(-1, 1.2)
  
  par(mar=c(5, 6.5, 3, 1), las=1, cex.axis=1.2)
  cexPlotTitle <- 1.4
  cexAxisLabel <- 1.3
  cexLegend <- 1.2
  
  plot(timeGrid, 1 - exp(logHR), xlab="", ylab="", type="n", xaxt="n", yaxt="n", xlim=c(0, xMax), ylim=yLim, bty="l")
  
  axis(1, at=seq(0, xMax, by=16 * 7), labels=seq(0, xMax, by=16 * 7) / 7)
  #axis(1, at=seq(8 * 7, xMax, by=16 * 7), labels=seq(8 * 7, xMax, by=16 * 7) / 7)
  axis(2, at=seq(-1, 1, by=0.25), labels=seq(-100, 100, by=25))
  
  mtext("Weeks since Enrollment", side=1, line=2.8, cex=cexAxisLabel)
  mtext("Instantaneous Hazard-Ratio PE (%)", side=2, line=3.3, cex=cexAxisLabel, las=0)
  mtext(title, side=3, line=0.5, cex=cexPlotTitle, las=0)
  
  abline(h=0, col="gray50", lwd=2)
  lines(timeGrid, 1 - exp(logHR), lwd=3.5)
  lines(timeGrid, 1 - exp(ptLB), lty="dashed", lwd=3)
  lines(timeGrid, 1 - exp(ptUB), lty="dashed", lwd=3)
  lines(timeGrid, 1 - exp(smLB), lty="dotted", lwd=3)
  lines(timeGrid, 1 - exp(smUB), lty="dotted", lwd=3)
  
  if (showLegend){
    legend(0.2 * xMax, 1.2, legend=c("95% pointwise CI", "95% simultaneous CI"), cex=cexLegend, lwd=3, lty=c("dashed", "dotted"), bty="n", x.intersp=0.5, y.intersp=0.7)
  }
}

# txNames = character vector; if specified, the first component is assumed to be the placebo group label
# loadFile = the output .RData from confidenceband.loghazardratio() to be loaded
# saveFile = the output .RData from confidenceband.loghazardratio() to be saved
plotSmoothHazPE3 <- function(data, txVar, txNames, title, xMax=NULL, showLegend=TRUE, nBoot=NULL, loadFile=NULL, saveFile=NULL, saveDir=NULL){
  source("t:/vaccine/p704/analysis/efficacy/code/confidenceband.loghazardratio.r")
  
  # standardize variable names
  data$tx <- data[, txVar]
  data$fTime <- data$fudayswk104
  data$fInd <- data$statuswk104
  
  # keep only data for the specified comparison
  data <- subset(data, tx %in% txNames)
  data$tx <- factor(data$tx, levels=txNames)
  
  # set input parameters
  N <- 10000
  
  # set u1=t1 to the maximum of the smallest treatment-specific failure times
  t1 <- max(tapply(data$fTime, data$tx, min))
  
  # set u2=t2 to the minimum of the 2nd largest treatment-specific failure times
  t2 <- min(tapply(data$fTime, data$tx, max2))
  
  if (is.null(nBoot)){
    band12 <- band11 <- (t2 - t1) / 2  
  } else {
    band12 <- band11 <- 0  
  }
  
  # default settings for the rest of the parameters
  band1 <- 0
  band2 <- 0
  ngrid <- 100
  biasadjust <- TRUE
  tailsl <- TRUE
  tailsu <- TRUE
  Nv <- 0
  
  # revise format of 'tx' as required by 'confidenceband.loghazardratio'
  levels(data$tx) <- 1:length(levels(data$tx))
  data$tx <- as.numeric(data$tx)
  
  if (is.null(loadFile)){
    # if there are 2 treatment arms
    # if (max(data$tx)==2){
    logHRfit <- confidenceband.loghazardratio(data$fTime, data$fInd, data$tx, t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nBoot, Nv)  
    
    # save the output from confidenceband.loghazardratio() which is the bottleneck of this code
    if (!(is.null(saveFile) && is.null(saveDir))){
      save(logHRfit, file=file.path(saveDir, saveFile))  
    }
  } else if (!is.null(saveDir)){
    load(file.path(saveDir, loadFile))
  } else {
    stop("'saveDir' is missing and required to locate 'loadFile'.")
  }
  
  # if there are 3 treatment arms
  # } else if (max(data$tx)==3){
  #   logHRfit <- list(confidenceband.loghazardratio(data$fTime[data$tx<3], data$fInd[data$tx<3], data$tx[data$tx<3], t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nboot, Nv))
  #   logHRfit[[2]] <- confidenceband.loghazardratio(data$fTime[data$tx!=2], data$fInd[data$tx!=2], as.numeric(data$tx[data$tx!=2]==3) + 1, t1, t2, t1, t2, N, band11, band1, band12, band2, ngrid, biasadjust, tailsl, tailsu, nboot, Nv)
  # }
  
  timeGrid <- logHRfit[, 1]
  if (is.null(xMax)){ xMax <- max(timeGrid) }
  
  logHR <- logHRfit[timeGrid<=xMax, 2]
  ptLB <- logHRfit[timeGrid<=xMax, 3]
  ptUB <- logHRfit[timeGrid<=xMax, 4]
  smLB <- logHRfit[timeGrid<=xMax, 5]
  smUB <- logHRfit[timeGrid<=xMax, 6]
  timeGrid <- timeGrid[timeGrid<=xMax]
  
  yLim <- c(-1, 1.2)
  
  par(mar=c(5, 6.5, 3, 1), las=1, cex.axis=1.2)
  cexPlotTitle <- 1.4
  cexAxisLabel <- 1.3
  cexLegend <- 1.2
  
  plot(timeGrid, 1 - exp(logHR), xlab="", ylab="", type="n", xaxt="n", yaxt="n", xlim=c(0, xMax), ylim=yLim, bty="l")
  
  axis(1, at=seq(0, xMax, by=16 * 7), labels=seq(0, xMax, by=16 * 7) / 7)
  axis(1, at=108 * 7, labels=108)
  axis(2, at=seq(-1, 1, by=0.25), labels=seq(-100, 100, by=25))
  
  mtext("Weeks since Enrollment", side=1, line=2.8, cex=cexAxisLabel)
  mtext("Instantaneous Hazard-Ratio PE (%)", side=2, line=3.3, cex=cexAxisLabel, las=0)
  mtext(title, side=3, line=0.5, cex=cexPlotTitle, las=0)
  
  abline(h=0, col="gray50", lwd=2)
  lines(timeGrid, 1 - exp(logHR), lwd=3.5)
  lines(timeGrid, 1 - exp(ptLB), lty="dashed", lwd=3)
  lines(timeGrid, 1 - exp(ptUB), lty="dashed", lwd=3)
  #lines(timeGrid, 1 - exp(smLB), lty="dotted", lwd=3)
  #lines(timeGrid, 1 - exp(smUB), lty="dotted", lwd=3)
  
  if (showLegend){
    legend(0.2 * xMax, 1.2, legend=c("95% pointwise CI", "95% simultaneous CI"), cex=cexLegend, lwd=3, lty=c("dashed", "dotted"), bty="n", x.intersp=0.5, y.intersp=0.7)
  }
}


# 'data' is the original 'survival_infusions.csv' data frame
excludeLowerConfCases <- function(data, LBpInfTimeWithin4wkInfu=0.2, UBpInfTimeWithin4wkInfu=0.8){
  # exclude non-MITT participants and duplicate randomizations
  data <- subset(data, efficacy_flag==1)
  
  # at this point, all NAs in 'infufudays' are due to
  #   - Unable to estimate date of infection (see 'flag' variable)
  #   - Estimated date of infection (median) is prior to enrollment (see 'median_tsic_dee_rounded' and 'infudt' variables)
  
  # get the total number of cases in 'data'
  nAllCases <- length(unique(subset(data, infustatus==1 | is.na(infufudays))$ptid))
  
  # get the number of low-confidence cases to be excluded
  nLowConfCases <- length(unique(subset(data, is.na(infufudays) | (infustatus==1 & p_1st4_overall>=LBpInfTimeWithin4wkInfu & p_1st4_overall<=UBpInfTimeWithin4wkInfu))$ptid))
  
  # output data set with excluded infusion intervals with >0 probability of covering the true infection time (i.e., p_infect=1) among cases with 'p_1st4_overall' in 
  # the ['LBpInfTimeWithin4wkInfu', 'UBpInfTimeWithin4wkInfu'] interval
  data <- subset(data, (!is.na(infufudays)) | !(p_1st4_overall>=LBpInfTimeWithin4wkInfu & p_1st4_overall<=UBpInfTimeWithin4wkInfu & p_infect==1))
  
  return(list(data=data, nAllCases=nAllCases, nLowConfCases=nLowConfCases, propLowConfCases=nLowConfCases / nAllCases))
}
