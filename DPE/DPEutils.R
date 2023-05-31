#### Modified DVE script to accommodate stratified baseline hazards
DVE.f <- function (eventTime,eventInd,mark,tx,stratVar=rep(0,length(eventTime))){
  # eventTime- right-censored failure time
  # eventInd- indicator of observed infection (0=no, 1=yes)
  # mark- indicator of failure type (NA if not infected, 0 if infected and observed type 0, 1 if infected and observed type 1)
  # tx- Vaccination status (0=placebo, 1=vaccine)
  
  # Standard Prentice et al. (1978) competing risks analysis
  
  # Analysis of type 0

  eventInd0 <- eventInd==1&mark==0
  fit0 <- survival::coxph(Surv(eventTime,eventInd0) ~ tx + strata(stratVar) )
  
  # Analysis of type 1
  eventInd1 <- eventInd==1&mark==1
  fit1 <- survival::coxph(Surv(eventTime,eventInd1) ~ tx + strata(stratVar))
  
  # Evaluate  H0: VE(type 0) = VE(type 1) via the Lunn and McNeil (1995, Biometrics) trick
  # convert mark values for non-primary endpoints through tau to NA
  mark2 <-  ifelse(eventInd==0, NA, mark)
  
  fit = lunnMcneilTestS(eventTime, eventInd, mark2 + 1, tx, stratVar=stratVar)
  return( list( Diffp = fit$coef[3, 5], fit0=fit0, fit1=fit1 ) );
}

#Hazard-based VE forest plots for feature sets with adjusted p-values
PEtable.forestplot.withSieveT <- function(PEtable.x, xlim.x, ticks_at.x, outputDir, outputFileName){
  
  header <- c("Env AA\nPosition", "Seq\nFeature","No. Cases (VRC01 vs. P)\n(Incidence per PYR)", 
              "PE (%) (95% CI)","mean", "lower", "upper","P-value","Unadj. ",
              "FWER-adj.","FDR-adj.")
  colnames(PEtable.x) = header
  dt <- PEtable.x[,1:7]
  dt$` ` <- paste(rep(" ", 20), collapse = " ")#this determines how wide the forestplot is
  dt <- cbind(dt, PEtable.x[,8:11])
  
  #library(forestploter)
  tm <- forest_theme(base_size = 15,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = 15,
                     ci_lty = 1.5,
                     ci_lwd = 2.3,
                     ci_Theight = 0.3, # Set an T end at the end of CI 
                     # Reference line width/type/color
                     refline_lwd = 1,
                     refline_lty = "solid",
                     refline_col = "grey60",
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = "solid",
                     vertline_col = "grey60")
  p <- forest(dt[,c(1:4, 8:12)], #omit the columns of mean, lower, upper
              est = dt$mean,
              lower = dt$lower, 
              upper = dt$upper,
              sizes = 3.5,
              ci_cols = rep(c("black","royalblue","darkred"), length(dt$mean)/3),
              ci_column = 5,
              ref_line = 0,
              xlim = xlim.x,
              xlab = "PE (%) (95% CI)",
              ticks_at = ticks_at.x,
              theme = tm)
  p<- insert_text(p,
                  text = ,"2-sided Differential PE P-values",
                  col = 7:9,
                  part = "header",
                  just = "center",
                  gp = gpar(fontface = "bold", cex = 1.25))
 
  p <- add_underline(p,row = 0, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p,row = 1, col = c(7:9), part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = 2, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = length(dt$mean), part = "body", gp = gpar(lwd=2))
  
  
  ggplot2::ggsave(filename = outputFileName, 
                  plot = p, 
                  path = outputDir,
                  dpi = 320,
                  width = 14, height = 0.25*(length(dt$mean))+3.5,units = "in")  
}

#PE table for founder type analyses
PEtable.forestplot.withSieveT2 <- function(PEtable.x, xlim.x, ticks_at.x, outputDir, outputFileName){
  
  header <- c("Founder Multiplicity Type","No. Cases (VRC01 vs. P)\n(Incidence per PYR)", 
              "PE (%) (95% CI)","mean", "lower", "upper","P-value","2-sided Differential PE\nUnadjusted P Value ")
  colnames(PEtable.x) = header
  dt <- PEtable.x[,1:6]
  dt$` ` <- paste(rep(" ", 20), collapse = " ")#this determines how wide the forestplot is
  dt <- cbind(dt, PEtable.x[,7:8])
  
  #library(forestploter)
  tm <- forest_theme(base_size = 15,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = 15,
                     ci_lty = 1.5,
                     ci_lwd = 2.3,
                     ci_Theight = 0.3, # Set an T end at the end of CI 
                     # Reference line width/type/color
                     refline_lwd = 1,
                     refline_lty = "solid",
                     refline_col = "grey60",
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = "solid",
                     vertline_col = "grey60")
  p <- forest(dt[,c(1:3, 7:9)], #omit the columns of mean, lower, upper
              est = dt$mean,
              lower = dt$lower, 
              upper = dt$upper,
              sizes = 3.5,
              ci_cols = rep(c("black","royalblue","darkred"), length(dt$mean)/3),
              ci_column = 4,
              ref_line = 0,
              xlim = xlim.x,
              xlab = "PE (%) (95% CI)",
              ticks_at = ticks_at.x,
              theme = tm)
 
  p <- add_underline(p,row = 0, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = 1, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = length(dt$mean), part = "body", gp = gpar(lwd=2))
  p <- add_underline(p, row = 3, part = "body", gp = gpar(lwd=2))
  
  
  ggplot2::ggsave(filename = outputFileName, 
                  plot = p, 
                  path = outputDir,
                  dpi = 320,
                  width = 15, height = 0.25*(length(dt$mean))+3.5,units = "in")  
}



PEtable.forestplot.withSieveT3 <- function(PEtable.x, xlim.x, ticks_at.x, outputDir, outputFileName){
  
  header <- c("No. Lineages","No. Cases (VRC01 vs. P)\n(Incidence per PYR)", 
              "PE (%) (95% CI)","mean", "lower", "upper","P-value","Unadj. ",
              "FWER-adj.","FDR-adj.")
  colnames(PEtable.x) = header
  dt <- PEtable.x[,1:6]
  dt$` ` <- paste(rep(" ", 20), collapse = " ")#this determines how wide the forestplot is
  dt <- cbind(dt, PEtable.x[,7:10])
  
  #library(forestploter)
  tm <- forest_theme(base_size = 15,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = 15,
                     ci_lty = 1.5,
                     ci_lwd = 2.3,
                     ci_Theight = 0.3, # Set an T end at the end of CI 
                     # Reference line width/type/color
                     refline_lwd = 1,
                     refline_lty = "solid",
                     refline_col = "grey60",
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = "solid",
                     vertline_col = "grey60")
  p <- forest(dt[,c(1:3, 7:11)], #omit the columns of mean, lower, upper
              est = dt$mean,
              lower = dt$lower, 
              upper = dt$upper,
              sizes = 3.5,
              ci_cols = rep(c("black","royalblue","darkred"), length(dt$mean)/3),
              ci_column = 4,
              ref_line = 0,
              xlim = xlim.x,
              xlab = "PE (%) (95% CI)",
              ticks_at = ticks_at.x,
              theme = tm)
  p<- insert_text(p,
                  text = ,"2-sided Differential PE P-values",
                  col = 6:8,
                  part = "header",
                  just = "center",
                  gp = gpar(fontface = "bold", cex = 1.25))
  
  p <- add_underline(p,row = 0, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p,row = 1, col = c(6:8), part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = 2, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = length(dt$mean), part = "body", gp = gpar(lwd=2))
  p <- add_underline(p, row = 3, part = "body", gp = gpar(lwd=2))
  
  
  ggplot2::ggsave(filename = outputFileName, 
                  plot = p, 
                  path = outputDir,
                  dpi = 320,
                  width = 15.4, height = 0.25*(length(dt$mean))+3.5,units = "in")  
}


PEtable.forestplot.withSieveT4 <- function(PEtable.x,xlim.x, ticks_at.x, outputDir, 
                                           outputFileName, header, add.line = NULL){
 #browser()
  colnames(PEtable.x) = header
  num.col = dim(PEtable.x)[2]
  mean.col <- seq(1,num.col,1)[colnames(PEtable.x) == "mean"]
  lower.col <- seq(1,num.col,1)[colnames(PEtable.x) == "lower"]
  upper.col <- seq(1,num.col,1)[colnames(PEtable.x) == "upper"]
  dt <- PEtable.x[,1:(mean.col-1)]
  dt$` ` <- paste(rep(" ", 20), collapse = " ")#this determines how wide the forestplot is
  dt <- cbind(dt, PEtable.x[,mean.col:num.col])
  
  #library(forestploter)
  tm <- forest_theme(base_size = 15,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = 15,
                     ci_lty = 1.5,
                     ci_lwd = 2.3,
                     ci_Theight = 0.3, # Set an T end at the end of CI 
                     # Reference line width/type/color
                     refline_lwd = 1,
                     refline_lty = "solid",
                     refline_col = "grey60",
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = "solid",
                     vertline_col = "grey60")
  p <- forest(dt[,c(1:(mean.col), (upper.col+2):(num.col+1))], #omit the columns of mean, lower, upper
              est = dt$mean,
              lower = dt$lower, 
              upper = dt$upper,
              sizes = 3.5,
              ci_cols = rep(c("black","royalblue","darkred"), length(dt$mean)/3),
              ci_column = mean.col,
              xlog = FALSE,
              ref_line = 0,
              xlim = xlim.x,
              xlab = "PE (%) (95% CI)",
              ticks_at = ticks_at.x,
              theme = tm)
  
  p <- add_underline(p,row = 0, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = 1, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = length(dt$mean), part = "body", gp = gpar(lwd=2))
  p <- add_underline(p, row = add.line, part = "body", gp = gpar(lwd=2))
  
  p<- insert_text(p,
                  text = "HVTN 704/HPTN 085",
                  col = 1:2,
                  row = 1,
                  part = "body",
                  just = "left",
                  gp = gpar(fontface = "plain", cex = 1.25))
  p<- insert_text(p,
                  text = "HVTN 703/HPTN 081",
                  col = 1:2,
                  row = 11,
                  part = "body",
                  just = "left",
                  gp = gpar(fontface = "plain", cex = 1.25))
  
  ggplot2::ggsave(filename = outputFileName, 
                  plot = p, 
                  path = outputDir,
                  dpi = 320,
                  width = 15, height = 0.25*(length(dt$mean))+3.5,units = "in")  
}

PEtable.forestplot.withSieveT5 <- function(PEtable.x, xlim.x, ticks_at.x, outputDir, 
                                           outputFileName, header, add.line = NULL){
  
  colnames(PEtable.x) = header
  num.col = dim(PEtable.x)[2]
  mean.col <- seq(1,num.col,1)[colnames(PEtable.x) == "mean"]
  lower.col <- seq(1,num.col,1)[colnames(PEtable.x) == "lower"]
  upper.col <- seq(1,num.col,1)[colnames(PEtable.x) == "upper"]
  dt <- PEtable.x[,1:(mean.col-1)]
  dt$` ` <- paste(rep(" ", 20), collapse = " ")#this determines how wide the forestplot is
  dt <- cbind(dt, PEtable.x[,mean.col:num.col])
  
  #library(forestploter)
  tm <- forest_theme(base_size = 15,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = 15,
                     ci_lty = 1.5,
                     ci_lwd = 2.3,
                     ci_Theight = 0.3, # Set an T end at the end of CI 
                     # Reference line width/type/color
                     refline_lwd = 1,
                     refline_lty = "solid",
                     refline_col = "grey60",
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = "solid",
                     vertline_col = "grey60")
 
  p <- forest(dt[,c(1:(mean.col), (upper.col+2):(num.col+1))], #omit the columns of mean, lower, upper
              est = dt$mean,
              lower = dt$lower, 
              upper = dt$upper,
              sizes = 3.5,
              ci_cols = rep(c("black","royalblue","darkred"), length(dt$mean)/3),
              ci_column = mean.col,
              xlog = FALSE,
              ref_line = 0,
              xlim = xlim.x,
              xlab = "PE (%) (95% CI)",
              ticks_at = ticks_at.x,
              theme = tm)
  
  p <- add_underline(p,row = 0, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = 1, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = length(dt$mean), part = "body", gp = gpar(lwd=2))
  
  ggplot2::ggsave(filename = outputFileName, 
                  plot = p, 
                  path = outputDir,
                  dpi = 320,
                  width = 15, height = 0.25*(length(dt$mean))+3.5,units = "in")  
}


forestplot.logisticReg<- function(result, xlog = FALSE, ref_line = 0, ticks_digits, xlim.x, ticks_at.x, ci_cols,
                                  outputDir, outputFileName, header){
  #browser()
  colnames(result) = header
  num.col = dim(result)[2]
  mean.col <- seq(1,num.col,1)[colnames(result) == "mean"]
  lower.col <- seq(1,num.col,1)[colnames(result) == "lower"]
  upper.col <- seq(1,num.col,1)[colnames(result) == "upper"]
  dt <- result[,1:(mean.col-1)]
  dt$` ` <- paste(rep(" ", 20), collapse = " ")#this determines how wide the forestplot is
  dt <- cbind(dt, result[,mean.col:num.col])
  
  #library(forestploter)
  tm <- forest_theme(base_size = 15,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = 15,
                     ci_lty = 1.5,
                     ci_lwd = 2.3,
                     ci_Theight = 0.3, # Set an T end at the end of CI 
                     # Reference line width/type/color
                     refline_lwd = 1,
                     refline_lty = "solid",
                     refline_col = "grey60",
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = "solid",
                     vertline_col = "grey60")
  p <- forest(dt[,c(1:(mean.col), (upper.col+2):(num.col+1))], #omit the columns of mean, lower, upper
              est = dt$mean,
              lower = dt$lower, 
              upper = dt$upper,
              sizes = 3.5,
              ci_cols = ci_cols,
              ci_column = mean.col,
              xlog = xlog,
              ref_line = ref_line,
              ticks_digits = ticks_digits,
              xlim = xlim.x,
              xlab = "Odds Ratio (95% CI)",
              ticks_at = ticks_at.x,
              theme = tm)
  
  p <- add_underline(p,row = 0, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = 1, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = length(dt$mean), part = "body", gp = gpar(lwd=2))
  
  
  ggplot2::ggsave(filename = outputFileName, 
                  plot = p, 
                  path = outputDir,
                  dpi = 320,
                  width = 15, height = 0.25*(length(dt$mean))+3.5,units = "in")  
}

PEtable.forestplot <- function(PEtable.x, xlim.x, ticks_at.x, outputDir, outputFileName,
                               header = c("Env AA\nPosition", "Seq\nFeature","No. Cases (VRC01 vs. P)\n(Incidence per PYR)",
                                           "PE (%) (95% CI)","mean", "lower", "upper","P-value"),
                               ncolors = 2,addLine=NULL){
  
  
  colnames(PEtable.x) = header
  dt <- PEtable.x[,colnames(PEtable.x) != "P-value"]
  ncol1 <- dim(dt)[2]
  dt$` ` <- paste(rep(" ", 20), collapse = " ")#this determines how wide the forestplot is
  dt <- cbind(dt, PEtable.x[,"P-value"])
  ncol2 <- dim(dt)[2]
  #library(forestploter)
  tm <- forest_theme(base_size = 10,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = 15,
                     ci_lty = 1,
                     ci_lwd = 1.7,
                     ci_Theight = 0.3, # Set an T end at the end of CI 
                     # Reference line width/type/color
                     refline_lwd = 1,
                     refline_lty = "solid",
                     refline_col = "grey60",
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = "solid",
                     vertline_col = "grey60")
  if(ncolors == 2){ci_cols.x = rep(c("royalblue","darkred"), length(dt$mean)/2)}
  else{ci_cols.x = "royalblue"}
  p <- forest(dt[,c(1:(ncol1-3), (ncol2-1):(ncol2))], #omit the columns of mean, lower, upper
              est = dt$mean,
              lower = dt$lower, 
              upper = dt$upper,
              sizes = 2.7,
              ci_cols = ci_cols.x,
              ci_column = ncol1-2,
              ref_line = 0,
              xlim = xlim.x,
              xlab = "PE (%) (95% CI)",
              ticks_at = ticks_at.x,
              theme = tm)
  p <- add_underline(p,row = 0, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = 1, part = "header", gp = gpar(lwd=2))
  p <- add_underline(p, row = length(dt$mean), part = "body", gp = gpar(lwd=2))
  if(!is_null(addLine)){
    p <- add_underline(p, row = addLine, part = "body", gp = gpar(lwd=2))
  }
  
  
  ggplot2::ggsave(filename = outputFileName, 
                  plot = p, 
                  path = outputDir,
                  dpi = 320,
                  width = 9, height = 0.23*(length(dt$mean)+ ifelse(length(dt$mean)<=20, 4, 3)),units = "in")  
}
format.p <- function(p, ndigits=2){
  pp <- NULL
  for(i in 1:length(p)){
    if(is.na(p[i])){
      pp[i] <- "--"
    }else{
      if(p[i]<0.001){pp[i] <- " < 0.001"}
      else if (p[i]==1){pp[i] <- "1"}
      else{pp[i] <-as.character(format(as.numeric(p[i]),digits=ndigits,nsmall=ndigits))
      if(pp[i] == "1.00"){pp[i]="1"}
      } 
    }
  }
  return (pp)
}

format.p2 <- function(p, ndigits=2){
  pp <- NULL
  for(i in 1:length(p)){
    if(is.na(p[i])){
      pp[i] <- "--"
    }else{
      if(p[i]<0.001){pp[i] <- " < 0.001"}
      else if (p[i]==1){pp[i] <- "= 1"}
      else{pp[i] <-paste0(" = ",as.character(format(as.numeric(p[i]),digits=ndigits,nsmall=ndigits)))
           if(pp[i]== " = 1.00") {pp[i]= " = 1"}} 
    }
  }
  return (pp)
}
format.PE <- function(PE){
  if(is.na(PE)){
    return ("--")
  }else if (PE == "Inf"){
    return (Inf)
  }else if (PE == "-Inf"|as.numeric(PE)<(-1000)){
    return (-Inf)
  }else{
    return(as.character(format(round(as.numeric(PE),1),nsmall=1)))
  }
 }

table.seqFeatLabel.tier1 <- function(feature){
  if(feature == "num.pngs.v5"){
    Haplotype1 = "> 1"
    Haplotype0 = "<= 1"
  }else if (feature == "num.cysteine.gp120"){
    Haplotype1 = "> 18"
    Haplotype0 = "<= 18"
  }else if (feature %in% c("hxb2.156.pngs", "hxb2.229.pngs", "hxb2.234.pngs", "hxb2.616.pngs", "hxb2.824.pngs","hxb2.230.pngs")){
    Haplotype1 = "PNGS"
    Haplotype0 = "Not PNGS"
  }else{
    Haplotype1 = strsplit(feature,"\\.")[[1]][[3]]
    Haplotype0 = paste("Not", Haplotype1)
  }
  return(list(Haplotype1 = Haplotype1, Haplotype0 = Haplotype0))
}

table.seqFeatLabel.tier2 <- function(feature){
  Haplotype1 = strsplit(feature,"\\.")[[1]][[4]]
  Haplotype0 = paste("Not", Haplotype1)
  return(list(Haplotype1 = Haplotype1, Haplotype0 = Haplotype0))
}

#obtain the dataset for a specific trial and dosage
trial_dose_data <- function(master, trial, dose){
  if(trial == "704and703"){
    master.trial <- master
  }else{
    master.trial <- filter(master, protocol == paste0("HVTN ", trial) )
  }
  if(dose %in% c("T1", "T2")){
    master.trial <- filter(master.trial, tx %in% c("C3", dose))#dose-specific analysis exclude data from the other dosage
    master.trial$txLabel <- 1*(master.trial$tx == dose)
  }else{
    master.trial$txLabel <- 1*(master.trial$tx_pool == "T1+T2")
  } 
  master.trial <- filter(master.trial, !is.na(hiv1fpday))  # two 703 ppts with missing sequences have also a missing time-to-event
  
  master.trial <- mutate(master.trial, stratVar=case_when(protocol=="HVTN 704" & southAmerica==1 ~ "704SAm",  
                                                          protocol=="HVTN 704" & southAmerica==0 ~ "704notSAm",
                                                          protocol=="HVTN 703" & southAfrica==1 ~ "703SAf",
                                                          protocol=="HVTN 703" & southAfrica==0 ~ "703notSAf"))
  
  
  return(master.trial)
}

#summarize fits that has fit1, fit0 and a differential VE test p-value
result.summary <- function(fit.list, dose, .x){
  if(dose == "T1+T2"){
    .DPE.p.value <- fit.list$Diffp #Differential VE test is only done for dose pooled vs. Placebo
  }else{
    .DPE.p.value <- NA
  }
  # Extract VE values
  .type.1.PE <- ( 1 - summary(fit.list$fit1 )$conf[ , -2 ] ) * 100;
  .type.1.PE[ c( 2, 3 ) ] <- .type.1.PE[ c( 3, 2 ) ];
  .type.1.PE <- c( .type.1.PE, summary(fit.list$fit1 )$sc[ "pvalue" ] );
  .type.0.PE <- ( 1 - summary(fit.list$fit0 )$conf[ , -2 ] ) * 100;
  .type.0.PE[ c( 2, 3 ) ] <- .type.0.PE[ c( 3, 2 ) ];
  .type.0.PE <- c( .type.0.PE, summary(fit.list$fit0 )$sc[ "pvalue" ] );
  
  # Extract numbers of events
  .n.events.type.1.trt <- sum( .x == 1 & master.trial$txLabel == 1 & master.trial$hiv1event==1, na.rm = TRUE)
  .n.events.type.0.trt <- sum( .x == 0 & master.trial$txLabel == 1 & master.trial$hiv1event==1, na.rm = TRUE)
  .n.events.type.1.placebo <- sum( .x == 1 & master.trial$txLabel == 0 & master.trial$hiv1event==1, na.rm = TRUE)
  .n.events.type.0.placebo <- sum( .x == 0 & master.trial$txLabel == 0 & master.trial$hiv1event==1, na.rm = TRUE)
  
  if((.n.events.type.1.trt ==0 & .n.events.type.1.placebo == 0)){.type.1.PE[4] <- NA}
  if((.n.events.type.0.trt ==0 & .n.events.type.0.placebo == 0)){.type.0.PE[4] <- NA}
  
  if((.n.events.type.1.trt ==0 & .n.events.type.1.placebo == 0) |(.n.events.type.0.trt ==0 & .n.events.type.0.placebo == 0)){
    .DPE.p.value <- NA
  }
  # Obtain incidence rates (number of cases per person-years)
  .inc.1.trt <- .n.events.type.1.trt/(sum(master.trial$hiv1fpday[.x == 1 & master.trial$txLabel == 1], na.rm = TRUE)/365.25)
  .inc.0.trt <- .n.events.type.0.trt/(sum(master.trial$hiv1fpday[.x == 0 & master.trial$txLabel == 1], na.rm = TRUE)/365.25)
  .inc.1.placebo <- .n.events.type.1.placebo/(sum(master.trial$hiv1fpday[.x == 1 & master.trial$txLabel == 0], na.rm = TRUE)/365.25)
  .inc.0.placebo <- .n.events.type.0.placebo/(sum(master.trial$hiv1fpday[.x == 0 & master.trial$txLabel == 0], na.rm = TRUE)/365.25) 
  #paste the number of events and incidence rates together
  paste.n.inc <- function(n, inc){paste0(n, " (", ifelse(is.na(inc), "--", format(inc, nsmall=2, digits=2)), ")")}
  .n.events.type.1.trt <- paste.n.inc (.n.events.type.1.trt , .inc.1.trt)
  .n.events.type.0.trt <- paste.n.inc (.n.events.type.0.trt , .inc.0.trt)
  .n.events.type.1.placebo <- paste.n.inc (.n.events.type.1.placebo , .inc.1.placebo)
  .n.events.type.0.placebo <- paste.n.inc (.n.events.type.0.placebo , .inc.0.placebo)
  return(c( .n.events.type.1.trt, 
            .n.events.type.0.trt, .n.events.type.1.placebo, .n.events.type.0.placebo, .DPE.p.value,.type.1.PE, .type.0.PE ))
}
