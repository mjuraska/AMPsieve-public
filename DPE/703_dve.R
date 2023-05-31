# Purpose: Estimation of hazard ratio-based PE by AA position-specific features; hypothesis testing
# Method:  Lunn & McNeil 
#          Custom script allowing for stratification of baseline hazards, see dve3() below
# Input:  Master file (.csv) with time to event data, treatment group assignment, region (south america vs. other) and binary marks to be analyzed
# Output: .Rda containing results table with estimates of mark-specific PE, hazard rate ratio, and associated confidence intervals   
# Author:  Danica Shao
# Date:    May 13, 2022


rm(list=ls(all=TRUE))

figuresOutDir <- "/Volumes/trials/vaccine/p703/analysis/sieve/figures/DVE"
tablesOutDir <- "/Volumes/trials/vaccine/p703/analysis/sieve/tables/DVE"

options( width = 150, stringsAsFactors = FALSE );
setwd( "/Volumes/trials/vaccine/p703/analysis/sieve/code/DVE")
library( "survival" );

#### Modified DVE script to accomodate stratified baseline hazards
dve3 <- function (flrtime,flrstatus,flrtype,Vx,str=rep(0,length(flrtime))){
  # flrtime- right-censored failure time
  # flrstatus- indicator of observed infection (0=no, 1=yes)
  # flrtype- indicator of failure type (NA if not infected, 0 if infected and observed type 0, 1 if infected and observed type 1)
  # Vx- Vaccination status (0=placebo, 1=vaccine)
  
  # Standard Prentice et al. (1978) competing risks analysis
  
  # Analysis of type 0
  flrstatus0 <- flrstatus==1&flrtype==0
  fit0 <- coxph(Surv(flrtime,flrstatus0) ~ Vx + strata(str) )
  
  # Analysis of type 1
  flrstatus1 <- flrstatus==1&flrtype==1
  fit1 <- coxph(Surv(flrtime,flrstatus1) ~ Vx + strata(str))
  
  # Evaluate  H0: VE(type 0) = VE(type 1) via the Lunn and McNeil (1995, Biometrics) trick
  
  dblfutime <- c(flrtime,flrtime)
  dbldelta <- rep(c(0,1),each=length(flrtime))
  dblfustat <- c(flrstatus0,flrstatus1)
  dbltreatment <- c(Vx,Vx)+0
  dblstr <- c(str,str)
  fit <- coxph(Surv(dblfutime,dblfustat) ~ strata(dbldelta) + strata(dblstr) + dbltreatment + dbltreatment:dbldelta  ) 
  
  return( list( fit=fit, fit0=fit0, fit1=fit1 ) );
}

### Read master file with covariates and time to event data
master <-  read.csv( "../../adata/v703_survival_wk80_tau_sieve_v5_cam.csv" ); 

### Artificially censor CASES without neutralization data (variables gmt****)
master$hiv1event[ is.na( master$gmt50ms ) ] <- 0

### As per email from Michal 5/20/21 - we are analyzing the ordered categorical marks # PNGS in V5 and # Cys in gp120 as binary because they have very limited observed distributions. For the former we are splitting into 0:1 vs 2:3 and for the latter, 18 vs. 19:20
master$num.pngs.v5.greater.than.1 <- I( master$num.pngs.v5.mf > 1 )
master$num.cysteine.gp120.greater.than.18 <- I( master$num.cysteine.gp120.mf > 18 )

### List of features to analyze
features.to.analyze <- c( "hxb2.60.A.mf", "hxb2.170.Q.mf", "hxb2.230.D.mf", "hxb2.279.N.mf", "hxb2.280.N.mf", "hxb2.317.F.mf",
                          "hxb2.365.S.mf", "hxb2.429.E.mf", "hxb2.456.R.mf", "hxb2.458.G.mf", "hxb2.459.G.mf", "hxb2.471.G.mf",
                          "hxb2.156.pngs.mf", "hxb2.229.pngs.mf", "hxb2.234.pngs.mf", "hxb2.616.pngs.mf", "hxb2.824.pngs.mf", 
                          "num.pngs.v5.greater.than.1", "num.cysteine.gp120.greater.than.18" )

### Make some more readable feature names
feature.names <- sapply( strsplit(features.to.analyze,"\\."), function(.feature){
  paste( toupper( .feature[3]), "at", .feature[2])
}) 
feature.names[ 18:19 ] <- c("# PNGS in v5", "# Cys in gp120")


### List of columns in output 
results.colnames <- c( "feature", "feature.name", paste( ".n.events.type", rep(1:0,2), rep(c("trt","control"),each=2), sep = "."), "DVE.p.value", paste( "HR.Ratio", c( "mismatch/match", "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), paste( "VE.vs.Match", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), paste( "VE.vs.Mismatch", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ;


# We'll store the results here:
results <- data.frame( feature = features.to.analyze, feature.name = feature.names );
results[ , 3:( length( results.colnames) ) ] <- NA;
colnames( results ) <- results.colnames

# Begin for loop 
.results.row.i <- 0;
for( .feature in features.to.analyze ){
  
  .results.row.i <- .results.row.i + 1;
  
  # Be verbose:
  cat( "Computing results for feature", results[ .results.row.i, 1 ], "\n" );
  
  # Grab the column containing feature values from the master table
  .x <- master[ , .feature]
  
  ### Run DVE
  .dve.results <- dve3( master$hiv1fpday, master$hiv1event, .x, I( master$tx_pool == "T1+T2" ) )
  .DVE.p.value <- summary( .dve.results$fit )$coef[ 2, "Pr(>|z|)" ];
  
  # Extract hazard ratio 
  .HRR.results.is.mismatch.over.match <- summary( .dve.results$fit )$conf[ 2, 2 ] > 1;
  if( is.na( .HRR.results.is.mismatch.over.match ) ) {
    .HRR.results <- c( NA, NA, NA, NA );
  } else {
    .HRR.results <- summary( .dve.results$fit )$conf[ 2, ifelse( .HRR.results.is.mismatch.over.match, -1, -2 ) ];
    if( .HRR.results.is.mismatch.over.match ) {
      .HRR.results[ c( 2, 3 ) ] <- exp( -log( .HRR.results[ c( 3, 2 ) ] ) );
    }
    .HRR.results <- c( .HRR.results, summary( .dve.results$fit )$coefficients[ 2, "Pr(>|z|)" ] );
    names( .HRR.results )[ 4 ] <- "Pr(>|z|)";
  }
  # Extract VE values
  .type.1.VE <- ( 1 - summary( .dve.results$fit1 )$conf[ , -2 ] ) * 100;
  .type.1.VE[ c( 2, 3 ) ] <- .type.1.VE[ c( 3, 2 ) ];
  .type.1.VE <- c( .type.1.VE, summary( .dve.results$fit1 )$sc[ "pvalue" ] );
  .type.0.VE <- ( 1 - summary( .dve.results$fit0 )$conf[ , -2 ] ) * 100;
  .type.0.VE[ c( 2, 3 ) ] <- .type.0.VE[ c( 3, 2 ) ];
  .type.0.VE <- c( .type.0.VE, summary( .dve.results$fit0 )$sc[ "pvalue" ] );
  
  # Extract numbers of events
  .n.events.type.1.trt <- sum( .x == 1 & master$tx_pool == "T1+T2" & master$hiv1event, na.rm = T)
  .n.events.type.0.trt <- sum( .x == 0 & master$tx_pool == "T1+T2" & master$hiv1event, na.rm = T)
  .n.events.type.1.control <- sum( .x == 1 & master$tx_pool != "T1+T2" & master$hiv1event, na.rm = T)
  .n.events.type.0.control <- sum( .x == 0 & master$tx_pool != "T1+T2" & master$hiv1event, na.rm = T)
  
  # Insert into results table 
  results[ .results.row.i,  c( paste( ".n.events.type", rep(1:0,2), rep(c("trt","control"),each=2), sep = "."), "DVE.p.value", paste( "HR.Ratio", c( "mismatch/match", "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), paste( "VE.vs.Match", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." ), paste( "VE.vs.Mismatch", c( "estimate", "CI.low", "CI.high", "p.value" ), sep = "." )) ] <- c( .n.events.type.1.trt, .n.events.type.0.trt, .n.events.type.1.control, .n.events.type.0.control, .DVE.p.value, .HRR.results.is.mismatch.over.match, .HRR.results, .type.1.VE, .type.0.VE );
}

dve_results <- results

save( dve_results, file = file.path( tablesOutDir, "HVTN703_DVE_results_by_feature.Rda" )) 

### Forest plot
pdf( file = file.path( figuresOutDir, "HVTN703_DVE_results_by_feature_forest_plot.pdf" ), width = 6, height = 6)

display.large.numbers.as.infinity <- function( x, threshold = 10^8 ){
  replace <- abs(x) > threshold
  replace[ is.na( replace ) ] <- F
  x[ replace ] <- sign( x[ replace ] ) * Inf
  x
}

prettyPrint <- function( x, digits = 2, na.mask = NULL, na.replacement = "---", ... ){
  if( is.null( na.mask ) ) na.mask <- is.na(x) 
  x <- format( round( x, digits), nsmall = digits, ... ) 
  x[ na.mask ] <- na.replacement
  # Re-justify after replacing NAs 
  max.width <- max(nchar(x))
  for(i in 1:length(x)){
    if( nchar(x[i]) < max.width ) x[i] <- paste( c( rep(" ", max.width - nchar(x[i]) ), x[i]), collapse = "" )
  }
  x
}

### Set up the plotting area

# Define the plot region , print column headings and draw horizontal lines
par(xpd = 0, mar = c(0,0,0,0)) 
plot( 0, 0, xaxt = "n", yaxt= "n", xlab= "", ylab= '', bty = 'n', type='n', xlim = c(-0.7,11.2), ylim = c(-1.8,nrow(results)*1.5+0.8)) ## Empty plot
text( c(0,1,2,3,4,5.25,6.5,11), nrow(results)*1.5+0, c("Feature","Value","VRC01","Control","PE (%)","95% CI","P-value","Differential PE\nP-value"),cex=0.5,adj=0.5,pos=3) # Column header row 1
text( 2.5, nrow(results)*1.5+0.8, "Number of events (%)",cex=0.5,adj=0.5,pos=3)
abline( h = nrow(results)*1.5+0)
abline(h=-1) 

# Draw PE axis at bottom of figure
lines( c(7.5, 10.5), rep(-1.33, 2)) 
text( c(7.5, 8, 8.5, 9, 9.5, 10, 10.5 ), rep( -1.5, 7), rep( "|", 7), cex = 0.3 )
text( c(7.5, 8, 8.5, 9, 9.5, 10, 10.5 ), rep( -2, 7), c(-50,-25,0,25,50,75,100), cex = 0.5) 
text( 9,-2.8, "PE (%)", cex =0.5)

### Fill in the table values

# Feature names
text( 0.35, 1.5*nrow(results):1 - 1.5, results$feature.name, cex = 0.5, adj = 1 )

# Feature values (just T/F for binary marks, and the binary bins for the two categorical marks, see above)
text( 1, rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 , c( rep("T",nrow(results)-2), "","", rep("F",nrow(results)-2), "","" ), cex = 0.5, adj = 0.5 )

text( 1, 1.5*1+0.33,  bquote("">1), cex = 0.5, adj = 0.5 )
text( 1, 1.5*1-0.33,  bquote(""<=1), cex = 0.5, adj = 0.5 )
text( 1, 1.5*0+0.33,  bquote("">18), cex = 0.5, adj = 0.5 )
text( 1, 1.5*0-0.33,  bquote(""<=18), cex = 0.5, adj = 0.5 )

# No. of trt events (with percentage)
text( 2, rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 , 
      paste0( c( results$.n.events.type.1.trt, results$.n.events.type.0.trt), " (", prettyPrint( c( results$.n.events.type.1.trt, results$.n.events.type.0.trt) / sum(master$tx_pool == "T1+T2") * 100, 1 )  ,")"), cex = 0.5, adj = 0.5 )
# No. of control events (with percentage)
text( 3, rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 , 
      paste0( c( results$.n.events.type.1.control, results$.n.events.type.0.control), " (", prettyPrint( c( results$.n.events.type.1.control, results$.n.events.type.0.control) / sum(master$tx_pool != "T1+T2") * 100, 1 ) ,")"), cex = 0.5, adj = 0.5 )

# Protection efficacy (% ranging from -Inf to 100)
point.estimates <- c(  results$VE.vs.Match.estimate, results$VE.vs.Mismatch.estimate )
text( 4, rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 , prettyPrint( display.large.numbers.as.infinity(     point.estimates ), 1 ) , cex = 0.5, adj = 0.5 )

# 95% CI for protection efficacy
confidence.intervals <- paste0( "(", prettyPrint( c( results$VE.vs.Match.CI.low, results$VE.vs.Mismatch.CI.low), 1), ", ", prettyPrint( c( results$VE.vs.Match.CI.high, results$VE.vs.Mismatch.CI.high), 1)  , ")" )
confidence.intervals[ is.na( point.estimates) ] <- "-----"
text( 5.25, rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 , confidence.intervals , cex = 0.5, adj = 0.5 )

# P-value for testing protection efficacy =/= 0
text( 6.5, rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 ,  ifelse( is.na(point.estimates), "---", prettyPrint( c( results$VE.vs.Match.p.value, results$VE.vs.Mismatch.p.value ) , 3 )), cex = 0.5, adj = 0.5 )

# P-values testing for differential PE (PE_0 =/= PE_1)
text( 11, 1.5*nrow(results):1 - 1.5, prettyPrint( results$HR.Ratio.p.value, 3), cex = 0.5, adj = 0.5 )


### Draw the forest plot

## Vertical dashed line at PE = 0
lines( rep(8.5, 2), c(-1, nrow(results)*1.5+0), lty = 2, col = "gray")

# Point estimates
lower.bounds <- c(  results$VE.vs.Match.CI.low, results$VE.vs.Mismatch.CI.low )
upper.bounds <- c(  results$VE.vs.Match.CI.high, results$VE.vs.Mismatch.CI.high )

points( 8.5 + 0.02*point.estimates , rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5, col = rep( c("red","blue"), each = nrow(results)), cex = 0.5, pch = ifelse( point.estimates < -60, NA, 19 ) )

# Here we draw the line from the point estimate to the right confidence bound. Note that each of these will have an "arrow" (actually a small line perpendicular to the CI) denoting a closed interval, since we know all confidence intervals end on the right below PE = 100% and thus will be shown on the figure.  
arrows( 8.5 + 0.02 * pmax( -60, point.estimates ), rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 , 8.5 + 0.02 *  pmin( 100, upper.bounds ), rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5, col = rep( c("red","blue"), each = nrow(results)), code = 2, cex = 0.5, length = 0.02, angle = 90 )


# If the left confidence bound is less than minimum VE value (default -60), we draw a line with an arrow indicating the interval goes beyond what is notated in this figure. If not or if the bound does not exist (NA), don't draw anything (set lty to 0)
arrows( 8.5 + 0.02 * pmax( -60, lower.bounds ), rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 , 8.5 + 0.02 *  pmax( -60, upper.bounds ), rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5, col = rep( c("red","blue"), each = nrow(results)), code = 1, lty= ifelse( is.na(point.estimates), 0, ifelse (c(  results$VE.vs.Match.CI.low, results$VE.vs.Mismatch.CI.low ) < -60 , 1, 0 ) ) , length = 0.05, angle = 30 )

# If the left confidence bound is greater than than minimum VE value (default -60) we draw another perpendicular line indicating a closed interval. If not or if the bound does not exist (NA), don't draw anything (set lty to 0)
arrows( 8.5 + 0.02 * pmax( -60, lower.bounds ), rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5 , 8.5 + 0.02 *  pmax( -60, upper.bounds ), rep( 1.5*nrow(results):1, 2) + rep( c(0.33,-0.33), each = nrow(results) ) - 1.5, col = rep( c("red","blue"), each = nrow(results)), code = 1, lty= ifelse( is.na(point.estimates), 0, ifelse (c(  results$VE.vs.Match.CI.low, results$VE.vs.Mismatch.CI.low ) < -60 , 0, 1 ) ), length = 0.02, angle = 90 )

dev.off() 




