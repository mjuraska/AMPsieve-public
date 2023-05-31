# Running just the 230 png result that failed for given seed to see whether and how
# seed effects the results.

library(SuperLearner)
library(drtmle)

MARKS = expand.grid(name = c('hxb2.230.pngs'),
                    seq =c('ls'))
MARKS$mark = paste(MARKS$name, MARKS$seq, sep = '.')


# Read in the mark data file
dat = read.csv('/Volumes/trials/vaccine/p704/analysis/sieve/adata/amp_sieve_pooled_marks_final_v9.csv')
ss = subset(dat, hiv1event==1)
stopifnot(all(MARKS$mark %in% names(dat)))

# Save DRTMLE fits for each mark by pooled and individual trial
set.seed(1) 

A = ifelse(ss$tx_pool=='T1+T2', 1, 0) # It's okay to use tx_pooled for high (or low) dose comparisons below because A will  
                                      # be subset below in dose loop.
W = data.frame( 
  intercept = rep(1, nrow(ss)),
  protocol = ifelse(ss$protocol=="HVTN 703", 1, 0), 
  southAfrica = ss$southAfrica, 
  southAmerica = ss$southAmerica)
colsPooled = c('protocol', 'southAfrica', 'southAmerica')
cols703    = c('intercept', 'southAfrica')
cols704    = c('intercept', 'southAmerica')

fits.drtmle = list()
for( trial in c('HVTN 703')) {
  fits.drtmle[[trial]] = list()

  # Covariate matrix, W, columns depend on trial (Pooled, 703, or 704)  
  # Rows to 'keep' depend on trial
  if( trial=='Pooled' ) {
    colsW = colsPooled
    keep.trial = rep(TRUE, nrow(ss)) 
  } else if( trial=='HVTN 703') {
    colsW = cols703
    keep.trial = ss$protocol == trial
  } else if( trial=='HVTN 704') {
    colsW = cols704
    keep.trial = ss$protocol==trial
  } else {
    stop(sprintf('Unexpected trial: %s', trial))
  }
  
  for( dose in c('Pooled') ) {
    
    # DRTMLE fits often result in error for the trial and dose specific analyses,
    # analyses are only run if either is dose or trials are Pooled
    if( dose != 'Pooled' && trial != 'Pooled') next
    
    # rows to keep depend on dose 
    if( dose == 'Pooled' ) {
      keep.dose = rep(TRUE, nrow(ss))
    } else if( dose == 'Low' ) {
      keep.dose = ss$tx %in% c('C3', 'T1')
    } else if( dose == 'High' ) {
      keep.dose = ss$tx %in% c('C3', 'T2')
    } else {
      stop(sprintf('Unexpected dose: %s', dose))
    }
    
    for( m in MARKS$mark) {
      Y = ss[[m]]  
      keep.mark = !is.na(Y)
      
      keep = keep.trial & keep.dose & keep.mark
      Wkeep = W[keep,colsW,drop=FALSE]
      Akeep = A[keep]
      Ykeep = Y[keep]
     
      done = 0
      fit = NULL
      while( done < 10 ) {
        if( grepl('hxb2.230.pngs', m) ) {
          fit = try( { drtmle(W = Wkeep, A = Akeep, Y = Ykeep, 
                              stratify = FALSE, 
                              SL_g = c("SL.mean", "SL.glm", "SL.glmnet"),
                              SL_Q = c("SL.mean", "SL.glm", "SL.glmnet"), 
                              glm_Qr = "gn",
                              glm_gr = "Qn", 
                              cvFolds = 5, 
                              SL_method = "method.NNLS",
                              a_0 = c(1, 0),
                              family = binomial()) } )
          
        } else {
          fit = try( { drtmle(W = Wkeep, A = Akeep, Y = Ykeep, 
                              stratify = FALSE, 
                              SL_g = c("SL.mean", "SL.glm", "SL.glmnet"),
                              SL_Q = c("SL.mean", "SL.glm", "SL.glmnet"), 
                              glm_Qr = "gn",
                              glm_gr = "Qn", 
                              cvFolds = 5, 
                              SL_method = "method.NNLS",
                              a_0 = c(1, 0),
                              family = gaussian()) } )
        }
        
        done = done + 1
        print(is(fit))
        print(done)
        fits.drtmle[[done]] = fit
        
      }
    }
  }
}


#q(save="no")


