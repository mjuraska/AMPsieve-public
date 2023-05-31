# See the README in T:\vaccine\p704\analysis\sieve\adata\amp_sieve_marks_final_README.txt


library(SuperLearner)
library(drtmle)

# Marks of interest
MARKS = expand.grid( name = c("parscore1", "parscore2",
        "hxb2.60.A",  "hxb2.170.Q", "hxb2.230.D", "hxb2.279.N", "hxb2.280.N", 
        "hxb2.317.F", "hxb2.365.S", "hxb2.429.E", "hxb2.456.R", "hxb2.458.G", 
        "hxb2.459.G", "hxb2.471.G", "hxb2.156.pngs", "hxb2.229.pngs", "hxb2.234.pngs", 
        "hxb2.616.pngs", "hxb2.824.pngs",
        "length.gp120", "length.v1v2", "length.v5", 
        "num.pngs.gp120", "num.pngs.v1v2", "num.pngs.v5", 
        "num.cysteine.gp120",
        "hdist.zspace.sites.preselect.all", "hdist.zspace.sites.binding.all", 
        "epitope.dist.b", "epitope.dist.c", "epitope.dist.any"),
        seq = 'ls')
MARKS$mark = paste(MARKS$name, MARKS$seq, sep = '.')

# Read in the mark data file
dat = read.csv('/Volumes/trials/vaccine/p704/analysis/sieve/adata/amp_sieve_pooled_marks_final_v9.csv')
ss = subset(dat, hiv1event==1)
stopifnot(all(MARKS$mark %in% names(dat)))

# Dichotomize "num.pngs.v5.ls" as <=1 vs. >1 and "num.cysteine.gp120.ls" as <=18 vs. >18 
ss$num.pngs.v5.ls = ifelse(ss$num.pngs.v5.ls <= 1, 0, 1)
ss$num.cysteine.gp120.ls = ifelse(ss$num.cysteine.gp120.ls <= 18, 0, 1)

# Detect binary marks (should be 'hxb2' position marks plus the two dichotomized marks num.pngs.v5.ls and num.cysteine.gp120.ls)
MARKS$binary = sapply(MARKS$mark, function(n) all(ss[[n]] %in% c(0,1,NA)))

# Save DRTMLE fits for each mark by pooled and individual trial

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

# Save resutls into fits.drtmle
fits.drtmle = list()
for( trial in c('Pooled', 'HVTN 703', 'HVTN 704')) {
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
  
  for( dose in c('Pooled', 'Low', 'High') ) {
    
    # DRTMLE fits often result in error for the trial and dose specific analyses,
    # analyses are only run if either is dose or trials are Pooled
    if( dose != 'Pooled' && trial != 'Pooled') next
    
    fits.drtmle[[trial]][[dose]] = list()
    
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
    
    for( i in 1:nrow(MARKS)) {
      set.seed(1) 
      
      m = MARKS$mark[i]
      
      Y = ss[[m]]  
      keep.mark = !is.na(Y)
      
      keep = keep.trial & keep.dose & keep.mark
      Wkeep = W[keep,colsW,drop=FALSE]
      Akeep = A[keep]
      Ykeep = Y[keep]
      
      # for binary marks only run analysis if >=6 primary endpoints for each feature level
      if( MARKS$binary[i] & ((sum(Ykeep==0) < 6) | (sum(Ykeep==1) < 6)) ) {
        fit = 'not done'
      } else {
        # If DRTMLE fails to fit try until 'done'; set number of attempts in while loop
        # if while( done < 1 ) only one attempt is made
        done = 0
        fit = NULL
        while( done < 1 ) {
          if( MARKS$binary[i] ) {
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
          if( is(fit)=='drtmle' ) break
          
        }
      }      
      fits.drtmle[[trial]][[dose]][[m]] = fit

    }
  }
}

# save results
fp.rda = file.path('./data', 'drtmle_fits.rda')
save(fits.drtmle, file = fp.rda)

q(save="no")


