

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
-- 
-- README:  AMPsieve-public / mkMarkFiles
--          https://github.com/mjuraska/AMPsieve-public/tree/main/mkMarkFiles
--
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------


This subdirectory contains the R script ("mk_mark_files_v8.R") used to generate 
the analysis datasets (a.k.a. the marks files) for the AMP sieve analysis.

The source data files required to successfully run this script are available on 
this public repository for the AMP sieve dataset, on the Atlas data portal:

https://atlas.scharp.org/cpas/project/HVTN%20Public%20Data/HVTN%20704%20HPTN%20085%20and%20HVTN%20703%20HPTN%20081%20AMP/begin.view

Before running, the script needs to be customized to your local environment by
editing the variable "path.home" (in line 28).  This variable will need to be 
changed to reflect the path on your local system of the "mkMarkFiles" 
subdirectory as downloaded from the "AMPsieve-public" GitHub repository, which
should also contain the subdirectories of data sources from the Atlas portal.

Running the script will create the mark files in the "out" subdirectory.  Three
marks files are created:  one for each study and one for both studies pooled.


-- 

  C.A. Magaret
  cmagaret@fredhutch.org
  06 Sept 2023


