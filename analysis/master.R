# current version 09.06.2019
rm(list = ls() )

# libraries
library(foreign)
library(Amelia)
library(texreg)
library(psych)
library(frm)
library(Matching)
library(MASS)
library(gbm)
library(dplyr)
library(xtable)

# setup
#_____________________________________________________________________________________________________________________
# change this to your directory (containing master.R) using forward slashes instead of Windows back slashes
# e.g.: main <- "C:/Power and Transparency in Political Negotiations"
main <- "C:/analysis"
#_____________________________________________________________________________________________________________________

#_____________________________________________________________________________________________________________________
# The DEU data will be imputed for EP & Commission positions and salience, outcomes, reference points
# Choose the number of imputation data sets that you want to estimate
# default for replication is: N_imputations <- 1000
N_imputations <- 1000
#_____________________________________________________________________________________________________________________

#_____________________________________________________________________________________________________________________
# Latex fonts ON or OFF (turning this ON might not work on Mac)
# default is: latex <- 0

# 1 = prints tables and plots using computer modern fonts
# 0 = prints everything in standard fonts
latex <- 0

# if latex fonts are turned on, the following must be installed:
# 1) Rtools
# 2) R packages: rjson, fontcm, extrafont
# 3) Ghostscript

# if latex = 1, ghostscript location MUST be specified (otherwise ignore)
# e.g.: ghostscript.file.path <- "C:/Program Files/gs/gs9.27/bin/gswin64c.exe"
if (latex) ghostscript.file.path <- "C:/Program Files/gs/gs9.27/bin/gswin64c.exe"
#_____________________________________________________________________________________________________________________

# sub-directories (do not change)
subdirs <- list(scripts = paste(main, "scripts", sep = "/"),
                org.DEUII = paste(main, "original data sets", "DEU2", sep = "/"),
                leg.cycles = paste(main, "original data sets", "Legislative Cycles", sep = "/"),
                council.voting = paste(main, "original data sets", "Council Voting Weights", sep = "/"),
                corrected.EAs = paste(main, "original data sets", "EAsCorrected", sep = "/"),
                same.country = paste(main, "original data sets", "ConflictOfInterest", sep = "/"),
                type.data = paste(main, "original data sets", "type of file", sep = "/"),
                ea.data = paste(main, "original data sets", "Early Agreements", sep = "/"),
                workdata = paste(main, "work data", sep = "/"),
                graphs = paste(main, "figures", sep = "/"),
                tables = paste(main, "tables", sep = "/"))

#_____________________________________________________________________________________________________________________
# Sequence 
#_____________________________________________________________________________________________________________________
set.seed(123)

# record the runtime for each script
time.box <- list(script1 = NA, script2 = NA, script3 = NA, script4 = NA, script5 = NA, script6 = NA, script7 = NA, 
                 script8 = NA, script9 = NA, script10 = NA, script11 = NA, script12 = NA)

# Step 1: merging DEU2 data with Proctype (provided by Robert Thomson upon request)
# last run on 10.06.2019
setwd(subdirs$scripts)
source("01 merging.R", echo = TRUE) 

# Step 2: estimating the position of the Council and adding it to the data set
# last run on 10.06.2019
setwd(subdirs$scripts)
source("02 council position.R", echo = TRUE)

# Step 3: imputation
# last run on 10.06.2019
setwd(subdirs$scripts)
source("03 imputation.R", echo = TRUE)

# Step 4: coding main DV and IVs
# last run on 10.06.2019
setwd(subdirs$scripts)
source("04 variable coding.R", echo = TRUE) 

# Step 5: matching, balance statistics, main model
# last run on 10.06.2019
setwd(subdirs$scripts)
source("05 main model.R", echo = TRUE)

# Step 6: matching, balance and model with non-parametric propensity score
# last run on 10.06.2019
setwd(subdirs$scripts)
source("06 main model with non parametric propensity score.R", echo = TRUE)

# Step 7: treatment only (matching but without controlling for propensity score)
# last run on 10.06.2019
setwd(subdirs$scripts)
source("07 treatment only.R", echo = TRUE)

# Step 8: placebo test
# last run on 10.06.2019
setwd(subdirs$scripts)
source("08 placebo test.R", echo = TRUE)

# Step 9: multiverse (Note: commented out latex fonts for replication)
# last run on 10.06.2019
setwd(subdirs$scripts)
source("09 multiverse.R", echo = TRUE)

# Step 10: Descriptes
# last run on 10.06.2019
setwd(subdirs$scripts)
source("10 descriptives.R", echo = TRUE)

# Step 11: appendix model (without matching)
# last run on 10.06.2019
setwd(subdirs$scripts)
source("11 appendix model.R", echo = TRUE)

# Step 12: appendix comparative statics
# last run on 10.06.2019
setwd(subdirs$scripts)
if(latex){
  source("12 comparative statics latex.R", echo = FALSE)  
} else{
  source("12 comparative statics.R", echo = FALSE)
}

# Step 13: plots
# last run on 10.06.2019
setwd(subdirs$scripts)
if(latex) {
  source("13 plots latex.R", echo = FALSE)
} else{
  source("13 plots standard fonts.R", echo = FALSE)
}

# Step 14: weighted salience model with parametric propensity (appendix)
# last run 10.06.2019
setwd(paste(subdirs$scripts,"/weighted Council salience models", sep=""))
source("01 matching and balance relative salience weighted.R", echo = TRUE)

# Step 15: weighted salience model with non-parametric propensity (appendix)
# last run 10.06.2019
setwd(paste(subdirs$scripts,"/weighted Council salience models", sep=""))
source("02 gradient boosting relative salience weighted.R", echo = TRUE)

# Step 16: weighted salience model with treatment only (appendix)
# last run 10.06.2019
setwd(paste(subdirs$scripts,"/weighted Council salience models", sep=""))
source("03 models weighted salience.R", echo = TRUE)
