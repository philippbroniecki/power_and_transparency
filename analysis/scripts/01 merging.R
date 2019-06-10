# current version 22.05.2019 (rewrote as functions for enclosures)
# 1) merges DEU2 data sets, i.e. 'DEU2.csv' and 'deu15_27_9Nov12_proc_type.dta'
# from proc_type data only post04 can be kept because it is the only variable that does not vary over issues
# the proc_type data does not have an issue level identifier consequently proccde and type are dropped
# 2) merges the DEU with the legislative cycles data
# 3) merges DEU with Conflict of interest
# 4) merges DEU with type of file

###############################################################################
# 1) DEU
###############################################################################

deu.merge <- function(){
  
  # directory of the original DEU2 data set
  setwd(subdirs$org.DEUII)
  
  # load and order by prnr identifier
  deu <- read.csv("DEU2.csv")
  deu <- deu[order(deu$prnr),] 
  
  # load and order proctype data set
  proctype <- read.dta("deu15_27_9Nov12_proc_type.dta")
  proctype <- proctype[order(proctype$prnr),]
  
  # creating a unique ID
  deu$id <- seq(from = 1, to = length(deu[,1]), by = 1)
  proctype$id <- seq(from = 1, to = length(proctype[,1]), by = 1)
  proctype$prnr <- NULL
  
  # merging
  df <- merge(deu,proctype, by = "id")
  df$id <- NULL
  df$proccde <- NULL
  df$type <- NULL
  df <- df[order(df$isnrnmc),]
  
  return(df)
}

##############################################################################
# 2) Legislative Cycles
##############################################################################
leg.cycles.merge <- function(){
  
  # run function 1 (merge DEU2 and proctype)
  df <- deu.merge()
  
  # load data set from Reh et al
  setwd(subdirs$leg.cycles)
  reh <- read.csv("LCycles.csv")
  
  # merging (keep observations in df)
  df <- merge(df, reh, by = "cod", all.x = T)
  
  # generate a binary codecision variable (1 = codecision)
  df$codecision <- ifelse(grepl("COD", df$cod), yes = 1, no = 0)
  
  # convert words from factor to numeric
  if(is.factor(df$words)) df$words <- as.numeric(as.character(df$words))
  
  # early agreements 
  # ea = informal compromise 
  # ea.first = informal compromise and conclusion at first reading
  # ea.second = informal compromise and conclusion at second reading
  df$ea <- ifelse(df$early_agree == 1, yes = 1, no = 0)
  df$early_agree <- NULL
  df$ea.first <- ifelse(df$ea == 1 & df$stage == 1, yes = 1, no = 0)
  df$ea.second <- ifelse(df$ea == 1 & df$stage == 2, yes = 1, no = 0)
  df$stage <- NULL
  
  return(df)
}


###############################################################################
# 3) Conflict of interest
###############################################################################
# for all codecision files in the DEUII data
# I have collected whether the rapporteur is from a party
# that is in the government coalition and whether
# the country he is from held the presidency

conflict.of.interest.merge <- function(){
  
  # run previous functions
  df <- leg.cycles.merge()
  
  # location of conflict of interest data
  setwd(subdirs$same.country)
  
  # load conflict of interest/ 2nd principal data
  df.conflict <- read.csv("ConflictOfInterest.csv")
  df.conflict <- dplyr::select(df.conflict, isnrnmc, con_of_int)
  
  # merging with full data set
  df <- merge(df, df.conflict, by = "isnrnmc", all.x = TRUE)
  
  return(df)
}


###############################################################################
# 4) File type properties such as Directive, Committees asked for Opinion, Recitals
###############################################################################

file.type.merge <- function(){
  
  # run previous functions
  df <- conflict.of.interest.merge()
  
  # location of file type dat
  setwd(subdirs$type.data)
  
  # load type data
  df.type <- read.csv("typedata.csv")
  df.type$cod <- NULL
  
  # merge with full data
  df <- merge(df, df.type, by = "isnrnmc", all.x = TRUE)
  
  # no differentiation between (the 1) decision case and regulations
  df$directive[which(df$directive==3)] <- 0
  
  # save data
  setwd(subdirs$workdata)
  save(df, file = "01 step merged.RData")
  
}

# run everything (file.type.merge calls previous functions)
s.time <- Sys.time()
file.type.merge()
e.time <- Sys.time()
e.time - s.time
time.box$script1 <- e.time - s.time

# remove functions from global environment
rm(conflict.of.interest.merge, deu.merge, file.type.merge, leg.cycles.merge)