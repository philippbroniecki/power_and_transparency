# Imputation

imputation.fun <- function(){
  
  # location of step 2 data
  setwd(subdirs$workdata)
  load("02 step pcou added.RData")
  
  # indicator varibles (which rows have missings on pcom, pep, sep, rp, out)
  df$miss.ind <- 0
  df$miss.ind[which(is.na(df$pcom) | is.na(df$pep) | is.na(df$sep) |
                      is.na(df$rp) | is.na(df$out))] <- 1
  
  df$pcom.miss <- 0
  df$pcom.miss[which(is.na(df$pcom))] <- 1
  
  df$pep.miss <- 0
  df$pep.miss[which(is.na(df$pep))] <- 1
  
  df$sep.miss <- 0
  df$sep.miss[which(is.na(df$sep))] <- 1
  
  df$rp.miss <- 0
  df$rp.miss[which(is.na(df$rp))] <- 1
  
  df$out.miss <- 0
  df$out.miss[which(is.na(df$out))] <- 1
  
  # ctte is the responsible committee
  df$LIBE <- ifelse(df$ctte == 19, yes = 1, no = 0)
  df$TRAN <- ifelse(df$ctte == 13, yes = 1, no = 0)
  df$AFET <- ifelse(df$ctte == 25 | df$ctte == 1, yes = 1, no = 0)
  df$ITRE <- ifelse(df$ctte == 11, yes = 1, no = 0)
  df$JURI <- ifelse(df$ctte == 24 | df$ctte == 18, yes = 1, no = 0)
  df$ENVI <- ifelse(df$ctte == 10, yes = 1, no = 0)
  df$CULT <- ifelse(df$ctte == 17, yes = 1, no = 0)
  df$ECON <- ifelse(df$ctte == 8, yes = 1, no = 0)
  df$AGRI <- ifelse(df$ctte == 15, yes = 1, no = 0)
  df$REGI <- ifelse(df$ctte == 14, yes = 1, no = 0)
  df$EMPL <- ifelse(df$ctte == 9, yes = 1, no = 0)
  df$IMCO <- ifelse(df$ctte == 12, yes = 1, no = 0)
  df$ctte <- NULL
  
  # removing unnecessary variables
  df$id <- NULL
  df$fin <- NULL
  df$yearfin <- NULL
  df$pepppty <- NULL
  df$ppespty  <- NULL
  df$paldepty <- NULL
  df$pgreenpty <- NULL
  df$pguepty <- NULL
  df$puenpty <- NULL
  
  ##############################################
  #--------------Imputation  ------------------#
  ##############################################
  # m is number of imputations
  # p2s is print to screen (0 - nothing, 1 - little, 2 - thorough)
  # ts is time series identifier
  # cs is cross section identifier
  # ords is to keep the ordinal scale
  # noms for nominal vars
  # idvars are kept in the data set but not used for imputation
  
  # where the Banzhaf power index was used
  df.b <- amelia(df, m = N_imputations, p2s=2, cs="isnrnmc", incheck = TRUE,
                 idvars = c("prnrnmc", "prnr", "cod", "prname", "isnr",  "shu", "slv", "slt", 
                            "smt", "spl", "sro", "ssi", "ssk", "pbu", "pcy", "pee", "phu", "plv", 
                            "plt", "pmt", "ppl", "pro", "psi", "psk", "pcz", "miss.ind",
                            "pcom.miss", "pep.miss", "sep.miss", "rp.miss", "out.miss",
                            "dur_months", "ctte_op", "words", "recitals", "ea", "ea.first",
                            "ea.second", "con_of_int", "pcou_ss", #"pcou_p", #"pcou_avg",  
                            "start", "yearstart", "LIBE", "TRAN", "AFET", "ITRE", "JURI", "ENVI", 
                            "CULT", "ECON", "AGRI", "REGI", "EMPL", "IMCO", "directive","s_cou_w"))
  
  # where the Shapley-Shubik power index was used
  df.ss <- amelia(df, m = N_imputations, p2s=2, cs="isnrnmc", incheck = TRUE,
                  idvars = c("prnrnmc", "prnr", "cod", "prname", "isnr", "shu", "slv", "slt", 
                             "smt", "spl", "sro", "ssi", "ssk", "pbu", "pcy", "pee", "phu", "plv", 
                             "plt", "pmt", "ppl", "pro", "psi", "psk", "pcz", "miss.ind",
                             "pcom.miss", "pep.miss", "sep.miss", "rp.miss", "out.miss",
                             "dur_months", "ctte_op", "words", "recitals", "ea", "ea.first",
                             "ea.second", "con_of_int", "pcou_b", #"pcou_p", #"pcou_avg", 
                             "start", "yearstart", "LIBE", "TRAN", "AFET", "ITRE", "JURI","ENVI", 
                             "CULT", "ECON", "AGRI", "REGI", "EMPL", "IMCO", "directive","s_cou_w"))
  
  # I do not need member state positions anymore so I am dropping them
  for (e in 1: N_imputations){
    df.b$imputations[[e]] <- df.b$imputations[[e]] %>%
      dplyr::select(prnrnmc, isnrnmc, prnr, cod, prname, isnr, pcou_b, pcom, pcom, pep, rp, out, scou,
                    scom, sep, post04, dur_months, ctte_op, recitals, words, ea, ea.first, ea.second,
                    codecision, con_of_int, miss.ind, pcom.miss, pep.miss, sep.miss, rp.miss, out.miss, start, 
                    yearstart, LIBE, TRAN, AFET, ITRE, JURI, ENVI, CULT, ECON, AGRI, REGI, EMPL, IMCO, directive, s_cou_w) %>%
      rename(pcou = pcou_b)
    
    df.ss$imputations[[e]] <- df.ss$imputations[[e]] %>%
      dplyr::select(prnrnmc, isnrnmc, prnr, cod, prname, isnr, pcou_ss, pcom, pcom, pep, rp, out, scou,
                    scom, sep, post04, dur_months, ctte_op, recitals, words, ea, ea.first, ea.second,
                    codecision, con_of_int, miss.ind, pcom.miss, pep.miss, sep.miss, rp.miss, out.miss, start,
                    yearstart, LIBE, TRAN, AFET, ITRE, JURI, ENVI, CULT, ECON, AGRI, REGI, EMPL, IMCO, directive, s_cou_w) %>%
      rename(pcou = pcou_ss)
  }
  
  df.full <- list(banzhaf = df.b, shapshu = df.ss) #, pivot = df.p) #, simpleaavg = df.avg, )
  
  # Save the imputed data sets
  save(df.full, file = "03 step imputation.RData")
}

# run everything
imputation.fun()
rm(imputation.fun)
