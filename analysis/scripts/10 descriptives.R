# SQ in between institutions function   
no.change.func <- function(pep, rp, pcou){
  no.change <- as.numeric((pep < rp & pcou > rp) | (pcou < rp & pep > rp)) 
  return(no.change)
}

# midpoint between both institutions
mid.func <- function(pep, pcou){
  mid <- (pep + pcou) / 2
  return(mid)
}

descriptives <- function(){

  # location of matched data
  setwd(subdirs$workdata)
  load("05 step analysis data with relative salience.RData")
    
  # counting frequencies of scenarios
  for (e in 1:N_imputations){
    
    df <- df.matched[[1]]$imputations[[e]]
    
    # droping observations outside 1999-2009 period
    
    # is the status quo in between both institutions
    nochange <- no.change.func(pep = df$pep, rp = df$rp, pcou = df$pcou) 
    # indifference of actor closest to sq
    indiff <- rep(NA, length(df$pep))
    # midpoint between both institutions
    midpoint <- mid.func(pep = df$pep, pcou = df$pcou)
    
    ###### Policy constellations #####################################
    ## In scenarios 1 the SQ in not in beween both institutions and 
    ## the EP is closer to the SQ than the Council
    # Scenario 1a: --SQ---EP-------Cou----
    scenario1a <- nochange == 0 & abs(df$pep - df$rp) < abs(df$pcou - df$rp) & df$pep < df$pcou
    indiff <- ifelse(scenario1a == TRUE, yes = df$pep + abs(df$rp - df$pep), no = indiff)
    s1a.constrained <- scenario1a == TRUE & indiff < midpoint
    s1a.free <- scenario1a == TRUE & indiff > midpoint
    if(sum(scenario1a) != sum(s1a.constrained) + sum(s1a.free)) stop("check scenarios")
    
    # Scenario 1b: --Cou---------EP----SQ-
    scenario1b <- nochange == 0 & abs(df$pep - df$rp) < abs(df$pcou - df$rp) & df$pep > df$pcou  
    indiff <- ifelse(scenario1b == TRUE, yes = df$pep - abs(df$rp - df$pep), no = indiff)
    s1b.constrained <- scenario1b == TRUE & indiff > midpoint
    s1b.free <- scenario1b == TRUE & indiff < midpoint
    
    ## In scenarios 2 the SQ in not in beween both institutions and 
    ## the Council is closer to the SQ than the EP
    # Scenario 2a: --SQ---Cou--------EP---
    scenario2a <- nochange == 0 & abs(df$pcou - df$rp) < abs(df$pep - df$rp) & df$pcou < df$pep
    indiff <- ifelse(scenario2a == TRUE, yes = df$pcou + abs(df$rp - df$pcou), no = indiff)
    s2a.constrained <- scenario2a == TRUE & indiff < midpoint
    s2a.free <- scenario2a == TRUE & indiff > midpoint
    
    # Scenario 2b --EP-----------Cou---SQ
    scenario2b <- nochange == 0 & abs(df$pcou - df$rp) < abs(df$pep - df$rp) & df$pcou > df$pep
    indiff <- ifelse(scenario2b == TRUE, yes = df$pcou - abs(df$rp - df$pcou), no = indiff)
    s2b.constrained <- scenario2b == TRUE & indiff > midpoint
    s2b.free <- scenario2b == TRUE & indiff < midpoint
    
    ## Counting
    # Council closer to status quo
    # unconstrained
    ifelse(e == 1, yes = s1 <- sum(table(s2a.free)[2], table(s2b.free)[2], na.rm = TRUE), 
           no = s1 <- s1 + sum(table(s2a.free)[2], table(s2b.free)[2], na.rm = TRUE))
    # constrained
    ifelse(e == 1, yes = s2 <- sum(table(s2a.constrained)[2], table(s2b.constrained)[2], na.rm = TRUE), 
           no = s2 <- s2 + sum(table(s2a.constrained)[2], table(s2b.constrained)[2], na.rm = TRUE))
    
    # EP closer to status quo
    # unconstrained
    ifelse(e == 1, yes = s3 <- sum(table(s1a.free)[2], table(s1b.free)[2], na.rm = TRUE), 
           no = s3 <- s3 + sum(table(s1a.free)[2], table(s1b.free)[2], na.rm = TRUE))
    # constrained
    ifelse(e == 1, yes = s4 <- sum(table(s1a.constrained)[2], table(s1b.constrained)[2], na.rm = TRUE), 
           no = s4 <- s4 + sum(table(s1a.constrained)[2], table(s1b.constrained)[2], na.rm = TRUE))
    
    # no change environment
    ifelse(e == 1, s5 <- table(nochange)[2], s5 <- s5 + table(nochange)[2])
    
    # check for errors 
    t <- NA
    t[1] <- sum(table(s2a.free)[2], table(s2b.free)[2], na.rm = TRUE)
    t[2] <- sum(table(s2a.constrained)[2], table(s2b.constrained)[2], na.rm = TRUE)
    t[3] <- sum(table(s1a.free)[2], table(s1b.free)[2], na.rm = TRUE)
    t[4] <- sum(table(s1a.constrained)[2], table(s1b.constrained)[2], na.rm = TRUE)
    t[5] <- table(nochange)[2]
    
    if(sum(t, na.rm = TRUE) != 92) stop("observation(s) unaccounted for")
    
    # influence ratio # standard
    ifelse(e == 1, mean.ep.inf <- mean(df$ep.inf), 
           mean.ep.inf <- mean.ep.inf + mean(df$ep.inf))
    
    # influence by first reading EAs
    tmp2 <- df %>%
      group_by(ea.first) %>%
      summarise(infl = mean(ep.inf))
    ifelse(e == 1, infl.by.ea.first <- tmp2, infl.by.ea.first <- infl.by.ea.first + tmp2)
    
  } # end of loop over all imputations
  
  ## averaging 
  s1 / N_imputations
  s2 / N_imputations
  s3 / N_imputations
  s4 / N_imputations
  s5 / N_imputations
  if( sum(s1, s2, s3, s4, s5) / 1000 != 92) stop("counting error")
  
  # Table 2: Frequency of Policy Environments in the Data
  table2 <- matrix(c(s1, s2, sum(s1,s2, na.rm = TRUE),
                     s3, s4, sum(s3,s4, na.rm = TRUE),
                     NA, s5, s5,
                     sum(s1, s3, na.rm = TRUE),
                     sum(s2, s4, s5, na.rm = TRUE),
                     sum(sum(s1,s2, na.rm = TRUE), sum(s3,s4, na.rm = TRUE), s5, na.rm = TRUE)),
                   nrow = 4, ncol = 3, byrow = TRUE )
  table2 <- round(table2/N_imputations)
  
  colnames(table2) <- c("Unconstrained", "Constrained", "total")
  rownames(table2) <- c("1) SQ--Council--EP", "2) SQ--EP--Council", "3) Chamber1--SQ--Chamber2", "total")
  table2 <- xtable(table2, digits = 0, align = c("l","c","c","c"))
  print.xtable(table2, type = "html", file = paste(subdirs$tables,"/table2.html",sep=""))
  
  #_____________________________________________________________________
  # explaining the table
  #_____________________________________________________________________
  #
  # S9 = no change environment
  # EP closer to SQ: S1 - S4
  # S1 and S3 are constrained
  # Council closer to SQ: S5 - S8
  # S5 and S7 are constrained
  #_____________________________________________________________________
  
  # mean influence ratio
  mean.infl.ratio <- round(mean.ep.inf / N_imputations, digits = 2)
  
  # influence by first reading EA
  infl.by.ea.first <- round(infl.by.ea.first / N_imputations, digits = 2)
  
  # influence table
  x <- matrix( c(infl.by.ea.first[,2],mean.infl.ratio), nrow = 3, ncol = 1 )
  colnames(x) <- c("Relative Influence EP")
  rownames(x) <- c("EP Transparent", "EP Opaque", "Overall")
  influence.table <- xtable(x, digits = 2)
  print.xtable(influence.table, type = "html", file = paste(subdirs$tables,"/data_overview_section.html",sep=""))
  
  ##########################################
  # Summary stats table
  for (e in 1: N_imputations){
    
    df <- df.matched[[1]]$imputations[[e]] %>%
      filter(yearstart > 1998) %>%
      #mutate(strong = ifelse(ENVI==1|TRAN==1|JURI==1|CULT==1|ECON==1, 1, 0)) %>%
      mutate(recitals = sqrt(recitals)) %>%
      mutate(ctte_op = sqrt(ctte_op)) %>%
      mutate(words = log(words)) %>%
      mutate(time.trend = as.numeric(time.trend))
    
    # relative influence
    ifelse(e == 1, ep.inf <- df$ep.inf, ep.inf <- ep.inf + df$ep.inf)
    # 1st reading
    ifelse(e == 1, ea.first <- summary(df$ea.first), ea.first <- ea.first + summary(df$ea.first))
    # 2nd reading
    ifelse(e == 1, ea.second <- summary(df$ea.second), ea.second <- ea.second + summary(df$ea.second))
    # EU Enlargement
    ifelse(e == 1, post04 <- summary(df$post04), post04 <- post04 + summary(df$post04))
    # Com. for EP
    ifelse(e == 1, Com.support.EP <- summary(df$Com.support.EP), Com.support.EP <- Com.support.EP + 
             summary(df$Com.support.EP))
    # Relative Salience
    ifelse(e == 1, salience.relative <- summary(df$salience.relative), salience.relative <- salience.relative + summary(df$salience.relative))
    # Salience Commission
    ifelse(e == 1, scom <- summary(df$scom), scom <- scom + summary(df$scom))
    # 2nd Principle
    ifelse(e == 1, con_of_int <- summary(df$con_of_int), con_of_int <- con_of_int + summary(df$con_of_int))
    # number of committees asked for opinion
    ifelse(e == 1, ctte_op <- summary(df$ctte_op), ctte_op <- ctte_op + summary(df$ctte_op))
    # length of proposal
    ifelse(e == 1, words <- summary(df$words), words <- words + summary(df$words))
    # number of recitals
    ifelse(e == 1, recitals <- summary(df$recitals), recitals <- recitals + summary(df$recitals))
    # directive
    ifelse(e == 1, directive <- summary(df$directive), directive <- directive + summary(df$directive))
    # time trend
    ifelse(e == 1, time.trend <- summary(df$time.trend), time.trend <- time.trend + summary(df$time.trend))
  }
  
  # averaging
  ep.inf <- ep.inf / N_imputations
  ea.first <- ea.first / N_imputations
  ea.second <- ea.second / N_imputations
  post04 <- post04 / N_imputations
  Com.support.EP <- Com.support.EP / N_imputations
  salience.relative <- salience.relative / N_imputations
  scom <- scom / N_imputations
  con_of_int <- con_of_int / N_imputations
  ctte_op <- ctte_op / N_imputations
  words <- words / N_imputations
  recitals <- recitals / N_imputations
  directive <- directive / N_imputations
  time.trend <- time.trend / N_imputations
  
  # descriptive stats output table
  X3 <- as.matrix(rbind(summary(ep.inf), ea.first, ea.second, post04, Com.support.EP,
                        salience.relative, scom, con_of_int, ctte_op, words, recitals,
                        directive, time.trend))
  X3 <- X3[,c("Min.","Max.","Mean","Median")]
  
  # Online Appendix Summary Stats Table
  colnames(X3) <- c("Min.", "Max", "Mean", "Median")
  rownames(X3) <- c("Relative Influence", "EP Position Opaque", "Informal & Known EP Position",
                    "EU Enlargement", "Commission for EP", "Relative Salience",
                    "Salience Commission", "2nd Principle", "Committees Asked",
                    "Length Proposal", "Recitals", "Directive", "Time Trend")
  summary.stats.table <- xtable(X3, digits = 2)
  print.xtable(summary.stats.table, type = "html", file = paste(subdirs$tables,"/appendix_table2.html",sep=""))
  
  #### missingness tables
  # percent missings on these variables
  pmissing <- cbind(length(which(df$pcom.miss==1)) / 92,
                    length(which(df$pep.miss==1)) / 92,
                    length(which(df$sep.miss==1)) / 92,
                    length(which(df$rp.miss==1)) / 92,
                    length(which(df$out.miss==1)) / 92)
  colnames(pmissing) <- c("pcom","pep","sep","rp","out")
  pmissing
  
  # absolute number of missings
  amissing <- cbind(length(which(df$pcom.miss==1)),
                    length(which(df$pep.miss==1)),
                    length(which(df$sep.miss==1)),
                    length(which(df$rp.miss==1)),
                    length(which(df$out.miss==1)))
  colnames(amissing) <- c("pcom","pep","sep","rp","out")
  amissing
  
  # mean of imputed values and standard deviations
  imputed.vals <- matrix(NA, nrow = N_imputations, ncol = 8)
  colnames(imputed.vals) <- c("pcom.mean", "pcom.sd", "pep.mean", "pep.sd", "rp.mean", "rp.sd", "out.mean", "out.sd")
  tmp <- df %>%
    filter(pcom.miss == 1 | pep.miss == 1 | rp.miss == 1 | out.miss == 1)
  nrow(tmp)
  
}

s.time <- Sys.time()
# run everything
descriptives()
e.time <- Sys.time()
e.time - s.time
time.box$script10 <- e.time - s.time

# remove descriptives function from global environment
rm(descriptives)