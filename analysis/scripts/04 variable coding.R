# define functions

# Commission prefers EP position over Council position
com.supp.func <- function(pep, pcou, pcom){
  com.sup <- as.numeric(abs(pcom - pep) < abs(pcom - pcou))
  return(com.sup)
}

# SQ in between institutions function   
no.change.func <- function(pep, rp, pcou){
  no.change <- as.numeric((pep < rp & pcou > rp) | (pcou < rp & pep > rp)) 
  return(no.change)
}

# midpoint between both insitutions
mid.func <- function(pep, pcou){
  mid <- (pep + pcou) / 2
  return(mid)
}

# calculate the counterfactual based on the policy constellation
counterfact.func <- function(pep, pcou, rp, rp.miss, use){
  
  # counterctual expectation
  cf <- rep(NA, length(pep))
  
  # is the status quo in between both institutions
  nochange <- no.change.func(pep = pep, rp = rp, pcou = pcou)
  
  # midpoint between both institutions
  midpoint <- mid.func(pep = pep, pcou = pcou)
  
  # indifference of actor closest to sq
  indiff <- rep(NA, length(pep))
  
  ###### Policy constellations #####################################
  ## In scenarios 1 the SQ in not in beween both institutions and 
  ## the EP is closer to the SQ than the Council
  # Scenario 1a: --SQ---EP-------Cou---- 
  scenario1a <- nochange == 0 & abs(pep - rp) < abs(pcou - rp) & pep < pcou
  indiff <- ifelse(scenario1a == TRUE, yes = pep + abs(rp - pep), no = indiff)
  s1a.constrained <- scenario1a == TRUE & indiff < midpoint
  s1a.free <- scenario1a == TRUE & indiff > midpoint 
  
  # Scenario 1b: --Cou---------EP----SQ-
  scenario1b <- nochange == 0 & abs(pep - rp) < abs(pcou - rp) & pep > pcou  
  indiff <- ifelse(scenario1b == TRUE, yes = pep - abs(rp - pep), no = indiff)
  s1b.constrained <- scenario1b == TRUE & indiff > midpoint
  s1b.free <- scenario1b == TRUE & indiff < midpoint
  
  ## In scenarios 2 the SQ in not in beween both institutions and 
  ## the Council is closer to the SQ than the EP
  # Scenario 2a: --SQ---Cou--------EP---
  scenario2a <- nochange == 0 & abs(pcou - rp) < abs(pep - rp) & pcou < pep
  indiff <- ifelse(scenario2a == TRUE, yes = pcou + abs(rp - pcou), no = indiff)
  s2a.constrained <- scenario2a == TRUE & indiff < midpoint
  s2a.free <- scenario2a == TRUE & indiff > midpoint
  
  # Scenario 2b --EP----------SQ---Cou-
  scenario2b <- nochange == 0 & abs(pcou - rp) < abs(pep - rp) & pcou > pep
  indiff <- ifelse(scenario2b == TRUE, yes = pcou - abs(rp - pcou), no = indiff)
  s2b.constrained <- scenario2b == TRUE & indiff > midpoint
  s2b.free <- scenario2b == TRUE & indiff < midpoint
  
  ## assign counterfactual
  cf <- ifelse(s1a.constrained | s1b.constrained | s2a.constrained | s2b.constrained,
               yes = indiff, no = cf)
  cf <- ifelse(s1a.free | s1b.free | s2a.free | s2b.free,
               yes = midpoint, no = cf)
  cf <- ifelse(nochange == 1, yes = rp, no = cf)
  cf <- ifelse(pcou == pep, yes = pcou, no = cf)
  
  # if imputed reference points are not used set the counterfactual to the midpoint
  cf <- ifelse(use=="rp.not.imputed" & rp.miss == 1, yes = midpoint, no = cf)
  
  # return the counterfactual
  return(cf)
} #endof function

# influence function
real.inf <- function(pcou, pep, out, cf){
  
  ## max gains and losses to the EP
  # max gain
  max.gain <- abs(pep - cf)
  # max loss
  dist.to.pole <- ifelse(abs(pep-100) > abs(pep-0), yes = abs(pep-100), no = abs(pep-0) )
  max.loss <- dist.to.pole - max.gain
  
  # actual
  actual <- abs(pep - cf) - abs(pep - out)
  # well behaved outcome, i.e. out at either actor or between them
  normal <- pep <= out & out <= pcou | pcou <= out & out <= pep
  # EP closer to the outcome than Council
  ep.closer <- abs(pep - out) < abs(pcou - out)
  max.loss <- ifelse(normal==FALSE & ep.closer == TRUE, yes = max.loss, no = abs(pcou-cf))
  
  # percentage loss/gain to total possible
  p.of.total <- ifelse(actual > 0, yes = actual / max.gain,
                       no = actual / max.loss)
  # this adjusts for the case when outcome is outside the EP-Cou interval, 
  # and the EP is further from the outcome
  #p.of.total <- ifelse(normal==FALSE & ep.closer == FALSE,  yes = max.loss / actual, no = p.of.total)
  p.of.total[which(is.na(p.of.total))] <- 0
  
  # change from the counterfactual rescaled to what is left
  p.change <- p.of.total * .5
  
  # EP's relative influence
  ep.inf <- ifelse(normal==TRUE, 
                   yes = 0.5 + p.change, 
                   no = ifelse(ep.closer == TRUE, yes = 0.5 + p.change, no = (0.5 + p.change) * -1))
  
  # stay within the unit interval (this is only a problem if imputation returned
  # positions outside the DEU scale)
  ep.inf <- ifelse(ep.inf < 0, yes = 0, no = ep.inf)
  ep.inf <- ifelse(ep.inf > 1, yes = 1, no = ep.inf)
  
  # normal adjustment for strange outcomes
  ep.inf <- ifelse(normal==FALSE & ep.closer==FALSE, yes = 0, no = ep.inf)
  ep.inf <- ifelse(normal==FALSE & ep.closer==TRUE, yes = 1, no = ep.inf)
  
  return(ep.inf)
}

# control influence function
cont.inf <- function(pcou, pep, out){
  ep.inf <- ifelse( abs(pcou - out) + abs(pep - out) != 0,
                    yes = abs(pcou - out) / (abs(pcou - out) + abs(pep - out)),
                    no = 0.5)
  return(ep.inf)
}


variable.coding <- function(){

  # location of imputation data
  setwd(subdirs$workdata)
  load("03 step imputation.RData")
  
  # loop over f impuation sets
  for (f in 1:length(df.full)){
    # loop over e impuations
    for (e in 1:N_imputations){ 
      
      # 1) Commission support - binary (1 = Commission is closer to EP than Council)
      df.full[[f]][[1]][[e]]$Com.support.EP <- com.supp.func(pep = df.full[[f]][[1]][[e]]$pep,
                                                             pcou = df.full[[f]][[1]][[e]]$pcou,
                                                             pcom = df.full[[f]][[1]][[e]]$pcom)
      
      # 2) years since EA introduction
      # converting startdate to R date format
      df.full[[f]][[1]][[e]]$start <- as.Date( df.full[[f]][[1]][[e]]$start, "%d/%m/%Y")
      # introduction of early agreements
      ea.introduction <- "01/05/1999"
      df.full[[f]][[1]][[e]]$Amsterdam <- rep(as.Date(ea.introduction, "%d/%m/%Y"))
      # time passed since introduction of early agreents
      df.full[[f]][[1]][[e]]$time.trend <- difftime( df.full[[f]][[1]][[e]]$start,
                                                     df.full[[f]][[1]][[e]]$Amsterdam, units = "weeks")
      
      # 2a) Counterfactual
      cf <- counterfact.func(pep = df.full[[f]][[1]][[e]]$pep,
                             rp = df.full[[f]][[1]][[e]]$rp,
                             pcou = df.full[[f]][[1]][[e]]$pcou,
                             rp.miss = df.full[[f]][[1]][[e]]$rp.miss,
                             use = "rp.imputed")
      
      # this counterfactual does not use imputed reference points
      cf2 <- counterfact.func(pep = df.full[[f]][[1]][[e]]$pep,
                              rp = df.full[[f]][[1]][[e]]$rp,
                              pcou = df.full[[f]][[1]][[e]]$pcou,
                              rp.miss = df.full[[f]][[1]][[e]]$rp.miss,
                              use = "rp.not.imputed")
      
      # this counterfactual is the standard approach, i.e. always the midpoint
      cf3 <- mid.func(pep = df.full[[f]][[1]][[e]]$pep, 
                      pcou = df.full[[f]][[1]][[e]]$pcou)
      
      #___________________________
      # relative salience: salience ep / (salience ep + salience council)
      df.full[[f]][[1]][[e]]$salience.relative <-  df.full[[f]][[1]][[e]]$sep / 
        (df.full[[f]][[1]][[e]]$sep + df.full[[f]][[1]][[e]]$scou)
      
      #################################################
      ### DV
      
      # real relative influence of the EP
      df.full[[f]][[1]][[e]]$ep.inf <- real.inf(pcou = df.full[[f]][[1]][[e]]$pcou,
                                                pep = df.full[[f]][[1]][[e]]$pep,
                                                out = df.full[[f]][[1]][[e]]$out,
                                                cf = cf)
      
      # real relative influence of the EP (imputed reference points not used)
      df.full[[f]][[1]][[e]]$ep.inf.alt <- real.inf(pcou = df.full[[f]][[1]][[e]]$pcou,
                                                    pep = df.full[[f]][[1]][[e]]$pep,
                                                    out = df.full[[f]][[1]][[e]]$out,
                                                    cf = cf2)
      
      # real relative influence of the EP (counterfactual always the midpoint)
      df.full[[f]][[1]][[e]]$ep.inf.std <- real.inf(pcou = df.full[[f]][[1]][[e]]$pcou,
                                                    pep = df.full[[f]][[1]][[e]]$pep,
                                                    out = df.full[[f]][[1]][[e]]$out,
                                                    cf = cf3)
      
      # the standard measure of relative influence
      df.full[[f]][[1]][[e]] <- df.full[[f]][[1]][[e]] %>%
        mutate(ep.inf.old = ifelse((abs(pep-out) + abs(pcou-out)) != 0,
                                   yes = (abs(pcou-out) / (abs(pep-out) + abs(pcou-out))),
                                   no = 0.5))
      
    } # end of loop over e impuations
  } # end of loop over f imputations sets
  
  # writing new file
  save(df.full, file="04 step main variables coded.RData")
}

# run
variable.coding()
rm(com.supp.func, counterfact.func, real.inf,
   mid.func, no.change.func, cont.inf, variable.coding)
