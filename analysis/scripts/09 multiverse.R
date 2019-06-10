#######################################################
# multiverse analysis (DOI: 10.1177/1745691616658637)
# current version 22.05.2019
#######################################################

multiverse <- function(){
  
  # load data
  setwd(subdirs$workdata)
  load("05 step analysis data with relative salience.RData")
  load("df.matched2 relative salience.RData")
  
  # counterfactuals
  DVs <- grep("ep.inf", names(df.matched[[1]][[1]][[1]]), value = TRUE)
  
  # number of data sets
  N.dfs <- length(df.matched) * length(DVs)
  N.dfs2 <- length(df.matched2) * length(DVs)
  N.dfs <- N.dfs + N.dfs2
  
  # stores effect of intransparency seperately
  transp.multi <- matrix(NA, nrow = 1000 * N_imputations, ncol = N.dfs)
  
  # actual models run counter
  models.run <- nrow(transp.multi)
  
  # counters
  counter <- 0
  counter_no_imps <- 0
  count_imps <- 0
  
  # loop over Council models
  # 1 Banzhaf weights
  # 2 Shapley-Shubrik weights
  for (a in 1:(length(df.matched) + length(df.matched2))) {
    # loop over 3 counterfactuals (midpoint; SQ constraint with imputed; SQ constraint without imputed)
    for (b in 1:length(DVs)){
      counter_no_imps <- counter_no_imps + 1
      # loop over imputations
      for (c in 1:N_imputations){
        counter <- counter + 1
        
        # current power weight
        if (a <= length(df.matched)){
          df <- df.matched[[a]]
        } else {
          df <- df.matched2[[a-2]]
        }
        
        # current imputation
        df <- df$imputations[[c]]
        
        # current counterfactual
        DV <- DVs[b]
        
        # current propensity score
        if (a <= length(df.matched)){
          ps <- df$ps
        } else {
          ps <- df$ps.nonparametric
        }
        
        # DV within the unit interval
        df[,DV][which(df[,DV] < 0)] <- 0
        
        # covariates 
        X.sep <- as.matrix(cbind(df[,DV],
                                 df$prnrnmc, 
                                 df$ea.first,
                                 df$ps,
                                 df$TRAN,
                                 df$ECON,
                                 df$LIBE,
                                 df$AFET,
                                 df$ITRE,
                                 #df$REGI,
                                 df$IMCO))
        
        X.sep <- na.omit(X.sep)
        dv.sep <- X.sep[, 1]
        cluster.sep <- X.sep[, 2]
        X.sep <- X.sep[, -c(1, 2)]
        colnames(X.sep) <- c("EP Position Opaque", 
                             "Controls",
                             "TRAN",
                             "ECON",
                             "LIBE",
                             "AFET",
                             "ITRE",
                             #"REGI",
                             "IMCO")
        
        
        try(m.sep <- frm(y = dv.sep, x = X.sep, linkfrac = "logit", var.cluster = cluster.sep, var.type = "cluster"),
            silent = TRUE)
        if (!exists("m.sep")){
          models.run <- models.run - 1
          next
        }
        if (class(m.sep)=="try-error") {
          models.run <- models.run - 1
          next
        }
        
        # Simulations
        # $p returns the coefficients
        # $p.var returns the covariance 
        # generate sampling distrubtuion of coefficients
        S <- mvrnorm(1000, m.sep$p, m.sep$p.var)
        
        # stacking simulations from the multiverse
        if (c == 1){ S.combined <- S}
        else { S.combined <- rbind(S.combined, S)}
        
      } # end of loop over imputatation data sets
      
      # seperately store the effect of intransparency with multiverse index
      transp.multi[,counter_no_imps] <- S.combined[,"EP Position Opaque"]
      
      # keep stacking for overall multiverse results
      if (counter_no_imps==1) {multiverse <- S.combined}
      else { multiverse <- rbind(multiverse, S.combined)}
      rm(S.combined)
      
    } # end of loop over counterfactuals
  } # end of loop over power weights
  
  # coefficients and confidence intervals 50%, 90% CI, 95% CI, 99% CI
  outs <- apply(multiverse, 2, quantile, probs = c(.5, .05, .95, .025, .975, .01, 0.99))
  # for the standard regression table
  sds <- apply(multiverse, 2, sd)
  
  ########################################################
  # result !!!
  outs
  # saving 
  setwd(subdirs$workdata)
  save(transp.multi, file = "datachoices.RData")
  #######################################################

  # multiple model specifications
  #############################################################################
  
  rm(df, multiverse, outs, S, transp.multi, X.sep, a, b, c, cluster.sep, 
     counter, count_imps, counter_no_imps, dv.sep, sds)
  
  # baseline model
  base.vars <- c("prnrnmc","ea.first", "ps", "TRAN", "ECON", "LIBE", "AFET",
                 "ITRE")
  # variables that will be added iteratively
  add.vars <- c("IMCO", "REGI")
  
  # number of specifications
  N.dfs <- length(df.matched) * length(DVs) * (length(add.vars)+1)
  
  # stores the effect of transparency
  transp.multi <- matrix(NA, nrow=1000*N_imputations, ncol = N.dfs)
  
  # counters
  counter <- 0
  counter_no_imps <- 0
  count_imps <- 0
  add.var.counter <- 0
  
  # loop over Council models
  # 1 compromise model with Banzhaf weights
  # 2 compromise model with Shapley-Shubrik weights
  # 3 simple average
  # 4 pivot model
  for (a in 1:length(df.matched)) {
    # loop over 3 counterfactuals (midpoint; SQ constraint with imputed; SQ constraint without imputed)
    for (b in 1:length(DVs)){
      # loop over model formulas
      for (d in 1:(length(add.vars)+1)){
        counter_no_imps <- counter_no_imps+1
        # loop over imputations
        for (c in 1:N_imputations){
          counter <- counter+1
          
          # current power weight
          df <- df.matched[[a]]
          # current imputation
          df <- df$imputations[[c]]
          # current counterfactual
          DV <- DVs[b]
          
          # data pre-processing
          df <- df %>%
            # only cases from 1999 onwards
            filter( yearstart > 1998 ) %>%
            # transformations for count variables
            mutate(recitals = sqrt(recitals)) %>%
            mutate(ctte_op = sqrt(ctte_op)) %>%
            mutate(words = log(words)) %>%
            mutate(dur_months = sqrt(dur_months)) %>%
            # strong and weak committees
            mutate(strong = ifelse(ENVI==1|TRAN==1|JURI==1|CULT==1|ECON==1, 1, 0)) %>%
            mutate(weak = ifelse(LIBE==1|AFET==1|ITRE==1|REGI==1|EMPL==1|IMCO==1, 1, 0)) %>%
            mutate(time.trend = as.numeric(time.trend))
          
          # DV within the unit interval
          df[,DV][which(df[,DV] < 0)] <- 0
          
          # model matrix with DV and cluster ID
          X.sep <- as.matrix(cbind(DV=df[,DV],
                                   df[,c(base.vars)]))
          if (add.var.counter>0) {
            X.sep <- as.matrix(cbind(X.sep, df[,c(add.vars[1:add.var.counter])]))
            colnames(X.sep) <- c("DV",base.vars, add.vars[1:add.var.counter])
          }
          X.sep <- na.omit(X.sep)
          dv.sep <- X.sep[, "DV"]
          cluster.sep <- X.sep[, "prnrnmc"]
          X.sep <- X.sep[, -c(1, 2)]
          
          # run model
          m.sep <- try(frm(y = dv.sep, x = X.sep, linkfrac = "logit", var.cluster = cluster.sep, var.type = "cluster"),
                       silent =  TRUE)
          if (class(m.sep) == "try-error"){
            next
          }
          
          # Simulations
          # $p returns the coefficients
          # $p.var returns the covariance 
          # generate sampling distrubtuion of coefficients
          S <- mvrnorm(1000, m.sep$p, m.sep$p.var)
          
          # stacking simulations from imputations
          if (c == 1){ S.combined <- S}
          else { S.combined <- rbind(S.combined, S)}
          
        } # end of loop over imputations
        add.var.counter <- add.var.counter+1
        
        # seperately store the effect of transparency for each specification
        transp.multi[1:length(S.combined[,base.vars[2]]),counter_no_imps] <- S.combined[,base.vars[2]]
        rm(S.combined)
        
      } # end of loop over model formulas
      add.var.counter <- 0
    } # end of loop over counterfactuals
  } # end of loop over council models
  
  # saving
  setwd(subdirs$workdata)
  save(transp.multi, file = "allvariantes.RData")  
  
}

# run everything
s.time <- Sys.time()
multiverse()
e.time <- Sys.time()
e.time - s.time
time.box$script9 <- e.time - s.time
rm(multiverse)