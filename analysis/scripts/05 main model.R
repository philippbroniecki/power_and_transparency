# covariate balance before and after matching
# set up data to re-run analysis with matched sample
# run the main model with the parametric propensity score
# get substantial effects

main.model <- function(){

  # location of step 4 data
  setwd(subdirs$workdata)
  load("04 step main variables coded.RData")
  
  # create a matched data copy
  df.matched <- df.full
  rm(df.full)
  
  # saves simulated results
  results.box <- vector("list", length(df.matched))
  
  # container for balance statistics
  balance.box <- vector("list", N_imputations)
  
  # fit statistics
  fit.box <- rep(NA, N_imputations)

  # number of observations  
  n.box <- rep(NA, N_imputations)
    
  # set loop.counter
  loop.counter <- 0
  
  # loop over imputation sets
  for (f in 1: length(df.matched)){
    
    # models run
    models.run <- N_imputations
    
    # loop over imputation data sets
    for (e in 1: N_imputations){
      
      # loop counter
      loop.counter <- loop.counter + 1
      
      # data pre-processing 
      df <- df.matched[[f]][[1]][[e]] %>%
        # drop observations before 1999
        filter( yearstart > 1998 ) %>%
        # transformations for count variables
        mutate(recitals = sqrt(recitals)) %>%
        mutate(ctte_op = sqrt(ctte_op)) %>%
        mutate(words = log(words))
      
      # propensity score
      ps <- glm(ea.first ~ post04 + Com.support.EP + 
                  salience.relative + scom + con_of_int + 
                  ctte_op + words + directive + time.trend, 
                family = binomial(link = "logit"), data = df)
      
      # attach linear predicted (not to compress the ps close to 1 and 0)
      df$ps <- fitted.values(ps)
      
      ### just to check
      table(prediction = ifelse(df$ps>.5,1,0), actual = df$ea.first)
      
      # matching
      set.seed(123)
      m.mat <- Match(Y=df$ep.inf, Tr=df$ea.first, X=df$ps, ties = FALSE, Weight = 2)
      # summary(m.mat)
      
      # save the matching balance
      if (loop.counter < 1001) balance.box[[loop.counter]] <- MatchBalance(ea.first ~ ps, match.out = m.mat, data = df)
      
      # covariates
      X.sep <- as.matrix(cbind(
        m.mat$mdata$Tr,
        df[c(m.mat$index.treated, m.mat$index.control), c("ps")],
        df[c(m.mat$index.treated, m.mat$index.control), c("TRAN", "ECON", "LIBE", "AFET", "ITRE",
                                                          "REGI", "IMCO")]
      ))
      
      # response
      dv.sep <- m.mat$mdata$Y
      # cluster id
      cluster.sep <- df$prnrnmc[c(m.mat$index.treated,m.mat$index.control)]
      # variable names
      colnames(X.sep) <- c("EP Position Opaque", 
                           "Controls",
                           #"ENVI", # baseline
                           "TRAN",
                           "ECON",
                           "LIBE",
                           "AFET",
                           "ITRE",
                           "REGI",
                           "IMCO")
      
      # storing matched data as analysis data set
      df.matched[[f]][[1]][[e]] <- df[c(m.mat$index.treated, m.mat$index.control),]
      
      # flogit with propensity scores and committee fixed effects
      m.sep <- try(frm(y = dv.sep, x = X.sep, linkfrac = "logit", var.cluster = cluster.sep, var.type = "cluster"),
                   silent = TRUE)
      if (class(m.sep)=="try-error") {
        models.run <- models.run - 1
        next
      }
      
      # R^2
      fit.box[loop.counter] <- cor(m.sep$yhat, dv.sep)^2
      
      # Number of observations
      n.box[loop.counter] <- length(m.sep$yhat)
      
      # what is the number of unique controls
      length(unique(m.mat$index.control))
      
      # Simulations
      # $p returns the coefficients
      # $p.var returns the covariance 
      # generate sampling distrubtuion of coefficients
      S <- mvrnorm(10000, m.sep$p, m.sep$p.var)
      
      # combining the simulations to incorporate imputation uncertainty
      if (!exists("S.combined")){ 
        S.combined <- S
      } else { 
        S.combined <- rbind(S.combined, S)
      }
    } # end of loop over all impuations    
    
    results.box[[f]] <- S.combined
  }
 
  # coefficients and confidence intervals 50%, 90% CI, 95% CI, 99% CI
  outs <- apply(results.box[[1]], 2, quantile, probs = c(.5, .01, 0.99, 0.001, 0.999))
  # for the standard regression table
  sds <- apply(results.box[[1]], 2, sd)
  # fit statistics
  rsqrt <- mean(fit.box, na.rm = TRUE)
  # number of observations
  obs <- mean(n.box, na.rm = TRUE)
  
  # saving
  setwd(subdirs$workdata)
  box <- list(pointestimates = outs, standarderrors = sds, r2 = rsqrt, N = obs)
  save(box, file = "main_model.Rdata")
  #save(S.combined, file = "main_model_simulations.Rdata")
  save(balance.box, file = "balance_statistics_with_relative_salience.RData")

  
  #######################################################
  # balance
  #######################################################
  
  for (a.idx in 1: length(balance.box)){
    
    # before matching
    m.Tr.before <- ifelse(a.idx==1, balance.box[[a.idx]]$BeforeMatching[[1]]$mean.Tr,
                          m.Tr.before + balance.box[[a.idx]]$BeforeMatching[[1]]$mean.Tr)  
    m.Co.before <- ifelse(a.idx==1, balance.box[[a.idx]]$BeforeMatching[[1]]$mean.Co,
                          m.Co.before + balance.box[[a.idx]]$BeforeMatching[[1]]$mean.Co)
    Std.bias.before <- ifelse(a.idx==1, balance.box[[a.idx]]$BeforeMatching[[1]]$sdiff,
                              Std.bias.before + balance.box[[a.idx]]$BeforeMatching[[1]]$sdiff)
    tt.pval.before <- ifelse(a.idx==1, balance.box[[a.idx]]$BeforeMatching[[1]]$p.value,
                             tt.pval.before + balance.box[[a.idx]]$BeforeMatching[[1]]$p.value)
    ks.before <- ifelse(a.idx==1, balance.box[[a.idx]]$BeforeMatching[[1]]$ks$ks.boot.pvalue,
                        ks.before + balance.box[[a.idx]]$BeforeMatching[[1]]$ks$ks.boot.pvalue)
    # after matching
    m.Tr.after <- ifelse(a.idx==1, balance.box[[a.idx]]$AfterMatching[[1]]$mean.Tr,
                         m.Tr.after + balance.box[[a.idx]]$AfterMatching[[1]]$mean.Tr)  
    m.Co.after <- ifelse(a.idx==1, balance.box[[a.idx]]$AfterMatching[[1]]$mean.Co,
                         m.Co.after + balance.box[[a.idx]]$AfterMatching[[1]]$mean.Co)
    Std.bias.after <- ifelse(a.idx==1, balance.box[[a.idx]]$AfterMatching[[1]]$sdiff,
                             Std.bias.after + balance.box[[a.idx]]$AfterMatching[[1]]$sdiff)
    tt.pval.after <- ifelse(a.idx==1, balance.box[[a.idx]]$AfterMatching[[1]]$p.value,
                            tt.pval.after + balance.box[[a.idx]]$AfterMatching[[1]]$p.value)
    ks.after <- ifelse(a.idx==1, balance.box[[a.idx]]$AfterMatching[[1]]$ks$ks.boot.pvalue,
                       ks.after + balance.box[[a.idx]]$AfterMatching[[1]]$ks$ks.boot.pvalue)
  }
  b.out <- matrix(c(m.Tr.before,
                    m.Co.before,
                    Std.bias.before,
                    tt.pval.before,
                    ks.before,
                    m.Tr.after,
                    m.Co.after,
                    Std.bias.after,
                    tt.pval.after,
                    ks.after), nrow = 5, ncol = 2)
  b.out <- b.out / a.idx  
  colnames(b.out) <- c("Before Matching", "After Matching")
  rownames(b.out) <- c("Mean Treatment", "Mean Control", "Std. Bias", 
                       "T-test p-value", "K-S p-value")
  
  round(b.out,2)
  b.out[4:5,]
  
  # table 1
  table1 <- xtable(b.out, align = c("l","c","c"), NA.string = " ", sanitize.rownames.function = TRUE, digits = 2)
  print.xtable(table1, type = "html", file = paste(subdirs$tables,"/table1.html",sep=""))
  
  # save
  save(df.matched, file = "05 step analysis data with relative salience.RData")
  
  ######################################
  ## Substantial Effect interpretation
  ######################################
  
  for (e in 1: length(df.matched[[1]][[1]]) ){
    
    # load data
    df <- df.matched[[1]]$imputations[[e]]
    
    ## 1st Reading EAs
    # EP not opaque
    X1 <- cbind(1, 
                0, 
                mean(df$ps),
                median(df$TRAN),
                median(df$ECON),
                median(df$LIBE),
                median(df$AFET),
                median(df$ITRE),
                median(df$REGI),
                median(df$IMCO))
    
    # EP position opaque
    X2 <- cbind(1, 
                1, 
                mean(df$ps),
                median(df$TRAN),
                median(df$ECON),
                median(df$LIBE),
                median(df$AFET),
                median(df$ITRE),
                median(df$REGI),
                median(df$IMCO))
    
    # estimate latent
    latent1 <- S.combined %*% t(X1)
    latent2 <- S.combined %*% t(X2)  
    
    # through the logit link
    pv1 <- 1/ ( 1 + exp( -latent1))
    pv2 <- 1/ ( 1 + exp( -latent2))
    
    # get the first difference 
    fd.ea <- pv2 - pv1
    #fd.ea <- quantile(fd.ea, probs = c(0.025, 0.5, 975))
    
    ifelse(e == 1, fist.diff <- fd.ea, fist.diff <- c(fist.diff+fd.ea))
    
    cat(paste("\n", round(e/N_imputations, digits = 2), " percent done", sep = "" ))
  }
  
  fist.diff <- fist.diff / N_imputations
  fist.diff <- quantile(fist.diff, probs = c(0.025, 0.5, 0.975))
  setwd(subdirs$workdata)
  save(fist.diff, file = "substantial_predictions.RData")
  
  
}

# run
main.model()
rm(main.model)
