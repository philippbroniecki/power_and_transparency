# gradient boosting to estimate propensity score

w.non.parametric <- function(){

  setwd(subdirs$workdata)
  load("04 step main variables coded.RData")
  
  #_______________________________________________________________
  # treatment effect with gradient boosted propensity score
  #_______________________________________________________________
  # create a matched data copy
  df.matched2 <- df.full
  rm(df.full)
  
  # saves simulated results
  results.box <- vector("list", 2)
  
  # container for balance stats
  balance.box <- vector("list", N_imputations)
  
  # set loop.counter
  loop.counter <- 0
  
  # models run
  models.run <- N_imputations
  
  # fit statistics
  fit.box <- rep(NA, N_imputations)
  
  # number of observations  
  n.box <- rep(NA, N_imputations)
  
  # loop over imputation sets
  for (f in 1: length(df.matched2)){
    # loop over imputation data sets
    for (e in 1: N_imputations){
      
      # loop counter
      loop.counter <- loop.counter + 1
      
      # data pre-processing 
      df <- df.matched2[[f]][[1]][[e]] %>%
        # drop observations before 1999
        filter( yearstart > 1998 ) %>%
        # transformations for count variables
        mutate(recitals = sqrt(recitals)) %>%
        mutate(ctte_op = sqrt(ctte_op)) %>%
        mutate(words = log(words))
      
      data.train <- df %>%
        dplyr::select(
          ea.first, ep.inf, post04, Com.support.EP, 
          s_cou_w, scom, con_of_int, ctte_op,
          words, directive, time.trend
        ) %>%
        dplyr::mutate(
          time.trend = as.numeric(time.trend),
          #ea.first = as.factor(ea.first),
          post04 = as.factor(post04),
          Com.support.EP = as.factor(Com.support.EP),
          con_of_int = as.factor(con_of_int),
          directive = as.factor(directive)
        )
      
      # run gbm with optimal params
      set.seed(1234)
      ensemble <- gbm(
        ea.first ~ . -ep.inf, 
        distribution = "bernoulli", 
        data = data.train,
        n.trees = 400, 
        interaction.depth = 1, 
        n.minobsinnode = 5,
        bag.fraction = .5,
        # cv.folds = nrow(data.train),
        shrinkage = 0.03, 
        train.fraction = .8, 
        n.cores = 1)
      
      df$ps.nonparametric <- predict(
        ensemble, newdata = data.train,
        n.trees = 100,
        type = "response",allow.new.levels = TRUE
      )  
      table(actual = data.train$ea.first, prediction = ifelse(df$ps.nonparametric>.5,1,0))
      
      #_______________________________________________________________
      # matching on propensity score
      #_______________________________________________________________
      
      # matching
      set.seed(123)
      m.mat <- Match(Y=df$ep.inf, Tr=df$ea.first, X=df$ps.nonparametric, ties = FALSE, Weight = 2)
      # summary(m.mat)
      
      # save the matching balance
      if (loop.counter < 1001) balance.box[[loop.counter]] <- MatchBalance(ea.first ~ ps.nonparametric, match.out = m.mat, data = df)
      
      #_______________________________________________________________
      # fractional logit
      #_______________________________________________________________
      
      # covariates
      X.sep <- as.matrix(cbind(
        m.mat$mdata$Tr,
        df[c(m.mat$index.treated, m.mat$index.control), c("ps.nonparametric")],
        df[c(m.mat$index.treated, m.mat$index.control), c("TRAN", "ECON", "LIBE", "AFET", "ITRE","IMCO")]
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
                           #"REGI",
                           "IMCO")
      
      # storing matched data as analysis data set
      df.matched2[[f]][[1]][[e]] <- df[c(m.mat$index.treated, m.mat$index.control),]
      
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
      
      # Simulations
      # $p returns the coefficients
      # $p.var returns the covariance 
      # generate sampling distrubtuion of coefficients
      S <- mvrnorm(10000, m.sep$p, m.sep$p.var)
      
      # combining the simulations to incorporate imputation uncertainty
      if (e == 1){ 
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
  # fit statistics for to show that the models explain something
  rsqrt <- mean(fit.box,na.rm = TRUE)
  # number of observations
  obs <- mean(n.box, na.rm = TRUE)
  
  # saving
  setwd(subdirs$workdata)
  box_non_parametric <- list(pointestimates = outs , sds = sds, R2 = rsqrt, N = obs)
  save(box_non_parametric, file = "matching_model_with_relative_salience3.Rdata")
  
}

# run
w.non.parametric()