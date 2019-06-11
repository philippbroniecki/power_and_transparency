# first version 30.01.2019
# latest version 09.062019
# placebo test

placebo.test <- function(){
  
  setwd(subdirs$workdata)
  load("05 step analysis data with relative salience.RData")
  
  #_______________________________________________________________
  # treatment effect for placebo
  #_______________________________________________________________
  
  # container for simulated results
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
  for (f in 1: length(df.matched[1:1])){
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
        mutate(words = log(words)) %>%
        mutate(dur_months = log(dur_months)) 
      
      # covariates
      X.sep <- as.matrix(cbind(
        df$ea.second,
        #df$ps,
        df[, c("TRAN", "ECON", "LIBE", "AFET", "ITRE", "REGI", "IMCO")]
      ))
      
      # response
      dv.sep <- df$ep.inf
      # cluster id
      cluster.sep <- df$prnrnmc
      # variable names
      colnames(X.sep) <- c("Informal but transparent", 
                           #"Controls",
                           #"ENVI", # baseline
                           "TRAN",
                           "ECON",
                           "LIBE",
                           "AFET",
                           "ITRE",
                           "REGI",
                           "IMCO")
      
      # flogit with propensity scores and committee fixed effects
      m.sep <- try(frm(y = dv.sep, x = X.sep, linkfrac = "logit", var.cluster = cluster.sep, var.type = "cluster"),
                   silent = TRUE)
      if (class(m.sep)=="try-error") {
        models.run <- models.run - 1
        next
      }
      
      # R^2
      fit.box[loop.counter] <- cor(m.sep$yhat, df$ep.inf)^2
      
      # Number of observations
      n.box[loop.counter] <- length(m.sep$yhat)
      
      # Simulations
      # $p returns the coefficients
      # $p.var returns the covariance 
      # generate sampling distrubtuion of coefficients
      S <- mvrnorm(10000, m.sep$p, m.sep$p.var)
      
      # combining the simulations to incorporate imputation uncertainty
      if (e == 1){ S.combined <- S}
      else { S.combined <- rbind(S.combined, S)}
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
  box_placebo_test <- list(outs, sds, rsqrt, N = obs)
  save(box_placebo_test, file = "placebo_test_model.RData")
  
}

table4 <- function(){
  # table 4
  setwd(subdirs$workdata)
  load(file = "placebo_test_model.RData")
  
  x <- cbind(
    # rownames
    cbind(
      rbind("Informal Transparent Negotiations", " ", "Constant", " ", "N", "R^2 (correlation squared)")
    ),
    # output
    rbind(
      round(box_placebo_test[[1]]["50%", "Informal but transparent"], 2),
      paste("(", round(box_placebo_test[[2]]["Informal but transparent"], 2),")",sep=""),
      round(box_placebo_test[[1]]["50%", "INTERCEPT"], 2),
      paste("(", round(box_placebo_test[[2]]["INTERCEPT"], 2),")",sep=""),
      box_placebo_test$N,
      round(box_placebo_test[[3]],2)
    ))
  x    
  
  colnames(x) <- c(" ","Treatment Only")
  table4 <- xtable(x, align = c("l","c", "c"), NA.string = " ", sanitize.rownames.function = TRUE)
  print.xtable(table4, include.rownames = FALSE, include.colnames = TRUE, 
               type = "html", file = paste(subdirs$tables,"/table4.html",sep=""))
}

# run
placebo.test()
table4()
rm(placebo.test, table4)
