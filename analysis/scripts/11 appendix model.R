appendix_model <- function(){
  # load data
  setwd(subdirs$workdata)
  load("04 step main variables coded.RData")
  
  results.box <- vector("list", 2)
  
  # fit statistics
  fit.box <- rep(NA, N_imputations)
  
  loop.counter <- 0
  
  # loop over imputation sets
  for (f in 1: length(df.full[1:1])){
    # loop over imputation data sets
    for (e in 1: N_imputations){
      
      loop.counter <- loop.counter +1
      
      # data pre-processing 
      df <- df.full[[f]][[1]][[e]] %>%
        # drop observations before 1999
        dplyr::filter( yearstart > 1998 ) %>%
        # transformations for count variables
        dplyr::mutate(recitals = sqrt(recitals)) %>%
        dplyr::mutate(ctte_op = sqrt(ctte_op)) %>%
        dplyr::mutate(words = log(words)) 
      
      N.obs <- nrow(df)
      
      # covariates 
      # separate EA categories
      X.sep <- as.matrix(cbind(df$ep.inf,
                               df$prnrnmc, 
                               df$ea.first,
                               df$post04,
                               df$Com.support.EP,
                               df$salience.relative,
                               df$scom,
                               df$con_of_int,
                               df$ctte_op,
                               df$words,
                               df$directive,
                               df$time.trend,
                               df$TRAN,
                               df$JURI,
                               df$CULT,
                               df$ECON,
                               df$LIBE,
                               df$AFET,
                               df$ITRE,
                               df$REGI,
                               df$EMPL,
                               df$IMCO))
      
      # dropping missings, serparating out DV, cluster ID, and IVs
      X.sep <- na.omit(X.sep)
      dv.sep <- X.sep[, 1]
      cluster.sep <- X.sep[, 2]
      X.sep <- X.sep[, -c(1, 2)]
      colnames(X.sep) <- c("EP Position Opaque", 
                           "2004 Enlargement",
                           "Com Supports EP",
                           "Relative Salience",
                           "Salience Commission",
                           "2nd Principal",
                           "Committees Asked",
                           "Length of Proposal",
                           "Directive",
                           "Time Trend",
                           "TRAN",
                           "JURI",
                           "CULT",
                           "ECON",
                           "LIBE",
                           "AFET",
                           "ITRE",
                           "REGI",
                           "EMPL",
                           "IMCO")
      
      m.sep <- frm(y = dv.sep, x = X.sep, linkfrac = "logit", var.cluster = cluster.sep, var.type = "cluster")
      
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
      if (e == 1){ S.combined <- S}
      else { S.combined <- rbind(S.combined, S)}
    } # end of loop over all impuations    
    
    results.box[[f]] <- S.combined
  }
  
  # coefficients and confidence intervals 50%, 90% CI, 95% CI, 99% CI
  outs <- apply(results.box[[1]], 2, quantile, probs = c(.5, .05, .95, .025, .975, .01, 0.99))
  # for the standard regression table
  sds <- apply(results.box[[1]], 2, sd)
  ########################################################
  # result !!!
  outs
  
  # fit statistics
  rsqrt <- mean(fit.box,na.rm = TRUE)
  
  # saving
  setwd(subdirs$workdata)
  box <- list(outs, sds, rsqrt)
  
  save(box, file = "appendix_model.RData")
  
  # table 3 for the appendix
  appendix_table3 <- rbind(
    cbind(
      box[[1]]["50%",2:11],
      box[[2]][2:11]
    ),
    cbind(
      box[[1]]["50%","INTERCEPT"],
      box[[2]]["INTERCEPT"]
    ),
    cbind(N.obs,NA),
    cbind(box[[3]],NA)
  )
  rownames(appendix_table3) <- c("EP Position Opaque", "EU Enlargement", "Commission Supports EP",
                                 "Relative Salience", "Salience Commission", "2nd Prinicpal",
                                 "Committees Asked", "Length of Proposal", "Directive", "Time Trend", 
                                 "Constant", "N", "R2")
  colnames(appendix_table3) <- c("beta", "se")
  
  appendix_table3 <- xtable(appendix_table3, digits = 2)
  print.xtable(appendix_table3, type = "html", file = paste(subdirs$tables,"/appendix_table3.html",sep=""))
  
}

# run
s.time <- Sys.time()
# run everything
appendix_model()
e.time <- Sys.time()
e.time - s.time
time.box$script11 <- e.time - s.time