# treatment effect only in weighted salience model

w.treatment.only <- function(){
  
  setwd(subdirs$workdata)
  load("05 step analysis data with relative salience2.RData")
  
  #_______________________________________________________________
  # treatment effect without propensity score
  #_______________________________________________________________
  
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
        mutate(words = log(words)) 
      
      # covariates
      X.sep <- as.matrix(cbind(
        df$ea.first,
        #df$ps,
        df[, c("TRAN", "ECON", "LIBE", "AFET", "ITRE", "REGI", "IMCO")]
      ))
      
      # response
      dv.sep <- df$ep.inf
      # cluster id
      cluster.sep <- df$prnrnmc
      # variable names
      colnames(X.sep) <- c("EP Position Opaque", 
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
  box_treatment_only <- list(pointestimates = outs , sds = sds, R2 = rsqrt, N = obs)
  save(box_treatment_only, file = "matching_model_with_relative_salience4.Rdata")

}

# to plot the appendix table 4
appendix_table4 <- function(){

  # load model objects
  setwd(subdirs$workdata)
  load(file = "matching_model_with_relative_salience2.Rdata")
  load(file = "matching_model_with_relative_salience3.Rdata")
  load(file = "matching_model_with_relative_salience4.Rdata")
  
  # output
  x <- cbind(
    # rownames
    cbind(
      rbind("EP Position Opaque", " ", "Propensity Score", " ", "Constant", " ", "N", "R^2 (correlation squared)")
    ),
    # model with treatment effect only
    rbind(
      round(box_treatment_only[[1]]["50%", "EP Position Opaque"], 2),
      paste("(", round(box_treatment_only[[2]]["EP Position Opaque"], 2),")",sep=""),
      " ",
      " ",
      round(box_treatment_only[[1]]["50%", "INTERCEPT"], 2),
      paste("(", round(box_treatment_only[[2]]["INTERCEPT"], 2),")",sep=""),
      box_treatment_only$N,
      round(box_treatment_only[[3]],2)
    ),
    # model with parametric propensity score
    rbind(
      round(box$pointestimates["50%", "EP Position Opaque"], 2),
      paste("(", round(box$sds["EP Position Opaque"], 2),")",sep=""),
      round(box$pointestimates["50%", "Controls"], 2),
      paste("(", round(box$sds["Controls"], 2),")",sep=""),
      round(box$pointestimates["50%", "INTERCEPT"], 2),
      paste("(", round(box$sds["INTERCEPT"], 2),")",sep=""),
      box$N,
      round(box$R2, 2)
    ),
    # model with non-parametric propensity score
    rbind(
      round(box_non_parametric[[1]]["50%", "EP Position Opaque"], 2),
      paste("(", round(box_non_parametric[[2]]["EP Position Opaque"], 2),")",sep=""),
      round(box_non_parametric[[1]]["50%", "Controls"], 2),
      paste("(", round(box_non_parametric[[2]]["Controls"], 2),")",sep=""),
      round(box_non_parametric[[1]]["50%", "INTERCEPT"], 2),
      paste("(", round(box_non_parametric[[2]]["INTERCEPT"], 2),")",sep=""),
      box_non_parametric$N,
      round(box_non_parametric[[3]], 2)
    )
  )
  colnames(x) <- c(" ","Treatment Only", "Parametric Propensity Score", "Non-Parametric Propensity Score")
  table4_appendix <- xtable(x, align = c("l","c","c","c","c"), NA.string = " ", sanitize.rownames.function = TRUE)
  print.xtable(table4_appendix, include.rownames = FALSE, include.colnames = TRUE, 
               type = "html", file = paste(subdirs$tables,"/appendix_table4.html",sep=""))
    
}

# run
w.treatment.only()
appendix_table4()
