# current version: 22 June 2019 (functions for enclosures)
# load voting weights
# calculate position of the Council of Ministers according to the compromise model
# using 1) Banzhaf power index; 2) Shapley-Shubik power index

council.position <- function(){
  
  # location of Council voting weights with power indices
  setwd(subdirs$council.voting)
  eu15 <- read.csv("EU15.csv", header=TRUE, sep=",")
  eu27 <- read.csv("EU27.csv", header=TRUE, sep=",")
  
  # location of dataset from step 1)
  setwd(subdirs$workdata)
  load("01 step merged.RData")
  
  # 1) Using the Banzhaf power index
  
  # define the council position variable
  df$pcou_b <- rep(NA, nrow(df))
  
  # define the council salience variable
  df$scou <- rep(NA, nrow(df))
  
  # loop over the data set
  for (a in 1:length(df[,1])){ # start of loop over each row in data
    
    # EU 15
    if (df$post04[a] == 0){ # start of EU 15 if condition
      print("EU15")
      position <- c(df$pat[a], df$pbe[a], df$pdk[a], df$pfi[a], df$pfr[a], df$pde[a], df$pel[a],
                    df$pie[a], df$pit[a], df$plu[a], df$pnl[a], df$ppt[a], df$pes[a], df$pse[a],
                    df$puk[a])
      salience <- c(df$sat[a], df$sbe[a], df$sdk[a], df$sfi[a], df$sfr[a], df$sde[a], df$sel[a],
                    df$sie[a], df$sit[a], df$slu[a], df$snl[a], df$spt[a], df$ses[a], df$sse[a],
                    df$suk[a])
      banzhaf <- c(eu15[,"banzhaf"]) 
      eu15.matrix <- as.matrix(cbind(position,salience,banzhaf), nrow = 15, ncol = 3)
      
      # eliminiating missings
      eu15.matrix <- na.omit(eu15.matrix)
      
      # upper product
      up.prod <- eu15.matrix[,1] * eu15.matrix[,2] * eu15.matrix[,3]
      eu15.matrix <- cbind(eu15.matrix, up.prod)
      
      # lower product
      low.prod <- eu15.matrix[,2] * eu15.matrix[,3]
      eu15.matrix <- cbind(eu15.matrix, low.prod)
      
      # council position
      df$pcou_b[a] <- apply(eu15.matrix, 2, sum)[4] / apply(eu15.matrix, 2, sum)[5]
      print(eu15.matrix)
      
      # council salience
      df$scou[a] <- apply(eu15.matrix, 2, sum)[2] / length(eu15.matrix[,2])
      
    } # end of EU 15 if condition
    
    # EU 27
    if (df$post04[a] == 1){ # start of EU 15 if condition
      print("EU27")
      position <- c(df$pat[a], df$pbe[a], df$pbu[a], df$pcy[a], df$pcz[a], df$pdk[a], df$pee[a],
                    df$pfi[a], df$pfr[a], df$pde[a], df$pel[a], df$phu[a], df$pie[a], df$pit[a],
                    df$plv[a], df$plt[a], df$plu[a], df$pmt[a], df$pnl[a], df$ppl[a], df$ppt[a],
                    df$pro[a], df$psi[a], df$psk[a], df$pes[a], df$pse[a], df$puk[a])
      salience <- c(df$sat[a], df$sbe[a], df$sbu[a], df$scy[a], df$scz[a], df$sdk[a], df$see[a],
                    df$sfi[a], df$sfr[a], df$sde[a], df$sel[a], df$shu[a], df$sie[a], df$sit[a],
                    df$slv[a], df$slt[a], df$slu[a], df$smt[a], df$snl[a], df$spl[a], df$spt[a],
                    df$sro[a], df$ssi[a], df$ssk[a], df$ses[a], df$sse[a], df$suk[a])
      banzhaf <- c(eu27[,"banzhaf"]) 
      eu27.matrix <- as.matrix(cbind(position,salience,banzhaf), nrow = 27, ncol = 3)
      
      # eliminiating missings
      eu27.matrix <- na.omit(eu27.matrix)
      
      # upper product
      up.prod <- eu27.matrix[,1] * eu27.matrix[,2] * eu27.matrix[,3]
      eu27.matrix <- cbind(eu27.matrix, up.prod)
      
      # lower product
      low.prod <- eu27.matrix[,2] * eu27.matrix[,3]
      eu27.matrix <- cbind(eu27.matrix, low.prod)
      
      # council position
      df$pcou_b[a] <- apply(eu27.matrix, 2, sum)[4] / apply(eu27.matrix, 2, sum)[5]
      print(eu27.matrix)
      
      # council salience
      df$scou[a] <- apply(eu27.matrix, 2, sum)[2] / length(eu27.matrix[,2])
      
    } # end of EU 27 if condition
  } # end of loop over each row in data
  
  # 2) Shapley Shubrik
  
  # define the salience weighted council position
  df$pcou_ss <- rep(NA,length(df[,1]))
  
  # loop over the data set
  for (a in 1:length(df[,1])){ # start of loop over each row in data
    
    # EU 15
    if (df$post04[a] == 0){ # start of EU 15 if condition
      print("EU15")
      position <- c(df$pat[a], df$pbe[a], df$pdk[a], df$pfi[a], df$pfr[a], df$pde[a], df$pel[a],
                    df$pie[a], df$pit[a], df$plu[a], df$pnl[a], df$ppt[a], df$pes[a], df$pse[a],
                    df$puk[a])
      salience <- c(df$sat[a], df$sbe[a], df$sdk[a], df$sfi[a], df$sfr[a], df$sde[a], df$sel[a],
                    df$sie[a], df$sit[a], df$slu[a], df$snl[a], df$spt[a], df$ses[a], df$sse[a],
                    df$suk[a])
      shapleyshubrik <- c(eu15[,"shapleyshubik"]) 
      eu15.matrix <- as.matrix(cbind(position,salience,shapleyshubrik), nrow = 15, ncol = 3)
      
      # eliminiating missings
      eu15.matrix <- na.omit(eu15.matrix)
      
      # upper product
      up.prod <- eu15.matrix[,1] * eu15.matrix[,2] * eu15.matrix[,3]
      eu15.matrix <- cbind(eu15.matrix, up.prod)
      
      # lower product
      low.prod <- eu15.matrix[,2] * eu15.matrix[,3]
      eu15.matrix <- cbind(eu15.matrix, low.prod)
      
      # council position
      df$pcou_ss[a] <- apply(eu15.matrix, 2, sum)[4] / apply(eu15.matrix, 2, sum)[5]
      print(eu15.matrix)
      
    } # end of EU 15 if condition
    
    # EU 27
    if (df$post04[a] == 1){ # start of EU 15 if condition
      print("EU27")
      position <- c(df$pat[a], df$pbe[a], df$pbu[a], df$pcy[a], df$pcz[a], df$pdk[a], df$pee[a],
                    df$pfi[a], df$pfr[a], df$pde[a], df$pel[a], df$phu[a], df$pie[a], df$pit[a],
                    df$plv[a], df$plt[a], df$plu[a], df$pmt[a], df$pnl[a], df$ppl[a], df$ppt[a],
                    df$pro[a], df$psi[a], df$psk[a], df$pes[a], df$pse[a], df$puk[a])
      salience <- c(df$sat[a], df$sbe[a], df$sbu[a], df$scy[a], df$scz[a], df$sdk[a], df$see[a],
                    df$sfi[a], df$sfr[a], df$sde[a], df$sel[a], df$shu[a], df$sie[a], df$sit[a],
                    df$slv[a], df$slt[a], df$slu[a], df$smt[a], df$snl[a], df$spl[a], df$spt[a],
                    df$sro[a], df$ssi[a], df$ssk[a], df$ses[a], df$sse[a], df$suk[a])
      shapleyshubrik <- c(eu27[,"shapleyshubik"]) 
      eu27.matrix <- as.matrix(cbind(position,salience,shapleyshubrik), nrow = 27, ncol = 3)
      
      # eliminiating missings
      eu27.matrix <- na.omit(eu27.matrix)
      
      # upper product
      up.prod <- eu27.matrix[,1] * eu27.matrix[,2] * eu27.matrix[,3]
      eu27.matrix <- cbind(eu27.matrix, up.prod)
      
      # lower product
      low.prod <- eu27.matrix[,2] * eu27.matrix[,3]
      eu27.matrix <- cbind(eu27.matrix, low.prod)
      
      # council position
      df$pcou_ss[a] <- apply(eu27.matrix, 2, sum)[4] / apply(eu27.matrix, 2, sum)[5]
      print(eu27.matrix)
      
    } # end of EU 27 if condition
  } # end of loop over each row in data
  
  # define Council salience variable
  df$s_cou_w <- rep(NA,length(df[,1]))
  
  # loop over the data set
  for (a in 1:length(df[,1])){ # start of loop over each row in data
    
    # EU 15
    if (df$post04[a] == 0){ # start of EU 15 if condition
      print("EU15")
      salience <- c(df$sat[a], df$sbe[a], df$sdk[a], df$sfi[a], df$sfr[a], df$sde[a], df$sel[a],
                    df$sie[a], df$sit[a], df$slu[a], df$snl[a], df$spt[a], df$ses[a], df$sse[a],
                    df$suk[a])
      shapleyshubrik <- c(eu15[,"shapleyshubik"]) 
      eu15.matrix <- as.matrix(cbind(salience,shapleyshubrik), nrow = 15, ncol = 2)
      
      # eliminiating missings
      eu15.matrix <- na.omit(eu15.matrix)
      
      # council position
      df$s_cou_w[a] <- eu15.matrix[,1] %*% eu15.matrix[,2]
      print(eu15.matrix)
      
    } # end of EU 15 if condition
    
    # EU 27
    if (df$post04[a] == 1){ # start of EU 15 if condition
      print("EU27")
      salience <- c(df$sat[a], df$sbe[a], df$sbu[a], df$scy[a], df$scz[a], df$sdk[a], df$see[a],
                    df$sfi[a], df$sfr[a], df$sde[a], df$sel[a], df$shu[a], df$sie[a], df$sit[a],
                    df$slv[a], df$slt[a], df$slu[a], df$smt[a], df$snl[a], df$spl[a], df$spt[a],
                    df$sro[a], df$ssi[a], df$ssk[a], df$ses[a], df$sse[a], df$suk[a])
      shapleyshubrik <- c(eu27[,"shapleyshubik"]) 
      eu27.matrix <- as.matrix(cbind(salience,shapleyshubrik), nrow = 27, ncol = 2)
      
      # eliminiating missings
      eu27.matrix <- na.omit(eu27.matrix)
      
      # council position
      df$s_cou_w[a] <- eu27.matrix[,1] %*% eu27.matrix[,2]
      print(eu27.matrix)
      
    } # end of EU 27 if condition
  } # end of loop over each row in data
  
  # saving
  save(df, file = "02 step pcou added.RData")
  
}

s.time <- Sys.time()
# run everything
council.position()
e.time <- Sys.time()
e.time - s.time
time.box$script2 <- e.time - s.time

# remove council function from global environment
rm(council.position)