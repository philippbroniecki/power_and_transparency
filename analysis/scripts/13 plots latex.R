# plots the resutls
# current version 22.05.2019

# to run this specify ghost script file location and install libraries
library(rjson)
library(fontcm) # install computer modern font
library(extrafont) # to use computer modern
font_install('fontcm')
loadfonts() # to use computer modern
# ghost script location (to imbed fonts in PDFs)
Sys.setenv(R_GSCMD = ghostscript.file.path) 

#_______________________________________________________________________________
# data choices plot from multiverse analysis
# saving in graphs directory as datachoices.pdf
#_______________________________________________________________________________
setwd(subdirs$workdata)
load("datachoices.RData")
load("05 step analysis data with relative salience.RData")
setwd(subdirs$graphs)
pdf("appendix_figure1.pdf", family="CM Roman", width = 9, height = 7)
par(mfrow = c(1,1), mai=c(0,0,0,0), xpd=TRUE)
# the effect of intransparency by data set
betas <- apply(transp.multi, 2, mean)
sigmas <- apply(transp.multi, 2, sd)
# counterfactuals
DVs <- grep("ep.inf", names(df.matched[[1]][[1]][[1]]), value = TRUE)
# number of data sets
N.dfs <- (length(df.matched)*2) * length(DVs)

# figure
plot(1,
     type = "n",
     xlim = c(0, 8),
     ylim = c(0 , N.dfs+3),
     axes = FALSE,
     xlab = "",
     ylab = "")

# no effect
segments(x0 = 0, y0 = 2, x1 = 0, y1 = N.dfs+3, lty="dashed", lwd = 1.5)

# adding the betas
points(x = betas, y = seq(3.5,N.dfs+2.5), pch = 16 )
# 95% confidence intervals
segments(x0 = betas - qnorm(.975)* sigmas, y0 = seq(3.5,N.dfs+2.5), 
         x1 = betas + qnorm(.975)* sigmas, y1 = seq(3.5,N.dfs+2.5), lwd = 1.5)

# models
vars <- c("(1) BH weights; RP (SQ constraints & imputed SQs); PS (parametric)",
          "(2) BH weights; RP (SQ constraints & no imputed SQs); PS (parametric)",
          "(3) BH weights; RP (always midpoint); PS (parametric)",
          "(4) BH weights; control influence measure; PS (parametric)",
          "(5) SS weights; RP (SQ constraints & imputed SQs); PS (parametric)",
          "(6) SS weights; RP (SQ constraints & no imputed SQs); PS (parametric)",
          "(7) SS weights; RP (always midpoint); PS (parametric)",
          "(8) SS weights; control influence measure; PS (parametric)",
          "(9) BH weights; RP (SQ constraints & imputed SQs); PS (non-parametric)",
          "(10) BH weights; RP (SQ constraints & no imputed SQs); PS (non-parametric)",
          "(11) BH weights; RP (always midpoint); PS (non-parametric)",
          "(12) BH weights; control influence measure; PS (non-parametric)",
          "(13) SS weights; RP (SQ constraints & imputed SQs); PS (non-parametric)",
          "(14) SS weights; RP (SQ constraints & no imputed SQs); PS (non-parametric)",
          "(15) SS weights; RP (always midpoint); PS (non-parametric)",
          "(16) SS weights; control influence measure; PS (non-parametric)")

y.pos <- seq(3.5, length(vars)+2.5)
text(x = 2.5, y = y.pos, labels = vars, pos = 4)

# x axis ticks
axis (side = 1, at = seq(.1, 2.4, 0.2), pos = 2)
text(x = 1.2, y = .3, "Effect of Opaque EP Position")
dev.off()
embed_fonts("appendix_figure1.pdf")

rm(df.matched, transp.multi, betas, DVs, N.dfs, sigmas, vars, y.pos)