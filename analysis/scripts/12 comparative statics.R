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
  
  # normal adjustment
  ep.inf <- ifelse(normal==FALSE & ep.closer==FALSE, yes = 0, no = ep.inf)
  ep.inf <- ifelse(normal==FALSE & ep.closer==TRUE, yes = 1, no = ep.inf)
  
  
  # control adjustment for strange outcomes
  ep.inf <- ifelse(normal==FALSE, yes = ifelse(abs(pcou - out) + abs(pep - out) != 0,
                                               yes = abs(pcou - out) / (abs(pcou - out) + abs(pep - out)),
                                               no = 0.5),no = ep.inf)
  
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out <- seq(0,100,1)
cou <- rep(25, length(out))
ep <- rep(75, length(out))

# influence according to my measure
inf1 <- real.inf(pcou = cou,
                 pep = ep,
                 out = out,
                 cf = ((ep+cou)/2))

# influence according to the control measure
inf2 <- cont.inf(pcou = cou,
                 pep = ep,
                 out = out)

# comparative statics
# figure
setwd(subdirs$graphs)
pdf("appendix_figure2.pdf",  width = 10, height = 7)
plot(1,
     type = "n",
     xlim = c(0, 100),
     ylim = c(0, 1),
     axes = FALSE,
     xlab = "Position of the Outcome",
     ylab = "Relative Influence of the EP",
     cex.lab = 1.5)
lines(x = out, y = inf1, col = "gray25", lwd = 4)
lines(x = out, y = inf2, col = "gray50", lwd = 4, lty = "dashed")

# borders of well behaved scenarios
abline(v = ep[1], lty = "dashed", lwd = 1.5)
abline(v = cou[1], lty = "dashed", lwd = 1.5)

# split influence line
abline(h = .5, lty = "dashed", lwd = 1.5)

# axis
axis (side = 1, at = seq(0, 100, 5), cex = 1.5)
axis (side = 2, at = seq(0, 1, 0.2), cex = 1.5)

# position markers of the actors
text(x = ep[1], y = -0.01, "EP", pos = 4, cex = 1.5)
text(x = cou[1], y = -0.01, "Council", pos = 4, cex = 1.5)
text(x = 47.2, y = -.01, "Cf", pos = 4, cex = 1.5)

# legend
legend(x = 26, y = 1, legend = c("Main Measure", "Control Measure"), lty = c("solid", "dashed"),
       col = c("gray25", "gray50"), cex = 1.5, bty = "n",
       lwd = 3)

dev.off()

rm(cou, ep, inf1, inf2, out, cont.inf, real.inf)
