# Function for MM vs PS (or GS) plots
# Takes as input a character of the module colour of interest (e.g. "blue"), and can optionally specify colour for points (useful for e.g. module "lightblue")
# Pre-calculated data needed: probeTraitSignificance and probeModuleMembership, as well as a Trait variable (1 column with colname)
# Plot top 15% with highest MM (useful for modules with n probes >10.000)

MMPSplot<-function(module, col=module){
  
  column = match(paste0("ME",module), modNames)
  moduleGenes = moduleColors==module
  
  # Full plot of MM vs PS
  plot(abs(probeTraitSignificance[moduleGenes, 1])~abs(probeModuleMembership[moduleGenes, column]),
       xlab = paste("Module Membership in", module, "module"),
       ylab = paste("Probe Significance for", names(Trait)),cex.lab = 1.2, cex.axis = 1.2, col = col)
  if (cor.test(abs(probeTraitSignificance[moduleGenes, 1]),abs(probeModuleMembership[moduleGenes, column]),method="pearson",use="pairwise.complete.obs")$p.value==0){
    title(main = paste0("Module Membership vs. Probe Significance\n", 
                        "cor=",
                        round(cor.test(abs(probeTraitSignificance[moduleGenes, 1]),abs(probeModuleMembership[moduleGenes, column]),method="pearson",use="pairwise.complete.obs")$estimate,3),
                        ", p<1e-200"),cex.main = 1.2)
  } else {      
    title(main = paste0("Module Membership vs. Probe Significance\n", 
                        "cor=",
                        round(cor.test(abs(probeTraitSignificance[moduleGenes, 1]),abs(probeModuleMembership[moduleGenes, column]),method="pearson",use="pairwise.complete.obs")$estimate,3),
                        ", p=",
                        scientific(cor.test(abs(probeTraitSignificance[moduleGenes, 1]),abs(probeModuleMembership[moduleGenes, column]),method="pearson",use="pairwise.complete.obs")$p.value,digits=3))
          ,cex.main = 1.2)
  }  
  if (col=="black"){
    abline(lm(abs(probeTraitSignificance[moduleGenes, 1])~abs(probeModuleMembership[moduleGenes, column])),col="red")
  } else{
    abline(lm(abs(probeTraitSignificance[moduleGenes, 1])~abs(probeModuleMembership[moduleGenes, column])))}
 
  
  # Plot top 15% with highest MM (useful for modules with n probes >10.000)
  t<- probeModuleMembership[moduleGenes, column,drop=F]
  x<- t[abs(t[,1])>quantile(abs(t[,1]),0.85),,drop=F]
  
  plot(abs(probeTraitSignificance[rownames(x),])~abs(x[,1]),
       xlab = paste("Module Membership in", module, "module"),
       ylab = paste("Probe Significance for", names(Trait)),cex.lab = 1.2, cex.axis = 1.2, col = col)
  if (cor.test(abs(probeTraitSignificance[rownames(x),]),abs(x[,1]),method="pearson",use="pairwise.complete.obs")$p.value==0)
  {
    title(main = paste0("Module Membership vs. Probe Significance\n", 
                        "cor=",
                        round(cor.test(abs(probeTraitSignificance[rownames(x),]),abs(x[,1]),method="pearson",use="pairwise.complete.obs")$estimate,3),
                        ", p<1e-200"),cex.main = 1.2)
  } else {      
    title(main = paste0("Module Membership vs. Probe Significance\n", 
                        "cor=",
                        round(cor.test(abs(probeTraitSignificance[rownames(x),]),abs(x[,1]),method="pearson",use="pairwise.complete.obs")$estimate,3),
                        ", p=",
                        scientific(cor.test(abs(probeTraitSignificance[rownames(x),]),abs(x[,1]),method="pearson",use="pairwise.complete.obs")$p.value,digits=3)),
          cex.main = 1.2)
  } 
  if (col=="black"){
    abline(lm(abs(probeTraitSignificance[rownames(x),])~abs(x[,1])), col="red")
    } else{
  abline(lm(abs(probeTraitSignificance[rownames(x),])~abs(x[,1])))}
  
  rm(t,x)
}



