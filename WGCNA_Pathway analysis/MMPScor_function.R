# Functions for calculating the MMPS correlations
# output = dataframe MMPScorrelation

MMPScorrelation<-data.frame(Estimate=integer(),Pvalue=integer())
MMPScor<-function(module){
  est<- cor.test(abs(probeModuleMembership[moduleColors==module, match(paste0("ME",module), modNames) ]),
                 abs(probeTraitSignificance[moduleColors==module, 1]),
                 method="pearson",use="pairwise.complete.obs",exact=T)$estimate
  pval<- cor.test(abs(probeModuleMembership[moduleColors==module, match(paste0("ME",module), modNames) ]),
                  abs(probeTraitSignificance[moduleColors==module, 1]),
                  method="pearson",use="pairwise.complete.obs",exact=T)$p.value
  MMPScorrelation<<-rbind(MMPScorrelation,data.frame(est,pval))
  rownames(MMPScorrelation)[nrow(MMPScorrelation)]<<-paste0(module,"-",colnames(Trait))
  rm(est,pval)
}

MMPScorTop<-function(module){
  t<- probeModuleMembership[moduleColors==module, match(paste0("ME",module), modNames),drop=F]
  x<- t[abs(t[,1])>quantile(abs(t[,1]),0.85),,drop=F]
  est<- cor.test(abs(x[,1]),
                 abs(probeTraitSignificance[rownames(x),]),
                 method="pearson",use="pairwise.complete.obs",exact=T)$estimate
  pval<- cor.test(abs(x[,1]),
                  abs(probeTraitSignificance[rownames(x),]),
                  method="pearson",use="pairwise.complete.obs",exact=T)$p.value
  MMPScorrelation<<-rbind(MMPScorrelation,data.frame(est,pval))
  rownames(MMPScorrelation)[nrow(MMPScorrelation)]<<-paste0(module,"-",colnames(Trait),"_top15")
  rm(est,pval,t,x)
}