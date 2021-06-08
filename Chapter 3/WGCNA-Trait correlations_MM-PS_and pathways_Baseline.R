############ Module-Trait Correlations, MM/PS ############
# Calculating MM by correlating betas corrected for covariates (from EWAS) to corrected MEs

setwd("")
library(WGCNA)
library(wateRmelon)
library(lm.beta)
library(clusterProfiler)
library(KEGGREST)
options(stringsAsFactors = FALSE)

######## Load files ########

load("") # WGCNA modules = bwnet.fcx
load("") # probes used to create WGCNA modules (filtered set = dat)
load("") # full dasen normalised DNA methylation data, also contains phenotype file
Anno<-read.csv("")
Anno<-Anno[match(colnames(dat), Anno$Name),]

# Also load betas corrected for covariates
load("")
data.dasen.regr <- RegData
rm(RegData)

identical(rownames(pheno), rownames(dat))
identical(colnames(data.dasen.regr), rownames(dat))

######## 1. Assign colours and calculate MEs ########
moduleColors = labels2colors(bwnet.fcx$colors)

MEs0= moduleEigengenes(dat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
rownames(MEs)<-rownames(dat)

# Disregard grey module from analysis,
# as it only contains probes that do not fit into any module
MEs<- MEs[,c(colnames(MEs)!="MEgrey")]
identical(rownames(pheno), rownames(MEs))


######## 2. Create Traits data ####

# Create new variables for comparing baseline diagnosis groups
pheno$CTLvsAD<- ifelse(pheno$AtVisit_Group3 == 2, 1, ifelse(pheno$AtVisit_Group3==0,0,NA))
pheno$CTLvsMCI<- ifelse(pheno$AtVisit_Group3 == 1, 1, ifelse(pheno$AtVisit_Group3==0,0,NA))
pheno$MCIvsAD<- ifelse(pheno$AtVisit_Group3 == 2, 1, ifelse(pheno$AtVisit_Group3==1,0,NA))

# Select only traits of interest
Traits<- pheno[ , c(90:92,54,52,17,69:77,13:14,63:68,79:82)]
Traits$Sex<-as.character(Traits$Sex)
Traits$Sex<-as.numeric(Traits$Sex)

colnames(Traits)<-c("CTLvsAD","CTLvsMCI","MCIvsAD" ,"Education_Yrs", "APOE #4", "MMSE", 
                    "REV", "LEV", "TEV", "MET", "VV", "LHV",
                    "RHV", "THV", "WBV","Sex", "Age", "CD8T", "CD4T", "NK", "BCell", "Mono", "Gran",
                    "Batch1", "Batch2", "Batch3", "Batch4")

identical(as.character(rownames(Traits)), as.character(rownames(MEs)))

######## 3. Regress out data #### 
# Correct the MEs for covariates as well (seeing as they were based on uncorrected, normalised data)
MET<-t(MEs)
identical(rownames(pheno), colnames(MET))

# Function to regress out covariates
RegrOut <- function( row, Age, Sex ,CD8T, CD4T, NK, Bcell, Mono,
                     Plate1 , Plate2 , Plate3){
  
  residuals <- try (
    resid(lm( row ~ Age+ Sex + 
                CD8T+ CD4T+ NK+ Bcell+ Mono+ Plate1 + Plate2 + Plate3 ,na.action=na.exclude)),
    silent=TRUE
  )
  residuals
} 

### Apply model
RegData  <- {
  Age    <- as.numeric(pheno[ colnames(MET), 'Age_At_Visit'])
  Sex <- as.numeric(pheno[ colnames(MET), 'Sex'])
  CD8T <- as.numeric(pheno[ colnames(MET), 'CD8T'])
  CD4T <- as.numeric(pheno[ colnames(MET), 'CD4T'])
  NK <- as.numeric(pheno[ colnames(MET), 'NK'])
  Bcell <- as.numeric(pheno[ colnames(MET), 'Bcell'])
  Mono <- as.numeric(pheno[ colnames(MET), 'Mono'])
  Plate1 <-as.numeric(pheno[ colnames(MET), 'Plate1'])
  Plate2 <-as.numeric(pheno[ colnames(MET), 'Plate2'])
  Plate3 <-as.numeric(pheno[ colnames(MET), 'Plate3'])
  
  t( apply(MET, 1, RegrOut, Age, Sex ,CD8T, CD4T, NK, Bcell, Mono,
           Plate1 , Plate2 , Plate3))
}


as.data.frame(t(RegData))->RegData.t


######## 4. Stripcharts ########
# Inspect values within modules for possible outliers

stripcharts<-function(module){
  par(mfrow=c(1,1))
  stripchart(RegData.t[ ,module], method="jitter",main=paste0(module, " Module Regressed Data"))
  abline(v=c(mean(RegData.t[ ,module]) +2*sd(RegData.t[ ,module])), col="orange")
  abline(v=c(mean(RegData.t[ ,module]) -2*sd(RegData.t[ ,module])), col="orange")
  abline(v=c(mean(RegData.t[ ,module]) +5*sd(RegData.t[ ,module])), col="red")
  abline(v=c(mean(RegData.t[ ,module]) -5*sd(RegData.t[ ,module])), col="red")
  
}
pdf("Plots/Stripcharts_Reg.pdf")
for (i in 1:ncol(RegData.t)){
  stripcharts(colnames(RegData.t)[i])
}
dev.off()

### Adds names of outliers to an Overview table ####
data.frame(matrix(nrow=5,ncol=1))->Overview
for (i in 1:ncol(MEs)){
  module=colnames(MEs)[i]
  MEs[MEs[,module] < mean(MEs[,module],na.rm=T)-5*sd(MEs[,module],na.rm=T),module]->tes
  cbind(c(rownames(MEs)[MEs[,module] %in% tes], rep("NA", length(rownames(Overview))-length(tes))))->Overview[, paste0(module,"Min")]
  MEs[MEs[,module] > mean(MEs[,module],na.rm=T)+5*sd(MEs[,module],na.rm=T),module]->tes
  cbind(c(rownames(MEs)[MEs[,module] %in% tes], rep("NA", length(rownames(Overview))-length(tes))))->Overview[, paste0(module,"Max")]
  
  RegData.t[RegData.t[,module] < mean(RegData.t[,module],na.rm=T)-5*sd(RegData.t[,module],na.rm=T),module]->tes
  cbind(c(rownames(RegData.t)[RegData.t[,module] %in% tes], rep("NA", length(rownames(Overview))-length(tes))))->Overview[, paste0(module,"Min_reg")]
  RegData.t[RegData.t[,module] > mean(RegData.t[,module],na.rm=T)+5*sd(RegData.t[,module],na.rm=T),module]->tes
  cbind(c(rownames(RegData.t)[RegData.t[,module] %in% tes], rep("NA", length(rownames(Overview))-length(tes))))->Overview[, paste0(module,"Max_reg")]
}
Overview[,1]<-NULL

# There are three extreme outliers in lightcyan, purple and midnightblue (one in each)
# two are from the same sample
# Extreme outliers before and after regression are the same
samples<-c("8622007102_R01C01","8622007094_R05C01")
pheno[samples, c(2,4:14,17,21,22,52,54,63:78,89)]->OutlierPheno
write.csv(OutlierPheno, "PhenoOutliers.csv")
# 1 CTL, 1 MCI, do not see a clear connection between the two
# Set to NA in specific modules
RegData.t["8622007102_R01C01", "MElightcyan"]<-NA
RegData.t["8622007094_R05C01", "MEpurple"]<-NA

pdf("Plots/Stripcharts_Reg_OutRem.pdf")
for (i in 1:ncol(RegData.t)){
  stripcharts(colnames(RegData.t)[i])
}
dev.off()


######## 5. Module-Trait Correlations ########

### Using regressed out data ####
# Remove covariates from Traits
Traits<- Traits[ ,1:15]
Traits.continuous<- Traits[ ,c(4,6:15)]
Traits.ordinal<- Traits[ ,c(1:3,5)]

# Correlate, use spearman for ordinal variables, pearson for continuous
# pairwise correlations are carried out as the regressed ME data contains NAs
for (i in 1:ncol(Traits.ordinal)){
  Cor<- cor(RegData.t, Traits.ordinal[,i],method="spearman", use="pairwise.complete.obs")
  if(i==1){
    Cor->Correlations.ord
  } else{ 
    Correlations.ord<-cbind(Correlations.ord,Cor)
  }
}

colnames(Correlations.ord)<-colnames(Traits.ordinal)

for (i in 1:ncol(Traits.continuous)){
  Cor<- cor(RegData.t, Traits.continuous[,i],method="pearson", use="pairwise.complete.obs")
  if(i==1){
    Cor->Correlations.con
  } else{ 
    Correlations.con<-cbind(Correlations.con,Cor)
  }
}

colnames(Correlations.con)<-colnames(Traits.continuous)


for(o in 1:ncol(RegData.t)){
  for (i in 1:ncol(Traits.continuous)){
    x<-cor.test(RegData.t[,o], Traits.continuous[,i],method="pearson",use="pairwise.complete.obs")$p.value
    if(i==1){
      x->y
    } else{ 
      y<-cbind(y,x)
    }
  }
  if(o==1){
    CorrelationsPval.con<-y
  }else{
    CorrelationsPval.con<-rbind(CorrelationsPval.con,y)
  }
}

for(o in 1:ncol(RegData.t)){
  for (i in 1:ncol(Traits.ordinal)){
    x<-cor.test(RegData.t[,o], Traits.ordinal[,i],method="spearman",use="pairwise.complete.obs",exact=F)$p.value
    if(i==1){
      x->y
    } else{ 
      y<-cbind(y,x)
    }
  }
  if(o==1){
    CorrelationsPval.ord<-y
  }else{
    CorrelationsPval.ord<-rbind(CorrelationsPval.ord,y)
  }
}

colnames(CorrelationsPval.ord)<-colnames(Traits.ordinal)
rownames(CorrelationsPval.ord)<-colnames(RegData.t)

colnames(CorrelationsPval.con)<-colnames(Traits.continuous)
rownames(CorrelationsPval.con)<-colnames(RegData.t)


cbind(Correlations.ord,Correlations.con)->Correlations
cbind(CorrelationsPval.ord,CorrelationsPval.con)->CorrelationsPval

# save results for tables
colnames(CorrelationsPval)<-paste0(colnames(CorrelationsPval),".p")
CorrelationsAll<- cbind(Correlations,CorrelationsPval)
write.csv(CorrelationsAll, "Correlations_Reg.csv")


pdf("Plots/Heatmap_Traits_Regressed out.pdf")
textMatrix = paste(signif(Correlations, 2), "\n(",
                   signif(CorrelationsPval, 1), ")", sep = "");
dim(textMatrix) = dim(Correlations)
par(mar = c(5,5,3,1));


labeledHeatmap(Matrix = Correlations,
               xLabels = colnames(Correlations),
               yLabels = colnames(RegData.t),
               ySymbols = colnames(RegData.t),
               colorLabels = T,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

######## Calculate and correlate MM/PS for modules associated to traits ####
# MM: correlate probe methylation to regressed out MEs (outliers removed), = pearson
# PS: correlate probe methylation to trait of interest, = pearson for continuous, spearman for ordinal
# Then do pearson correlations between MM and PS
# For modules with more than 10.000 probes, look at top 15% of probes in module (based on MM)

### 1. Load Functions ####
# 1.1 MMPS plots function
# the function verboseScatterplot (for MMPS plots) is supposed to always calculate the pearson correlations,
# but I noticed it actually showed spearman correlations. Have created new functions to create MMPS (full and top 15% MM) plots
source("MMPSplot_function.R")

# 1.2. MMPS correlations function
source("MMPScor_function.R")

### 2. Calculate Module Memberships ####
# need to use regressed out data for probes

dat.original <- dat
data.dasen.regr <- t(data.dasen.regr)
dat <- data.dasen.regr[,colnames(dat.original)]

modNames = rownames(RegData)
probeModuleMembership = as.data.frame(cor(dat, RegData.t, method="pearson",use = "pairwise.complete.obs"))
# Wrote a new loop to calculate MM p-value, as original code corrected for total samples, not pairwise
 for(o in 1:ncol(RegData.t)){
   for (i in 1:ncol(dat)){
     x<-cor.test(RegData.t[,o], dat[,i],method="pearson",use="pairwise.complete.obs")$p.value
     if(i==1){
       x->Pval.MM
     } else{
       Pval.MM<-cbind(Pval.MM,x)
     }
   }
   if(o==1){
     MMPvalue<-Pval.MM
   }else{
     MMPvalue<-rbind(MMPvalue,Pval.MM)
   }
 }
 colnames(MMPvalue)<-colnames(dat)
 rownames(MMPvalue)<-colnames(RegData.t)
 MMPvalue<- t(MMPvalue)
 
names(probeModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


### 3. CTLvsMCI ####
# Significant modules: lightcyan, yellow, purple, brown 
Trait = as.data.frame(pheno$CTLvsMCI);
names(Trait) = "CTLvsMCI"

# Calculate Probe Significance (correlation between the corrected DNA methylation data and the trait of interest
probeTraitSignificance = as.data.frame(cor(dat, Trait, method="spearman", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="spearman",
                                                           use="pairwise.complete.obs",exact=F)$p.value))
colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="");


# 3.1 Create Supplementary Table 
probes = colnames(dat)
probes2annot = match(probes, Anno$IlmnID)
# Create the starting data frame
geneInfo0 = data.frame(ProbeID = probes,
                       UCSCGene = Anno$UCSC_RefGene_Name[probes2annot],
                       moduleColor = moduleColors,
                       probeTraitSignificance,
                       PSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(RegData.t, Trait, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(probeModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, probeModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# 3.2 Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("lightcyan", col="blue")
MMPSplot("yellow", "goldenrod4")
MMPSplot("purple")
MMPSplot("brown")
dev.off()

# 3.3 Correlate MM to PS
MMPScor("lightcyan")
MMPScor("yellow")
MMPScorTop("yellow")
MMPScor("purple")
MMPScor("brown")
MMPScorTop("brown")


### 4. APOE ####
# Significant module: cyan
Trait = as.data.frame(pheno$APOE_Number4);
names(Trait) = "APOE#4"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="spearman", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="spearman",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 4.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 4.2. Add MM-PS correlations to table
MMPScor(module="cyan")


# 4.3. Plot 
pdf(paste0("Plots/MM-PS plot ",names(Trait),".pdf"))
MMPSplot("cyan","darkcyan")
dev.off()


### 5. Edu.yrs ####
# Significant module: brown
Trait = as.data.frame(pheno$Demographics.Fulltime_Education_Years);
names(Trait) = "Edu.yrs"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="spearman", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="spearman",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 5.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 5.2. Add MM-PS correlations to table
MMPScor(module="brown")
MMPScorTop(module="brown")

# 5.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("brown")
dev.off()

### 6. REV ####
# Sign: purple
Trait = as.data.frame(pheno$REV);
names(Trait) = "REV"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="pearson", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="pearson",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 6.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 6.2. Add MM-PS correlations to table
MMPScor(module="purple")

# 6.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("purple")
dev.off()

### 7. TEV ####
# Sign: brown
Trait = as.data.frame(pheno$TEV);
names(Trait) = "TEV"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="pearson", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="pearson",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 7.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 7.2. Add MM-PS correlations to table
MMPScor(module="purple")

# 7.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("purple")
dev.off()

### 8. MET ####
# Sign: purple,yellow
Trait = as.data.frame(pheno$MET);
names(Trait) = "MET"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="pearson", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="pearson",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 8.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 8.2. Add MM-PS correlations to table
MMPScor("purple")
MMPScor("yellow")
MMPScorTop("yellow")

# 8.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("purple")
MMPSplot("yellow")
dev.off()

### 9. WBV ####
# Sign: purple
Trait = as.data.frame(pheno$WBV);
names(Trait) = "WBV"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="pearson", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="pearson",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 9.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 9.2. Add MM-PS correlations to table
MMPScor("purple")

# 9.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("purple")
dev.off()

### 10. VV ####
# Sign: purple
Trait = as.data.frame(pheno$VV);
names(Trait) = "VV"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="pearson", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="pearson",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 10.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 10.2. Add MM-PS correlations to table
MMPScor("purple")

# 10.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("purple")
dev.off()

### 11. LHV ####
# Sign: purple
Trait = as.data.frame(pheno$LHV);
names(Trait) = "LHV"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="pearson", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="pearson",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

#1.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 11.2. Add MM-PS correlations to table
MMPScor("purple")

# 11.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("purple")
dev.off()

### 12. THV ####
# Sign: purple
Trait = as.data.frame(pheno$THV);
names(Trait) = "THV"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="pearson", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="pearson",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 12.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 12.2. Add MM-PS correlations to table
MMPScor("purple")

# 12.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("purple")
dev.off()

### 13. RHV ####
# Sign: purple
Trait = as.data.frame(pheno$RHV);
names(Trait) = "RHV"

probeTraitSignificance = as.data.frame(cor(dat, Trait, method="pearson", use = "pairwise.complete.obs"));
PSPvalue <-as.data.frame(apply(dat,2, function(r) cor.test(r,Trait[,1],method="pearson",
                                                           use="pairwise.complete.obs",exact=F)$p.value))

colnames(PSPvalue)<-colnames(Trait)
names(probeTraitSignificance) = paste("PS.", names(Trait), sep="");
names(PSPvalue) = paste("p.PS.", names(Trait), sep="")

# 12.1. Add PS values to table
geneInfo0<- merge(geneInfo0, probeTraitSignificance,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID
geneInfo0<- merge(geneInfo0, PSPvalue,by="row.names")
rownames(geneInfo0)<-geneInfo0$ProbeID

# 12.2. Add MM-PS correlations to table
MMPScor("purple")

# 12.3. Plot 
pdf(paste0("Plots/MM-PS plot ", names(Trait),".pdf"))
MMPSplot("purple")
dev.off()


##### Save Supplementary file ####
# Order the genes in the geneInfo variable first by module color, then by probeTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$PS.CTLvsMCI));
geneInfo = geneInfo0[geneOrder, ]
geneInfo[ ,1:20]<-NULL

save(geneInfo,geneInfo0, file="geneInfo files.RData")
write.csv(geneInfo, file="geneInfo.csv")

write.csv(MMPScorrelation,"MMPScorrelations.csv")


######### Run GO/KEGG pathway analysis for significant modules #######

dat<-t(dat)
names(moduleColors)<-rownames(dat)

# Specify script version (prevent overwriting previous runs)
scriptVersion <- "-v1"
## Load required packages / functions
source("crystalmeth_fun.R")

allcpPS <- rownames(data.dasen) # list of all probes

#### Yellow module ####
colour<- "yellow"
# Specify project name to be appended to file names
projectName <- paste0(colour," module CvM_MET")

sigcpPS <- names(which(moduleColors == colour)) # list of probes in the module


### GO gene set analysis ====================================================###

## Run GO gene set analysis
go<- crystalmeth(sig.cpg = sigcpPS, 
                 all.cpg = allcpPS,
                 array.type = "450K")
head(go[ , -7], 10)
sum(go$FDR < 0.05) # 0

# Save GO analysis results

write.csv(go, file=paste0("Pathways/",colour, " GO Gene Analysis Output WGCNA_fullbg.csv"))
save(go, file = paste0("Pathways/",colour, " GO Gene Analysis Output WGCNA_fullbg.Rdata"))

### KEGG gene set analysis ==================================================###

kegg <- crystalmeth(sig.cpg = sigcpPS, 
                    all.cpg = allcpPS,
                    collection = "KEGG",
                    array.type = "450K")
head(kegg, 10)
sum(kegg$FDR < 0.05) # 1

# Save results

write.csv(kegg, file=paste0("Pathways/",colour, " KEGG Gene Analysis Output WGCNA_fullbg.csv"))
save(kegg, file = paste0("Pathways/",colour, " KEGG Gene Analysis Output WGCNA_fullbg.Rdata"))


#### purple module ####
colour<- "purple"

# Specify project name to be appended to file names
projectName <- paste0(colour," module CvM")

sigcpPS <- names(which(moduleColors == colour)) # list of probes in the module

### GO gene set analysis ====================================================###

## Run GO gene set analysis
go<- crystalmeth(sig.cpg = sigcpPS, 
                 all.cpg = allcpPS,
                 array.type = "450K")
head(go[ , -7], 10)
sum(go$FDR < 0.05) # 0

# Save GO analysis results
write.csv(go, file=paste0("Pathways/",colour, " GO Gene Analysis Output WGCNA_fullbg.csv"))
save(go, file = paste0("Pathways/",colour, " GO Gene Analysis Output WGCNA_fullbg.Rdata"))

### KEGG gene set analysis ==================================================###

kegg <- crystalmeth(sig.cpg = sigcpPS, 
                    all.cpg = allcpPS,
                    collection = "KEGG",
                    array.type = "450K")
head(kegg, 10)
sum(kegg$FDR < 0.05) # 1

# Save results
write.csv(kegg, file=paste0("Pathways/",colour, " KEGG Gene Analysis Output WGCNA_fullbg.csv"))
save(kegg, file = paste0("Pathways/",colour, " KEGG Gene Analysis Output WGCNA_fullbg.Rdata"))



#### brown module ####
colour<- "brown"

# Specify project name to be appended to file names
projectName <- paste0(colour," module CvM_Edu.yrs")

sigcpPS <- names(which(moduleColors == colour)) # list of probes in the module


### GO gene set analysis ====================================================###

## Run GO gene set analysis
go<- crystalmeth(sig.cpg = sigcpPS, 
                 all.cpg = allcpPS,
                 array.type = "450K")
head(go[ , -7], 10)
sum(go$FDR < 0.05) # 0

# Save GO analysis results
write.csv(go, file=paste0("Pathways/",colour, " GO Gene Analysis Output WGCNA_fullbg.csv"))
save(go, file = paste0("Pathways/",colour, " GO Gene Analysis Output WGCNA_fullbg.Rdata"))

### KEGG gene set analysis ==================================================###

kegg <- crystalmeth(sig.cpg = sigcpPS, 
                    all.cpg = allcpPS,
                    collection = "KEGG",
                    array.type = "450K")
head(kegg, 10)
sum(kegg$FDR < 0.05) # 1

# Save results
write.csv(kegg, file=paste0("Pathways/",colour, " KEGG Gene Analysis Output WGCNA_fullbg.csv"))
save(kegg, file = paste0("Pathways/",colour, " KEGG Gene Analysis Output WGCNA_fullbg.Rdata"))



#### cyan module ####
colour<- "cyan"

# Specify project name to be appended to file names
projectName <- paste0(colour," module APOE#4")

sigcpPS <- names(which(moduleColors == colour)) # list of probes in the module

### GO gene set analysis ====================================================###

## Run GO gene set analysis
go<- crystalmeth(sig.cpg = sigcpPS, 
                 all.cpg = allcpPS,
                 array.type = "450K")
head(go[ , -7], 10)
sum(go$FDR < 0.05) # 0

# Save GO analysis results
write.csv(go, file=paste0("Pathways/",colour, " GO Gene Analysis Output WGCNA_fullbg.csv"))
save(go, file = paste0("Pathways/",colour, " GO Gene Analysis Output WGCNA_fullbg.Rdata"))

### KEGG gene set analysis ==================================================###

kegg <- crystalmeth(sig.cpg = sigcpPS, 
                    all.cpg = allcpPS,
                    collection = "KEGG",
                    array.type = "450K")
head(kegg, 10)
sum(kegg$FDR < 0.05) # 1

# Save results
write.csv(kegg, file=paste0("Pathways/",colour, " KEGG Gene Analysis Output WGCNA_fullbg.csv"))
save(kegg, file = paste0("Pathways/",colour, " KEGG Gene Analysis Output WGCNA_fullbg.Rdata"))



#### Yellow core probes ####
# extract probes with top 15% MM as this is a large module
colour="yellow"

moduleGenes = moduleColors==colour
geneInfo<-geneInfo[order(rownames(geneInfo)),]
identical(rownames(geneInfo), names(moduleColors))
t<- geneInfo[moduleGenes,"MM.MEyellow",drop=F]
x<- t[abs(t[,1])>quantile(abs(t[,1]),0.85),,drop=F]

sigcpPS <- rownames(x) # probes with top 15% highest MM

# Specify project name to be appended to file names
projectName <- paste0(colour," module CvM_MET_core probes")


### GO gene set analysis ====================================================###

## Run GO gene set analysis
go<- crystalmeth(sig.cpg = sigcpPS, 
                 all.cpg = allcpPS,
                 array.type = "450K")
head(go[ , -7], 10)
sum(go$FDR < 0.05) # 0

# Save GO analysis results
write.csv(go, file=paste0("Pathways/Core",colour, " Core GO Gene Analysis Output WGCNA_fullbg.csv"))
save(go, file = paste0("Pathways/Core",colour, " Core GO Gene Analysis Output WGCNA_fullbg.Rdata"))

### KEGG gene set analysis ==================================================###

kegg <- crystalmeth(sig.cpg = sigcpPS, 
                    all.cpg = allcpPS,
                    collection = "KEGG",
                    array.type = "450K")
head(kegg, 10)
sum(kegg$FDR < 0.05) # 1

# Save results
write.csv(kegg, file=paste0("Pathways/Core",colour, "Core KEGG Gene Analysis Output WGCNA_fullbg.csv"))
save(kegg, file = paste0("Pathways/Core",colour, "Core KEGG Gene Analysis Output WGCNA_fullbg.Rdata"))


#### Brown core probes ####
# extract probes with top 15% MM
colour="brown"

moduleGenes = moduleColors==colour
geneInfo<-geneInfo[order(rownames(geneInfo)),]
identical(rownames(geneInfo), names(moduleColors))
t<- geneInfo[moduleGenes,"MM.MEbrown",drop=F]
x<- t[abs(t[,1])>quantile(abs(t[,1]),0.85),,drop=F]

sigcpPS <- rownames(x)

# Specify project name to be appended to file names
projectName <- paste0(colour," module CvM_Edu.yrs_core probes")


### GO gene set analysis ====================================================###

## Run GO gene set analysis
go<- crystalmeth(sig.cpg = sigcpPS, 
                 all.cpg = allcpPS,
                 array.type = "450K")
head(go[ , -7], 10)
sum(go$FDR < 0.05) # 0

# Save GO analysis results
write.csv(go, file=paste0("Pathways/Core",colour, " Core GO Gene Analysis Output WGCNA_fullbg.csv"))
save(go, file = paste0("Pathways/Core",colour, " Core GO Gene Analysis Output WGCNA_fullbg.Rdata"))

### KEGG gene set analysis ==================================================###

kegg <- crystalmeth(sig.cpg = sigcpPS, 
                    all.cpg = allcpPS,
                    collection = "KEGG",
                    array.type = "450K")
head(kegg, 10)
sum(kegg$FDR < 0.05) # 1

# Save results
write.csv(kegg, file=paste0("Pathways/Core",colour, " Core KEGG Gene Analysis Output WGCNA_fullbg.csv"))
save(kegg, file = paste0("Pathways/Core",colour, " Core KEGG Gene Analysis Output WGCNA_fullbg.Rdata"))

