######## Linear regression of WGCNA modules & Conversion to AD ########
# calculating MM as correlation between corrected MEs and corrected data from EWAS

setwd("")
load("") # DNA methylation data of MCI sample subset, corrected for covariates
load("") # dat object used to generate modules
load("") # WGCNA modules generated for the conversion analysis
load("") # normalised DNA methylation data for the MCI sample subset, also contains pheno file
Anno<-read.csv("")
Anno<-Anno[match(colnames(dat), Anno$Name),]

data.dasen.regr<-RegData
rm(RegData)

library(WGCNA)
library(wateRmelon)
library(lm.beta)
options(stringsAsFactors = FALSE)

identical(as.character(rownames(dat)), as.character(rownames(pheno)))
identical(as.character(colnames(data.dasen.regr)), as.character(rownames(pheno)))
identical(as.character(colnames(data.dasen)), as.character(rownames(pheno)))


######## 1. Assign colours and calculate MEs ########
moduleColors = labels2colors(bwnet.fcx$colors)

MEs0= moduleEigengenes(dat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
rownames(MEs)<-rownames(dat)

# Disregard grey module from analysis,
# as it only contains probes that do not fit into any module
MEs<- MEs[,c(colnames(MEs)!="MEgrey")]
identical(rownames(pheno), rownames(MEs))


######## Regress out covariates ########
MET<-t(MEs)
identical(rownames(pheno), colnames(MET))

# Function to regress out covariates
RegrOut <- function( row, Age, Sex ,CD8T, CD4T, NK, Bcell, Mono,
                     Plate1 , Plate2 , Plate3,MMSE){
  
  residuals <- try (
    resid(lm( row ~ Age+ Sex + 
                CD8T+ CD4T+ NK+ Bcell+ Mono+ Plate1 + Plate2 + Plate3 +MMSE ,na.action=na.exclude)),
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
  MMSE <- as.numeric(pheno[colnames(MET), 'MMSE.MMSE_Total'])
  
  t( apply(MET, 1, RegrOut, Age, Sex ,CD8T, CD4T, NK, Bcell, Mono,
           Plate1 , Plate2 , Plate3,MMSE))
}
save(RegData,pheno, file="Regressed out data.RData")

######## Stripcharts: identify outliers ########
# Inspect values within modules for possible outliers

RegData.t<-t(RegData)
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

# Extreme outliers in:
# black, salmon, darkturquoise, cyan, orange, grey60, darkorange, lightyellow, saddlebrown, white, lightgreen
# steelblue, midnightblue, darkgreen,skyblue, purple,darkgrey,darkred

data.frame(matrix(nrow=5,ncol=1))->Overview
for (i in 1:ncol(MEs)){
  module=colnames(MEs)[i]
  
  RegData.t[RegData.t[,module] < mean(RegData.t[,module],na.rm=T)-5*sd(RegData.t[,module],na.rm=T),module]->tes
  cbind(c(rownames(RegData.t)[RegData.t[,module] %in% tes], rep("NA", length(rownames(Overview))-length(tes))))->Overview[, paste0(module,"Min_reg")]
  RegData.t[RegData.t[,module] > mean(RegData.t[,module],na.rm=T)+5*sd(RegData.t[,module],na.rm=T),module]->tes
  cbind(c(rownames(RegData.t)[RegData.t[,module] %in% tes], rep("NA", length(rownames(Overview))-length(tes))))->Overview[, paste0(module,"Max_reg")]
}
Overview[,1]<-NULL


# There are 16 extreme outliers in 16 modules (one in each)
# Look at pheno for these samples to see if they have anything in common
samples<-as.character(Overview[1,])[-which(as.character(Overview[1,])=="NA")]
pheno[samples, c(2,4:14,17,21,22,52,54,63:78,89)]->OutlierPheno
write.csv(OutlierPheno, "PhenoOutliers.csv")
# no clear commonalities; different genders, centres, chips, ages etc. all are MCI, but some conv/some non conv
# Set to NA in specific modules
RegData.t["8667045168_R06C01", "MEroyalblue"]<-NA
RegData.t["8667053124_R02C02", "MEdarkred"]<-NA
RegData.t["8667053073_R02C01", "MEsalmon"]<-NA
RegData.t["8622007175_R02C01", "MElightgreen"]<-NA
RegData.t["8667053165_R03C02", "MEsaddlebrown"]<-NA

RegData.t["8622007248_R05C01", "MEdarkturquoise"]<-NA
RegData.t["8667045168_R06C01", "MEgreenyellow"]<-NA
RegData.t["8622007102_R01C01", "MEblack"]<-NA
RegData.t["8667053148_R04C01", "MElightyellow"]<-NA
RegData.t["8667053163_R01C01", "MEpaleturquoise"]<-NA

RegData.t["8667053161_R03C02", "MEskyblue"]<-NA
RegData.t["8667053013_R03C02", "MEdarkgrey"]<-NA
RegData.t["8795207086_R02C02", "MEcyan"]<-NA
RegData.t["8667045168_R04C02", "MEdarkgreen"]<-NA
RegData.t["8622007175_R06C01", "MEgrey60"]<-NA
RegData.t["8667053163_R04C02", "MEsteelblue"]<-NA


pdf("Plots/Stripcharts_Reg_OutRem.pdf")
for (i in 1:ncol(RegData.t)){
  stripcharts(colnames(RegData.t)[i])
}
dev.off()

save(RegData.t, file="Regressed data_OutRem.RData")


# Linear regression of MEs####
t(RegData.t)->RegData
model <- function( row, Diagnosis ){
  
  fit <- try (
    lm( row ~ Diagnosis),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA,4))
  as.numeric(summary(fit)$coeff[,c(1,4)])
  if(length(as.numeric(c(summary(fit)$coeff[ , c(1, 4)]))) != 4) return(rep(NA, 4))
  as.numeric(c(summary(fit)$coeff[ , c(1, 4)]))
} 



### Apply model
model  <- {
  Diagnosis   <- as.numeric(pheno$Converters)
  
  
  t(apply(RegData, 1, model,  Diagnosis))
}


colnames(model)<-c("InterceptEst", "DiagnosisEst", "Intercept.p", "Diagnosis.p")

write.csv(model, "lm_MEs.csv")

# The orange module is significantly associated with future conversion to AD

######## MM/GS plot #################
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
t(MMPvalue)->MMPvalue

names(probeModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


### 3. MCI-MCI vs MCI-AD ####
# Significant module: orange 
Trait = as.data.frame(pheno$Converters);
names(Trait) = "MCI-MCIvsMCI-AD"

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
MMPSplot("orange")
dev.off()

# 3.3 Correlate MM to PS
MMPScor("orange")
MMPScorTop("orange")

write.csv(MMPScorrelation,"MMPScorrelation.csv")
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$PS.MCI.MCIvsMCI.AD));
geneInfo = geneInfo0[geneOrder, ]


save(geneInfo,geneInfo0, file="geneInfo files.RData")
write.csv(geneInfo, file="geneInfo.csv")


######### Run GO/KEGG pathway analysis for significant modules #######

dat<-t(dat)
names(moduleColors)<-rownames(dat)
allcpgs <- rownames(data.dasen) # create list of all probes

# Specify script version (prevent overwriting previous runs)
scriptVersion <- "-v1"
## Load required packages / functions
source("crystalmeth_fun.R")

#### Orange module pathway analysis ####
colour<- "orange"
sigcpgs <- names(which(moduleColors == colour))

# Specify project name to be appended to file names
projectName <- paste0(colour," module Converters")

### GO gene set analysis ====================================================###

## Run GO gene set analysis
go<- crystalmeth(sig.cpg = sigcpgs, 
                 all.cpg = allcpgs, array.type="450K")
head(go[ , -7], 10)
sum(go$FDR < 0.05) # 0

# Save GO analysis results
write.csv(go, file=paste0(colour, " GO Gene Analysis Output WGCNA_fullbg.csv"))
save(go, file = paste0(colour, " GO Gene Analysis Output WGCNA_fullbg.Rdata"))

### KEGG gene set analysis ==================================================###

kegg <- crystalmeth(sig.cpg = sigcpgs, 
                    all.cpg = allcpgs,
                    collection = "KEGG", array.type="450K")
head(kegg, 10)
sum(kegg$FDR < 0.05) # 1

# Save results
write.csv(kegg, file=paste0(colour, " KEGG Gene Analysis Output WGCNA_fullbg.csv"))
save(kegg, file = paste0(colour, " KEGG Gene Analysis Output WGCNA_fullbg.Rdata"))


