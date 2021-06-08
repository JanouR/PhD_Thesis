###################### QC and Normalisation 450K data ###################### 
# AddNeuroMed cohort

wd<-""
setwd(wd)
library(methylumi)
library(wateRmelon)


###################### load data  ######################
# data required: 450K IDATS, phenotype data, probe annotation information,
# and cross-hybridising/SNP probes lists


idatPath1<-""
idatPath2<-""
 
batch1<-readEPIC(idatPath1, barcodes=NULL, pdat=NULL,parallel=F,n=T,oob=F,force=F)
batch2<-readEPIC(idatPath2, barcodes=NULL, pdat=NULL,parallel=F,n=T,oob=F,force=F)
data<-combo(batch1, batch2)
 
dim(data)
#Features  Samples
# 485577      312

Anno<-read.csv("")
pheno<-read.csv("")
cross<- read.table("")

if(!dir.exists("Plots")){
  dir.create("Plots")
}

rm(idatPath1,idatPath2,batch1,batch2)

###################### Recoding pheno variables ######################

pheno$EU.code.1 <-NULL

### correct CDR variable (factor to numeric)

as.character(pheno$CDR.sum.of.boxes)->pheno$CDRsum
as.numeric(pheno$CDRsum)->pheno$CDRsum
pheno$CDR.sum.of.boxes<-NULL

### Create dummy variables for later analyses

pheno$Plate1<- ifelse(pheno$Plate == "Plate KATIE 1", 1, 0)
pheno$Plate2<- ifelse(pheno$Plate == "Plate KATIE 2", 1, 0)
pheno$Plate3<- ifelse(pheno$Plate == "Plate KATIE 3", 1, 0)
pheno$Plate4<- ifelse(pheno$Plate == "Plate KATIE 4", 1, 0)

pheno$Kuopio<- ifelse(pheno$Source == 0, 1, 0)
pheno$Lodz<- ifelse(pheno$Source == 1, 1, 0)
pheno$London<- ifelse(pheno$Source == 2, 1, 0)
pheno$Perugia<- ifelse(pheno$Source == 3, 1, 0)
pheno$Thessaloniki<- ifelse(pheno$Source == 4, 1, 0)
pheno$Toulouse<- ifelse(pheno$Source == 5, 1, 0)

###################### Add Houseman cell type proportions ######################

cells<-read.csv("BetasforClock_Filtered.output.csv")
cells$Sentrix<-gsub("X", "", cells$SampleID )
rownames(cells)<-cells$Sentrix
cells<- cells[, c(109, 67:72, 2,56,57)]
identical(rownames(pheno), rownames(cells))

pheno<-cbind(pheno,cells)
identical(pheno[,4], pheno[,126])
pheno[,126]<-NULL


###################### Data QC  ######################


###### QC 1. Calculate bisulfite conversion efficiency for BS only ######
# There are several fully methylated control probes included which can be used to generate a score
# to estimate the success of the bisulfite conversion step. As they are fully methylated these should have DNA methylation
# values of ~1. The bisulfite conversion score is essentially the median of these probes, and value < 80 is taken as a failure. 

bs<-bscon(data)
summary(bs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 84.22   90.19   91.44   91.15   92.29   94.59 

pdf("Plots/Plot_BSconversion.pdf")
hist(bs, xlab = "Median % BS conversion", main = "", col="steelblue")
dev.off()

rm(bs)

###### QC 2. Remove Cross-hybridising and SNP probes ######

crosslist <- cross[,1]
length(crosslist) #72067
data1 <-data[ ! rownames(data) %in% crosslist, ]
dim(data1) #  413510      307

identical(as.character(rownames(pheno)),as.character(colnames(data1)))

rm(crosslist, cross)

###### QC 3. Multidimensional Scaling of Sex Chomosomes ######

# this is a check that gender is correct (i.e. for sample mix up), 
# by comparing phenotype sex to profiles on X and Y chromosome

betas <- betas(data1)
Anno <- Anno[match(rownames(betas), Anno$Name),]

pheno$Sex<- as.factor(pheno$Sex)
l.bs<-pheno$Sadman.Code

### Testing the sex labels and profile based on bisulfite data
betas.x<-betas[as.character(Anno[which(Anno$CHR == "X"),]$IlmnID),]
fit.x<-data.frame(cmdscale(dist(t(betas.x)),2))
betas.y<-betas[as.character(Anno[which(Anno$CHR == "Y"),]$IlmnID),]
fit.y<-data.frame(cmdscale(dist(t(betas.y)),2))

sex_palette<-cbind(c("0", "1"), c("darkred", "darkcyan"))
sex_col<-sex_palette[match(pheno$Sex, sex_palette[,1]),2]

### Create MDS plot 
pdf("Plots/SexChr_gendercheck_noLab.pdf")

par(mfrow = c(1,2), xpd=NA, mar = par()$mar + c(2, 0, 0, 0))
plot(fit.x[,1], fit.x[,2], xlab = "MDS Coordinate 1", ylab = "MDS Coordinate 2", pch = 16, col = sex_col, main = "X Chromosome")
legend("bottom",inset=c(0,-0.255), horiz = TRUE, c("M", "F"), col = c("darkred", "darkcyan"), pch = 16)

plot(fit.y[,1], fit.y[,2], xlab = "MDS Coordinate 1", ylab = "MDS Coordinate 2", pch = 16, col = sex_col, main = "Y Chromosome")
legend("bottom", inset=c(0,-0.255),horiz = TRUE, c("M", "F"), col = c("darkred", "darkcyan"), pch = 16)

dev.off()

rm(betas.x, betas.y, fit.x,fit.y, sex_col,sex_palette,l.bs)

###### QC 4. Check genetically identical samples correlate across SNP probes #######

# This is a useful step which compares all pairs of samples across the 59 SNP probes and has two opposing goals 
# depending on your study design. Either you are checking that matched samples are genetically identical 
# (E.g. longitudinal study, twin design etc.) or you have an unrelated sample and you are checking that your samples are
# in fact not genetically identical. In the latter case, genetically identical smaples or duplicate samples may indicating
# a sample collection error or a experimental processing error and should be investigated. Note this approach is only able to 
# identify samples that are 100% genetically identical and no lower proportion of genetic relatedness (i.e. siblings).
#
# First, produce a histogram of the maximum correlation for each sample with all other samples (except for itself). 
# Any samples for which this is > 0.95 is suggestive that they are not genetically unique. 

# Use the data object that still contains SNP probes
betas <- betas(data)
pheno1<-pheno[match(colnames(betas), pheno$SentrixID),]

betas.rs<-betas[grep("rs", rownames(betas)),]
snpCor<-cor(betas.rs)
names(snpCor)<-pheno$SentrixID
for(i in 1:ncol(betas.rs)){
  snpCor[i,i]<-NA
}
corMax<-apply(snpCor, 1, max, na.rm = TRUE) ## calculates the maximum correlation for each sample with all other samples (except for itself)

pdf("Plots/SNPCorrelations.pdf")
hist(corMax, xlab = "Max. correlation with all other samples", main = "", col="steelblue")
dev.off()

###### QC 5. Log intensity boxplots of methylated and unmethylated raw data values #######

### extract sample intensities and plot
m_intensities<-methylated(data1)
u_intensities<-unmethylated(data1)
betas<-betas(data1)

M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

pdf("Plots/Histograms+Scatterplot_SampleIntensities_raw.pdf")
par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity",col="steelblue")
hist(U.median, xlab = "Median U intensity",col="steelblue")
par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity",ylim=c(1500*1.04, 4500*1.04),col="steelblue")
dev.off()


### Plot all samples as individual lines, + density plot of averages
pdf("Plots/Raw betas+averages.pdf")
par(mfrow=c(2,1))
densityPlot(betas, main = "Raw Betas")
plot(density(betas,na.rm=T), ylim=c(0,4), col='black', main = "Averages of Raw Data")
dev.off()

rm(m_intensities, u_intensities, M.median, U.median,colours)


###################### Pfilter processed beta values ######################

data.pf<- pfilter(data1)
# 1 samples having 1 % of sites with a detection p-value greater than 0.05 were removed 
# Samples removed: 8667045168_R04C01 
# 201 sites were removed as beadcount <3 in 5 % of samples 
# 2426 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 

### remove pfilter removed sample from pheno
rownames(pheno)-> samples
pheno.pf<- pheno[samples!= "8667045168_R04C01", ]
identical(as.character(rownames(pheno.pf)), as.character(colnames(data.pf)))
# [1] TRUE
rm(keep,samples)


###################### Remove outliers ######################
# Apply function outlix to detect outliers

outliers <- outlyx(betas)
# No samples are labelled as outliers

###################### Dasen normalisation ######################

data.dasen <-dasen(data.pf)

###### Intensity plots normalised data #####
# Extract sample intensities and plot
m_intensities<-methylated(data.dasen)
u_intensities<-unmethylated(data.dasen)
betasd<-betas(data.dasen)

M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

pdf("Plots/Histograms+Scatterplot_SampleIntensities_Normalised.pdf")
par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity",col="steelblue")
hist(U.median, xlab = "Median U intensity",col="steelblue")
par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity",col="steelblue")
dev.off()


# Median by Chip
pdf("Plots/Boxplot Median Intensities by Chip_Normalised.pdf")
par(mfrow = c(1,2))
par(mar = c(7, 5, 1, 1))
nCol<-length(unique(pheno.pf$Chip))
boxplot(M.median ~ pheno.pf$Chip, ann=F, xlab = "Chip", las = 2, col = rainbow(nCol), xaxt="n", cex.axis=0.9) 
title(ylab="Median M intensity", line=3.8)
title(xlab="Chip", line=5)
axis(side=1, las=2,labels=unique(pheno.pf$Chip),at=c(1:26), cex.axis=0.7)
boxplot(U.median ~ pheno.pf$Chip, ann=F, las = 2, col = rainbow(nCol), xaxt="n", cex.axis=0.9)
title(ylab="Median U intensity", line=3.8)
title(xlab="Chip", line=5)
axis(side=1, las=2,labels=unique(pheno.pf$Chip),at=c(1:26), cex.axis=0.7)
dev.off()

### Plot by plate
pdf("Plots/Boxplot Median Intensities by Plate_Normalised.pdf")
par(mfrow = c(1,2))
par(mar = c(5, 5, 1, 1))
nCol<-length(unique(pheno.pf$Plate))
boxplot(M.median ~ pheno.pf$Plate,ann=F, las = 2, col = rainbow(nCol), xaxt="n") 
title(ylab="Median M intensity", line=3.8)
title(xlab="Plate", line=3)
axis(side=1, las=0)
boxplot(U.median ~ pheno.pf$Plate, ann=F, las = 2, col = rainbow(nCol), xaxt="n")
title(ylab="Median U intensity", line=3.8)
title(xlab="Plate", line=3)
axis(side=1, las=0)
dev.off()

###### Density plots for pfiltered and normalised data, for all probes and per probe type #####

# Calculate betas of pfiltered data
betas.pf   <- betas(data.pf)
betas.pfI  <- betas(data.pf[which(data.pf@featureData$DESIGN == "I"),])
betas.pfII <- betas(data.pf[which(data.pf@featureData$DESIGN == "II"),])

# Calculate betas of dasen data
betas.dasen   <- betas(data.dasen)
betas.dasenI  <- betas(data.dasen[which(data.dasen@featureData$DESIGN == "I"),])
betas.dasenII <- betas(data.dasen[which(data.dasen@featureData$DESIGN == "II"),])


# Generate density plots of pfiltered betas and normalised betas
pdf("Plots/Data.density_Average_p-filtered+normalised.pdf")

plot(density(betas.pf, na.rm = TRUE), main = "P-filtered Data", 
     xlim = c(0, 1), ylim = c(0, 12), col = "black", lty = 1)
lines(density(betas.pfI, na.rm = TRUE), main = "P-filtered Data", 
      xlim = c(0, 1), ylim = c(0, 12), col = "blue", lty = 1)
lines(density(betas.pfII, na.rm = TRUE), main = "P-filtered Data", 
      xlim = c(0, 1), ylim = c(0, 12), col = "red", lty = 1)
legend("topright", legend = c("Combined", "Type I", "Type II"), lty = 1, lwd = 2,
       col = c("black", "blue", "red"), merge = TRUE, bty = "n")

plot(density(betas.dasen, na.rm = TRUE), main = "Dasen Normalised Data", 
     xlim = c(0, 1), ylim = c(0, 12), col = "black", lty = 1)
lines(density(betas.dasenI, na.rm = TRUE), main = "Dasen Normalised Data", 
      xlim = c(0, 1), ylim = c(0, 12), col = "blue", lty = 1)
lines(density(betas.dasenII, na.rm = TRUE), main = "Dasen Normalised Data", 
      xlim = c(0, 1), ylim = c(0, 12), col = "red", lty = 1)
legend("topright", legend = c("Combined", "Type I", "Type II"), lty = 1, lwd = 2,
       col = c("black", "blue", "red"), merge = TRUE, bty = "n")

dev.off()

###### DMRSE, Genki  ######

### Calculate DMRSE of normal preprocessed data.
dasendmrse.dasen<- dmrse(betas.dasen)
dasendmrse.dasenI <- dmrse(betas.dasenI)
dasendmrse.dasenII <- dmrse(betas.dasenII)
print(dasendmrse.dasen)
# [1] 0.007437589
print(dasendmrse.dasenI)
# [1] 0.007335743
print(dasendmrse.dasenII)
# [1] 0.007300075

rm(dasendmrse.dasen,dasendmrse.dasenI,dasendmrse.dasenII)

### Genki
# Uses 65 SNP probes that were previously removed from normalised data, still included in data object
pfilter(data)->d
dasen(d)->dd

# Calculate betas of dasen data
betas.dd   <- betas(dd)
betas.ddI  <- betas(dd[which(dd@featureData$DESIGN == "I"),])
betas.ddII <- betas(dd[which(dd@featureData$DESIGN == "II"),])
genki.dasen <- genki(betas.dd)
genki.dasenI <- genki(betas.ddI)
genki.dasenII <- genki(betas.ddII)

print(genki.dasen)
#[1] 3.473558e-05 4.206714e-05 2.057160e-05
print(genki.dasenI)
#[1] 5.095002e-05 4.534114e-05 2.957786e-05
print(genki.dasenII)
#[1] 2.487705e-05 4.000249e-05 1.486633e-05


rm(betas.dd,betas.ddI,betas.ddII,genki.dasen,genki.dasenII,genki.dasenI,d,dd)

###################### Calculate smoking scores ######################

smokingScore<-function(betas){
  
  load("SmokingScoreRefData.rda")
  #contains Illig_data, Illig_data_up and Illig_data_down
  
  #subsetting own data 
  #SABRE_data_down
  select <- rownames(betas) %in% Illig_data_down$cpgs
  A_down <- subset(betas, select =="TRUE")
  
  #SABRE_data_up
  select <- rownames(betas) %in% Illig_data_up$cpgs
  A_up <- subset(betas, select =="TRUE")
  
  #sort SABRE data by Cpg name
  A_up <- A_up[order(rownames(A_up)),]
  A_down <- A_down[order(rownames(A_down)),]
  
  #match Illig data by by Cpg name
  Illig_data_up<-Illig_data_up[match(rownames(A_up), Illig_data_up$cpgs),]
  Illig_data_down<-Illig_data_down[match(rownames(A_down), Illig_data_down$cpgs),]
  
  #as some outliers have been removed and replaced with NAs need to handle missing values
  matrix_up_A<- matrix(nrow=nrow(Illig_data_up), ncol=ncol(A_up))
  for (i in 1:ncol(A_up)){
    matrix_up_A[,i]<- (A_up[,i])-(Illig_data_up$reference_never_median_beta_all)}
  colnames(matrix_up_A)<- colnames(A_up)
  rownames(matrix_up_A)<- Illig_data_up$cpgs
  
  #Calculate scores - ##UP##
  #calculate scores for each individual in the dataset
  scores_up_A<- as.numeric(rep(NA,ncol(A_up)))
  for (i in 1:ncol(A_up)){
    scores_up_A[i]<-sum(matrix_up_A[,i]*Illig_data_up$weights)}
  
  
  #Calculate diffs between SABRE beta values and the reference for each site - ##DOWN###
  matrix_down_A<- matrix(nrow=nrow(Illig_data_down), ncol=ncol(A_down))
  for (i in 1:ncol(A_down)){
    matrix_down_A[,i]<- (Illig_data_down$reference_never_median_beta_all)-(A_down[,i])}
  colnames(matrix_down_A)<- colnames(A_down)
  rownames(matrix_down_A)<- Illig_data_down$cpgs
  
  #Calculate scores - ##DOWN##
  #calculate scores for each individual in the dataset
  scores_down_A<- as.numeric(rep(NA,ncol(A_down)))
  for (i in 1:ncol(A_down)){
    scores_down_A[i]<-sum(matrix_down_A[,i]*Illig_data_down$weights)}
  
  ##combine scores
  scores_combined_A <- scores_up_A + scores_down_A
  
  return(scores_combined_A)
}
smokingScore(betas.dasen)->pheno.pf$Smoking

pdf("Plots/SmokingDistribution.pdf")
plot(density(pheno.pf$Smoking),main="Distribution of Smoking Scores")
dev.off()


###################### Create final data sets ######################

save(data.dasen,pheno.pf,pheno ,file="data.dasen_pheno.Rdata") 

### Remove XY probes
# Create dataset to be used in Chapters 3 and 4
sexchrX <- Anno[which(Anno$CHR== "X"),]
as.character(sexchrX$IlmnID)->X
sexchrY <- Anno[which(Anno$CHR== "Y"),]
as.character(sexchrY$IlmnID)->Y

data.dasen <- data.dasen[ ! rownames(data.dasen) %in% X, ]
data.dasen <- data.dasen[ ! rownames(data.dasen) %in% Y, ]
dim(data.dasen)
# removed 9676 probes
# Features  Samples 
# 401266      300

identical(rownames(pheno.pf), colnames(data))
pheno <- pheno.pf

### Remove samples <65 years
pheno<- pheno[pheno$Age_At_Visit>64, ]
dasen<- data.dasen[,rownames(pheno1)]
dim(data.dasen)
# Features  Samples 
# 401266      284
identical(rownames(pheno), colnames(data.dasen))

save(data.dasen,pheno, file="data.dasen_pheno_XYrem.RData")

### Converters-only data

# Keep only MCI and MCI converters
phenoC <- pheno[c(pheno$AtVisit_Group3==1),]
dataC <- data.dasen[,rownames(phenoC)]

#recode MCI (0) and converters (1)
phenoC$Converters<- ifelse(phenoC$AtVisit_Group4 == 1, 0, 1)

# Remove samples with unknown (3 toulouse samples) or 1 late conversion date (LNDMCI027) 
samples<-rownames(phenoC)
keep<- samples!="8667053013_R04C02" & samples!="8622007177_R02C02" & samples!="8667053127_R04C02" & samples!="8667053088_R02C01"
pheno<- phenoC[keep,]
data.dasen<-dataC[,keep]

identical(rownames(pheno), colnames(betas(data.dasen)))

save(data.dasen,pheno, file="data.dasen_pheno_Converters.RData")

rm(betas,X,Y, sexchrY,sexchrX,data1, keep, samples)

###################### Check for differences in baseline MMSE converters ######################
var.test(pheno$MMSE.MMSE_Total~pheno$Converters)

# F test to compare two variances
# 
# data:  pheno$MMSE.MMSE_Total by pheno$Converters
# F = 0.55894, num df = 66, denom df = 37, p-value = 0.03922
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.3067026 0.9716525
# sample estimates:
#   ratio of variances 
# 0.5589388 

t.test(pheno$MMSE.MMSE_Total~pheno$Converters, alternative="greater",var.equal=F)

# Welch Two Sample t-test
# 
# data:  pheno$MMSE.MMSE_Total by pheno$Converters
# t = 2.6142, df = 60.754, p-value = 0.00563
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.3915478       Inf
# sample estimates:
#   mean in group 0 mean in group 1 
# 27.26866        26.18421 

# Baseline MMSE is not equal between converters and non-converters



###################### Session info ###################### 
sessionInfo()
# R version 3.6.0 (2019-04-26)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] wateRmelon_1.30.0                                  illuminaio_0.28.0                                 
# [3] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0 ROC_1.62.0                                        
# [5] lumi_2.38.0                                        limma_3.42.2                                      
# [7] methylumi_2.32.0                                   minfi_1.32.0                                      
# [9] bumphunter_1.28.0                                  locfit_1.5-9.4                                    
# [11] iterators_1.0.12                                   foreach_1.5.0                                     
# [13] Biostrings_2.54.0                                  XVector_0.26.0                                    
# [15] SummarizedExperiment_1.16.1                        DelayedArray_0.12.2                               
# [17] BiocParallel_1.20.1                                FDb.InfiniumMethylation.hg19_2.2.0                
# [19] org.Hs.eg.db_3.10.0                                TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2           
# [21] GenomicFeatures_1.38.2                             AnnotationDbi_1.48.0                              
# [23] GenomicRanges_1.38.0                               GenomeInfoDb_1.22.1                               
# [25] IRanges_2.20.2                                     S4Vectors_0.24.3                                  
# [27] matrixStats_0.56.0                                 ggplot2_3.3.0                                     
# [29] reshape2_1.4.3                                     scales_1.1.0                                      
# [31] Biobase_2.46.0                                     BiocGenerics_0.32.0                               
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.4-1         ellipsis_0.3.0           siggenes_1.60.0          mclust_5.4.5             base64_2.0              
# [6] rstudioapi_0.11          affyio_1.56.0            bit64_0.9-7              fansi_0.4.1              xml2_1.2.5              
# [11] codetools_0.2-16         splines_3.6.0            scrime_1.3.5             knitr_1.28               Rsamtools_2.2.3         
# [16] annotate_1.64.0          dbplyr_1.4.2             HDF5Array_1.14.3         BiocManager_1.30.10      readr_1.3.1             
# [21] compiler_3.6.0           httr_1.4.1               assertthat_0.2.1         Matrix_1.2-18            cli_2.0.2               
# [26] prettyunits_1.1.1        tools_3.6.0              affy_1.64.0              gtable_0.3.0             glue_1.3.2              
# [31] GenomeInfoDbData_1.2.2   dplyr_0.8.5              rappdirs_0.3.1           doRNG_1.8.2              Rcpp_1.0.2              
# [36] vctrs_0.2.4              multtest_2.42.0          preprocessCore_1.48.0    nlme_3.1-145             rtracklayer_1.46.0      
# [41] DelayedMatrixStats_1.8.0 xfun_0.12                stringr_1.4.0            lifecycle_0.2.0          rngtools_1.5            
# [46] XML_3.99-0.3             beanplot_1.2             nleqslv_3.3.2            zlibbioc_1.32.0          MASS_7.3-51.5           
# [51] hms_0.5.3                rhdf5_2.30.1             GEOquery_2.54.1          RColorBrewer_1.1-2       curl_4.3                
# [56] memoise_1.1.0            biomaRt_2.42.1           reshape_0.8.8            stringi_1.4.3            RSQLite_2.2.0           
# [61] genefilter_1.68.0        rlang_0.4.5              pkgconfig_2.0.3          bitops_1.0-6             nor1mix_1.3-0           
# [66] lattice_0.20-40          purrr_0.3.3              Rhdf5lib_1.8.0           GenomicAlignments_1.22.1 bit_1.1-15.2            
# [71] tidyselect_1.0.0         plyr_1.8.4               magrittr_1.5             R6_2.4.1                 DBI_1.1.0               
# [76] mgcv_1.8-31              pillar_1.4.3             withr_2.1.2              survival_3.1-11          RCurl_1.98-1.1          
# [81] tibble_3.0.0             crayon_1.3.4             KernSmooth_2.23-16       BiocFileCache_1.10.2     progress_1.2.2          
# [86] grid_3.6.0               data.table_1.12.8        blob_1.2.1               digest_0.6.25            xtable_1.8-4            
# [91] tidyr_1.0.2              openssl_1.4.1            munsell_0.5.0            askpass_1.1              quadprog_1.5-8