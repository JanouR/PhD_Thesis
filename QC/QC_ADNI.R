###################### QC and Normalisation 450K data ###################### 
# ADNI cohort

wd<-""

setwd(wd)
library(methylumi)
library(wateRmelon)
library(dplyr)


###################### load data ######################
# data required: EPIC IDATS, phenotype data, probe annotation information,
# and cross-hybridising/SNP/bad quality probes lists

idatPath<-""
data<-readEPIC(idatPath, barcodes=pheno.ADNI$Sentrix, pdat=NULL,parallel=F,n=T,oob=F,force=F)

dim(data)
# # Features  Samples 
# # 866895      443 

pheno<-read.csv("")

Anno<-read.csv("", skip = 7)
rownames(Anno)<-Anno$IlmnID

cross<-read.table("", stringsAsFactors = FALSE)
crosslist<-cross[,1]
length(crosslist)
#[1]44210

snpProbes<-read.table("", stringsAsFactors = FALSE, header = TRUE)
dim(snpProbes)
# 340327     18
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]
dim(snpProbes)
# 10888    18
SNPlist<-snpProbes[,1]

badProbes<-read.csv("", stringsAsFactors = FALSE, header = TRUE)
dim(badProbes)
badlist<-badProbes[,1]
# [1] 977  48

if(!dir.exists("Plots")){
  dir.create("Plots")
}

rm(idatPath)

###################### Organise data ######################

pheno$Sentrix<-as.character(pheno$Sentrix)
rownames(pheno)<-pheno$Sentrix
dim(pheno)
# 443  124

identical(colnames(data), rownames(pheno))

pheno$Sex<-ifelse(pheno$PTGENDER=="Male", 0, ifelse(pheno$PTGENDER=="Female",1, "NA"))

###################### Add Houseman cell type proportions ######################

cells<-read.csv("BetasforClock_Filtered.output.csv")
cells$Sentrix<-gsub("X", "", cells$SampleID )
rownames(cells)<-cells$Sentrix
cells<- cells[, c(109, 67:72, 2,56,57)]
identical(rownames(pheno), rownames(cells))

pheno<-cbind(pheno,cells)
identical(pheno[,4], pheno[,126])
pheno[,126]<-NULL


###################### DATA QC  ######################


###### QC 1. Calculate bisulfite conversion efficiency for BS only ######
# There are several fully methylated control probes included which can be used to generate a score
# to estimate the success of the bisulfite conversion step. As they are fully methylated these should have DNA methylation
# values of ~1. The bisulfite conversion score is essentially the median of these probes, and value < 80 is taken as a failure. 

bs<-bscon(data)
summary(bs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 83.95   88.42   89.66   89.54   90.73   93.67 

pdf("Plots/Plot_BSconversion.pdf")
hist(bs, xlab = "Median % BS conversion", main = "", col="steelblue")
dev.off()

rm(bs)

###### QC 2. Remove Cross-hybridising and SNP probes ######

data1 <-data[ ! rownames(data) %in% crosslist, ]
dim(data1) #  822685      443 

data1 <-data1 [! rownames(data1) %in% SNPlist,]
dim(data1) 	# 812771     443

data1 <- data1[ ! rownames(data1) %in% badlist,]
dim(data1) # 811909 443

# There are 65 SNP probes on the 450K array, of which 59 are also on the EPIC array
# these probes are not included in the lists above, so remove manually
names<-rownames(data1)[-grep("rs.", rownames(data1))]

data1<-data1[rownames(data1) %in% names,]

identical(as.character(rownames(pheno)),as.character(colnames(data1)))

###### QC 3. Multidimensional Scaling of Sex Chomosomes ######

# this is a check that gender is correct (i.e. for sample mix up), 
# by comparing phenotype sex to profiles on X and Y chromosome
# This method is different from the 450K method, as MDScaling doesn't work well for EPIC data

betas <- betas(data1)
Anno <- Anno[match(rownames(betas), Anno$Name),]
dim(Anno)	

targets<-pheno
targets$Basename<-targets$Sentrix
RGSet <- read.metharray.exp(base= "", targets = targets, recursive=T, verbose=T)

MSet<-preprocessRaw(RGSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)

predictedSex1 <- getSex(GRset, cutoff = -2)

Sexplot_data <- as.data.frame(predictedSex1)
identical(rownames(pheno), rownames(Sexplot_data))
predictedSex <-Sexplot_data$predictedSex
pheno <- cbind(pheno,predictedSex)

#here we add the reported sex data 
Sexplot_data <- cbind(Sexplot_data, pheno$Sex)
colnames(Sexplot_data)[colnames(Sexplot_data)=="pheno$Sex"] <- "Reported Sex"

#replace blanks with NA
Sexplot_data$`Reported Sex`[Sexplot_data$`Reported Sex` == ""] <- "NA"
levels(Sexplot_data$`Reported Sex`)<-c("Male", "Female")

pdf("Plots/SexCheck.pdf")
Sexplot <- ggplot(Sexplot_data, aes(Sexplot_data$xMed,Sexplot_data$yMed, 
                                    colour = Sexplot_data$`Reported Sex`)) +
  geom_point() + theme_bw()+ scale_color_manual(values=c("darkred", "darkcyan"))+ 
  labs(x= "X Chr, median total intensity (log2)", 
       y ="Y Chr, median total intensity (log2)", 
       colour = "Reported Sex")
dev.off()

pheno$predictedSex2<-ifelse(pheno$predictedSex=="F", 1, ifelse(pheno$predictedSex=="M",0,"NA"))
identical(pheno$Sex,pheno$predictedSex2)
# True, all samples are correct

rm(MSet,RSet, Sexplot_data,targets,GRset, predictedSex, predictedSex1)


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
data<- data[,colnames(data1)]
betas <- betas(data) 
pheno1<-pheno[match(colnames(betas), pheno$Sentrix),]

betas.rs<-betas[grep("rs", rownames(betas)),]
snpCor<-cor(betas.rs)
names(snpCor)<-pheno$Sentrix 
for(i in 1:ncol(betas.rs)){
  snpCor[i,i]<-NA
}
corMax<-apply(snpCor, 1, max, na.rm = TRUE) ## calculates the maximum correlation for each sample with all other samples (except for itself)

pdf("Plots/SNPCorrelations.pdf")
hist(corMax, xlab = "Max. correlation with all other samples", main = "", col="steelblue")
dev.off()


###### QC 5. Log intensity boxplots of methylated and unmethylated raw data values #######

### extract sample intensities and plot
m_intensities<-methylated(data3)
u_intensities<-unmethylated(data3)
betas<-betas(data3)

M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)
# 
pdf("Plots/Histograms+Scatterplot_SampleIntensities_raw.pdf")
par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity",col="steelblue")
hist(U.median, xlab = "Median U intensity",col="steelblue")
  par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", xlim=c(1500,8100),ylab = "Median U intensity",col="steelblue")
dev.off()

### Plot all samples as individual lines, + density plot of averages
pdf("Plots/Raw betas+averages.pdf")
par(mfrow=c(2,1))
densityPlot(betas, main = "Raw Betas")
plot(density(betas,na.rm=T),  col='black', main = "Averages of Raw Data")
dev.off()

rm(m_intensities, u_intensities, M.median, U.median,colours)


###################### Pfilter processed beta values ######################

data.pf<- pfilter(data1)
# 0 samples having 1 % of sites with a detection p-value greater than 0.05 were removed 
# Samples removed:  
#   162 sites were removed as beadcount <3 in 5 % of samples 
# 1605 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 

identical(as.character(rownames(pheno)), as.character(colnames(data.pf)))
# [1] TRUE

###################### Remove outliers ######################
# Apply function outlix to detect outliers

outliers <- outlyx(betas)
# no outliers identified

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
hist(M.median, xlab = "Median M intensity")
hist(U.median, xlab = "Median U intensity")
par(mfrow = c(1,1))
plot(M.median, U.median, col="steelblue", pch = 16, xlab = "Median M intensity", ylab = "Median U intensity")
dev.off()


# Median by Chip
pdf("Plots/Boxplot Median Intensities by Chip_Normalised.pdf", width = 10)
nCol<-length(unique(pheno$Chip))
boxplot(M.median ~ pheno$Chip, ylab = "Median M intensity", xlab = "Chip", las = 2, col = rainbow(nCol), xaxt="n", yaxt="n")
axis(2,cex.axis=0.7, las =1)
axis(side=1, las=2,at=c(1:224), labels=c(rep("",224)),cex.axis=1)
boxplot(U.median ~ pheno$Chip, ylab = "Median U intensity", xlab = "Chip", las = 2, col = rainbow(nCol), xaxt="n", yaxt="n")
axis(2,cex.axis=0.7, las =1)
axis(side=1, las=2,at=c(1:224), labels=c(rep("",224)),cex.axis=1)
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
# [1]  0.005552258
print(dasendmrse.dasenI)
# [1] 0.005925479
print(dasendmrse.dasenII)
# [1] 0.004922134

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
#[1] 3.187153e-05 4.603833e-05 1.728323e-05
print(genki.dasenI)
#[1] 5.197269e-05 5.130380e-05 2.613607e-05
print(genki.dasenII)
#[1] 2.179637e-05 4.321114e-05 1.170169e-05


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
smokingScore(betas.dasen)->pheno$Smoking

pdf("Plots/SmokingDistribution.pdf")
plot(density(pheno$Smoking),main="Distribution of Smoking Scores")
dev.off()

###################### Create final data sets ######################

save(data.dasen,pheno ,file="data.dasen_pheno.Rdata") 

### Remove XY probes
# For replication in Chapter 4
sexchrX <- Anno[which(Anno$CHR== "X"),]
as.character(sexchrX$IlmnID)->X
sexchrY <- Anno[which(Anno$CHR== "Y"),]
as.character(sexchrY$IlmnID)->Y

data.dasen <- data.dasen[ ! rownames(data.dasen) %in% X, ]
data.dasen <- data.dasen[ ! rownames(data.dasen) %in% Y, ]
dim(data.dasen)
# how many were removed
# Features  Samples 
# 792366      443 
identical(rownames(pheno), colnames(data.dasen))

save(data.dasen,pheno, file="data.dasen_pheno_XYrem.RData")

###################### Cohort comparison ######################
### Compare pheno ADNI and AddNeuroMed ####

# T test of age between CTL:
var.test(pheno[pheno$Diagnosis==0, "Age"], pheno.ANM[pheno.ANM$AtVisit_Group3==0, 'Age_At_Visit'])

# F test to compare two variances
# 
# data:  pheno[pheno$Diagnosis == 0, "Age"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 0, "Age_At_Visit"]
# F = 1.2287, num df = 139, denom df = 88, p-value = 0.2969
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.8336631 1.7814413
# sample estimates:
#   ratio of variances 
# 1.22868

t.test(pheno[pheno$Diagnosis==0, "Age"], pheno.ANM[pheno.ANM$AtVisit_Group3==0, 'Age_At_Visit'], var.equal = T)
# Two Sample t-test
# 
# data:  pheno[pheno$Diagnosis == 0, "Age"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 0, "Age_At_Visit"]
# t = 3.3167, df = 227, p-value = 0.001061
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   1.039058 4.080801
# sample estimates:
#   mean of x mean of y 
# 76.32397  73.76404 


# T test of age between MCIs:
var.test(pheno[pheno$Diagnosis==1, "Age"], pheno.ANM[pheno.ANM$AtVisit_Group3==1, 'Age_At_Visit'])
# 
# F test to compare two variances
# 
# data:  pheno[pheno$Diagnosis == 1, "Age"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 1, "Age_At_Visit"]
# F = 1.1693, num df = 215, denom df = 108, p-value = 0.3633
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.8344454 1.6090545
# sample estimates:
#   ratio of variances 
# 1.169251

t.test(pheno[pheno$Diagnosis==1, "Age"], pheno.ANM[pheno.ANM$AtVisit_Group3==1, 'Age_At_Visit'], var.equal = T)
# Two Sample t-test
# 
# data:  pheno[pheno$Diagnosis == 1, "Age"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 1, "Age_At_Visit"]
# t = 0.50266, df = 323, p-value = 0.6155
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.9983715  1.6836239
# sample estimates:
#   mean of x mean of y 
# 75.88391  75.54128 


# T test of age between ADs:
var.test(pheno[pheno$Diagnosis==2, "Age"], pheno.ANM[pheno.ANM$AtVisit_Group3==2, 'Age_At_Visit'])

# F test to compare two variances
# 
# data:  pheno[pheno$Diagnosis == 2, "Age"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 2, "Age_At_Visit"]
# F = 1.1591, num df = 86, denom df = 85, p-value = 0.4966
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.756242 1.775642
# sample estimates:
#   ratio of variances 
# 1.159111 


t.test(pheno[pheno$Diagnosis==2, "Age"], pheno.ANM[pheno.ANM$AtVisit_Group3==2, 'Age_At_Visit'], var.equal = T)
# Two Sample t-test

# data:  pheno[pheno$Diagnosis == 2, "Age"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 2, "Age_At_Visit"]
# t = 1.8512, df = 171, p-value = 0.06586
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1093739  3.4092795
# sample estimates:
#   mean of x mean of y 
# 78.40577  76.75581


# test MMSE #
# T test of MMSE between CTL:

var.test(pheno[pheno$Diagnosis==0, "MMSE"], pheno.ANM[pheno.ANM$AtVisit_Group3==0, 'MMSE.MMSE_Total'])
# 
# F test to compare two variances
# 
# data:  pheno[pheno$Diagnosis == 0, "MMSE"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 0, "MMSE.MMSE_Total"]
# F = 1.0195, num df = 139, denom df = 88, p-value = 0.9318
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.6917366 1.4781609
# sample estimates:
#   ratio of variances 
# 1.019504 
t.test(pheno[pheno$Diagnosis==0, "MMSE"], pheno.ANM[pheno.ANM$AtVisit_Group3==0, 'MMSE.MMSE_Total'], var.equal = T)
# Two Sample t-test
# data:  pheno[pheno$Diagnosis == 0, "MMSE"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 0, "MMSE.MMSE_Total"]
# t = 0.58378, df = 227, p-value = 0.5599
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2350559  0.4329693
# sample estimates:
#   mean of x mean of y 
# 29.12143  29.02247 

# T test of MMSE between MCIs:

var.test(pheno[pheno$Diagnosis==1, "MMSE"], pheno.ANM[pheno.ANM$AtVisit_Group3==1, 'MMSE.MMSE_Total'])
# 
# F test to compare two variances
# 
# data:  pheno[pheno$Diagnosis == 1, "MMSE"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 1, "MMSE.MMSE_Total"]
# F = 1.1399, num df = 215, denom df = 108, p-value = 0.4475
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.8135192 1.5687027
# sample estimates:
#   ratio of variances 
# 1.139928 

t.test(pheno[pheno$Diagnosis==1, "MMSE"], pheno.ANM[pheno.ANM$AtVisit_Group3==1, 'MMSE.MMSE_Total'],var.equal = T)
# Two Sample t-test
# 
# data:  pheno[pheno$Diagnosis == 1, "MMSE"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 1, "MMSE.MMSE_Total"]
# t = 1.9669, df = 323, p-value = 0.05006
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.0001153967  0.9428608945
# sample estimates:
#   mean of x mean of y 
# 27.37963  26.90826  

# T test of MMSE between ADs:
var.test(pheno[pheno$Diagnosis==2, "MMSE"], pheno.ANM[pheno.ANM$AtVisit_Group3==2, 'MMSE.MMSE_Total'])
# 
# F test to compare two variances
# 
# data:  pheno[pheno$Diagnosis == 2, "MMSE"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 2, "MMSE.MMSE_Total"]
# F = 1.2301, num df = 85, denom df = 85, p-value = 0.3417
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.8017299 1.8872067
# sample estimates:
#   ratio of variances 
# 1.230053 

t.test(pheno[pheno$Diagnosis==2, "MMSE"], pheno.ANM[pheno.ANM$AtVisit_Group3==2, 'MMSE.MMSE_Total'],var.equal = T)
# Two Sample t-test
# 
# data:  pheno[pheno$Diagnosis == 2, "MMSE"] and pheno.ANM[pheno.ANM$AtVisit_Group3 == 2, "MMSE.MMSE_Total"]
# t = -1.0657, df = 170, p-value = 0.2881
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.2221224  0.6639829
# sample estimates:
#   mean of x mean of y 
# 20.05814  20.83721 


# APOE
pheno.ANM$Diag<- ifelse(pheno.ANM$AtVisit_Group3==0, "CTL", ifelse(pheno.ANM$AtVisit_Group3==1, "MCI", ifelse(pheno.ANM$AtVisit_Group3==2, "AD","NA")))
pheno$Diag<- ifelse(pheno$Diagnosis==0, "CTL", ifelse(pheno$Diagnosis==1, "MCI", ifelse(pheno$Diagnosis==2, "AD","NA")))
tab.ANM<- table(pheno.ANM$Diag, pheno.ANM$APOE_Number4)
tab.ADNI<-table(pheno$Diag, pheno$APOE4)


#CTL

matrix(nrow=2, ncol=3, c(tab.ADNI['CTL',1], tab.ANM['CTL',1], 
                         tab.ADNI['CTL',2], tab.ANM['CTL',2], 
                         tab.ADNI['CTL',3], tab.ANM['CTL',3]), dimnames= list(c("ADNI", "ANM"),c("0","1","2")))->apoe.f
fisher.test(t(apoe.f))

# 
# Fisher's Exact Test for Count Data
# 
# data:  apoe.f
# p-value = 0.2436
# alternative hypothesis: two.sided

#MCI
matrix(nrow=2, ncol=3, c(tab.ADNI['MCI',1], tab.ANM['MCI',1], 
                         tab.ADNI['MCI',2], tab.ANM['MCI',2], 
                         tab.ADNI['MCI',3], tab.ANM['MCI',3]), dimnames= list(c("ADNI", "ANM"),c("0","1","2")))->apoe.f
fisher.test(apoe.f)

# Fisher's Exact Test for Count Data
# 
# data:  apoe.f
# p-value =  0.4702
# alternative hypothesis: two.sided

#AD
matrix(nrow=2, ncol=3, c(tab.ADNI['AD',1], tab.ANM['AD',1], 
                         tab.ADNI['AD',2], tab.ANM['AD',2], 
                         tab.ADNI['AD',3], tab.ANM['AD',3]), dimnames= list(c("ADNI", "ANM"),c("0","1","2")))->apoe.f
fisher.test(apoe.f)
# Fisher's Exact Test for Count Data
# 
# data:  apoe.f
# p-value = 0.07129
# alternative hypothesis: two.sided


### Converters-only data ####

# Keep only MCI and MCI converters
pheno <- pheno[c(pheno$Diagnosis==1),]
data.dasen <- data.dasen[,rownames(pheno)]

#recode MCI (0) and converters (1)
pheno$Converters<- ifelse(pheno$Diagnosis4 == 1, 0, 1)

save(data.dasen,pheno, file="data.dasen_pheno_Converters.RData")

rm(betas,X,Y, sexchrY,sexchrX,data1)

###################### Check for differences in baseline MMSE converters ######################
var.test(pheno$MMSE~pheno$Converters)

# F test to compare two variances
# 
# data:  pheno$MMSE by pheno$Converters
# F = 0.68489, num df = 131, denom df = 83, p-value = 0.05242
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.4591802 1.0039596
# sample estimates:
#   ratio of variances 
# 0.6848869  

t.test(pheno$MMSE~pheno$Converters, alternative="greater",var.equal=T)

# Two Sample t-test
# 
# data:  pheno$MMSE by pheno$Converters
# t = 6.5028, df = 214, p-value = 2.738e-10
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   1.291696      Inf
# sample estimates:
#   mean in group 0 mean in group 1 
# 28.05303        26.32143 

# Baseline MMSE is not equal between converters and non-converters

###################### Session info ###################### 
sessionInfo()
# R version 3.5.2 (2018-12-20)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0 IlluminaHumanMethylationEPICmanifest_0.3.0          dplyr_1.0.1                                        
# [4] wateRmelon_1.26.0                                   illuminaio_0.24.0                                   IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0 
# [7] ROC_1.58.0                                          lumi_2.34.0                                         limma_3.38.3                                       
# [10] methylumi_2.28.0                                    minfi_1.28.4                                        bumphunter_1.24.5                                  
# [13] locfit_1.5-9.4                                      iterators_1.0.12                                    foreach_1.5.0                                      
# [16] Biostrings_2.50.2                                   XVector_0.22.0                                      SummarizedExperiment_1.12.0                        
# [19] DelayedArray_0.8.0                                  BiocParallel_1.16.6                                 FDb.InfiniumMethylation.hg19_2.2.0                 
# [22] org.Hs.eg.db_3.7.0                                  TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2             GenomicFeatures_1.34.8                             
# [25] AnnotationDbi_1.44.0                                GenomicRanges_1.34.0                                GenomeInfoDb_1.18.2                                
# [28] IRanges_2.16.0                                      S4Vectors_0.20.1                                    matrixStats_0.56.0                                 
# [31] ggplot2_3.3.2                                       reshape2_1.4.4                                      scales_1.1.1                                       
# [34] Biobase_2.42.0                                      BiocGenerics_0.28.0                                
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.4-1         ellipsis_0.3.1           siggenes_1.56.0          mclust_5.4.5             base64_2.0               rstudioapi_0.11         
# [7] farver_2.0.3             affyio_1.52.0            bit64_4.0.2              fansi_0.4.1              xml2_1.3.2               codetools_0.2-16        
# [13] splines_3.5.2            Rsamtools_1.34.1         annotate_1.60.1          HDF5Array_1.10.1         BiocManager_1.30.10      readr_1.3.1             
# [19] compiler_3.5.2           httr_1.4.2               assertthat_0.2.1         Matrix_1.2-18            cli_2.0.2                prettyunits_1.1.1       
# [25] tools_3.5.2              gtable_0.3.0             glue_1.4.2               GenomeInfoDbData_1.2.0   affy_1.60.0              doRNG_1.8.2             
# [31] Rcpp_1.0.5               vctrs_0.3.2              multtest_2.38.0          preprocessCore_1.44.0    nlme_3.1-142             rtracklayer_1.42.2      
# [37] DelayedMatrixStats_1.4.0 stringr_1.4.0            lifecycle_0.2.0          rngtools_1.5             XML_3.99-0.3             beanplot_1.2            
# [43] nleqslv_3.3.2            zlibbioc_1.28.0          MASS_7.3-51.6            hms_0.5.3                rhdf5_2.26.2             GEOquery_2.50.5         
# [49] RColorBrewer_1.1-2       yaml_2.2.1               memoise_1.1.0            biomaRt_2.38.0           reshape_0.8.8            stringi_1.5.3           
# [55] RSQLite_2.2.0            genefilter_1.64.0        rlang_0.4.7              pkgconfig_2.0.3          bitops_1.0-6             nor1mix_1.3-0           
# [61] lattice_0.20-41          purrr_0.3.4              Rhdf5lib_1.4.3           labeling_0.3             GenomicAlignments_1.18.1 bit_4.0.4               
# [67] tidyselect_1.1.0         plyr_1.8.6               magrittr_1.5             R6_2.4.1                 generics_0.0.2           DBI_1.1.0               
# [73] pillar_1.4.6             withr_2.2.0              mgcv_1.8-31              survival_3.2-3           RCurl_1.98-1.2           tibble_3.0.3            
# [79] crayon_1.3.4             utf8_1.1.4               KernSmooth_2.23-16       progress_1.2.2           grid_3.5.2               data.table_1.13.0       
# [85] blob_1.2.1               digest_0.6.25            xtable_1.8-4             tidyr_1.1.0              openssl_1.4.2            munsell_0.5.0           
# [91] askpass_1.1              quadprog_1.5-8 
