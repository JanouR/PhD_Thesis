######## WGCNA: Create Modules ########
# for baseline analysis of CTL, MCI, and AD: based on all samples in the AddNeuroMed cohort
# for future conversion analysis of MCI-MCI and MCI-AD: based on the MCI subset of samples in the AddNeuroMed cohort
# The only adjustment is the loaded data and the softThreshold set to 9 for baseline and 8 for conversion analysis

library(cluster)
library(wateRmelon)
library(lm.beta)
library(WGCNA)
options(stringsAsFactors = FALSE)


### load data
# Required: normalised DNA methylation data and phenotype file
load("")
dat<-as.data.frame(betas(data.dasen))


#### Create modules ####

# 1.1. Gene filtering improves interpretability, reduces bias and saves computation time ####
var	<- (apply(dat, 1,function(x){var(x, na.rm=TRUE)})) #calculate the variance
med <- median(var)

pdf("varplot.pdf")
plot(density((var)), xlim=c(0, 0.02)) 
dev.off()


#remove non nariable probes (probes with variance <=0.0005092303 
non_var_probes<-which(var<=med)
betas_filtered<-dat[-non_var_probes,]
dim(betas_filtered)
#[1] 200633    284
rm(med,non_var_probes,var)

# 1.2. Create modules ####

dat<-t(betas_filtered)

#### 2. Check for extreme outliers ####

# Cluster samples on Euclidean distance
sampleTree = hclust(dist(dat), method = "average");

# Clustering dendrogram:
pdf(file = "")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",labels=F, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2,cex=0.4)
dev.off()
# no extreme outliers observed

##Set powers to investigate
enableWGCNAThreads(16)
powers = 1:20
##Call the network topology analysis function, set to unsigned network as this allows for correlations which are negative and positive.
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5, networkType = "unsigned")


##plot results
pdf('')
## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, Unsigned R^2", type = "n", main = paste("Scale independence"),ylim=c(0.2,1))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels= powers, cex = 0.8, col = "red")
##this line corresponds to using an R^2 cut off of h
abline(h = 0.9, col = "red")
abline(h = 0.8, col = "orange")
##Mean connectivity as
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex=0.8, col = "red")
dev.off()

### Generate modules ####
softPower = 9; # set to 8 for converters data
enableWGCNAThreads(24)
bwnet.fcx = blockwiseModules(dat, maxBlockSize = 10000, power = softPower, TOMType = "unsigned", minModuleSize = 100, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs= FALSE, verbose = 3)
save(bwnet.fcx, file = "")


