###################### Stratified anaylysis of the X and Y chromosomes  ######################
# ANOVA: comparison of CTL to MCI to AD
# Replication analysis in the ADNI cohort

setwd("")

# load packages
library(wateRmelon)
library(methylumi)
library(dplyr)
library(qqman)


if(!dir.exists(paste0(wd,"Plots"))){
  dir.create(paste0(wd,"Plots"))
}

if(!dir.exists(paste0(wd,"Combp"))){
  dir.create(paste0(wd,"Combp"))
}

######## Load in normalised data and annotation file ########

load("") # Load normalised data still containing the X and Y chromosome data
Anno<-read.csv("", skip=7)
data<-betas(data.dasen)
Anno<-Anno[match(rownames(data), Anno$Name),]

identical(as.character(rownames(pheno)), as.character(colnames(data)))


# split by sex ####

phenoF<-pheno[pheno$Sex==1, ]
phenoM<-pheno[pheno$Sex==0, ]

datF<-data.dasen[, colnames(data.dasen) %in% rownames(phenoF)]
datM<-data.dasen[, colnames(data.dasen) %in% rownames(phenoM)]

identical(rownames(phenoF), colnames(datF))
identical(rownames(phenoM), colnames(datM))

# Select X and Y only
sexchrX <- Anno[which(Anno$CHR== "X"),]
as.character(sexchrX$IlmnID)->X
sexchrY <- Anno[which(Anno$CHR== "Y"),]
as.character(sexchrY$IlmnID)->Y
XY<-c(X,Y)

datF <- datF[ rownames(datF) %in% X, ]
datM <- datM[ rownames(datM) %in% XY, ]
dim(datF)
dim(datM)


dataF<-as.data.frame(betas(datF))
dataM<-as.data.frame(betas(datM))

identical(rownames(phenoF), colnames(dataF))
identical(rownames(phenoM), colnames(dataM))

######## Regress out covariates ########

# Function
RegrOut <- function( row, Age, CD8T, CD4T, NK, Bcell, Mono){
  
  residuals <- try (
    resid(lm( row ~ Age+ 
                CD8T+ CD4T+ NK+ Bcell+ Mono ,na.action=na.exclude)),
    silent=TRUE
  )
  residuals
} 


# Set cores for parallel processing 
cores<-32
cl<- makeCluster(cores)

### Apply model
RegDataF  <- {
  Age    <- as.numeric(phenoF[ colnames(dataF), 'Age'])
  CD8T <- as.numeric(phenoF[ colnames(dataF), 'CD8T'])
  CD4T <- as.numeric(phenoF[ colnames(dataF), 'CD4T'])
  NK <- as.numeric(phenoF[ colnames(dataF), 'NK'])
  Bcell <- as.numeric(phenoF[ colnames(dataF), 'Bcell'])
  Mono <- as.numeric(phenoF[ colnames(dataF), 'Mono'])
  
  t( parApply(cl, dataF, 1, RegrOut, Age,CD8T, CD4T, NK, Bcell, Mono))
}

RegDataM  <- {
  Age    <- as.numeric(phenoM[ colnames(dataM), 'Age'])
  CD8T <- as.numeric(phenoM[ colnames(dataM), 'CD8T'])
  CD4T <- as.numeric(phenoM[ colnames(dataM), 'CD4T'])
  NK <- as.numeric(phenoM[ colnames(dataM), 'NK'])
  Bcell <- as.numeric(phenoM[ colnames(dataM), 'Bcell'])
  Mono <- as.numeric(phenoM[ colnames(dataM), 'Mono'])
  
  t( parApply(cl, dataM, 1, RegrOut, Age,CD8T, CD4T, NK, Bcell, Mono))
}

stopCluster(cl)

identical(rownames(phenoF), colnames(RegDataF))
identical(rownames(phenoM), colnames(RegDataM))

save(RegDataF,RegDataM,phenoF,phenoM, file="Regressed out data.RData")


######## Run ANOVA & Tukey's HSD ########

# Function
tukeys <- function( row, Diag){     
  
  anova.result <- aov( row ~ factor(Diag) )
  
  tukey.result <- TukeyHSD(anova.result)
  
  result <-  t(tukey.result[[1]])
  
  c(summary(anova.result)[[1]][,4][1], summary(anova.result)[[1]][,5][1], result[1:4,1],result[1:4,2],result[1:4,3]
  )
  
}

# Set cores for parallel processing
cores<-32
cl<- makeCluster(cores)

### Apply model

DiagtukeyF  <- {
  
  Diag     <- phenoF[ colnames(RegDataF), "Diagnosis"]
  
  t(parApply(cl, RegDataF, 1, tukeys, Diag))
}

DiagtukeyM  <- {
  
  Diag     <- phenoM[ colnames(RegDataM), "Diagnosis"]
  
  t(parApply(cl, RegDataM, 1, tukeys, Diag))
}

stopCluster(cl)

######## Annotate results ########

rownames(Anno) <-Anno$IlmnID
Anno.f<-Anno[,c(1,2,12,13,16:20)]
rownames(Anno.f)<-Anno.f$IlmnID

merge(DiagtukeyF, Anno.f[rownames(DiagtukeyF), ], by="row.names") -> TukeyResultsF
merge(DiagtukeyM, Anno.f[rownames(DiagtukeyM), ], by="row.names") -> TukeyResultsM

### (Tukey p-values have been corrected for multiple comparisons)
colnames(TukeyResultsF)[c(7,11,15)]<-c("pCvM","pCvA","pMvA")
colnames(TukeyResultsM)[c(7,11,15)]<-c("pCvM","pCvA","pMvA")

p.adjust(TukeyResultsF$pCvM, method="fdr")->TukeyResultsF$FDR_CTLvsMCI
p.adjust(TukeyResultsF$pCvA, method="fdr")->TukeyResultsF$FDR_CTLvsAD
p.adjust(TukeyResultsF$pMvA, method="fdr")->TukeyResultsF$FDR_MCIvsAD
p.adjust(TukeyResultsF$V2,method="fdr")->TukeyResultsF$FDR_ANOVA

p.adjust(TukeyResultsM$pCvM, method="fdr")->TukeyResultsM$FDR_CTLvsMCI
p.adjust(TukeyResultsM$pCvA, method="fdr")->TukeyResultsM$FDR_CTLvsAD
p.adjust(TukeyResultsM$pMvA, method="fdr")->TukeyResultsM$FDR_MCIvsAD
p.adjust(TukeyResultsM$V2,method="fdr")->TukeyResultsM$FDR_ANOVA

save(TukeyResultsM,TukeyResultsF, DiagtukeyM, DiagtukeyF,  file="Results AnovaTukey.RData")


######## Tables for paper: ########


TukeyResultsF$Position<- paste0("chr ",TukeyResultsF$CHR,": ", TukeyResultsF$MAPINFO)
TukeyRedF<-TukeyResultsF[,c(1,29,2,3,28,4:7,25,8:11,26,12:15,27,20,22,24)]
colnames(TukeyRedF)<-c("ProbeID","Position", "F", "p (F)","FDR.Anova","Difference", "CI lower", "CI upper", "pCvM","FDR",
                       "Difference", "CI lower", "CI upper", "pCvA","FDR",
                       "Difference", "CI lower", "CI upper", "pMvA","FDR",
                       "UCSC Gene", "UCSC Gene Group", "UCSC Relation to CpG Island")
rownames(TukeyRedF) <- TukeyRedF$ProbeID

TukeyResultsM$Position<- paste0("chr ",TukeyResultsM$CHR,": ", TukeyResultsM$MAPINFO)
TukeyRedM<-TukeyResultsM[,c(1,29,2,3,28,4:7,25,8:11,26,12:15,27,20,22,24)]
colnames(TukeyRedM)<-c("ProbeID","Position", "F", "p (F)","FDR.Anova","Difference", "CI lower", "CI upper", "pCvM","FDR",
                       "Difference", "CI lower", "CI upper", "pCvA","FDR",
                       "Difference", "CI lower", "CI upper", "pMvA","FDR",
                       "UCSC Gene", "UCSC Gene Group", "UCSC Relation to CpG Island")
rownames(TukeyRedM) <- TukeyRedM$ProbeID

# Top 1000 ANOVA
TukeyRedF <- TukeyRedF[order(TukeyRedF$`p (F)`),]
ANOVAtopF <- TukeyRedF[1:1000, ]
write.csv(ANOVAtopF, file="Top 1000 ANOVA_F.csv",row.names = F)
ANOVA.F <- rownames(TukeyRedF[1:10, ]) # For top 10 plots

TukeyRedM <- TukeyRedM[order(TukeyRedM$`p (F)`),]
ANOVAtopM <- TukeyRedM[1:1000, ]
write.csv(ANOVAtopM, file="Top 1000 ANOVA_M.csv",row.names = F)
ANOVA.M <- rownames(TukeyRedM[1:10, ]) # For top 10 plots


# Top 1000 CTL vs MCI
TukeyRedF <- TukeyRedF[order(TukeyRedF$pCvM),]
CvMtopF <- TukeyRedF[1:1000, ]
write.csv(CvMtopF, file="Top 1000 CTL vs MCI_F.csv",row.names = F)
CTLvsMCI.F <- rownames(TukeyRedF[1:10, ])

TukeyRedM <- TukeyRedM[order(TukeyRedM$pCvM),]
CvMtopM <- TukeyRedM[1:1000, ]
write.csv(CvMtopM, file="Top 1000 CTL vs MCI_M.csv",row.names = F)
CTLvsMCI.M <- rownames(TukeyRedM[1:10, ])

# Top 1000 CTL vs AD
TukeyRedF <- TukeyRedF[order(TukeyRedF$pCvA),]
CvAtopF <- TukeyRedF[1:1000, ]
write.csv(CvAtopF, file="Top 1000 CTL vs AD_F.csv",row.names = F)
CTLvsAD.F <- rownames(TukeyRedF[1:10, ])

TukeyRedM <- TukeyRedM[order(TukeyRedM$pCvA),]
CvAtopM <- TukeyRedM[1:1000, ]
write.csv(CvAtopM, file="Top 1000 CTL vs AD_M.csv",row.names = F)
CTLvsAD.M <- rownames(TukeyRedM[1:10, ])

# Top 1000 MCI vs AD
TukeyRedF <- TukeyRedF[order(TukeyRedF$pMvA),]
MvAtopF <- TukeyRedF[1:1000, ]
write.csv(MvAtopF, file="Top 1000 MCI vs AD_F.csv",row.names = F)
MCIvsAD.F <- rownames(TukeyRedF[1:10, ])

TukeyRedM <- TukeyRedM[order(TukeyRedM$pMvA),]
MvAtopM <- TukeyRedM[1:1000, ]
write.csv(MvAtopM, file="Top 1000 MCI vs AD_M.csv",row.names = F)
MCIvsAD.M <- rownames(TukeyRedM[1:10, ])

######## Create files for comb-p ########

DMRinput <- function(results,  chr = "CHR", 
                     mapinfo = "MAPINFO", p = "P"){
  dmr <- data.frame("chrom" = paste0("chr", as.character(results[ , chr])),
                    "start" = results[, mapinfo],
                    "end" = results[, mapinfo] + 1,
                    "pvalue" = results[ , p])
  dmr <- summarise(group_by(dmr, chrom, start), end = end, pvalue = pvalue)
  colnames(dmr) <- c("#chrom", "start", "end", "pvalue")
  return(dmr)
}

CM.f<-DMRinput(TukeyResultsF, p="pCvM")
CA.f<-DMRinput(TukeyResultsF, p="pCvA")
MA.f<-DMRinput(TukeyResultsF, p="pMvA")
ANOVA2.f<-DMRinput(TukeyResultsF, p="V2")

CM.m<-DMRinput(TukeyResultsM, p="pCvM")
CA.m<-DMRinput(TukeyResultsM, p="pCvA")
MA.m<-DMRinput(TukeyResultsM, p="pMvA")
ANOVA2.m<-DMRinput(TukeyResultsM, p="V2")

write.table(CM.f, file = "Combp/CTLvsMCI/CTLvsMCI_DMRinput_F.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(CA.f, file = "Combp/CTLvsAD/CTLvsAD_DMRinput_F.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(MA.f, file = "Combp/MCIvsAD/MCIvsAD_DMRinput_F.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(ANOVA2.f, file = "Combp/ANOVA/ANOVA_DMRinput_F.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)


write.table(CM.m, file = "Combp/CTLvsMCI/CTLvsMCI_DMRinput_M.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(CA.m, file = "Combp/CTLvsAD/CTLvsAD_DMRinput_M.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(MA.m, file = "Combp/MCIvsAD/MCIvsAD_DMRinput_M.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(ANOVA2.m, file = "Combp/ANOVA/ANOVA_DMRinput_M.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)

######## Create figures for publication ########

#### QQ plot ####

pdf("Plots/QQplot.pdf")
qq(TukeyResultsF$V2, main = expression(paste('Q-Q plot of ',italic('p'), '-values in females')),
   sub = paste0("Lambda: ", round( median(qchisq(1-TukeyResultsF$V2,1))/qchisq(0.5,1), 3)))
qq(TukeyResultsM$V2, main = expression(paste('Q-Q plot of ',italic('p'), '-values in males')),
   sub = paste0("Lambda: ", round( median(qchisq(1-TukeyResultsM$V2,1))/qchisq(0.5,1), 3)))

dev.off()

#### Manhattan Plot ####

# Recode CHR variable from factor to numeric
as.character(TukeyResultsF$CHR)->TukeyResultsF$Chromosome
TukeyResultsF$Chromosome<- ifelse(TukeyResultsF$Chromosome=="X", 23, "NA")
as.character(TukeyResultsM$CHR)->TukeyResultsM$Chromosome
TukeyResultsM$Chromosome<- ifelse(TukeyResultsM$Chromosome=="X", 23, 24)

# CTL vs MCI
pdf("Plots/Manhattan CTL vs MCI.pdf")
manhattan(TukeyResultsF,main="CTL vs MCI in females", chr="Chromosome", bp="MAPINFO",p="pCvM", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(9e-8, base = 10),ylim=c(0,8), chrlabs=c("X"), xlab="Chromosome X position")
manhattan(TukeyResultsM,main="CTL vs MCI in males", chr="Chromosome", bp="MAPINFO",p="pCvM", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(9e-8, base = 10),ylim=c(0,8), chrlabs=c("X","Y"))
dev.off()

# CTL vs AD
pdf("Plots/Manhattan CTL vs AD.pdf")
manhattan(TukeyResultsF,main="CTL vs AD in females", chr="Chromosome", bp="MAPINFO",p="pCvA", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(9e-8, base = 10),ylim=c(0,8), chrlabs=c("X"), xlab="Chromosome X position")
manhattan(TukeyResultsM,main="CTL vs AD in males", chr="Chromosome", bp="MAPINFO",p="pCvA", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(9e-8, base = 10),ylim=c(0,8), chrlabs=c("X","Y"))
dev.off()

# MCI vs AD
pdf("Plots/Manhattan MCI vs AD.pdf")
manhattan(TukeyResultsF,main="MCI vs AD in females", chr="Chromosome", bp="MAPINFO",p="pMvA", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(9e-8, base = 10),ylim=c(0,8), chrlabs=c("X"), xlab="Chromosome X position")
manhattan(TukeyResultsM,main="MCI vs AD in males", chr="Chromosome", bp="MAPINFO",p="pMvA", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(9e-8, base = 10),ylim=c(0,8), chrlabs=c("X","Y"))
dev.off()

# ANOVA
pdf("Plots/Manhattan ANOVA.pdf")
manhattan(TukeyResultsF,main=expression(paste('ANOVA ',italic('P'), '-values in Females')), chr="Chromosome", bp="MAPINFO",p="V2", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(9e-8, base = 10),ylim=c(0,8), chrlabs=c("X"), xlab="Chromosome X position")
manhattan(TukeyResultsM,main=expression(paste('ANOVA ',italic('P'), '-values in Males')), chr="Chromosome", bp="MAPINFO",p="V2", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(9e-8, base = 10),ylim=c(0,8), chrlabs=c("X","Y"))
dev.off()


######## Session Info ########

sessionInfo()
# R version 3.5.2 (2018-12-20)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] qqman_0.1.4                                        dplyr_1.0.1                                        wateRmelon_1.26.0                                 
# [4] illuminaio_0.24.0                                  IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0 ROC_1.58.0                                        
# [7] lumi_2.34.0                                        methylumi_2.28.0                                   minfi_1.28.4                                      
# [10] bumphunter_1.24.5                                  locfit_1.5-9.4                                     iterators_1.0.12                                  
# [13] foreach_1.5.0                                      Biostrings_2.50.2                                  XVector_0.22.0                                    
# [16] SummarizedExperiment_1.12.0                        DelayedArray_0.8.0                                 BiocParallel_1.16.6                               
# [19] FDb.InfiniumMethylation.hg19_2.2.0                 org.Hs.eg.db_3.7.0                                 TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2           
# [22] GenomicFeatures_1.34.8                             AnnotationDbi_1.44.0                               GenomicRanges_1.34.0                              
# [25] GenomeInfoDb_1.18.2                                IRanges_2.16.0                                     S4Vectors_0.20.1                                  
# [28] ggplot2_3.3.2                                      reshape2_1.4.4                                     scales_1.1.1                                      
# [31] matrixStats_0.56.0                                 limma_3.38.3                                       Biobase_2.42.0                                    
# [34] BiocGenerics_0.28.0                               
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.4-1         ellipsis_0.3.1           siggenes_1.56.0          mclust_5.4.5             base64_2.0               rstudioapi_0.11         
# [7] affyio_1.52.0            bit64_4.0.2              fansi_0.4.1              xml2_1.3.2               codetools_0.2-16         splines_3.5.2           
# [13] Rsamtools_1.34.1         annotate_1.60.1          HDF5Array_1.10.1         BiocManager_1.30.10      readr_1.3.1              compiler_3.5.2          
# [19] httr_1.4.2               assertthat_0.2.1         Matrix_1.2-18            cli_2.0.2                prettyunits_1.1.1        tools_3.5.2             
# [25] gtable_0.3.0             glue_1.4.2               GenomeInfoDbData_1.2.0   affy_1.60.0              doRNG_1.8.2              Rcpp_1.0.5              
# [31] vctrs_0.3.2              multtest_2.38.0          preprocessCore_1.44.0    nlme_3.1-142             rtracklayer_1.42.2       DelayedMatrixStats_1.4.0
# [37] stringr_1.4.0            lifecycle_0.2.0          rngtools_1.5             XML_3.99-0.3             beanplot_1.2             nleqslv_3.3.2           
# [43] zlibbioc_1.28.0          MASS_7.3-51.6            hms_0.5.3                rhdf5_2.26.2             GEOquery_2.50.5          RColorBrewer_1.1-2      
# [49] yaml_2.2.1               memoise_1.1.0            biomaRt_2.38.0           calibrate_1.7.7          reshape_0.8.8            stringi_1.5.3           
# [55] RSQLite_2.2.0            genefilter_1.64.0        rlang_0.4.7              pkgconfig_2.0.3          bitops_1.0-6             nor1mix_1.3-0           
# [61] lattice_0.20-41          purrr_0.3.4              Rhdf5lib_1.4.3           GenomicAlignments_1.18.1 bit_4.0.4                tidyselect_1.1.0        
# [67] plyr_1.8.6               magrittr_1.5             R6_2.4.1                 generics_0.0.2           DBI_1.1.0                pillar_1.4.6            
# [73] withr_2.2.0              mgcv_1.8-31              survival_3.2-3           RCurl_1.98-1.2           tibble_3.0.3             crayon_1.3.4            
# [79] KernSmooth_2.23-16       progress_1.2.2           grid_3.5.2               data.table_1.13.0        blob_1.2.1               digest_0.6.25           
# [85] xtable_1.8-4             tidyr_1.1.0              openssl_1.4.2            munsell_0.5.0            askpass_1.1              quadprog_1.5-8  