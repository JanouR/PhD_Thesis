###################### Analysis of the Xi chromosome  ######################
# ANOVA: comparison of CTL to MCI to AD
# AddNeuroMed cohort

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


######## Load data and annotation file ########

load("") # dasen normalised data (including X and Y chromosomes)
Anno<-read.csv("")
data<-betas(data.dasen)
Anno<-Anno[match(rownames(data), Anno$Name),]

identical(as.character(rownames(pheno.pf)), as.character(colnames(data)))

### Remove samples >65 years (were originally filtered out after removal of X and Y chr)
pheno <- pheno.pf[pheno.pf$Age_At_Visit>64, ] 
data.dasen<- data.dasen[,rownames(pheno)]
dim(data.dasen)
# Features  Samples 
# 410942      284


######## split by sex ######## 

phenoF<-pheno[pheno$Sex==1, ]
phenoM<-pheno[pheno$Sex==0, ]

datF<-data.dasen[, colnames(data.dasen) %in% rownames(phenoF)]
datM<-data.dasen[, colnames(data.dasen) %in% rownames(phenoM)]

identical(rownames(phenoF), colnames(datF))
identical(rownames(phenoM), colnames(datM))

# Select Xonly
sexchrX <- Anno[which(Anno$CHR== "X"),]
as.character(sexchrX$IlmnID)->X


datF <- datF[ rownames(datF) %in% X, ]
datM <- datM[ rownames(datM) %in% X, ]
dim(datF)
dim(datM)


dataF<-as.data.frame(betas(datF))
dataM<-as.data.frame(betas(datM))

identical(rownames(phenoF), colnames(dataF))
identical(rownames(phenoM), colnames(dataM))


######## Regress out covariates ########

# Function
RegrOut <- function( row, Age, CD8T, CD4T, NK, Bcell, Mono,
                     Plate1 , Plate2 , Plate3){
  
  residuals <- try (
    resid(lm( row ~ Age+ 
                CD8T+ CD4T+ NK+ Bcell+ Mono+ Plate1 + Plate2 + Plate3 ,na.action=na.exclude)),
    silent=TRUE
  )
  residuals
} 


# Set cores for parallel processing 
cores<-32
cl<- makeCluster(cores)

### Apply model
RegDataF  <- {
  Age    <- as.numeric(phenoF[ colnames(dataF), 'Age_At_Visit'])
  CD8T <- as.numeric(phenoF[ colnames(dataF), 'CD8T'])
  CD4T <- as.numeric(phenoF[ colnames(dataF), 'CD4T'])
  NK <- as.numeric(phenoF[ colnames(dataF), 'NK'])
  Bcell <- as.numeric(phenoF[ colnames(dataF), 'Bcell'])
  Mono <- as.numeric(phenoF[ colnames(dataF), 'Mono'])
  Plate1 <-as.numeric(phenoF[ colnames(dataF), 'Plate1'])
  Plate2 <-as.numeric(phenoF[ colnames(dataF), 'Plate2'])
  Plate3 <-as.numeric(phenoF[ colnames(dataF), 'Plate3'])
  
  t( parApply(cl, dataF, 1, RegrOut, Age,CD8T, CD4T, NK, Bcell, Mono,
              Plate1 , Plate2 , Plate3))
}

RegDataM  <- {
  Age    <- as.numeric(phenoM[ colnames(dataM), 'Age_At_Visit'])
  CD8T <- as.numeric(phenoM[ colnames(dataM), 'CD8T'])
  CD4T <- as.numeric(phenoM[ colnames(dataM), 'CD4T'])
  NK <- as.numeric(phenoM[ colnames(dataM), 'NK'])
  Bcell <- as.numeric(phenoM[ colnames(dataM), 'Bcell'])
  Mono <- as.numeric(phenoM[ colnames(dataM), 'Mono'])
  Plate1 <-as.numeric(phenoM[ colnames(dataM), 'Plate1'])
  Plate2 <-as.numeric(phenoM[ colnames(dataM), 'Plate2'])
  Plate3 <-as.numeric(phenoM[ colnames(dataM), 'Plate3'])
  
  t( parApply(cl, dataM, 1, RegrOut, Age,CD8T, CD4T, NK, Bcell, Mono,
              Plate1 , Plate2 , Plate3))
}

stopCluster(cl)

identical(rownames(phenoF), colnames(RegDataF))
identical(rownames(phenoM), colnames(RegDataM))


#### Calculate Xi by: ######## 
# calculate average methylation per site per group M
# subtract that average from methylation per diagnostic group F

# Split by group
M.AD <- as.data.frame(RegDataM[ ,colnames(RegDataM) %in% phenoM[phenoM$AtVisit_Group3==2,'SentrixID'] ])
F.AD <- as.data.frame(RegDataF[ ,colnames(RegDataF) %in% phenoF[phenoF$AtVisit_Group3==2,'SentrixID'] ])

M.MCI <- as.data.frame(RegDataM[ ,colnames(RegDataM) %in% phenoM[phenoM$AtVisit_Group3==1,'SentrixID'] ])
F.MCI <- as.data.frame(RegDataF[ ,colnames(RegDataF) %in% phenoF[phenoF$AtVisit_Group3==1,'SentrixID'] ])

M.CTL <- as.data.frame(RegDataM[ ,colnames(RegDataM) %in% phenoM[phenoM$AtVisit_Group3==0,'SentrixID'] ])
F.CTL <- as.data.frame(RegDataF[ ,colnames(RegDataF) %in% phenoF[phenoF$AtVisit_Group3==0,'SentrixID'] ])

# Calculate averages by site
M.AD$Mean<-rowMeans(M.AD)
M.MCI$Mean<-rowMeans(M.MCI)
M.CTL$Mean<-rowMeans(M.CTL)

# subtract
F.ADs <- F.AD-M.AD$Mean
F.MCIs <- F.MCI-M.MCI$Mean
F.CTLs <- F.CTL-M.CTL$Mean

# combine
identical(rownames(F.ADs), rownames(F.MCIs))
identical(rownames(F.ADs), rownames(F.CTLs))

F.s<-cbind(F.CTLs,F.MCIs,F.ADs)
F.s<-F.s[, order(colnames(F.s))]
identical(as.character(colnames(F.s)), as.character(rownames(phenoF)))

save(F.s,phenoF,file="Xi data.Rdata")

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
  
  Diag     <- phenoF[ colnames(F.s), "AtVisit_Group3"]
  
  t(parApply(cl, F.s, 1, tukeys, Diag))
}

stopCluster(cl)

######## Annotate results ########

rownames(Anno) <-Anno$IlmnID

merge(DiagtukeyF, Anno[rownames(DiagtukeyF), ], by="row.names") -> TukeyResultsF

### Apply multiple testing correction 
### (Tukey p-values have been corrected for multiple comparisons)
colnames(TukeyResultsF)[c(7,11,15)]<-c("pCvM","pCvA","pMvA")

p.adjust(TukeyResultsF$pCvM, method="fdr")->TukeyResultsF$FDR_CTLvsMCI
p.adjust(TukeyResultsF$pCvA, method="fdr")->TukeyResultsF$FDR_CTLvsAD
p.adjust(TukeyResultsF$pMvA, method="fdr")->TukeyResultsF$FDR_MCIvsAD
p.adjust(TukeyResultsF$V2,method="fdr")->TukeyResultsF$FDR_ANOVA

save(TukeyResultsF, DiagtukeyF,  file="Results AnovaTukey.RData")


######## Tables for paper: ########

TukeyResultsF$Position<- paste0("chr ",TukeyResultsF$CHR,": ", TukeyResultsF$MAPINFO)
TukeyResultsF$GREAT1<- paste0(TukeyResultsF$GREAT_anno1," (", TukeyResultsF$GREAT_anno1_site,")")
TukeyResultsF$GREAT2<- paste0(TukeyResultsF$GREAT_anno2," (", TukeyResultsF$GREAT_anno2_site,")")
TukeyRedF<-TukeyResultsF[,c(1,91,2,3,90,4:7,87,8:11,88,12:15,89,37,39,41,67,66,92,93)]
colnames(TukeyRedF)<-c("ProbeID","Position", "F", "p (F)","FDR.Anova","Difference", "CI lower", "CI upper", "pCvM","FDR",
                       "Difference", "CI lower", "CI upper", "pCvA","FDR",
                       "Difference", "CI lower", "CI upper", "pMvA","FDR",
                       "UCSC Gene", "UCSC Gene Group", "UCSC Relation to CpG Island", "Gene Closest TSS","Distance Closest TSS", "GREAT Annotation 1",
                       "GREAT Annotation 2")
rownames(TukeyRedF) <- TukeyRedF$ProbeID


# Top 1000 ANOVA
TukeyRedF <- TukeyRedF[order(TukeyRedF$`p (F)`),]
ANOVAtopF <- TukeyRedF[1:1000, ]
write.csv(ANOVAtopF, file="Top 1000 ANOVA_F.csv",row.names = F)
ANOVA.F <- rownames(TukeyRedF[1:10, ]) # For top 10 plots


# Top 1000 CTL vs MCI
TukeyRedF <- TukeyRedF[order(TukeyRedF$pCvM),]
CvMtopF <- TukeyRedF[1:1000, ]
write.csv(CvMtopF, file="Top 1000 CTL vs MCI_F.csv",row.names = F)
CTLvsMCI.F <- rownames(TukeyRedF[1:10, ])

# Top 1000 CTL vs AD
TukeyRedF <- TukeyRedF[order(TukeyRedF$pCvA),]
CvAtopF <- TukeyRedF[1:1000, ]
write.csv(CvAtopF, file="Top 1000 CTL vs AD_F.csv",row.names = F)
CTLvsAD.F <- rownames(TukeyRedF[1:10, ])


# Top 1000 MCI vs AD
TukeyRedF <- TukeyRedF[order(TukeyRedF$pMvA),]
MvAtopF <- TukeyRedF[1:1000, ]
write.csv(MvAtopF, file="Top 1000 MCI vs AD_F.csv",row.names = F)
MCIvsAD.F <- rownames(TukeyRedF[1:10, ])


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


write.table(CM.f, file = "Combp/CTLvsMCI/CTLvsMCI_DMRinput_F.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(CA.f, file = "Combp/CTLvsAD/CTLvsAD_DMRinput_F.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(MA.f, file = "Combp/MCIvsAD/MCIvsAD_DMRinput_F.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(ANOVA2.f, file = "Combp/ANOVA/ANOVA_DMRinput_F.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)

######## Create figures for publication ########

#### QQ plot ####

pdf("Plots/QQplot.pdf")
qq(TukeyResultsF$V2, main = expression(paste('Q-Q plot of ',italic('p'), '-values in females')),
   sub = paste0("Lambda: ", round( median(qchisq(1-TukeyResultsF$V2,1))/qchisq(0.5,1), 3)))
dev.off()

#### Manhattan Plot ####

# Recode CHR variable from factor to numeric
as.character(TukeyResultsF$CHR)->TukeyResultsF$Chromosome
TukeyResultsF$Chromosome<- ifelse(TukeyResultsF$Chromosome=="X", 23, "NA")
#slightly adjusted qqman function to allow adjusted font size in plot annotations
source("qqman-manhattan.R")
library(calibrate)

# CTL vs MCI
pdf("Plots/Manhattan CTL vs MCI.pdf")
manhattan(TukeyResultsF,main="CTL vs MCI", chr="Chromosome", bp="MAPINFO",p="pCvM", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(5.18e-6, base = 10),ylim=c(0,15), xlab="Chromosome X position",snp="Row.names",annotatePval = 2.4e-7,annotateTop = F)
dev.off()

# CTL vs AD
pdf("Plots/Manhattan CTL vs AD.pdf")
manhattan(TukeyResultsF,main="CTL vs AD", chr="Chromosome", bp="MAPINFO",p="pCvA", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(5.18e-6, base = 10),ylim=c(0,14), xlab="Chromosome X position",snp="Row.names",annotatePval = 2.4e-7,annotateTop = F)

dev.off()

# MCI vs AD
pdf("Plots/Manhattan MCI vs AD.pdf")
manhattan(TukeyResultsF,main="MCI vs AD", chr="Chromosome", bp="MAPINFO",p="pMvA", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(5.18e-6, base = 10),ylim=c(0,14), xlab="Chromosome X position",snp="Row.names",annotatePval = 2.4e-7,annotateTop = F)

dev.off()

# ANOVA
# only annotating top 10
pdf("Plots/Manhattan ANOVA.pdf")
par(mar=c(5.1,4.1,4.1,3.1))
manhattan(TukeyResultsF,main=expression(paste('ANOVA ',italic('P'), '-values')), chr="Chromosome", bp="MAPINFO",p="V2", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(5.18e-6, base = 10),ylim=c(0,21), xlab="Chromosome X position",snp="Row.names",annotatePval = 4.2e-12,annotateTop = F)

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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] qqman_0.1.4                                        dplyr_1.0.1                                       
# [3] wateRmelon_1.26.0                                  illuminaio_0.24.0                                 
# [5] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0 ROC_1.58.0                                        
# [7] lumi_2.34.0                                        methylumi_2.28.0                                  
# [9] minfi_1.28.4                                       bumphunter_1.24.5                                 
# [11] locfit_1.5-9.4                                     iterators_1.0.12                                  
# [13] foreach_1.5.0                                      Biostrings_2.50.2                                 
# [15] XVector_0.22.0                                     SummarizedExperiment_1.12.0                       
# [17] DelayedArray_0.8.0                                 BiocParallel_1.16.6                               
# [19] FDb.InfiniumMethylation.hg19_2.2.0                 org.Hs.eg.db_3.7.0                                
# [21] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2            GenomicFeatures_1.34.8                            
# [23] AnnotationDbi_1.44.0                               GenomicRanges_1.34.0                              
# [25] GenomeInfoDb_1.18.2                                IRanges_2.16.0                                    
# [27] S4Vectors_0.20.1                                   ggplot2_3.3.2                                     
# [29] reshape2_1.4.4                                     scales_1.1.1                                      
# [31] matrixStats_0.56.0                                 limma_3.38.3                                      
# [33] Biobase_2.42.0                                     BiocGenerics_0.28.0                               
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.4-1         ellipsis_0.3.1           siggenes_1.56.0          mclust_5.4.5            
# [5] base64_2.0               rstudioapi_0.11          affyio_1.52.0            bit64_4.0.2             
# [9] fansi_0.4.1              xml2_1.3.2               codetools_0.2-16         splines_3.5.2           
# [13] Rsamtools_1.34.1         annotate_1.60.1          HDF5Array_1.10.1         BiocManager_1.30.10     
# [17] readr_1.3.1              compiler_3.5.2           httr_1.4.2               assertthat_0.2.1        
# [21] Matrix_1.2-18            cli_2.0.2                prettyunits_1.1.1        tools_3.5.2             
# [25] gtable_0.3.0             glue_1.4.2               GenomeInfoDbData_1.2.0   affy_1.60.0             
# [29] doRNG_1.8.2              Rcpp_1.0.5               vctrs_0.3.2              multtest_2.38.0         
# [33] preprocessCore_1.44.0    nlme_3.1-142             rtracklayer_1.42.2       DelayedMatrixStats_1.4.0
# [37] stringr_1.4.0            lifecycle_0.2.0          rngtools_1.5             XML_3.99-0.3            
# [41] beanplot_1.2             nleqslv_3.3.2            zlibbioc_1.28.0          MASS_7.3-51.6           
# [45] hms_0.5.3                rhdf5_2.26.2             GEOquery_2.50.5          RColorBrewer_1.1-2      
# [49] yaml_2.2.1               memoise_1.1.0            biomaRt_2.38.0           calibrate_1.7.7         
# [53] reshape_0.8.8            stringi_1.5.3            RSQLite_2.2.0            genefilter_1.64.0       
# [57] rlang_0.4.7              pkgconfig_2.0.3          bitops_1.0-6             nor1mix_1.3-0           
# [61] lattice_0.20-41          purrr_0.3.4              Rhdf5lib_1.4.3           GenomicAlignments_1.18.1
# [65] bit_4.0.4                tidyselect_1.1.0         plyr_1.8.6               magrittr_1.5            
# [69] R6_2.4.1                 generics_0.0.2           DBI_1.1.0                pillar_1.4.6            
# [73] withr_2.2.0              mgcv_1.8-31              survival_3.2-3           RCurl_1.98-1.2          
# [77] tibble_3.0.3             crayon_1.3.4             KernSmooth_2.23-16       progress_1.2.2          
# [81] grid_3.5.2               data.table_1.13.0        blob_1.2.1               digest_0.6.25           
# [85] xtable_1.8-4             tidyr_1.1.0              openssl_1.4.2            munsell_0.5.0           
# [89] askpass_1.1              quadprog_1.5-8      