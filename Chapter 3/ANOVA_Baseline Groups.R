###################### ANOVA: CTL vs MCI vs AD comparison  ######################
# AddNeuroMed 

setwd("")

# load packages
library(wateRmelon)
library(methylumi)
library(dplyr)
library(qqman)

# Create a folder for storing output plots
if(!dir.exists(paste0(wd,"Plots"))){
  dir.create(paste0(wd,"Plots"))
}

if(!dir.exists(paste0(wd,"Combp"))){
  dir.create(paste0(wd,"Combp"))
}

######## Load in normalised data and annotation file ########

load("") # Normalised DNA methylation data + pheno file
Anno<-read.csv("") # Annotation data
data<-betas(data.dasen)
Anno<-Anno[match(rownames(data), Anno$Name),]


######## Regress out covariates ########

# Function
RegrOut <- function( row, Age, Sex ,CD8T, CD4T, NK, Bcell, Mono,
                      Plate1 , Plate2 , Plate3){
  
  residuals <- try (
    resid(lm( row ~ Age+ Sex + 
                CD8T+ CD4T+ NK+ Bcell+ Mono+ Plate1 + Plate2 + Plate3 ,na.action=na.exclude)),
    silent=TRUE
  )
  residuals
} 


# Set cores for parallel processing 
cores<-32
cl<- makeCluster(cores)

### Apply model
RegData  <- {
  Age    <- as.numeric(pheno[ colnames(data), 'Age_At_Visit'])
  Sex <- as.numeric(pheno[ colnames(data), 'Sex'])
  CD8T <- as.numeric(pheno[ colnames(data), 'CD8T'])
  CD4T <- as.numeric(pheno[ colnames(data), 'CD4T'])
  NK <- as.numeric(pheno[ colnames(data), 'NK'])
  Bcell <- as.numeric(pheno[ colnames(data), 'Bcell'])
  Mono <- as.numeric(pheno[ colnames(data), 'Mono'])
  Plate1 <-as.numeric(pheno[ colnames(data), 'Plate1'])
  Plate2 <-as.numeric(pheno[ colnames(data), 'Plate2'])
  Plate3 <-as.numeric(pheno[ colnames(data), 'Plate3'])
  
  t( parApply(cl, data, 1, RegrOut, Age, Sex ,CD8T, CD4T, NK, Bcell, Mono,
              Plate1 , Plate2 , Plate3))
}

stopCluster(cl)

identical(rownames(pheno), colnames(RegData))

save(RegData, file="Regressed out data.RData")


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
cores<-16
cl<- makeCluster(cores)

### Apply model

Diagtukey  <- {
  
  Diag     <- pheno[ colnames(RegData), "AtVisit_Group3"]
  
  t(parApply(cl, RegData, 1, tukeys, Diag))
}

stopCluster(cl)

######## Annotate results ########

rownames(Anno) <-Anno$IlmnID

merge(Diagtukey, Anno[rownames(Diagtukey), ], by="row.names") -> TukeyResults

### Apply multiple testing correction 
### (Tukey p-values have been corrected for multiple comparisons)
colnames(TukeyResults)[c(7,11,15)]<-c("pCvM","pCvA","pMvA")

p.adjust(TukeyResults$pCvM, method="fdr")->TukeyResults$FDR_CTLvsMCI
p.adjust(TukeyResults$pCvA, method="fdr")->TukeyResults$FDR_CTLvsAD
p.adjust(TukeyResults$pMvA, method="fdr")->TukeyResults$FDR_MCIvsAD
p.adjust(TukeyResults$V2,method="fdr")->TukeyResults$FDR_ANOVA

save(TukeyResults, Diagtukey,  file="Results AnovaTukey.RData")

######## Tables for paper: ########

TukeyResults$Position<- paste0("chr ",TukeyResults$CHR,": ", TukeyResults$MAPINFO)
TukeyResults$GREAT1<- paste0(TukeyResults$GREAT_anno1," (", TukeyResults$GREAT_anno1_site,")")
TukeyResults$GREAT2<- paste0(TukeyResults$GREAT_anno2," (", TukeyResults$GREAT_anno2_site,")")
TukeyRed<-TukeyResults[,c(1,91,2,3,90,4:7,87,8:11,88,12:15,89,37,39,41,67,66,92,93)]
colnames(TukeyRed)<-c("ProbeID","Position", "F", "p (F)","FDR.Anova","Difference", "CI lower", "CI upper", "pCvM","FDR",
                      "Difference", "CI lower", "CI upper", "pCvA","FDR",
                      "Difference", "CI lower", "CI upper", "pMvA","FDR",
                      "UCSC Gene", "UCSC Gene Group", "UCSC Relation to CpG Island", "Gene Closest TSS","Distance Closest TSS", "GREAT Annotation 1",
                      "GREAT Annotation 2")
rownames(TukeyRed) <- TukeyRed$ProbeID

# Top 1000 ANOVA
TukeyRed <- TukeyRed[order(TukeyRed$`p (F)`),]
ANOVAtop <- TukeyRed[1:1000, ]
write.csv(ANOVAtop, file="Top 1000 ANOVA.csv",row.names = F)
ANOVA <- rownames(TukeyRed[1:10, ]) # For top 10 plots

# Top 1000 CTL vs MCI
TukeyRed <- TukeyRed[order(TukeyRed$pCvM),]
CvMtop <- TukeyRed[1:1000, ]
write.csv(CvMtop, file="Top 1000 CTL vs MCI.csv",row.names = F)
CTLvsMCI <- rownames(TukeyRed[1:10, ])

# Top 1000 CTL vs AD
TukeyRed <- TukeyRed[order(TukeyRed$pCvA),]
CvAtop <- TukeyRed[1:1000, ]
write.csv(CvAtop, file="Top 1000 CTL vs AD.csv",row.names = F)
CTLvsAD <- rownames(TukeyRed[1:10, ])

# Top 1000 MCI vs AD
TukeyRed <- TukeyRed[order(TukeyRed$pMvA),]
MvAtop <- TukeyRed[1:1000, ]
write.csv(MvAtop, file="Top 1000 MCI vs AD.csv",row.names = F)
MCIvsAD <- rownames(TukeyRed[1:10, ])

# for sharing:
# save(ANOVAtop,CvMtop,CvAtop, MvAtop, file="Top 1000 tables.RData")
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

CM<-DMRinput(TukeyResults, p="pCvM")
CA<-DMRinput(TukeyResults, p="pCvA")
MA<-DMRinput(TukeyResults, p="pMvA")
ANOVA2<-DMRinput(TukeyResults, p="V2")

write.table(CM, file = "/Combp/CTLvsMCI_DMRinput.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(CA, file = "/Combp/CTLvsAD_DMRinput.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(MA, file = "/Combp/MCIvsAD_DMRinput.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
write.table(ANOVA2, file = "/Combp/ANOVA_DMRinput.bed", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)

######## Create figures for publication ########

#### QQ plot ####

pdf("Plots/QQplot.pdf")
qq(TukeyResults$V2, main = expression(paste('Q-Q plot of ',italic('p'), '-values')),
   sub = paste0("Lambda: ", round( median(qchisq(1-TukeyResults$V2,1))/qchisq(0.5,1), 3)))
dev.off()


#### Manhattan Plot ####

# Recode CHR variable from factor to numeric
as.character(TukeyResults$CHR)->TukeyResults$Chromosome
as.numeric(TukeyResults$Chromosome)->TukeyResults$Chromosome

# CTL vs MCI
pdf("Plots/Manhattan CTL vs MCI.pdf")
manhattan(TukeyResults,main="CTL vs MCI", chr="Chromosome", bp="MAPINFO",p="pCvM", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(2.4e-7, base = 10),ylim=c(0,8))
dev.off()

# CTL vs AD
pdf("Plots/Manhattan CTL vs AD.pdf")
manhattan(TukeyResults,main="CTL vs AD", chr="Chromosome", bp="MAPINFO",p="pCvA", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(2.4e-7, base = 10),ylim=c(0,8))
dev.off()

# MCI vs AD
pdf("Plots/Manhattan MCI vs AD.pdf")
manhattan(TukeyResults,main="MCI vs AD", chr="Chromosome", bp="MAPINFO",p="pMvA", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(2.4e-7, base = 10),ylim=c(0,8))
dev.off()

# ANOVA
pdf("Plots/Manhattan ANOVA.pdf")
manhattan(TukeyResults,main="ANOVA p-values", chr="Chromosome", bp="MAPINFO",p="V2", col = c("#b2182b","#053061"),suggestiveline = FALSE, 
          genomewideline = -log(2.4e-7, base = 10),ylim=c(0,8))
dev.off()


######## Session Info ########

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
#   [1] qqman_0.1.4                                        dplyr_0.8.5                                       
# [3] wateRmelon_1.30.0                                  illuminaio_0.28.0                                 
# [5] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0 ROC_1.62.0                                        
# [7] lumi_2.38.0                                        limma_3.42.2                                      
# [9] methylumi_2.32.0                                   minfi_1.32.0                                      
# [11] bumphunter_1.28.0                                  locfit_1.5-9.4                                    
# [13] iterators_1.0.12                                   foreach_1.5.0                                     
# [15] Biostrings_2.54.0                                  XVector_0.26.0                                    
# [17] SummarizedExperiment_1.16.1                        DelayedArray_0.12.2                               
# [19] BiocParallel_1.20.1                                FDb.InfiniumMethylation.hg19_2.2.0                
# [21] org.Hs.eg.db_3.10.0                                TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2           
# [23] GenomicFeatures_1.38.2                             AnnotationDbi_1.48.0                              
# [25] GenomicRanges_1.38.0                               GenomeInfoDb_1.22.1                               
# [27] IRanges_2.20.2                                     S4Vectors_0.24.3                                  
# [29] matrixStats_0.56.0                                 ggplot2_3.3.0                                     
# [31] reshape2_1.4.3                                     scales_1.1.0                                      
# [33] Biobase_2.46.0                                     BiocGenerics_0.32.0                               
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.4-1         ellipsis_0.3.0           siggenes_1.60.0          mclust_5.4.5             base64_2.0              
# [6] rstudioapi_0.11          affyio_1.56.0            bit64_0.9-7              fansi_0.4.1              xml2_1.2.5              
# [11] codetools_0.2-16         splines_3.6.0            scrime_1.3.5             knitr_1.28               Rsamtools_2.2.3         
# [16] annotate_1.64.0          dbplyr_1.4.2             HDF5Array_1.14.3         BiocManager_1.30.10      readr_1.3.1             
# [21] compiler_3.6.0           httr_1.4.1               assertthat_0.2.1         Matrix_1.2-18            cli_2.0.2               
# [26] prettyunits_1.1.1        tools_3.6.0              affy_1.64.0              gtable_0.3.0             glue_1.3.2              
# [31] GenomeInfoDbData_1.2.2   rappdirs_0.3.1           doRNG_1.8.2              Rcpp_1.0.2               vctrs_0.2.4             
# [36] multtest_2.42.0          preprocessCore_1.48.0    nlme_3.1-145             rtracklayer_1.46.0       DelayedMatrixStats_1.8.0
# [41] xfun_0.12                stringr_1.4.0            lifecycle_0.2.0          rngtools_1.5             XML_3.99-0.3            
# [46] beanplot_1.2             nleqslv_3.3.2            zlibbioc_1.32.0          MASS_7.3-51.5            hms_0.5.3               
# [51] rhdf5_2.30.1             GEOquery_2.54.1          RColorBrewer_1.1-2       curl_4.3                 memoise_1.1.0           
# [56] biomaRt_2.42.1           calibrate_1.7.5          reshape_0.8.8            stringi_1.4.3            RSQLite_2.2.0           
# [61] genefilter_1.68.0        rlang_0.4.5              pkgconfig_2.0.3          bitops_1.0-6             nor1mix_1.3-0           
# [66] lattice_0.20-40          purrr_0.3.3              Rhdf5lib_1.8.0           GenomicAlignments_1.22.1 bit_1.1-15.2            
# [71] tidyselect_1.0.0         plyr_1.8.4               magrittr_1.5             R6_2.4.1                 DBI_1.1.0               
# [76] mgcv_1.8-31              pillar_1.4.3             withr_2.1.2              survival_3.1-11          RCurl_1.98-1.1          
# [81] tibble_3.0.0             crayon_1.3.4             KernSmooth_2.23-16       BiocFileCache_1.10.2     progress_1.2.2          
# [86] grid_3.6.0               data.table_1.12.8        blob_1.2.1               digest_0.6.25            xtable_1.8-4            
# [91] tidyr_1.0.2              openssl_1.4.1            munsell_0.5.0            askpass_1.1              quadprog_1.5-8 