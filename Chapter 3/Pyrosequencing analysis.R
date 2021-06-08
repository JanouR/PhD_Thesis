###################### Analysis of HOXB6 pyrosequencing data ###################### 
# AddNeuroMed cohort

setwd("")

library(wateRmelon)
library(qqman)
library(yarrr)


# Load Data
pdata <- read.csv("") # Data from pyrosequencing (post QC)

load("") # Load arrat methylation data, for pheno file
data <- betas(data.dasen)

# Filter out NTC and Meth+ve, and match samples
pdat <- pdata[grep("0", pdata$Sample.ID),]
rownames(pdat) <- pdat$Sample.ID
rownames(pheno) <- pheno$Sadman.Code
pheno.s <- pheno[,c(2,13,14,21,22,63:68)]
overlap <- intersect(rownames(pheno), rownames(pdat))
pheno.s <- pheno.s[overlap ,]
pdat <- pdat[overlap,]
pdat <- pdat[order(pdat$Sample.ID),]
pheno.s <- pheno.s[order(pheno.s$Sadman.Code),]
identical(rownames(pheno.s), rownames(pdat))
#264 samples overlapping

pheno.s$Plates<-pdat$Plate

pheno.s$Plate1<-ifelse(pheno.s$Plates=="Plate 1",1,0)
pheno.s$Plate2<-ifelse(pheno.s$Plates=="Plate 2",1,0)
pheno.s$Plate3<-ifelse(pheno.s$Plates=="Plate 3",1,0)

t(pdat[,c(4:8)])->pdat

# Set blanks and - to NA
pdat[pdat==""]<- NA
pdat[pdat=="-"]<- NA


######## Regress out covariates ########

# Function
RegrOut <- function( row, Age, Sex ,CD8T, CD4T, NK, Bcell, Mono,
                     Plate1 , Plate2 ){
  
  residuals <- try (
    resid(lm( row ~ Age+ Sex + 
                CD8T+ CD4T+ NK+ Bcell+ Mono+ Plate1 + Plate2  ,na.action=na.exclude)),
    silent=TRUE
  )
  residuals
} 



### Apply model
RegData  <- {
  Age    <- as.numeric(pheno.s[ colnames(pdat), 'Age_At_Visit'])
  Sex <- as.numeric(pheno.s[ colnames(pdat), 'Sex'])
  CD8T <- as.numeric(pheno.s[ colnames(pdat), 'CD8T'])
  CD4T <- as.numeric(pheno.s[ colnames(pdat), 'CD4T'])
  NK <- as.numeric(pheno.s[ colnames(pdat), 'NK'])
  Bcell <- as.numeric(pheno.s[ colnames(pdat), 'Bcell'])
  Mono <- as.numeric(pheno.s[ colnames(pdat), 'Mono'])
  Plate1 <-as.numeric(pheno.s[ colnames(pdat), 'Plate1'])
  Plate2 <-as.numeric(pheno.s[ colnames(pdat), 'Plate2'])
  
  
  t( apply( pdat, 1, RegrOut, Age, Sex ,CD8T, CD4T, NK, Bcell, Mono,
              Plate1 , Plate2 ))
}


identical(rownames(pheno.s), colnames(RegData))

save(RegData, file="Regressed out data.RData")

#### Make row with average methylation
RegData <- as.data.frame(RegData)

for (i in 1:ncol(RegData)){
  RegData[,i]<-as.character(RegData[,i])
  RegData[,i]<-as.numeric(RegData[,i])
}

RegData <- rbind(RegData, colMeans(RegData,na.rm=T ))
rownames(RegData)[6]<-"Mean"

######## Run ANOVA & Tukey's HSD ########

# Function
tukeys <- function( row, Diag){     
  
  anova.result <- aov( row ~ factor(Diag) )
  
  tukey.result <- TukeyHSD(anova.result)
  
  result <-  t(tukey.result[[1]])
  
  c(summary(anova.result)[[1]][,4][1], summary(anova.result)[[1]][,5][1], result[1:4,1],result[1:4,2],result[1:4,3]
  )
  
}


### Apply model

Diagtukey  <- {
  
  Diag     <- pheno.s[ colnames(RegData), "AtVisit_Group3"]
  
  t(apply( RegData, 1, tukeys, Diag))
}

colnames(Diagtukey)[c(1,2,6,10,14)]<-c("ANOVA.F", "ANOVA.p","pCvM","pCvA","pMvA")
as.data.frame(Diagtukey)->Diagtukey

save(pdat,pheno.s,Diagtukey, file="Raw data and Results.RData")
write.csv(Diagtukey,"Results pyro ANOVA.csv")

