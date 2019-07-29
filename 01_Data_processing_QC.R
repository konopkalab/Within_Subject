### Stefano Berto, 10/2018 - 
### Code for SME - gene expression correlation analysis
### Perform expression filtering, regression, correlation and cross-validations]

# Load libraries
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(sqldf))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(made4))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(matrixStats))
source("Utils.R")

# Create Output directories
dir.create("processing_qc")

# Load the input expression and filter for more than 0 for all samples
load(here("rawdata","expData.RData"))

logCPM <- exp[,grep("Lega",names(exp))]
perc <- 90
vec=round((ncol(logCPM) * perc)/100)
notAllZero = (rowSums(logCPM>0)>=vec)
logCPM=logCPM[notAllZero,]

# Load Covariates
pheno$Class <- NULL # remove the different data class
pheno$Pmi <- NULL # same PMI for resected tissues
pheno <- pheno[match(colnames(logCPM),rownames(pheno)),]
pheno <- droplevels(pheno)

# Distribution before quantile
pdf("processing_qc/Distribution_logCPM_LEGA.pdf",width=6,height=5)
colramp = colorRampPalette(c(3,"white",2))(ncol(logCPM))
plot(density(logCPM[,1]),col=colramp[1],lwd=3,ylim=c(0,.3))
for(i in 2:ncol(logCPM)){lines(density(logCPM[,i]),lwd=3,col=colramp[i])}
dev.off()

# Quantile normalization
p <- normalize.quantiles(as.matrix(logCPM))
rownames(p) <- rownames(logCPM)
colnames(p) <- colnames(logCPM)

# Variance explained
covariates <- pheno[-c(9)]
var <- VarExp(p,covariates,5,FALSE)
pdf("processing_qc/Variance_Explained_Lega.pdf",width=8,height=6)
plotVarExp(var,"Variance Explained")
dev.off()

# Adjust expression of all resected tissue
avebeta.lm<-lapply(1:nrow(p), function(x){
  lm(unlist(p[x,])~., data=covariates)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
boot_data<-residuals+matrix(apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(boot_data)<-rownames(p)

# Load SME data and filter for matched data
load(here("rawdata","SME_Values.RData"))
pheno.ws <- na.omit(pheno)
pheno.ws <- pheno.ws[match(colnames(Lateral_resected),pheno.ws$UT.code),]
within_data <- boot_data[,colnames(boot_data) %in% rownames(pheno.ws)]
within_data <- within_data[,match(rownames(pheno.ws),colnames(within_data))]



# Adjust expression of all resected tissue + frozed brain
phenoAll <- pheno_all
phenoAll$EpDur[is.na(phenoAll$EpDur)]=0
logCPM_ALL <- exp
logCPM_ALL <- logCPM_ALL[,match(rownames(phenoAll),colnames(logCPM_ALL))]
logCPM_ALL <- logCPM_ALL[rownames(logCPM_ALL)%in%rownames(p),]
pAll <- normalize.quantiles(as.matrix(logCPM_ALL))
rownames(pAll) <- rownames(logCPM_ALL)
colnames(pAll) <- colnames(logCPM_ALL)

# get confounder for lm regression
covariatesAll <- phenoAll[-c(1,11)]
covariatesAll <- droplevels(covariatesAll)

var <- VarExp(pAll,covariatesAll,5,FALSE)
pdf("processing_qc/Variance_Explained_All.pdf",width=8,height=6)
plotVarExp(var,"Variance Explained")
dev.off()

# Adjust expression
avebeta.lm<-lapply(1:nrow(pAll), function(x){
  lm(unlist(pAll[x,])~., data=covariatesAll)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
boot_data_all<-residuals+matrix(apply(pAll, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(boot_data_all)<-rownames(pAll)

# Resave the file
save(within_data,pheno.ws,boot_data,pheno,boot_data_all,phenoAll,file = "rawdata/expAdj.RData")

# sessionInfo
sessionInfo()
