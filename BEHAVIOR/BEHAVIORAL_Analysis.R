### Stefano Berto, 10/2018 - 
### Code for Behavioral Task - gene expression correlation analysis
### Perform expression filterin and regression, correlation and cross-validations]

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
source("Utils.R")

# Create Output directories
folder_names <- c("OUTPUTS", "PLOTS")
sapply(folder_names, dir.create)

# Load the input expression and filter for more than 0.5 for all samples
load("RAW_DATA/expData.RData")
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

# Quantile normalization
p <- normalize.quantiles(as.matrix(logCPM))
rownames(p) <- rownames(logCPM)
colnames(p) <- colnames(logCPM)

# Adjust expression
covariates <- pheno[-c(9)]
avebeta.lm<-lapply(1:nrow(p), function(x){
  lm(unlist(p[x,])~., data=covariates)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
adj.residuals<-residuals+matrix(apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(adj.residuals)<-rownames(p)
rownames(residuals) <- rownames(p)

# Load SME data 
behav <- read.table("RAW_DATA/behavioral_data.txt",header=T,sep="\t",row.names=1)
behavT <- t(behav)

# Check the UTcodes and subset data for that
pheno <- na.omit(pheno)
pheno.tmp <- pheno[match(colnames(behavT),pheno$UT.code),]
within_data <- adj.residuals[,colnames(adj.residuals) %in% rownames(pheno.tmp)]
within_data <- within_data[,match(rownames(pheno.tmp),colnames(within_data))]

# Dummy data for bootstrap
boot_data <- adj.residuals

# Correlation ANALYSIS
res <- list()
for(i in 1:nrow(behavT))
{
res[[i]] <- FastCor(within_data,behavT[i,],method="spearman",alternative="two.sided",cores=10,override=TRUE)
}

SME_cor_BEHAVIOR <- data.frame()
for(i in 1:length(res)){
  res[[i]] <- as.data.frame(res[[i]])
  res[[i]]$Gene <- rownames(res[[i]])
  rownames(res[[i]]) <- NULL
  res[[i]]$Waves <- rep(paste(rownames(behavT)[i]),nrow(res[[i]]))
}

SME_cor_BEHAVIOR <- do.call(rbind,res)
SME_cor_BEHAVIOR <- SME_cor_BEHAVIOR[c(3,4,1,2)]
save(SME_cor_BEHAVIOR,file = "OUTPUTS/BEHAVIOR_Genes.RData")

###########################
# Bootstrap Resected data #
###########################
# Select 16 subject randomly 100 times
# Recalculate correlation
# Compare with the observed one

B=100  ## select number of bootstrap resamples
index.b <- list()
Y.b <- list()
for (i in 1:B){
  set.seed(i*B+1)
  print(i)
      index.b[[i]]=sample(x=1:ncol(boot_data), size=16, replace=TRUE) # Sampling 16 columns at time.
      Y.b[[i]]=boot_data[,as.numeric(index.b[[i]])]
}

# Design lists for each frequency
PercRec <- list()
EncWord <- list()
Sessions <- list()

# Correlation analysis for each of the bootstrap. 
for(i in 1:B) {
    PercRec[[i]] <- FastCor(Y.b[[i]],mriT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    EncWord[[i]] <- FastCor(Y.b[[i]],mriT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    Sessions[[i]] <- FastCor(Y.b[[i]],mriT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultPercRec_B <- do.call(cbind,PercRec)
resultPercRec_B <- resultPercRec_B[,grep("Rho",colnames(resultPercRec_B))]

resultEncWord_B <- do.call(cbind,EncWord)
resultEncWord_B <- resultEncWord_B[,grep("Rho",colnames(resultEncWord_B))]

resultSessions_B <- do.call(cbind,Sessions)
resultSessions_B <- resultSessions_B[,grep("Rho",colnames(resultSessions_B))]

resultBoot <- rbind(resultPercRec_B,resultEncWord_B,resultSessions_B)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultBoot)){
SME_cor_BEHAVIOR$BootP[i] <- sum(abs(resultBoot[i,]) >= abs(SME_cor_BEHAVIOR$Rho[i]))/B
}

# Resave the file
save(SME_cor_BEHAVIOR,file = "OUTPUTS/BEHAVIOR_Genes.RData")

#############################
# Permutation Resected data #
#############################
# Permute the expression data 100 times
# Recalculate correlation
# Compare with the observed one

P=100  ## select number of permutation resamples
Y.b <- mclapply(1:P,mc.cores =12, function(i){apply(within_data, 2, function(col){sample(col)})}) #Randomize 100 times the within_data

# Design lists for each frequency
PercRec <- list()
EncWord <- list()
Sessions <- list()

# Correlation analysis for each of the permutation 
for(i in 1:P) {
    PercRec[[i]] <- FastCor(Y.b[[i]],mriT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    EncWord[[i]] <- FastCor(Y.b[[i]],mriT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    Sessions[[i]] <- FastCor(Y.b[[i]],mriT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultPercRec_P <- do.call(cbind,PercRec)
resultPercRec_P <- resultPercRec_P[,grep("Rho",colnames(resultPercRec_P))]

resultEncWord_P <- do.call(cbind,EncWord)
resultEncWord_P <- resultEncWord_P[,grep("Rho",colnames(resultEncWord_P))]

resultSessions_P <- do.call(cbind,Sessions)
resultSessions_P <- resultSessions_P[,grep("Rho",colnames(resultSessions_P))]

resultPerm <- rbind(resultPercRec_P,resultEncWord_P,resultSessions_P)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultPerm)){
SME_cor_BEHAVIOR$PermP[i] <- sum(abs(resultPerm[i,]) >= abs(SME_cor_BEHAVIOR$Rho[i]))/P
}

# Resave the file
save(SME_cor_BEHAVIOR,file = "OUTPUTS/BEHAVIOR_Genes.RData")

####################################
# Bootstrap Resected + Frozen Data #
####################################
# Select 16 subject randomly 100 times
# Recalculate correlation
# Compare with the observed one

phenoAll <- pheno_all
phenoAll$EpDur[is.na(phenoAll$EpDur)]=0
logCPM_ALL <- exp
logCPM_ALL <- logCPM_ALL[rownames(logCPM_ALL)%in%rownames(p),]
pAll <- normalize.quantiles(as.matrix(logCPM_ALL))
rownames(pAll) <- rownames(logCPM_ALL)
colnames(pAll) <- colnames(logCPM_ALL)

# get confounder for lm regression
covariatesAll <- phenoAll[-c(1,11)]
covariatesAll <- droplevels(covariatesAll)

# Adjust expression
avebeta.lm<-lapply(1:nrow(pAll), function(x){
  lm(unlist(pAll[x,])~., data=covariatesAll)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
adj.residuals_All<-residuals+matrix(apply(pAll, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(adj.residuals_All)<-rownames(pAll)
rownames(residuals) <- rownames(pAll)
boot_data_all <- adj.residuals_All

B=100  ## select number of bootstrap resamples
index.b <- list()
Y.b <- list()
for (i in 1:B){
  set.seed(i*B+1)
  print(i)
      index.b[[i]]=sample(x=1:ncol(boot_data_all), size=16, replace=TRUE) # Sampling 16 columns at time.
      Y.b[[i]]=boot_data_all[,as.numeric(index.b[[i]])]
}

# Design lists for each frequency
PercRec <- list()
EncWord <- list()
Sessions <- list()

# Correlation analysis for each of the bootstrap. 
for(i in 1:B) {
    PercRec[[i]] <- FastCor(Y.b[[i]],mriT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    EncWord[[i]] <- FastCor(Y.b[[i]],mriT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    Sessions[[i]] <- FastCor(Y.b[[i]],mriT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultPercRec_A <- do.call(cbind,PercRec)
resultPercRec_A <- resultPercRec_A[,grep("Rho",colnames(resultPercRec_A))]

resultEncWord_A <- do.call(cbind,EncWord)
resultEncWord_A <- resultEncWord_A[,grep("Rho",colnames(resultEncWord_A))]

resultSessions_A <- do.call(cbind,Sessions)
resultSessions_A <- resultSessions_A[,grep("Rho",colnames(resultSessions_A))]

resultBootAll <- rbind(resultPercRec_A,resultEncWord_A,resultSessions_A)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultBootAll)){
SME_cor_BEHAVIOR$BootP_All[i] <- sum(abs(resultBootAll[i,]) >= abs(SME_cor_BEHAVIOR$Rho[i]))/B
}

# Resave the file
save(SME_cor_BEHAVIOR,file = "OUTPUTS/BEHAVIOR_Genes.RData")

# Filter and save the significant genes. 
BEHAVIOR_Significant_Genes <- SME_cor_BEHAVIOR[SME_cor_BEHAVIOR$Pval < 0.05 & SME_cor_BEHAVIOR$BootP < 0.05 & SME_cor_BEHAVIOR$PermP < 0.05 & SME_cor_BEHAVIOR$BootP_All < 0.05,]
write.table(BEHAVIOR_Significant_Genes,"OUTPUTS/BEHAVIOR_Significant_Genes.txt",sep="\t",quote=F)


