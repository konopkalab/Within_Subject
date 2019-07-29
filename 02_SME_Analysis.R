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
dir.create("processing_memory")

# Load inputs
load(here("rawdata","expAdj.RData"))
load(here("rawdata","SME_Values.RData"))

# Correlation ANALYSIS
res <- list()
for(i in 1:nrow(Lateral_resected))
{
res[[i]] <- FastCor(within_data,Lateral_resected[i,],method="spearman",alternative="two.sided",cores=12,override=TRUE)%>%
                      as.data.frame() %>%
                      rownames_to_column('Gene') 
}

SME_cor_LR <- data.frame()
for(i in 1:length(res)){
	res[[i]]$Waves <- rep(paste(rownames(Lateral_resected)[i]),nrow(res[[i]]))
}

SME_cor_LR <- do.call(rbind,res)
SME_cor_LR <- SME_cor_LR[c(1,4,2,3)]
save(SME_cor_LR,file = "processing_memory/SME_Genes.RData")

# Bootstrap Resected Data
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
deltaBoot <- list()
thetaBoot <- list()
alphaBoot <- list()
betaBoot <- list()
gammaBoot <- list()
highgammaBoot <- list()

# Correlation analysis for each of the bootstrap. 
for(i in 1:B) {
    deltaBoot[[i]] <- FastCor(Y.b[[i]],Lateral_resected[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
   	thetaBoot[[i]] <- FastCor(Y.b[[i]],Lateral_resected[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    alphaBoot[[i]] <- FastCor(Y.b[[i]],Lateral_resected[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    betaBoot[[i]] <- FastCor(Y.b[[i]],Lateral_resected[4,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    gammaBoot[[i]] <- FastCor(Y.b[[i]],Lateral_resected[5,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    highgammaBoot[[i]] <- FastCor(Y.b[[i]],Lateral_resected[6,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultDeltaB <- do.call(cbind,deltaBoot)
resultDeltaB <- resultDeltaB[,grep("Rho",colnames(resultDeltaB))]

resultThetaB <- do.call(cbind,thetaBoot)
resultThetaB <- resultThetaB[,grep("Rho",colnames(resultThetaB))]

resultAlphaB <- do.call(cbind,alphaBoot)
resultAlphaB <- resultAlphaB[,grep("Rho",colnames(resultAlphaB))]

resultBetaB <- do.call(cbind,betaBoot)
resultBetaB <- resultBetaB[,grep("Rho",colnames(resultBetaB))]

resultGammaB <- do.call(cbind,gammaBoot)
resultGammaB <- resultGammaB[,grep("Rho",colnames(resultGammaB))]

resultHighGammaB <- do.call(cbind,highgammaBoot)
resultHighGammaB <- resultHighGammaB[,grep("Rho",colnames(resultHighGammaB))]

resultBoot <- rbind(resultDeltaB,resultThetaB,resultAlphaB,resultBetaB,resultGammaB,resultHighGammaB)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultBoot)){
SME_cor_LR$BootP[i] <- sum(abs(resultBoot[i,]) >= abs(SME_cor_LR$Rho[i]))/B
}

# Resave the file
save(SME_cor_LR,file = "processing_memory/SME_Genes.RData")

# Permutation Within Data
# Permute the expression data 100 times
# Recalculate correlation
# Compare with the observed one

P=100  ## select number of permutation resamples
Y.b <- mclapply(1:P,mc.cores =12, function(i){apply(within_data, 2, function(col){sample(col)})}) #Randomize 100 times the within_data

# Design lists for each frequency
deltaPerm <- list()
thetaPerm <- list()
alphaPerm <- list()
betaPerm <- list()
gammaPerm <- list()
highgammaPerm <- list()

# Correlation analysis for each of the permutation 
for(i in 1:P) {
    deltaPerm[[i]] <- FastCor(Y.b[[i]],Lateral_resected[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
   	thetaPerm[[i]] <- FastCor(Y.b[[i]],Lateral_resected[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    alphaPerm[[i]] <- FastCor(Y.b[[i]],Lateral_resected[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    betaPerm[[i]] <- FastCor(Y.b[[i]],Lateral_resected[4,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    gammaPerm[[i]] <- FastCor(Y.b[[i]],Lateral_resected[5,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    highgammaPerm[[i]] <- FastCor(Y.b[[i]],Lateral_resected[6,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultDeltaP <- do.call(cbind,deltaPerm)
resultDeltaP <- resultDeltaP[,grep("Rho",colnames(resultDeltaP))]

resultThetaP <- do.call(cbind,thetaPerm)
resultThetaP <- resultThetaP[,grep("Rho",colnames(resultThetaP))]

resultAlphaP <- do.call(cbind,alphaPerm)
resultAlphaP <- resultAlphaP[,grep("Rho",colnames(resultAlphaP))]

resultBetaP <- do.call(cbind,betaPerm)
resultBetaP <- resultBetaP[,grep("Rho",colnames(resultBetaP))]

resultGammaP <- do.call(cbind,gammaPerm)
resultGammaP <- resultGammaP[,grep("Rho",colnames(resultGammaP))]

resultHighGammaP <- do.call(cbind,highgammaPerm)
resultHighGammaP <- resultHighGammaP[,grep("Rho",colnames(resultHighGammaP))]

resultPerm <- rbind(resultDeltaP,resultThetaP,resultAlphaP,resultBetaP,resultGammaP,resultHighGammaP)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultPerm)){
SME_cor_LR$PermP[i] <- sum(abs(resultPerm[i,]) >= abs(SME_cor_LR$Rho[i]))/P
}

# Resave the file
save(SME_cor_LR,file = "processing_memory/SME_Genes.RData")

# Bootstrap Resected + Frozen Data
# Select 16 subject randomly 100 times
# Recalculate correlation
# Compare with the observed one

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
deltaBootAll <- list()
thetaBootAll <- list()
alphaBootAll <- list()
betaBootAll <- list()
gammaBootAll <- list()
highgammaBootAll <- list()

# Correlation analysis for each of the bootstrap. 
for(i in 1:B) {
    deltaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
   	thetaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    alphaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    betaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[4,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    gammaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[5,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    highgammaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[6,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultDeltaBall <- do.call(cbind,deltaBootAll)
resultDeltaBall <- resultDeltaBall[,grep("Rho",colnames(resultDeltaBall))]

resultThetaBall <- do.call(cbind,thetaBootAll)
resultThetaBall <- resultThetaBall[,grep("Rho",colnames(resultThetaBall))]

resultAlphaBall <- do.call(cbind,alphaBootAll)
resultAlphaBall <- resultAlphaBall[,grep("Rho",colnames(resultAlphaBall))]

resultBetaBall <- do.call(cbind,betaBootAll)
resultBetaBall <- resultBetaBall[,grep("Rho",colnames(resultBetaBall))]

resultGammaBall <- do.call(cbind,gammaBootAll)
resultGammaBall <- resultGammaBall[,grep("Rho",colnames(resultGammaBall))]

resultHighGammaBall <- do.call(cbind,highgammaBootAll)
resultHighGammaBall <- resultHighGammaBall[,grep("Rho",colnames(resultHighGammaBall))]

resultBootAll <- rbind(resultDeltaBall,resultThetaBall,resultAlphaBall,resultBetaBall,resultGammaBall,resultHighGammaBall)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultBootAll)){
SME_cor_LR$BootP_All[i] <- sum(abs(resultBootAll[i,]) >= abs(SME_cor_LR$Rho[i]))/B
}

# Resave the file
save(SME_cor_LR,file = "processing_memory/SME_Genes.RData")

# Filter and save the significant genes. 
SME_Significant_Genes <- SME_cor_LR %>%
                              filter(
                                Pval < 0.05,
                                BootP < 0.05,
                                PermP < 0.05, 
                                BootP_All < 0.05)
                              
write.table(SME_Significant_Genes,"processing_memory/SME_Significant_Genes.txt",sep="\t",quote=F)

# sessionInfo
sessionInfo()
