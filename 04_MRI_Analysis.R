## Correlation Analysis
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

dir.create("processing_mri")

# Load the input expression and filter for more than 0.5 for all samples
load(here("rawdata","expAdj.RData"))
load(here("rawdata","MRI_Values.RData"))

# Correlation ANALYSIS
res <- list()
for(i in 1:nrow(mriT))
{
res[[i]] <- FastCor(within_data,mriT[i,],method="spearman",alternative="two.sided",cores=10,override=TRUE) %>%
                      as.data.frame() %>%
                      rownames_to_column('Gene')
}


MRI_cor_LR <- data.frame()
for(i in 1:length(res)){
  res[[i]]$Waves <- rep(paste(rownames(mriT)[i]),nrow(res[[i]]))
}

MRI_cor_LR <- do.call(rbind,res)
MRI_cor_LR <- MRI_cor_LR[c(1,4,2,3)]
save(MRI_cor_LR,file = "processing_mri/MRI_cor_LR_within.RData")

###############################################################
###############################################################
################### Bootstrap Resected Data ###################
###############################################################
###############################################################
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
STG_GM_Thickness <- list()
BrainSeg <- list()
IpsHipVol <- list()

# Correlation analysis for each of the bootstrap. 
for(i in 1:B) {
    STG_GM_Thickness[[i]] <- FastCor(Y.b[[i]],mriT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
   	BrainSeg[[i]] <- FastCor(Y.b[[i]],mriT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    IpsHipVol[[i]] <- FastCor(Y.b[[i]],mriT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultSTG_B <- do.call(cbind,STG_GM_Thickness)
resultSTG_B <- resultSTG_B[,grep("Rho",colnames(resultSTG_B))]

resultBrainSeg_B <- do.call(cbind,BrainSeg)
resultBrainSeg_B <- resultBrainSeg_B[,grep("Rho",colnames(resultBrainSeg_B))]

resultIps_B <- do.call(cbind,IpsHipVol)
resultIps_B <- resultIps_B[,grep("Rho",colnames(resultIps_B))]

resultBoot <- rbind(resultSTG_B,resultBrainSeg_B,resultIps_B)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultBoot)){
MRI_cor_LR$BootP[i] <- sum(abs(resultBoot[i,]) >= abs(MRI_cor_LR$Rho[i]))/B
}

# Resave the file
save(MRI_cor_LR,file = "processing_mri/MRI_cor_LR_within.RData")

###############################################################
###############################################################
################### Permutation Within Data ###################
###############################################################
###############################################################
P=100  ## select number of permutation resamples
Y.b <- mclapply(1:P,mc.cores =12, function(i){apply(within_data, 2, function(col){sample(col)})}) #Randomize 100 times the within_data

# Design lists for each frequency
STG_GM_Thickness <- list()
BrainSeg <- list()
IpsHipVol <- list()

# Correlation analysis for each of the permutation 
for(i in 1:P) {
    STG_GM_Thickness[[i]] <- FastCor(Y.b[[i]],mriT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
   	BrainSeg[[i]] <- FastCor(Y.b[[i]],mriT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    IpsHipVol[[i]] <- FastCor(Y.b[[i]],mriT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultSTG_P <- do.call(cbind,STG_GM_Thickness)
resultSTG_P <- resultSTG_P[,grep("Rho",colnames(resultSTG_P))]

resultBrainSeg_P <- do.call(cbind,BrainSeg)
resultBrainSeg_P <- resultBrainSeg_P[,grep("Rho",colnames(resultBrainSeg_P))]

resultIps_P <- do.call(cbind,IpsHipVol)
resultIps_P <- resultIps_P[,grep("Rho",colnames(resultIps_P))]

resultPerm <- rbind(resultSTG_P,resultBrainSeg_P,resultIps_P)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultPerm)){
MRI_cor_LR$PermP[i] <- sum(abs(resultPerm[i,]) >= abs(MRI_cor_LR$Rho[i]))/P
}

# Resave the file
save(MRI_cor_LR,file = "processing_mri/MRI_cor_LR_within.RData")

###############################################################
###############################################################
#################### Bootstrap Frozen Data ####################
###############################################################
###############################################################

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
STG_GM_Thickness <- list()
BrainSeg <- list()
IpsHipVol <- list()

# Correlation analysis for each of the bootstrap. 
for(i in 1:B) {
    STG_GM_Thickness[[i]] <- FastCor(Y.b[[i]],mriT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    BrainSeg[[i]] <- FastCor(Y.b[[i]],mriT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    IpsHipVol[[i]] <- FastCor(Y.b[[i]],mriT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
}

# Define and combine the bootstrapped data
resultSTG_A <- do.call(cbind,STG_GM_Thickness)
resultSTG_A <- resultSTG_A[,grep("Rho",colnames(resultSTG_A))]

resultBrainSeg_A <- do.call(cbind,BrainSeg)
resultBrainSeg_A <- resultBrainSeg_A[,grep("Rho",colnames(resultBrainSeg_A))]

resultIps_A <- do.call(cbind,IpsHipVol)
resultIps_A <- resultIps_A[,grep("Rho",colnames(resultIps_A))]

resultBootAll <- rbind(resultSTG_A,resultBrainSeg_A,resultIps_A)

# Calculate the P from bootstrapped data and attach to Observed one
for (i in 1:nrow(resultBootAll)){
MRI_cor_LR$BootP_All[i] <- sum(abs(resultBootAll[i,]) >= abs(MRI_cor_LR$Rho[i]))/B
}

# Resave the file
save(MRI_cor_LR,file = "processing_mri/MRI_cor_LR_within.RData")

# Filter and save the significant genes. 
MRI_Significant_Genes <- MRI_cor_LR %>%
                              filter(
                                Pval < 0.05,
                                BootP < 0.05,
                                PermP < 0.05, 
                                BootP_All < 0.05)

write.table(MRI_Significant_Genes,"processing_mri/MRI_Significant_Genes.txt",sep="\t",quote=F)

# sessionInfo
sessionInfo()

