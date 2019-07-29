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

dir.create("processing_behavior")

# Load the input expression and filter for more than 0.5 for all samples
load(here("rawdata","expAdj.RData"))
load(here("rawdata","BEHAVIORAL_Values.RData"))

# Correlation ANALYSIS
res <- list()
for(i in 1:nrow(behavT))
{
res[[i]] <- FastCor(within_data,behavT[i,],method="spearman",alternative="two.sided",cores=10,override=TRUE) %>%
                      as.data.frame() %>%
                      rownames_to_column('Gene')
}

BEHAVIORAL_cor_LR <- data.frame()
for(i in 1:length(res)){
  res[[i]]$Waves <- rep(paste(rownames(behavT)[i]),nrow(res[[i]]))
}

BEHAVIORAL_cor_LR <- do.call(rbind,res)
BEHAVIORAL_cor_LR <- BEHAVIORAL_cor_LR[c(3,4,1,2)]
save(BEHAVIORAL_cor_LR,file = "processing_behavior/BEHAVIOR_Genes.RData")

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
    PercRec[[i]] <- FastCor(Y.b[[i]],behavT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    EncWord[[i]] <- FastCor(Y.b[[i]],behavT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    Sessions[[i]] <- FastCor(Y.b[[i]],behavT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
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
BEHAVIORAL_cor_LR$BootP[i] <- sum(abs(resultBoot[i,]) >= abs(BEHAVIORAL_cor_LR$Rho[i]))/B
}

# Resave the file
save(BEHAVIORAL_cor_LR,file = "processing_behavior/BEHAVIOR_Genes.RData")

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
    PercRec[[i]] <- FastCor(Y.b[[i]],behavT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    EncWord[[i]] <- FastCor(Y.b[[i]],behavT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    Sessions[[i]] <- FastCor(Y.b[[i]],behavT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
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
BEHAVIORAL_cor_LR$PermP[i] <- sum(abs(resultPerm[i,]) >= abs(BEHAVIORAL_cor_LR$Rho[i]))/P
}

# Resave the file
save(BEHAVIORAL_cor_LR,file = "processing_behavior/BEHAVIOR_Genes.RData")

####################################
# Bootstrap Resected + Frozen Data #
####################################
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
PercRec <- list()
EncWord <- list()
Sessions <- list()

# Correlation analysis for each of the bootstrap. 
for(i in 1:B) {
    PercRec[[i]] <- FastCor(Y.b[[i]],behavT[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    EncWord[[i]] <- FastCor(Y.b[[i]],behavT[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
    Sessions[[i]] <- FastCor(Y.b[[i]],behavT[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame()
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
BEHAVIORAL_cor_LR$BootP_All[i] <- sum(abs(resultBootAll[i,]) >= abs(BEHAVIORAL_cor_LR$Rho[i]))/B
}

# Resave the file
save(BEHAVIORAL_cor_LR,file = "processing_behavior/BEHAVIOR_Genes.RData")

# Filter and save the significant genes. 
BEHAVIOR_Significant_Genes <- BEHAVIORAL_cor_LR %>%
                              filter(
                                Pval < 0.05,
                                BootP < 0.05,
                                PermP < 0.05, 
                                BootP_All < 0.05)
                              
write.table(BEHAVIOR_Significant_Genes,"processing_behavior/BEHAVIOR_Significant_Genes.txt",sep="\t",quote=F)


