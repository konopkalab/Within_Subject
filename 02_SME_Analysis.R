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

##########################
## Subsampling Analysis ##
##########################

# Create Output directories
dir.create("processing_memory/PermAnalysis")

# Load inputs
load(here("rawdata","expAdj.RData"))
load(here("rawdata","SME_Values.RData"))

# Bootstrap Resected + Frozen Data
B=100  ## select number of bootstrap resamples
index.b <- list()
Y.b <- list()
for (i in 1:B){
  set.seed(i*B+1)
  print(i)
      boot_data_all <- boot_data_all[,!(colnames(boot_data_all)%in%colnames(within_data))]
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
    deltaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[1,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame() %>% rownames_to_column('Gene')
    thetaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[2,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame() %>% rownames_to_column('Gene')
    alphaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[3,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame() %>% rownames_to_column('Gene')
    betaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[4,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame() %>% rownames_to_column('Gene')
    gammaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[5,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame() %>% rownames_to_column('Gene')
    highgammaBootAll[[i]] <- FastCor(Y.b[[i]],Lateral_resected[6,],method="spearman",alternative="two.sided",cores=12,override=TRUE) %>% as.data.frame() %>% rownames_to_column('Gene')
}

save(deltaBootAll,thetaBootAll,alphaBootAll,betaBootAll,gammaBootAll,highgammaBootAll,file="processing_memory/PermAnalysis/temporary_Perms.RData")


# Compare with WS results
load(here("processing_memory","SME_Genes.RData"))

# Delta genes compared with permuted
delta_genes <- SME_cor_LR %>% filter(Waves == "Delta", Pval < 0.05)
intercept <- nrow(delta_genes)
tmp <- map(deltaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmp, nrow))
a <- gghistogram(df, 
      x = "Perm", 
      color = "orange",
      add = "mean", rug = FALSE,
      bins = 30)+
      ylim(0,20)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      xlab("# of Genes")

# Delta genes compared with permuted
theta_genes <- SME_cor_LR %>% filter(Waves == "Theta", Pval < 0.05)
intercept <- nrow(theta_genes)
tmp <- map(thetaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmp, nrow))
b <- gghistogram(df, 
      x = "Perm", 
      color = "grey60",
      add = "mean", rug = FALSE,
      bins = 30)+
      ylim(0,20)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      xlab("# of Genes")

# Alpha genes compared with permuted
alpha_genes <- SME_cor_LR %>% filter(Waves == "Alpha", Pval < 0.05)
intercept <- nrow(alpha_genes)
tmp <- map(alphaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmp, nrow))
c <- gghistogram(df, 
      x = "Perm", 
      color = "blue",
      add = "mean", rug = FALSE,
      bins = 30)+
      ylim(0,20)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      xlab("# of Genes")

# Beta genes compared with permuted
Beta_genes <- SME_cor_LR %>% filter(Waves == "Beta", Pval < 0.05)
intercept <- nrow(Beta_genes)
tmp <- map(betaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmp, nrow))
d <- gghistogram(df, 
      x = "Perm", 
      color = "red",
      add = "mean", rug = FALSE,
      bins = 30)+
      ylim(0,20)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      xlab("# of Genes")

# Gamma genes compared with permuted
Gamma_genes <- SME_cor_LR %>% filter(Waves == "Gamma", Pval < 0.05)
intercept <- nrow(Gamma_genes)
tmp <- map(gammaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmp, nrow))
e <- gghistogram(df, 
      x = "Perm", 
      color = "magenta",
      add = "mean", rug = FALSE,
      bins = 30)+
      ylim(0,20)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      xlab("# of Genes")


# High Gamma genes compared with permuted
HGamma_genes <- SME_cor_LR %>% filter(Waves == "High.Gamma", Pval < 0.05)
intercept <- nrow(HGamma_genes)
tmp <- map(highgammaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmp, nrow))
f <- gghistogram(df, 
      x = "Perm", 
      color = "green",
      add = "mean", rug = FALSE,
      bins = 30)+
      ylim(0,20)+
      geom_vline(xintercept=intercept, colour="black",size=0.5) +
      xlab("# of Genes")


plot3by3 <- cowplot::plot_grid(a,b,c,d,e,f, ncol = 3)
cowplot::save_plot("processing_memory/PermAnalysis/Permutation.pdf", plot3by3, ncol = 2,base_height=4,base_width=4)


# Percentages
# Delta genes compared with permuted
delta_genes <- SME_cor_LR %>% filter(Waves == "Delta", Pval < 0.05)
tmp <- map(deltaBootAll, ~ (.x %>% filter(Gene %in% delta_genes$Gene,Pval < 0.05)))
tmpA <- map(deltaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmpA, nrow), Saved = sapply(tmp, nrow))
df$Perc <- (df$Saved*100)/df$Perm
a <- gghistogram(df, 
      x = "Perc", 
      color = "orange",
      rug = FALSE,
      add = "mean",
      bins = 60)+  
      xlim(-1,50)+
      ylim(0,50)+
      xlab("Percentage")

# Theta genes compared with permuted
theta_genes <- SME_cor_LR %>% filter(Waves == "Theta", Pval < 0.05)
tmp <- map(thetaBootAll, ~ (.x %>% filter(Gene %in% theta_genes$Gene,Pval < 0.05)))
tmpA <- map(thetaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmpA, nrow), Saved = sapply(tmp, nrow))
df$Perc <- (df$Saved*100)/df$Perm
b <- gghistogram(df, 
      x = "Perc", 
      color = "grey60",
      rug = FALSE,
      add = "mean",
      bins = 60)+  
      xlim(-1,50)+
      ylim(0,50)+
      xlab("Percentage")

# Alpha genes compared with permuted
alpha_genes <- SME_cor_LR %>% filter(Waves == "Alpha", Pval < 0.05)
tmp <- map(alphaBootAll, ~ (.x %>% filter(Gene %in% alpha_genes$Gene,Pval < 0.05)))
tmpA <- map(alphaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmpA, nrow), Saved = sapply(tmp, nrow))
df$Perc <- (df$Saved*100)/df$Perm
c <- gghistogram(df, 
      x = "Perc", 
      color = "blue",
      rug = FALSE,
      add = "mean",
      bins = 60)+  
      xlim(-1,50)+
      ylim(0,50)+
      xlab("Percentage")

# Beta genes compared with permuted
Beta_genes <- SME_cor_LR %>% filter(Waves == "Beta", Pval < 0.05)
tmp <- map(betaBootAll, ~ (.x %>% filter(Gene %in% Beta_genes$Gene,Pval < 0.05)))
tmpA <- map(betaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmpA, nrow), Saved = sapply(tmp, nrow))
df$Perc <- (df$Saved*100)/df$Perm
d <- gghistogram(df, 
      x = "Perc", 
      color = "red",
      rug = FALSE,
      add = "mean",
      bins = 60)+  
      xlim(-1,50)+
      ylim(0,50)+
      xlab("Percentage")

# Gamma genes compared with permuted
Gamma_genes <- SME_cor_LR %>% filter(Waves == "Gamma", Pval < 0.05)
tmp <- map(gammaBootAll, ~ (.x %>% filter(Gene %in% Gamma_genes$Gene,Pval < 0.05)))
tmpA <- map(gammaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmpA, nrow), Saved = sapply(tmp, nrow))
df$Perc <- (df$Saved*100)/df$Perm
e <- gghistogram(df, 
      x = "Perc", 
      color = "magenta",
      rug = FALSE,
      add = "mean",
      bins = 60)+  
      xlim(-1,50)+
      ylim(0,50)+
      xlab("Percentage")


# High Gamma genes compared with permuted
HGamma_genes <- SME_cor_LR %>% filter(Waves == "High.Gamma", Pval < 0.05)
tmp <- map(highgammaBootAll, ~ (.x %>% filter(Gene %in% HGamma_genes$Gene,Pval < 0.05)))
tmpA <- map(highgammaBootAll, ~ (.x %>% filter(Pval < 0.05)))
df <- data.frame(Perm = sapply(tmpA, nrow), Saved = sapply(tmp, nrow))
df$Perc <- (df$Saved*100)/df$Perm
f <- gghistogram(df, 
      x = "Perc", 
      color = "green",
      rug = FALSE,
      add = "mean",
      bins = 60)+  
      xlim(-1,50)+
      ylim(0,50)+
      xlab("Percentage")

plot3by3 <- cowplot::plot_grid(a,b,c,d,e,f, ncol = 3)
cowplot::save_plot("processing_memory/PermAnalysis/Permutation_Percentage.pdf", plot3by3, ncol = 2,base_height=4,base_width=3)

# sessionInfo
sessionInfo()
