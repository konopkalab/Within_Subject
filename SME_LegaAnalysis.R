
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
source("BertoEtal_SME_Functions.R")

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

# Distribution before quantile
pdf("PLOTS/Distribution_logCPM_LEGA.pdf",width=6,height=5)
colramp = colorRampPalette(c(3,"white",2))(ncol(logCPM))
plot(density(logCPM[,1]),col=colramp[1],lwd=3,ylim=c(0,.3))
for(i in 2:ncol(logCPM)){lines(density(logCPM[,i]),lwd=3,col=colramp[i])}
dev.off()

# Quantile normalization
p <- normalize.quantiles(as.matrix(logCPM))
rownames(p) <- rownames(logCPM)
colnames(p) <- colnames(logCPM)
write.table(p,"OUTPUTS/LogCPM_LEGA_QuantNorm.txt",sep="\t",quote=F)

# get confounder for lm regression
covariates <- pheno[-c(9)]

# Variance explained
var <- VarExp(p,covariates,5,FALSE)
pdf("PLOTS/Variance_Explained_Lega.pdf",width=8,height=6)
plotVarExp(var,"Variance Explained")
dev.off()

# Adjust expression
avebeta.lm<-lapply(1:nrow(p), function(x){
  lm(unlist(p[x,])~., data=covariates)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
adj.residuals<-residuals+matrix(apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(adj.residuals)<-rownames(p)
rownames(residuals) <- rownames(p)
write.table(adj.residuals, "OUTPUTS/logCPM_LEGA_QuantNorm_LMreg.txt",sep="\t",quote=F)

# Load SME data 
load("RAW_DATA/SME_Values.RData")

# Check the UTcodes and subset data for that
pheno <- na.omit(pheno)
pheno.tmp <- pheno[match(colnames(Lateral_resected),pheno$UT.code),]
within_data <- adj.residuals[,colnames(adj.residuals) %in% rownames(pheno.tmp)]
within_data <- within_data[,match(rownames(pheno.tmp),colnames(within_data))]
write.table(within_data, "OUTPUTS/logCPM_LEGA_QuantNorm_LMreg_WithinData.txt",sep="\t",quote=F)
save(within_data, Lateral_resected,pheno.tmp,file="OUTPUTS/INPUT_DATA_SME.RData")

# Dummy data for bootstrap
boot_data <- adj.residuals

# Correlation ANALYSIS
res <- list()
for(i in 1:nrow(Lateral_resected))
{
res[[i]] <- FastCor(within_data,Lateral_resected[i,],method="spearman",alternative="two.sided",cores=10,override=TRUE)
}

SME_cor_LR <- data.frame()
for(i in 1:length(res)){
	res[[i]] <- as.data.frame(res[[i]])
	res[[i]]$Gene <- rownames(res[[i]])
	rownames(res[[i]]) <- NULL
	res[[i]]$Waves <- rep(paste(rownames(Lateral_resected)[i]),nrow(res[[i]]))
}

SME_cor_LR <- do.call(rbind,res)
SME_cor_LR <- SME_cor_LR[c(3,4,1,2)]
save(SME_cor_LR,file = "OUTPUTS/SME_cor_LR_within.RData")

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
save(SME_cor_LR,file = "OUTPUTS/SME_cor_LR_within.RData")

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
save(SME_cor_LR,file = "OUTPUTS/SME_cor_LR_within.RData")

# Bootstrap Resected + Frozen Data
# Select 16 subject randomly 100 times
# Recalculate correlation
# Compare with the observed one

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

# Adjust expression
avebeta.lm<-lapply(1:nrow(pAll), function(x){
  lm(unlist(pAll[x,])~., data=covariatesAll)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
adj.residuals_All<-residuals+matrix(apply(pAll, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(adj.residuals_All)<-rownames(pAll)
rownames(residuals) <- rownames(pAll)
write.table(adj.residuals_All, "OUTPUTS/logCPM_ALL_QuantNorm_LMreg.txt",sep="\t",quote=F)
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
save(SME_cor_LR,file = "OUTPUTS/SME_cor_LR_within.RData")

# Filter and save the significant genes. 
SME_cor_LR_sign <- SME_cor_LR[SME_cor_LR$Pval < 0.05 & SME_cor_LR$BootP < 0.05 & SME_cor_LR$PermP < 0.05 & SME_cor_LR$BootP_All < 0.05,]
write.table(SME_cor_LR_sign,"OUTPUTS/SME_cor_LR_BothCor_P005.txt",sep="\t",quote=F)


