# Load libraries
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(quantreg))
source("UTILS/Utils.R")

# Create Output directories
dir.create("processing_memory_QuantReg")

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
pheno$Batch <- as.factor(pheno$Batch) # same PMI for resected tissues
pheno <- pheno[match(colnames(logCPM),rownames(pheno)),]
pheno <- droplevels(pheno)

# Quantile normalization
p <- normalize.quantiles(as.matrix(logCPM))
rownames(p) <- rownames(logCPM)
colnames(p) <- colnames(logCPM)

# Load SME data and filter for matched data
load(here("rawdata","SME_Values.RData"))
pheno.ws <- na.omit(pheno)
pheno.ws <- pheno.ws[match(colnames(Lateral_resected),pheno.ws$UT.code),]
within_data <- p[,colnames(p) %in% rownames(pheno.ws)]
within_data <- within_data[,match(rownames(pheno.ws),colnames(within_data))]
pheno.ws$UT.code <- NULL


## Quantile Regression
plan("multiprocess", workers=3)

# Create the models
model <- 'geneExpr ~ Wave + Age + Sex + Race + Ethnicity + EpDur + RIN + Hemisphere + Batch' # Full model
modelNull <- 'geneExpr ~ Age + Sex + Race + Ethnicity + EpDur + RIN + Hemisphere + Batch' # Null model

# Function for fitting the two models with/without the "diagnosis".
#rho <- function(u,tau=.5)u*(tau - (u < 0))

fit_rq <- function(vectorizedExpression) {
  tmpMetaData <- cbind(metadata, data.frame(geneExpr = unname(vectorizedExpression)))
  residuals1 <- rq(model,data=tmpMetaData, tau=0.5)
  residuals2 <- rq(modelNull,data=tmpMetaData, tau=0.5)
  #Vrq <- 1 - residuals2$rho/residuals1$rho
  pval <- lrtest(residuals1, residuals2)$`Pr(>Chisq)`[2]
  effect_size <- as.numeric(summary(residuals1)$coefficients[2])
  c(EffSize_QR = effect_size, Pval_QR = pval)
}

expression_fit_rq <- function(vectorizedExpression) {
  tryCatch(fit_rq(vectorizedExpression), error = function(e) c(EffSize_QR = NA,Pval_QR = NA))
}

# 
# obtaining a p-value below the chosen significance level from these two tests indicates sufficient 
# evidence in favor of rejecting the null hypothesis claiming the two models are equivalent.
res <- list()
for(i in 1:nrow(Lateral_resected))
{ 
metadata <- pheno.ws
metadata$Wave <- as.numeric(Lateral_resected[i,])

res[[i]] <- future_apply(within_data, 1, expression_fit_rq) %>%
              t() %>%
              as.data.frame() %>% 
              rownames_to_column("Gene") %>%
              mutate(FDR_QR = p.adjust(Pval_QR,method="BH"), BONF_QR = p.adjust(Pval_QR,method="bonferroni")) %>%
              mutate(Waves = rep(paste(rownames(Lateral_resected)[i]), nrow(within_data)))
}


SME_rq_LR <- do.call(rbind, res)

SME_rq_LR <- SME_rq_LR %>% 
              select(Gene, Waves,EffSize_QR, Pval_QR, FDR_QR,BONF_QR)

# Combine the data with the spearman analysis. 
load(here("processing_memory","SME_Genes.RData"))

SME_Cor_QuantReg <- Reduce(dplyr::full_join, list(SME_cor_LR, SME_rq_LR))

openxlsx::write.xlsx(SME_Cor_QuantReg, 
                     file = "processing_memory_QuantReg/SME_Cor_QuantReg.xlsx", 
                     colNames = TRUE, 
                     borders = "columns",
                     sheetName="Full Table")

pdf(file="processing_memory_QuantReg/GeneCorrelations_ScatterCor_QuantReg.pdf",width=8,height=6)
ggscatter(SME_Cor_QuantReg, x = "Rho", y = "EffSize_QR",
   color = "lightgray", size = 0.5, 
   add = "reg.line",  
   add.params = list(color = "red", fill = "blue"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.sep = "\n"),
   ellipse = FALSE)+ 
xlab("Original Rho") + ylab("Effect Size Rq") +
geom_hline(yintercept = 0, linetype="dotted", color = "black", size=1)+ 
geom_vline(xintercept = 0, linetype="dotted", color = "black", size=1)+
facet_wrap(~Waves,scales="free")
dev.off()

SME_Cor_QuantReg_Sign <- SME_Cor_QuantReg %>% 
                          filter(Pval < 0.05,BootP < 0.05, PermP < 0.05, BootP_All < 0.05, FDR_QR < 0.05)


openxlsx::write.xlsx(SME_Cor_QuantReg_Sign, 
                     file = "processing_memory_QuantReg/SME_Cor_QuantReg_FDR.xlsx", 
                     colNames = TRUE, 
                     borders = "columns",
                     sheetName="Sign Table")

# Check the concordance of the estimates
pdf(file="processing_memory_QuantReg/GeneCorrelations_SMEgenes_ScatterCor_QuantReg.pdf",width=8,height=6)
ggscatter(SME_Cor_QuantReg_Sign, x = "Rho", y = "EffSize_QR",
   color = "lightgray", size = 0.5, 
   add = "reg.line",  
   add.params = list(color = "red", fill = "blue"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.sep = "\n"),
   ellipse = FALSE)+ 
xlab("Original Rho") + ylab("Effect Size Rq") +
geom_hline(yintercept = 0, linetype="dotted", color = "black", size=1)+ 
geom_vline(xintercept = 0, linetype="dotted", color = "black", size=1)+
facet_wrap(~Waves,scales="free")
dev.off()

# Save the data
save(SME_rq_LR,SME_Cor_QuantReg,SME_Cor_QuantReg_Sign, file = "processing_memory_QuantReg/SME_QuantReg.RData")


# New barplot
tmp <- SME_Cor_QuantReg %>% 
              filter(Pval < 0.05,BootP < 0.05, PermP < 0.05, BootP_All < 0.05)

tmp <- tmp %>%
        mutate(QR_Saved = case_when(FDR_QR < 0.05 ~ "Sign", FDR_QR > 0.05  ~ "NotSign"))

pdf("processing_memory_QuantReg/Barplot_SME_QR_Saved.pdf",width=4,height=3,useDingbats=FALSE)
df <- reshape2::melt(table(tmp$Waves,tmp$QR_Saved))
df$Var1 <- factor(df$Var1,levels = c("Delta","Theta","Alpha","Beta","Gamma","High.Gamma"))
ggbarplot(df, 
  "Var1", 
  "value",
     fill = "Var2", 
     color = "Var2", 
     palette = c("red","black"),
     position = position_dodge()) + 
  xlab("") + 
  ylab("# of Observed SME genes")+
  theme_classic()+
  theme(legend.position="right")+
  rotate_x_text(angle = 45)
dev.off()


# New barplot 2 
sign <- tmp %>% filter(FDR_QR < 0.05)
df <- as.data.frame(table(sign$Waves))
df$Var1 <- factor(df$Var1,levels = c("Delta","Theta","Alpha","Beta","Gamma","High.Gamma"))

pdf("processing_memory_QuantReg/Barplot_SME_QR_FDR005.pdf",width=4,height=4,useDingbats=FALSE)
ggbarplot(df, x = "Var1", y = "Freq",
          fill = "Var1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c("orange","grey60","blue","red","magenta","green"),
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45,           # Rotate vertically x axis texts
          label = TRUE, label.pos = "out")+
ggtitle("Genes FDR < 0.05")+
xlab("")+ 
ylab("# of Genes")+
theme(legend.position="none")+
ylim(0,1000)
dev.off()

# New barplot 3
sign <- tmp %>% filter(FDR_QR < 0.01)
df <- as.data.frame(table(sign$Waves))
df$Var1 <- factor(df$Var1,levels = c("Delta","Theta","Alpha","Beta","Gamma","High.Gamma"))

pdf("processing_memory_QuantReg/Barplot_SME_QR_FDR001.pdf",width=4,height=4,useDingbats=FALSE)
ggbarplot(df, x = "Var1", y = "Freq",
          fill = "Var1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c("orange","grey60","blue","red","magenta","green"),
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45,           # Rotate vertically x axis texts
          label = TRUE, label.pos = "out")+
ggtitle("Genes FDR < 0.05")+
xlab("")+ 
ylab("# of Genes")+
theme(legend.position="none")+
ylim(0,1000)
dev.off()

# sessionInfo
sessionInfo()
