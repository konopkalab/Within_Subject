suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
source("UTILS/Utils.R")
enableWGCNAThreads()

dir.create("enrichments_SME")

##############################
# Behavioral task enrichment #
##############################
load("processing_memory_QuantReg/SME_QuantReg.RData")
tab <- SME_Cor_QuantReg_Sign
tab <- tab %>%
        select(Gene,Waves) %>%
        rename(DEFINITION = Waves) %>%
        droplevels()
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists

GeneSets <- split(tab,tab$DEFINITION)
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)
for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq  
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file="enrichments_SME/SMEvsSME_Cross-Enrichment.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,7]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='BH') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherAdj[FisherAdj == 0] <- 2.520179e-176
FisherOR[!(is.finite(FisherOR))] <- 200
FisherOR[FisherOR < 1]=0

pdf("enrichments_SME/SMEvsSME_Cross-Enrichment.pdf",width=6,height=6)
vec1 <- c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta")
FisherAdj <- FisherAdj[match(vec1,rownames(FisherAdj)),match(rev(vec1),colnames(FisherAdj)),drop=FALSE]
FisherOR <- FisherOR[match(vec1,rownames(FisherOR)),match(rev(vec1),colnames(FisherOR)),drop=FALSE]
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 3), "\n(",signif(FisherAdj, 2), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "red"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5)
dev.off()

FisherOR[FisherOR < 1]=0
pdf("enrichments_SME/SMEvsBehavior_Cross-Enrichment.pdf",width=3,height=3)
vec1 <- c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta")
FisherAdj <- FisherAdj[match(vec1,rownames(FisherAdj)),,drop=FALSE]
FisherOR <- FisherOR[match(vec1,rownames(FisherOR)),,drop=FALSE]
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "red"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5)
dev.off()









##############################
# Behavioral task enrichment #
##############################
load("processing_memory_QuantReg/SME_QuantReg.RData")
tab <- SME_Cor_QuantReg_Sign
tab <- tab %>%
        select(Gene,Waves) %>%
        rename(DEFINITION = Waves) %>%
        droplevels()
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
tmp <- read.table("additional_data/BEHAVIOR_Significant_Genes.txt",sep="\t",header=T)
tmp <- tmp %>%
		select(Gene,Waves) %>%
		filter(Waves == "Perc_recall") %>%
		droplevels()

GeneSets <- split(tmp,tmp$Waves)
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)
for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq  
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file="enrichments_SME/SMEvsBehavior_Cross-Enrichment.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,7]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='BH') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0
pdf("enrichments_SME/SMEvsBehavior_Cross-Enrichment.pdf",width=3,height=3)
vec1 <- c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta")
FisherAdj <- FisherAdj[match(vec1,rownames(FisherAdj)),,drop=FALSE]
FisherOR <- FisherOR[match(vec1,rownames(FisherOR)),,drop=FALSE]
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "red"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5)
dev.off()

########################
# MATH task enrichment #
########################
load("processing_memory_QuantReg/SME_QuantReg.RData")
tab <- SME_Cor_QuantReg_Sign
tab <- tab %>%
        select(Gene,Waves) %>%
        rename(DEFINITION = Waves) %>%
        droplevels()
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes=as.data.frame(table(tab$DEFINITION))


# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
tmp <- read.table("additional_data/MATH_Significant_Genes.txt",sep="\t",header=T)
tmp <- tmp %>%
    select(Gene,Waves) 
tmp$Waves <- paste(tmp$Waves,"_MATH",sep="")
GeneSets <- split(tmp,tmp$Waves)
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)
for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq  
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file="enrichments_SME/SMEvsMath_Cross-Enrichment.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,7]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='BH') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0
pdf("enrichments_SME/SMEvsMath_Cross-Enrichment.pdf",width=4,height=4)
vec1 <- c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta")
vec2 <- c("High.Gamma_MATH","Gamma_MATH","Beta_MATH","Alpha_MATH","Theta_MATH","Delta_MATH")
FisherAdj <- FisherAdj[match(vec1,rownames(FisherAdj)),match(rev(vec2),colnames(FisherAdj))]
FisherOR <- FisherOR[match(vec1,rownames(FisherOR)),match(rev(vec2),colnames(FisherOR))]
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "red"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5)
dev.off()

############################
# MRI Thickness Enrichment #
############################
load("processing_memory_QuantReg/SME_QuantReg.RData")
tab <- SME_Cor_QuantReg_Sign
tab <- tab %>%
        select(Gene,Waves) %>%
        rename(DEFINITION = Waves) %>%
        droplevels()
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes=as.data.frame(table(tab$DEFINITION))


# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
tmp <- read.table("additional_data/MRI_Significant_Genes.txt",sep="\t",header=T)
tmp <- tmp %>%
    select(Gene,Waves) %>%
    filter(Waves == "STG_GM_Thickness_mm") %>%
    droplevels()

GeneSets <- split(tmp,tmp$Waves)
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)
for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq  
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file="enrichments_SME/SMEvsMRI_Cross-Enrichment.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,7]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='BH') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0
pdf("enrichments_SME/SMEvsMRI_Cross-Enrichment.pdf",width=3,height=3)
vec1 <- c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta")
FisherAdj <- FisherAdj[match(vec1,rownames(FisherAdj)),,drop=FALSE]
FisherOR <- FisherOR[match(vec1,rownames(FisherOR)),,drop=FALSE]
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "red"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5)
dev.off()
