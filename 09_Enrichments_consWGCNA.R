suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
source("Utils.R")
enableWGCNAThreads()

dir.create("enrichments_wgcna")

#######################
# Cer Cort enrichment #
#######################
tab <- read.table("processing_wgcna/ModuleOutput_WITHIN.txt",sep="\t",header=T)
tab <- tab %>%
        select(Gene,ModuleName) %>%
        rename(DEFINITION = ModuleName)
tab$DEFINITION <- factor(tab$DEFINITION,levels=paste("WM",1:length(unique(tab$DEFINITION)),sep=""))
Genes=as.data.frame(table(tab$DEFINITION))

# Cer Cort genes
GeneSets <- load(here("rawdata","geneset","GeneSets_CerCort.RData")) %>% get()
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
save(FisherMat,file="enrichments_wgcna/FisherMat_BertoEtal2017.RData")

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


p <- as.data.frame(FisherAdj)
cor <- as.data.frame(FisherOR)


p $Module <- rownames(p )
cor$Module <- rownames(cor)

p <- melt(p)
cor <- melt(cor)
p$cor <- cor$value
p$log <- -log10(p$value)
p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))

p$log[p$log < 1.3] <- NA
p$abs <- abs(p$cor)
p$OR<- ifelse(is.na(p$log), p$log, p$abs)

pdf("enrichments_wgcna/Modules_CerCort_Bubble.pdf",width=8,height=2)
ggscatter(p, 
            x = "variable",
            y = "Module",
            size="OR",
            color="log",
            alpha = 0.8,
            xlab = "",
            ylab = "",) +
            theme_minimal() +
            rotate_x_text(angle = 45)+
        coord_flip()+
        scale_size(range = c(2, 10))+ 
        gradient_color(c("green","darkgreen"))
dev.off()


###############################
# psychEncode DEGs enrichment #
###############################
tab <- read.table("processing_wgcna/ModuleOutput_WITHIN.txt",sep="\t",header=T)
tab <- tab %>%
        select(Gene,ModuleName) %>%
        rename(DEFINITION = ModuleName)
tab$DEFINITION <- factor(tab$DEFINITION,levels=paste("WM",1:length(unique(tab$DEFINITION)),sep=""))
Genes=as.data.frame(table(tab$DEFINITION))

# Cer Cort genes
GeneSets <- load(here("rawdata","geneset","PsychENCODE_DEGs.RData")) %>% get()
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
save(FisherMat,file="enrichments_wgcna/psychENCODE_EnrichMat.RData")

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


p <- as.data.frame(FisherAdj)
cor <- as.data.frame(FisherOR)


p $Module <- rownames(p )
cor$Module <- rownames(cor)

p <- melt(p)
cor <- melt(cor)
p$cor <- cor$value
p$log <- -log10(p$value)
p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))

p$log[p$log < 1.3] <- NA
p$abs <- abs(p$cor)
p$OR<- ifelse(is.na(p$log), p$log, p$abs)

pdf("enrichments_wgcna/Modules_psychENCODE_Bubble.pdf",width=8,height=2)
ggscatter(p, 
            x = "variable",
            y = "Module",
            size="OR",
            color="log",
            alpha = 0.8,
            xlab = "",
            ylab = "",) +
            theme_minimal() +
            rotate_x_text(angle = 45)+
        coord_flip()+
        scale_size(range = c(2, 10))+ 
        gradient_color(c("steelblue","darkblue"))
dev.off()


###############################
# psychEncode DEGs enrichment #
###############################
tab <- read.table("processing_wgcna/ModuleOutput_WITHIN.txt",sep="\t",header=T)
tab <- tab %>%
        select(Gene,ModuleName) %>%
        rename(DEFINITION = ModuleName)
tab$DEFINITION <- factor(tab$DEFINITION,levels=paste("WM",1:length(unique(tab$DEFINITION)),sep=""))
Genes=as.data.frame(table(tab$DEFINITION))

# Cer Cort genes
GeneSets <- load(here("rawdata","geneset","PsychENCODE_DEGs.RData")) %>% get()
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
save(FisherMat,file="enrichments_wgcna/psychENCODE_EnrichMat.RData")

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


p <- as.data.frame(FisherAdj)
cor <- as.data.frame(FisherOR)


p $Module <- rownames(p )
cor$Module <- rownames(cor)

p <- melt(p)
cor <- melt(cor)
p$cor <- cor$value
p$log <- -log10(p$value)
p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))

p$log[p$log < 1.3] <- NA
p$abs <- abs(p$cor)
p$OR<- ifelse(is.na(p$log), p$log, p$abs)

pdf("enrichments_wgcna/Modules_psychENCODE_Bubble.pdf",width=8,height=2)
ggscatter(p, 
            x = "variable",
            y = "Module",
            size="OR",
            color="log",
            alpha = 0.8,
            xlab = "",
            ylab = "",) +
            theme_minimal() +
            rotate_x_text(angle = 45)+
        coord_flip()+
        scale_size(range = c(2, 10))+ 
        gradient_color(c("steelblue","darkblue"))
dev.off()


##################################
# psychEncode Modules enrichment #
##################################
tab <- read.table("processing_wgcna/ModuleOutput_WITHIN.txt",sep="\t",header=T)
tab <- tab %>%
        select(Gene,ModuleName) %>%
        rename(DEFINITION = ModuleName)
tab$DEFINITION <- factor(tab$DEFINITION,levels=paste("WM",1:length(unique(tab$DEFINITION)),sep=""))
Genes=as.data.frame(table(tab$DEFINITION))

# Cer Cort genes
GeneSets <- load(here("rawdata","geneset","PsychEncode_Modules.RData")) %>% get()
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
save(FisherMat,file="enrichments_wgcna/psychENCODE_Modules_EnrichMat.RData")

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


vec <- paste("geneM",0:34,sep="")
FisherAdj <- FisherAdj[,match(vec,colnames(FisherAdj))]
FisherP <- FisherP[,match(vec,colnames(FisherP))]
FisherOR <- FisherOR[,match(vec,colnames(FisherOR))]

FisherP[FisherP>0.05]=1
FisherOR[FisherOR < 1]=0

p <- as.data.frame(FisherAdj)
cor <- as.data.frame(FisherOR)


p $Module <- rownames(p )
cor$Module <- rownames(cor)

p <- melt(p)
cor <- melt(cor)
p$cor <- cor$value
p$log <- -log10(p$value)
p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))

p$log[p$log < 1.3] <- NA
p$abs <- abs(p$cor)
p$OR<- ifelse(is.na(p$log), p$log, p$abs)

pdf("enrichments_wgcna/Modules_psychENCODE_Modules_Bubble.pdf",width=8,height=5)
ggscatter(p, 
            x = "variable",
            y = "Module",
            size="OR",
            color="log",
            alpha = 0.8,
            xlab = "",
            ylab = "",) +
            theme_minimal() +
            rotate_x_text(angle = 45)+
        coord_flip()+
        scale_size(range = c(2, 10))+ 
        gradient_color(c("magenta","purple"))
dev.off()
