library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
enableWGCNAThreads()

# Brain Exp
rm(list=ls())
wgcna = list.files(pattern = 'ModuleOutput')
tab=read.table(wgcna,sep="\t",header=T)
tab <- tab[c(1,4)]
colnames(tab)=c("Gene","DEFINITION")
tab$DEFINITION <- factor(tab$DEFINITION,levels=paste("WM",1:length(unique(tab$DEFINITION)),sep=""))

Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists

files1=list.files(pattern="SME_cor_BEHAVIOR_BothCor_P005")
tmp=as.data.frame(lapply(files1,read.table,sep="\t")[[1]])
tmp <- tmp[c(1,2)]
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

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #Protein Coding bkg = 19901 Ortho bkg = 14373 Brain Expressed bkg = 15585
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
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
save(FisherMat,file="GeneSet_RECALL.RData")

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
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
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

pdf("GeneSet_RECALL_Genes.pdf",width=4,height=7)
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


FisherP[FisherP>0.05]=1
pdf("GeneSet_RECALL_Genes_NomilanP.pdf",width=4,height=7)
df=-log10(FisherP)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherP, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
xLabels = colnames(df), 
yLabels = rownames(df), 
colorLabels =FALSE,
colors=colorRampPalette(c("white", "blue"))(50),
textMatrix=LabelMat, 
setStdMargins = FALSE, 
cex.text = 0.5)
dev.off()


