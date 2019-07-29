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
tmp <- list.files(pattern = 'MATH_cor')
tmp=read.table(tmp,sep="\t",header=T)
tmp <- tmp[c(1,2)]
tmp$Waves <- paste(tmp$Waves,"_MATH",sep="")
GeneSets <- split(tmp,tmp$Waves)
#for(i in 1:length(GeneSets))
#{
#	GeneSets[[i]] <- GeneSets[[i]][GeneSets[[i]]$Gene %in% tab$Gene,]
#}
#GeneSets <- GeneSets[1:8]
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
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585
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
save(FisherMat,file="Module_Enrich_MATH.RData")

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
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0
pdf("Heatmap_WM_MATH.pdf",width=4,height=7)
vec2 <- c("Delta_MATH","Theta_MATH","Alpha_MATH","Beta_MATH","Gamma_MATH","High.Gamma_MATH")
FisherAdj <- FisherAdj[,match(vec2,colnames(FisherAdj))]
FisherOR <- FisherOR[,match(vec2,colnames(FisherOR))]
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
pdf("Heatmap_WM_MATH_NomilanP.pdf",width=4,height=7)
FisherP <- FisherP[,match(vec2,colnames(FisherP))]
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

