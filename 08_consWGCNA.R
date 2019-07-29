# WGCNA

suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors = FALSE);
enableWGCNAThreads()


dir.create("processing_wgcna")
# Load tables 
load(here("rawdata","expAdj.RData"))
load(here("rawdata","SME_Values.RData"))

datExpr <- as.data.frame(t(within_data));
names(datExpr) <- rownames(within_data);
rownames(datExpr) <- names(within_data);

datAll <- as.data.frame(t(boot_data));
names(datAll) <- rownames(boot_data);
rownames(datAll) <- names(boot_data);

datFull <- as.data.frame(t(boot_data_all));
names(datFull) <- rownames(boot_data_all);
rownames(datFull) <- names(boot_data_all);


datTraits=as.data.frame(t(Lateral_resected))
colnames(datTraits) = paste(colnames(datTraits),"_LR",sep="")


## Powers analysis
powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr,powerVector=powers,verbose = 5, blockSize= 16000, networkType = "signed",RsquaredCut = 0.85) 
pdf("processing_wgcna/SoftThresholdingPower_signed.pdf")
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

######################################################################################################################
############################################ MODULE construction #####################################################
############################################    signed network   #####################################################
PWR=sft$powerEstimate

nSets = 3
multiExpr<-vector(mode="list",length=nSets)
multiExpr=list(Within=list(data=datExpr), All=list(data=datAll),Full=list(data=datFull))
net = blockwiseConsensusModules(multiExpr,
corType="bicor",
maxBlockSize = 16000,
networkType="signed",
minCoreKME = 0.5, 
minKMEtoStay = 0.6,
power=PWR, 
checkMissingData = TRUE,
minModuleSize=75,
nThreads=15,
TOMType = "signed",
TOMDenom = "mean",
deepSplit=4,
verbose=1,
mergeCutHeight=0.15,
reassignThreshold = 1e-10,
numericLabels=TRUE,
saveIndividualTOMs=TRUE,
saveConsensusTOMs=TRUE)
save(net,file="processing_wgcna/WITHIN_ConsNet.RData")


load("processing_wgcna/consensusTOM-block.1.RData")
consensusTOM_final <- consTomDS 
geneTree = hclust(1-as.dist(consensusTOM_final), method="average")
consMEs = net$multiMEs;
consTree = net$dendrograms[[1]];
moduleLabels = net$colors;
moduleColors = labels2colors(moduleLabels)
write.table(moduleColors, "WITHIN_ConsColors.txt",sep="\t",quote=F)

pdf("processing_wgcna/NetworkDendrogram_WITHIN.pdf",width=10,height=3)
#plotDendroAndColors(geneTree, moduleColors, colorHeight = 0.1, dendroLabels=FALSE,hang = 0.03, addGuide = TRUE, guideHang = 0.05 )
#dev.off()
x=cor(datTraits,datExpr,method="pearson")
x[x > 0.25]=6
x[x < -0.25]=30
x[ x > -0.25 & x < 0.25] = 27
listLabels2 <- t(labels2colors(x))
rownames(listLabels2) = colnames(x) 
colnames(listLabels2) = rownames(x)
plotColors=as.data.frame(cbind(moduleColors,listLabels2));
colnames(plotColors)=c("Module",names(datTraits)) 
plotDendroAndColors(geneTree, plotColors,colorHeight = 0.1, dendroLabels=FALSE,hang = 0.03, addGuide = TRUE, guideHang = 0.05 )
dev.off()


#KMEs
KMEs<-consensusKME(multiExpr, moduleColors,signed = TRUE,)
kme=data.frame(rownames(within_data), moduleColors, KMEs)
colnames(kme)[1]="Symbol"
rownames(kme)=NULL
write.table(kme,"processing_wgcna/ConsKME_WITHIN.txt",sep="\t",quote=F)

# Dat Trait Heatmap
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColors,softPower = PWR,impute = TRUE)$eigengenes
MEs0$MEgrey=NULL
modTraitCor= cor(MEs0, datTraits,method="spearman")
write.table(modTraitCor,"modTraitCor_WITHIN.txt",sep="\t",quote=F)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
write.table(modTraitP,"modTraitP_WITHIN.txt",sep="\t",quote=F)
textMatrix = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor) 
par(mar = c(6, 8.5, 3, 3))
pdf("processing_wgcna/Heatmap_DatTraits.pdf",width=7,height=7)
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), yLabels = names(MEs0), ySymbols = names(MEs0), colorLabels =FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix, setStdMargins = FALSE, cex.text = 0.4, zlim = c(-1,1), main = paste("Module Association"))
dev.off()

#EigenGeneNet
MET=orderMEs(MEs0)
pdf("processing_wgcna/EigengeneNetworks.pdf")
plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.5,xLabelsAngle=90)
dev.off()

#GGplotInput
MEs0$Rows=colnames(within_data)
write.table(MEs0, "processing_wgcna/Matrix_module_correlation.txt",sep="\t",quote=F)

# Adjacency matrix
Adj = adjacency(datExpr, power = PWR,type="signed",corFnc = "bicor")
moduleOutput <- data.frame(rownames(within_data))
moduleOutput[,2]<- as.factor(moduleColors)
intraCon <- intramodularConnectivity(Adj, moduleColors)
moduleOutput[,3]<-intraCon$kWithin
colnames(moduleOutput) <- c("Gene", "ModuleColor", "kWithin")
moduleOutput <- droplevels(moduleOutput[!(moduleOutput$ModuleColor=="grey"),])
moduleOutput$ModuleName <- plyr::mapvalues(moduleOutput$ModuleColor, as.character(sort(unique(moduleOutput$ModuleColor))), paste("W","M",1:length(levels(moduleOutput$ModuleColor)),sep=""))
moduleOutput <- moduleOutput[order(moduleOutput$ModuleColor),]
write.table(moduleOutput, "processing_wgcna/ModuleOutput_WITHIN.txt", sep="\t", quote=F,row.names=FALSE)

# TO connectivity (you need the table as single gene list column in the directory)
TOM = TOMsimilarityFromExpr(datExpr, power= PWR,corType = "bicor",networkType="signed",TOMType="signed",TOMDenom = "mean",nThreads = 15,verbose = 5, indent = 0)
colnames(TOM)=rownames(TOM)=colnames(datExpr)
save(TOM,file="processing_wgcna/TOM_SIGNED.RData")
Connectivity=apply(TOM,1,sum)
save(Connectivity,file="processing_wgcna/Connectivity.RData")

# CytoScape output
dir.create("processing_wgcna/Cyto")
for(module in unique(moduleColors)){
inModule <- is.finite(match(moduleColors, module))
modTOM <- TOM[inModule, inModule]
cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste("processing_wgcna/Cyto/CytoEdge",paste(module,collapse="-"),".txt",sep=""), nodeFile=paste("processing_wgcna/Cyto/CytoNode",paste(module,collapse="-"),".txt",sep=""), weighted = TRUE, threshold = 0, nodeAttr = moduleColors[inModule], nodeNames = names(datExpr)[inModule])
}


# Making bubbly Traits
p <- read.table("processing_wgcna/modTraitP_WITHIN.txt")
rownames(p) <- gsub("ME","",rownames(p))
cor <- read.table("processing_wgcna/modTraitCor_WITHIN.txt")
rownames(cor) <- gsub("ME","",rownames(cor))

rownames(p) <- paste("W","M",1:nrow(p),sep="")
rownames(cor) <- paste("W","M",1:nrow(cor),sep="")
p$Module <- rownames(p)
cor$Module <- rownames(cor)

p <- melt(p)
cor <- melt(cor)
p$cor <- cor$value
p$log <- -log10(p$value)
p$Direction <- ifelse(p$cor > 0,"Pos","Neg")
p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))
p$variable <- gsub("_LR","",p$variable)

p$log[p$log < 1.3] <- NA
p$abs <- abs(p$cor)
p$abs[p$abs < 0.5] <- NA

pdf("processing_wgcna/Module_Traits_Bubble.pdf",width=8,height=1.8)
ggscatter(p, 
  			x = "variable",
  			y = "Module",
   			size="abs",
   			color="Direction",
   			palette=c("blue","red"),
   			alpha = 0.8,
   			xlab = "",
            ylab = "",) +
   			theme_minimal() +
   			rotate_x_text(angle = 45)+
        coord_flip()
dev.off()


# Module Viz
dir.create("processing_wgcna/MODULES/")
exp <- within_data
mod <- read.table("processing_wgcna/ModuleOutput_WITHIN.txt",header=T)

selection <- "WM4"

tmp <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- tmp %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("processing_wgcna/MODULES/WM4_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()

selection <- "WM11"
tmp <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- tmp %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("processing_wgcna/MODULES/WM11_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()


# WM11
selection <- "WM12"
tmp <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- tmp %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("processing_wgcna/MODULES/WM12_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()

# WM19
selection <- "WM19"
tmp <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- tmp %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("processing_wgcna/MODULES/WM19_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()

# WM21
selection <- "WM21"
tmp <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- tmp %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("processing_wgcna/MODULES/WM21_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()

# WM22
selection <- "WM22"
tmp <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- tmp %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)

pdf("processing_wgcna/MODULES/WM22_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()
