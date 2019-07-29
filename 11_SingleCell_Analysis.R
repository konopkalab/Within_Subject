# Single Cell Analysis
rm(list = ls())
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))

source("Utils.R")

# Create directories for outputs
folder_names <- c("processing_scRNAseq","processing_scRNAseq/OUTPUTS","processing_scRNAseq/PLOTS")
sapply(folder_names, dir.create)

## QC ##
#load("SME2_FA25X_FA09X_FA26X_CLUST.RData")
#data <- as.data.frame(as.matrix(sme2@raw.data))
#colnames(data) <- gsub("FA926X_","",colnames(data))
#seuObjSME <- seuratify(data, "ALL_DATA")
#plotbinhist(seuObjSME)
#seuObjSME <- plotqcmito(seuObjSME)
#save(seuObjSME, file = "seuObjSME.RData")

load(here("rawdata","seuObjSME.RData"))
# seuObjDark
##-------------------------------------------------------
## nUMI HISTOGRAM
##-------------------------------------------------------
## Generate histograms for UMI data
all.data <- as.data.frame(as.matrix(seuObjSME@assays$RNA@counts))
dim(all.data)
# 16651  4690
dataSME <- as.data.frame(colSums(all.data))
colnames(dataSME) <- "UMI_Count"

## Plot histograms UMI Distribution (auto-scale)
histUMI <- ggplot(dataSME, aes(x = log10(UMI_Count+10^-10))) + 
geom_histogram(bins = 100, colour="black", fill="grey") + 
labs(title = "BA38", x = "Log10(Number of UMIs)", y = "Number of Cells")
ggsave("processing_scRNAseq/PLOTS/SME_DATA_HIST.pdf", plot = histUMI, width = 6, height = 4, units = "in", dpi = 150)

##-------------------------------------------------------
## DATA FILTERING
##-------------------------------------------------------
nUlo <- -Inf
nUhi <- round((max(seuObjSME@meta.data$nCount_RNA)*85)/100)
nGlo <- -Inf
nGhi <- round((max(seuObjSME@meta.data$nFeature_RNA)*85)/100)
pMlo <- -Inf
pMhi <- 0.1 #(10%)

plotSME_UM <- ggplot(seuObjSME@meta.data, aes(x = nCount_RNA, y = pMito)) + 
				geom_point(size = 1, alpha = 0.5) + 
				labs(title = paste("SME", round(cor(seuObjSME@meta.data$nCount_RNA, seuObjSME@meta.data$pMito), 2), sep = ", "), x = "Number of UMIs", y = "Percent Mito") + 
				annotate("rect", xmin = nUlo, xmax = nUhi, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") + 
				annotate("rect", xmin = -Inf, xmax = Inf, ymin = pMlo, ymax = pMhi, alpha = 0.1, fill = "green")

plotSME_UG <- ggplot(seuObjSME@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) + 
				geom_point(size = 1, alpha = 0.5) + 
				labs(title = paste("SME", round(cor(seuObjSME@meta.data$nCount_RNA, seuObjSME@meta.data$nFeature_RNA), 2), sep = ", "), x = "Number of UMIs", y = "Number of Genes") + 
				annotate("rect", xmin = nUlo, xmax = nUhi, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") + 
				annotate("rect", xmin = -Inf, xmax = Inf, ymin = nGlo, ymax = nGhi, alpha = 0.1, fill = "green")

plotSME_GM <- ggplot(seuObjSME@meta.data, aes(x = nFeature_RNA, y = pMito)) + 
				geom_point(size = 1, alpha = 0.5) + 
				labs(title = paste("SME", round(cor(seuObjSME@meta.data$nFeature_RNA, seuObjSME@meta.data$pMito), 2), sep = ", "), x = "Number of Genes", y = "Percent Mito") + 
				annotate("rect", xmin = nGlo, xmax = nGhi, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") + 
				annotate("rect", xmin = -Inf, xmax = Inf, ymin = pMlo, ymax = pMhi, alpha = 0.1, fill = "red")

plot2by2 <- plot_grid(plotSME_UM, plotSME_UG, plotSME_GM,NULL,labels=c("A", "B", "C", "D"), ncol = 2)
save_plot("processing_scRNAseq/PLOTS/BA38_QC.pdf", plot2by2, ncol = 2,base_height=5,base_width=4)


### Filtering Stuff 

tmp <- SubsetData(object = seuObjSME, subset.name = "nCount_RNA", low.threshold = nUlo, high.threshold = nUhi)
tmp2 <- SubsetData(object = tmp, subset.name = "nFeature_RNA", low.threshold = nGlo, high.threshold = nGhi)
seuObjSME_Filt <- SubsetData(object = tmp2, subset.name = "pMito", low.threshold = pMlo, high.threshold = pMhi)
dim(seuObjSME_Filt@assays$RNA@data)
# 16651  4536
save(seuObjSME_Filt, nUlo, nUhi, nGlo, nGhi, pMlo, pMhi, file = "processing_scRNAseq/OUTPUTS/seuObjSME_Filt.RData")

#####################
## Remove MT genes ##
#####################

all.data.filt <- as.data.frame(as.matrix(seuObjSME_Filt@assays$RNA@counts))
all.data.filt <- all.data.filt[-grep("MT-",rownames(all.data.filt)),]

# Initialize the Seurat object with the raw (non-normalized data).
seuObjSME_final <- CreateSeuratObject(counts = all.data.filt, project = "SME filter")

## Add pMito info from meta data for all cells before filtering
metaAll <- as.data.frame(seuObjSME_Filt@meta.data)


seuObjSME_final <- AddMetaData(object = seuObjSME_final, metadata = as.data.frame(seuObjSME_Filt@meta.data))
#table(seuObjSME_final@meta.data$orig.ident)
#    FA09X FA25X FA26X 
# 1702   881  1953 

save(seuObjSME_final, file = "processing_scRNAseq/OUTPUTS/seuObjSME_Final.RData")

# LogNormalize the data using scaling factor
seuObjSME_Norm <- NormalizeData(object = seuObjSME_final, normalization.method = "LogNormalize", scale.factor = 10000)

# nfeatures = 2000
seuObjSME_Norm <- FindVariableFeatures(object = seuObjSME_Norm,
					mean.function = ExpMean, 
					dispersion.function = LogVMR, 
					x.low.cutoff = 0.0125, 
					x.high.cutoff = 3, 
					y.cutoff = 0.5,
					nfeatures = 2000)

#length(x = seuObjSME_Norm@assays$RNA@var.features)

# Scale the data by regressing out for un-intended variation (nUMI per cell, batch effects, cell cycle, mito.genes, etc.)
# regressing to nUMI only as all mitochondrial genes have been removed during filtering
# also cells with high mitichondrial contents (> 0.3%) have been removed

demoMeta <- read.table("rawdata/demoMeta.txt",header=T)
metaAll <- as.data.frame(seuObjSME_Norm@meta.data)
metaAll$Rows <- rownames(metaAll)

newmeta <- merge(metaAll,demoMeta,by.x="orig.ident",by.y = "Sample",all=T)
rownames(newmeta) <- newmeta$Rows
newmeta$Rows <- NULL

newmeta <- newmeta[match(rownames(metaAll),rownames(newmeta)),]
seuObjSME_Norm@meta.data <- newmeta

seuObjSME_Norm <- ScaleData(object = seuObjSME_Norm, 
					vars.to.regress = c("nCount_RNA", "pMito","Age","EpDur"), 
					model.use = "linear",
					do.par=TRUE,
					num.cores=4)

seuObjSME_Norm_PC <- RunPCA(object = seuObjSME_Norm, 
					features=NULL, 
					weight.by.var = TRUE, 
					ndims.print = 1:5, 
					nfeatures.print = 30, 
					npcs = 50, 
					reduction.name = "pca")

# Examine and visualize PCA results a few different ways
pdf("processing_scRNAseq/PLOTS/BA38_Norm_PC_PCA1.pdf", width = 9, height = 21)
VizDimLoadings(object = seuObjSME_Norm_PC, pcs.use = 1:30)
dev.off()

pdf("processing_scRNAseq/PLOTS/seuObjSME_Norm_PC_2.pdf", width = 9, height = 6)
DimPlot(object = seuObjSME_Norm_PC, dim.1 = 1, dim.2 = 2)
dev.off()

seuObjSME_Norm_PC <- ProjectDim(object = seuObjSME_Norm_PC, reduction = "pca")

pdf("processing_scRNAseq/PLOTS/seuObjSME_Norm_PC_3.pdf", width = 9, height = 21)
DimHeatmap(object = seuObjSME_Norm_PC, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

pdf("processing_scRNAseq/PLOTS/seuObjSME_Norm_PC_4.pdf", width = 6, height = 4)
ElbowPlot(object = seuObjSME_Norm_PC, ndims = 50, reduction = "pca")
dev.off()


# Determine statistically significant principal components, takes long time
seuObjSME_Norm_PC <- JackStraw(object = seuObjSME_Norm_PC, 
						reduction = "pca", 
						dims = 50, 
						num.replicate = 100,
						maxit = 1000)

seuObjSME_Norm_PC <- ScoreJackStraw(object = seuObjSME_Norm_PC, 
						dims = 1:50)

pdf("processing_scRNAseq/PLOTS/seuObjSME_Norm_PC_5.pdf", width = 12, height = 8)
JackStrawPlot(object = seuObjSME_Norm_PC, dims = 1:50, reduction = "pca")
dev.off()


# Cluster the cells - graph based (resolution: 0.6 - 1.2)

numPCs <- 25

seuObjSME_Norm_PC_Clust <- FindNeighbors(object = seuObjSME_Norm_PC, 
							reduction = "pca", 
							dims = 1:numPCs)

seuObjSME_Norm_PC_Clust <- FindClusters(object = seuObjSME_Norm_PC_Clust, 
							resolution = 0.8, 
							algorithm = 1,
							n.iter = 100,
							save.SNN = TRUE)

##-------------------------------------------------------
## CLUSTERING using TSNE
##-------------------------------------------------------

# Run Non-linear dimensional reduction (tSNE)
seuObjSME_Norm_PC_Clust <- RunTSNE(object = seuObjSME_Norm_PC_Clust, 
							dims = 1:numPCs, 
							reduction = "pca")

pdf("processing_scRNAseq/PLOTS/seuObjSME_Norm_PC_Clust_TSNE.pdf", width = 7, height = 6)
DimPlot(object = seuObjSME_Norm_PC_Clust, reduction = "tsne", label = TRUE, pt.size = 0.5)
dev.off()

# Add UMAP
# sudo pip --proxy=https://S157784@proxy.swmed.edu:3128 install umap-learn
seuObjSME_Norm_PC_Clust <- RunUMAP(object = seuObjSME_Norm_PC_Clust, 
							reduction = "pca", 
							dims = 1:numPCs)

pdf("processing_scRNAseq/PLOTS/seuObjSME_Norm_PC_Clust_UMAP.pdf", width = 7, height = 6)
DimPlot(object = seuObjSME_Norm_PC_Clust, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

# Save Object
save(seuObjSME_Norm_PC_Clust, file = "processing_scRNAseq/OUTPUTS/seuObjSME_Norm_PC_Clust_tSNE_UMAP.RData")


##-------------------------------------------------------
## DEG
##-------------------------------------------------------
seuObjSME.markers <- FindAllMarkers(object = seuObjSME_Norm_PC_Clust,
						only.pos = TRUE, min.pct = 0.25, thresh.use = 0.2)

write.table(seuObjSME.markers, "processing_scRNAseq/OUTPUTS/seuObjSME.markers_OriginalClusters.txt", row.names = T, col.names = T, quote = F, sep = "\t")

# Relabel the cluster based on allen data

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new.cluster.ids <- c("Astro",
					"Exc_1",
					"Exc_2",
					"Exc_3",
					"Inh_1",
					"Inh_2",
					"Exc_4",
					"Inh_3",
					"Oligo",
					"Exc_5",
					"Inh_4",
					"Exc_6",
					"Exc_7",
					"OPC",
					"Exc_8")

seuObjSME_Norm_PC_Clust@active.ident <- plyr::mapvalues(x = seuObjSME_Norm_PC_Clust@active.ident, from = current.cluster.ids, to = new.cluster.ids)

pdf("processing_scRNAseq/PLOTS/seuObjSME_Norm_PC_Clust_UMAP_mapped.pdf", width = 7, height = 6)
DimPlot(object = seuObjSME_Norm_PC_Clust, reduction = "umap", label = TRUE, pt.size = 0.5)+
theme(legend.position="none")
dev.off()

save(seuObjSME_Norm_PC_Clust, file = "processing_scRNAseq/OUTPUTS/seuObjSME_Norm_PC_Clust_tSNE_UMAP_ReMap.RData")

seuObjSME.markers_map <- FindAllMarkers(object = seuObjSME_Norm_PC_Clust,
						only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

seuObjSME.markers_map2 <- seuObjSME.markers_map %>%
							filter(p_val_adj < 0.05)

write.table(seuObjSME.markers_map, "processing_scRNAseq/OUTPUTS/seuObjSME.markers_NewMap.txt", row.names = T, col.names = T, quote = F, sep = "\t")


# Home made violin
seuObjSME_Norm_PC_Clust@active.ident <- factor(seuObjSME_Norm_PC_Clust@active.ident, levels = c("Astro",
																								"Oligo",
																								"OPC",
																								"Exc_1",
																								"Exc_2",
																								"Exc_3",
																								"Exc_4",
																								"Exc_5",
																								"Exc_6",
																								"Exc_7",
																								"Exc_8",
																								"Inh_1",
																								"Inh_2",
																								"Inh_3",
																								"Inh_4"
																								))

# Plot Inhibitory Markers
genes <- c(
"CCK",
"SST",
"GAD1",
"GAD2",
"FGF13",
"SV2C",
"ANK1",
"RELN"
)

plots <- VlnPlot(object = seuObjSME_Norm_PC_Clust, features = genes, slot = "data",pt.size = 0, same.y.lims = TRUE,log=TRUE,ncol=1,combine=FALSE)

n <- length(genes)

for(i in 1:n) {
  plots[[i]] <- plots[[i]] + 
  theme_bw() + 
  theme(legend.position = 'none')+ 
  labs(title="", x ="", y = paste(genes[[i]])) + 
  theme(axis.title.y = element_text(angle = 90,hjust=1)) + #axis.text.x=element_blank(),axis.ticks.x=element_blank(),
  #theme(plot.margin = unit(c(0,1,0.5,1), "cm")) + # c(bottom, left, top, right)
  rotate_x_text(angle = 45)
}

#cowplot::plot_grid(plotlist = plots, ncol = 2,align="v")
pdf("processing_scRNAseq/PLOTS/Markers_Cells_Inhibitory.pdf",width=6,height=5)
cowplot::plot_grid(plotlist = plots, ncol = 2,align="v")
dev.off()


# Plot Exchitatory Markers
genes <- c(
"RORB",
"FOXP2",
"SLC17A7",
"COL22A1",
"CBLN2",
"LRRK1",
"HTR2C",
"RARB"
)

plots <- VlnPlot(object = seuObjSME_Norm_PC_Clust, features = genes, slot = "data",pt.size = 0, same.y.lims = TRUE,log=TRUE,ncol=1,combine=FALSE)

n <- length(genes)

for(i in 1:n) {
  plots[[i]] <- plots[[i]] + 
  theme_bw() + 
  theme(legend.position = 'none')+ 
  labs(title="", x ="", y = paste(genes[[i]])) + 
  theme(axis.title.y = element_text(angle = 90,hjust=1)) + #axis.text.x=element_blank(),axis.ticks.x=element_blank(),
  #theme(plot.margin = unit(c(0,1,0.5,1), "cm")) + # c(bottom, left, top, right)
  rotate_x_text(angle = 45)
}

#cowplot::plot_grid(plotlist = plots, ncol = 2,align="v")
pdf("processing_scRNAseq/PLOTS/Markers_Cells_Excitatory.pdf",width=6,height=5)
cowplot::plot_grid(plotlist = plots, ncol = 2,align="v")
dev.off()


# Plot glia Markers
genes <- c(
"FGFR3",
"CABLES1",
"SLC14A1",
"MOBP",
"ST18",
"RNF220",
"VCAN",
"MEGF11"
)

plots <- VlnPlot(object = seuObjSME_Norm_PC_Clust, features = genes, slot = "data",pt.size = 0, same.y.lims = TRUE,log=TRUE,ncol=1,combine=FALSE)

n <- length(genes)

for(i in 1:n) {
  plots[[i]] <- plots[[i]] + 
  theme_bw() + 
  theme(legend.position = 'none')+ 
  labs(title="", x ="", y = paste(genes[[i]])) + 
  theme(axis.title.y = element_text(angle = 90,hjust=1)) + #axis.text.x=element_blank(),axis.ticks.x=element_blank(),
  #theme(plot.margin = unit(c(0,1,0.5,1), "cm")) + # c(bottom, left, top, right)
  rotate_x_text(angle = 45)
}

#cowplot::plot_grid(plotlist = plots, ncol = 2,align="v")
pdf("processing_scRNAseq/PLOTS/Markers_Cells_Glias.pdf",width=6,height=5)
cowplot::plot_grid(plotlist = plots, ncol = 2,align="v")
dev.off()


# Cool Genes Markers
genes <- c(
"IL1RAPL2",
"CCDC85A",
"SMAD3",
"NEDD4L"
)

plots <- VlnPlot(object = seuObjSME_Norm_PC_Clust, features = genes, slot = "data",pt.size = 0, same.y.lims = TRUE,log=TRUE,ncol=1,combine=FALSE)

n <- length(genes)

for(i in 1:n) {
  plots[[i]] <- plots[[i]] + 
  theme_bw() + 
  theme(legend.position = 'none')+ 
  labs(title="", x ="", y = paste(genes[[i]])) + 
  theme(axis.title.y = element_text(angle = 90,hjust=1)) + #axis.text.x=element_blank(),axis.ticks.x=element_blank(),
  #theme(plot.margin = unit(c(0,1,0.5,1), "cm")) + # c(bottom, left, top, right)
  rotate_x_text(angle = 45)
}

#cowplot::plot_grid(plotlist = plots, ncol = 2,align="v")
pdf("processing_scRNAseq/PLOTS/Markers_Cells_CoolGenes.pdf",width=6,height=3)
cowplot::plot_grid(plotlist = plots, ncol = 2,align="v")
dev.off()



pdf("processing_scRNAseq/PLOTS/Markers_Cells_CoolGenes_DotPlot.pdf",width=5,height=6)
DotPlot(seuObjSME_Norm_PC_Clust, features = genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis()
dev.off()


