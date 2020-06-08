# Single Cell Analysis
rm(list = ls())
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
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
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(future))
source("UTILS/Utils.R")

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 5000 * 1024^2)

# Create directories for outputs
folder_names <- c("integration_scRNAseq","integration_scRNAseq/OUTPUTS","integration_scRNAseq/PLOTS")
sapply(folder_names, dir.create)

# Load first data
dt1 <- read.table(url("https://www.dropbox.com/s/d2prk166dl4nx92/Lega_31_Counts.txt?dl=1"),sep="\t",header=T)
dt2 <- read.table(url("https://www.dropbox.com/s/43ve3ewc37a00q0/Lega_33_Counts.txt?dl=1"),sep="\t",header=T)
dt3 <- read.table(url("https://www.dropbox.com/s/0267zft4s2z0lal/Lega_41_Counts.txt?dl=1"),sep="\t",header=T)

mat <- cbind(dt1,dt2,dt3)

meta <- read.table(url("https://www.dropbox.com/s/o4wyfm45kt2ozk6/metadata_integration.txt?dl=1"),sep="\t",header=T)


demo_2 <- data.frame(ID = colnames(mat), 
				   Subject = c(rep("Lega31",ncol(dt1)), rep("Lega33",ncol(dt2)),rep("Lega41",ncol(dt3))))

demo_2 <- merge(demo_2,meta,by.x = "Subject",by.y="Sample",all=F)

rownames(demo_2) <- paste(demo_2$Subject,demo_2$ID, sep="_")

colnames(mat) <- rownames(demo_2)

mat_batch2 <- mat
demo_batch2 <- demo_2

# Load first batch
load(url("https://www.dropbox.com/s/y01v4jy3l9w1wuk/seuObjSME.RData?dl=1"))

mat_batch1 <- as.data.frame(as.matrix(seuObjSME@assays$RNA@counts))

colnames(mat_batch1) <- gsub("FA09X","Lega42",colnames(mat_batch1))
colnames(mat_batch1) <- gsub("FA25X","Lega15",colnames(mat_batch1))
colnames(mat_batch1) <- gsub("FA26X","Lega45",colnames(mat_batch1))

demo_1 <- data.frame(ID = do.call(rbind,strsplit(colnames(mat_batch1),"_"))[,2], 
				   Subject = do.call(rbind,strsplit(colnames(mat_batch1),"_"))[,1])

demo_1 <- merge(demo_1,meta,by.x = "Subject",by.y="Sample",all=F)

rownames(demo_1) <- paste(demo_1$Subject,demo_1$ID, sep="_")

demo_batch1 <- demo_1

# Matching the rows
genes <- intersect(rownames(mat_batch1),rownames(mat_batch2))

mat_batch1 <- mat_batch1[rownames(mat_batch1) %in% genes,]
mat_batch2 <- mat_batch2[rownames(mat_batch2) %in% genes,]

mat_batch1 <- mat_batch1[match(genes,rownames(mat_batch1)),]
mat_batch2 <- mat_batch2[match(genes,rownames(mat_batch2)),]

full_mat <- cbind(mat_batch1,mat_batch2)

full_demo <- rbind(demo_batch1,demo_batch2)

## Start the analysis
dataSME <- as.data.frame(colSums(full_mat))
colnames(dataSME) <- "UMI_Count"
histUMI <- ggplot(dataSME, aes(x = log10(UMI_Count+10^-10))) + 
geom_histogram(bins = 100, colour="black", fill="grey") + 
labs(title = "BA38", x = "Log10(Number of UMIs)", y = "Number of Cells")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_HIST.pdf", plot = histUMI, width = 6, height = 4, units = "in", dpi = 150)

# N of cells
ncells <- table(full_demo$Subject) %>% 
	as.data.frame() %>% 
	mutate(Batch = c("Batch1","Batch1","Batch1","Batch2","Batch2","Batch2")) %>% 
	ggbarplot(
  	 	"Var1", 
  	 	"Freq", 
  	 	color = "Batch",
   		palette = c("#E7B800", "#FC4E07")) + 
  		theme_classic() +
		rotate_x_text(angle = 45) + 
		xlab("") + 
		ylab("# of Cells")+
		theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_NofCells.pdf", plot = ncells, width = 4, height = 4, units = "in", dpi = 150)

# Seurat Object
sme2 <- CreateSeuratObject(counts = full_mat, project = "Sme_New_Data")

sme2[["pMito"]] <- PercentageFeatureSet(sme2, pattern = "^MT-")

sme2@meta.data

sme2@meta.data <- sme2@meta.data %>%
			        dplyr::rename(nUMI = nCount_RNA,
			                      nGene = nFeature_RNA)

sme2@meta.data$Age <- full_demo$Age
sme2@meta.data$EpDur <- full_demo$EpDur
sme2@meta.data$RIN <- full_demo$RIN
sme2@meta.data$Batch <- full_demo$Batch

pdf("integration_scRNAseq/PLOTS/SME_Integration_QC_Plot.pdf", width=5,height=4)
VlnPlot(sme2, features = c("nGene", "nUMI", "pMito"), ncol = 3,pt.size=0)
dev.off()

ngenes <- sme2@meta.data %>% 
		ggboxplot( 
			x = "orig.ident", 
			y = "nGene",
 			color = "Batch",
 			palette =c("#E7B800", "#FC4E07")) + 
  		theme_classic() +
		rotate_x_text(angle = 45) + 
		xlab("") + 
		ylab("# of Genes")+
		theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_nGENES.pdf", plot = ngenes, width = 4, height = 4, units = "in", dpi = 150)


numis <- sme2@meta.data %>% 
		ggboxplot( 
			x = "orig.ident", 
			y = "nUMI",
 			color = "Batch",
 			palette =c("#E7B800", "#FC4E07")) + 
  		theme_classic() +
		rotate_x_text(angle = 45) + 
		xlab("") + 
		ylab("# of UMI")+
		theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_nUMI.pdf", plot = numis, width = 4, height = 4, units = "in", dpi = 150)

numislog <- sme2@meta.data %>% 
		ggboxplot( 
			x = "orig.ident", 
			y = "nUMI",
 			color = "Batch",
 			palette =c("#E7B800", "#FC4E07")) + 
  		theme_classic() +
  		scale_y_log10() +
  		annotation_logticks(sides = "l") +
		rotate_x_text(angle = 45) + 
		xlab("") + 
		ylab("log10(UMI)")+
		theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_log10_nUMI.pdf", plot = numislog, width = 4, height = 4, units = "in", dpi = 150)

nmito <- sme2@meta.data %>% 
		ggboxplot( 
			x = "orig.ident", 
			y = "pMito",
 			color = "Batch",
 			palette =c("#E7B800", "#FC4E07")) + 
  		theme_classic() +
		rotate_x_text(angle = 45) + 
		xlab("") + 
		ylab("% Mito")+
		theme(legend.position="none") + 
		ylim(c(0,10))
ggsave("integration_scRNAseq/PLOTS/SME_Integration_pMito.pdf", plot = nmito, width = 4, height = 4, units = "in", dpi = 150)

pdf("integration_scRNAseq/PLOTS/SME_Integrated_UMIvsGene.pdf", width=5,height=4)
ggplot(sme2@meta.data, aes(x=nUMI, y=nGene, color=orig.ident)) +
ggrastr::geom_point_rast(size=0.5) + 
theme_classic()
dev.off()

pdf("integration_scRNAseq/PLOTS/SME_Integrated_pMitoVsUMI.pdf", width=5,height=4)
ggplot(sme2@meta.data, aes(x=nUMI, y=pMito, color=orig.ident)) +
ggrastr::geom_point_rast(size=0.5) + 
theme_classic()
dev.off()


# Visualize the correlation between genes detected and number of UMIs. D
# Determine whether strong presence of cells with low numbers of genes/UMIs
scatter <- sme2@meta.data %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=pMito)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~orig.ident)
ggsave("integration_scRNAseq/PLOTS/SME_Integration_GeneXUMI.pdf", plot = scatter, width = 6, height = 4, units = "in", dpi = 150)


##-------------------------------------------------------
## DATA FILTERING
##-------------------------------------------------------
SME_Integration_Filt <- subset(x = sme2, subset = nUMI < 10000 & pMito < 5)

#####################
## Remove MT genes ##
#####################

all.data.filt <- as.data.frame(as.matrix(SME_Integration_Filt@assays$RNA@counts))
all.data.filt <- all.data.filt[-grep("MT-",rownames(all.data.filt)),]

# Initialize the Seurat object with the raw (non-normalized data).
SME_Integration_Final <- CreateSeuratObject(counts = all.data.filt, project = "SME filter")

## Add pMito info from meta data for all cells before filtering
metaAll <- as.data.frame(SME_Integration_Filt@meta.data)
SME_Integration_Final <- AddMetaData(object = SME_Integration_Final, metadata = as.data.frame(SME_Integration_Filt@meta.data))
SME_Integration_Final@meta.data$nCount_RNA <- NULL
SME_Integration_Final@meta.data$nFeature_RNA <- NULL

# N cell surviving
ncells <- table(SME_Integration_Final@meta.data$orig.ident) %>% 
	as.data.frame() %>% 
	mutate(Batch = c("Batch1","Batch1","Batch1","Batch2","Batch2","Batch2")) %>% 
	ggbarplot(
  	 	"Var1", 
  	 	"Freq", 
  	 	color = "Batch",
   		palette = c("#E7B800", "#FC4E07")) + 
  		theme_classic() +
		rotate_x_text(angle = 45) + 
		xlab("") + 
		ylab("# of Cells")+
		theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_NofCells_AfterFiltering.pdf", plot = ncells, width = 4, height = 4, units = "in", dpi = 150)

# Save Filtered Data
save(SME_Integration_Final, file = "integration_scRNAseq/OUTPUTS/SME_Integration_Final.RData")

#####################
# Downstream analysis
#####################

# Data integration by SCtransform!
load("integration_scRNAseq/OUTPUTS/SME_Integration_Final.RData")

# Integration
SME_Integration_Split <- SplitObject(SME_Integration_Final, split.by = "Batch")

SME_Integration_Split <- SME_Integration_Split[c("Batch1", "Batch2")]

for (i in 1:length(SME_Integration_Split)) {
    SME_Integration_Split[[i]] <- SCTransform(SME_Integration_Split[[i]], 
				    						vars.to.regress = c("nUMI","pMito","Age","EpDur","RIN"), 
											verbose = FALSE)
    }

integ_features <- SelectIntegrationFeatures(object.list = SME_Integration_Split, 
											nfeatures = 3000) 

SME_Integration_Split <- PrepSCTIntegration(object.list = SME_Integration_Split, 
											anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = SME_Integration_Split, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

SME_Data_Integrated <- IntegrateData(
								anchorset = integ_anchors,
								new.assay.name = "integrated",
								normalization.method = "SCT",
								dims = 1:30,
								k.weight = 100,
								sd.weight = 1,
								do.cpp = TRUE,
								eps = 0.5,
								verbose = TRUE
								)

SME_Data_Integrated <- RunPCA(object = SME_Data_Integrated, 
								features=NULL, 
								weight.by.var = TRUE, 
								ndims.print = 1:5, 
								nfeatures.print = 30, 
								npcs = 30, 
								reduction.name = "pca")

SME_Data_Integrated <- FindNeighbors(object = SME_Data_Integrated, 
										reduction = "pca", 
										dims = 1:30, 
										nn.eps = 0.5)

SME_Data_Integrated <- FindClusters(object = SME_Data_Integrated, 
										resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4), 
										algorithm = 1,
										n.iter = 1000)

pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_Clustree.pdf", width = 12, height = 6)
clustree(SME_Data_Integrated@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Select resolution and run UMAP
Idents(object = SME_Data_Integrated) <- "integrated_snn_res.0.8"

SME_Data_Integrated <- RunUMAP(object = SME_Data_Integrated, 
										reduction = "pca", 
										dims = 1:30)

pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP.pdf", width = 7, height = 6)
DimPlot(object = SME_Data_Integrated, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

# Select the RNA counts slot to be the default assay
DefaultAssay(SME_Data_Integrated) <- "RNA"
SME_Data_Integrated <- NormalizeData(object = SME_Data_Integrated, 
						normalization.method = "LogNormalize", 
						scale.factor = 10000)

save(SME_Data_Integrated, file = "integration_scRNAseq/OUTPUTS/SME_Data_Integrated.RData")

# Differential Expression
Integrated.Markers <- FindAllMarkers(object = SME_Data_Integrated,
						min.cells.gene = 5,
						only.pos = TRUE, 
						logfc.threshold = 0.3)

Integrated.Markers$cluster <- paste("Cluster",Integrated.Markers$cluster,sep="_")

Integrated.Markers.Sign <- Integrated.Markers %>%
								filter(p_val_adj < 0.05,pct.1 > 0.5,pct.2 < 0.6)
openxlsx::write.xlsx(Integrated.Markers.Sign, file = "integration_scRNAseq/OUTPUTS/SME_Integrated_DGE_Stats.xlsx", colNames = TRUE, borders = "columns")

Integrated.Top.Markers <- tbl_df(Integrated.Markers.Sign) %>% 
                  group_by(cluster) %>% 
                  top_n(n = 5, wt = avg_logFC) %>%
                  as.data.frame()

save(Integrated.Markers,Integrated.Markers.Sign, Integrated.Markers.Sign, Integrated.Top.Markers,
	file = "integration_scRNAseq/OUTPUTS/SME_Integrated_scRNAseq_DGEdata.RData")

df <- Integrated.Markers.Sign %>% select(gene,cluster) %>% rename(Gene = gene)
write.table(df,"integration_scRNAseq/OUTPUTS/SME_Integrated_DGE.txt",sep="\t",quote=F,row.names=F)



# Run enrichment for markers of MTG (Hodge et al. 2019)
# Output will be a new cluter label that can be use to label the current data 
system("R CMD BATCH --vanilla Fisher_enrich_Markers.R")

# Relabel
labels <- read.table("Labels_Clusters.txt",header=T,sep="\t")

new <- labels %>% 
		separate(variable, c("Cell", "Layer","Marker1","Marker2"),"_",remove = FALSE) %>%
		separate(Rows, c("Class", "Cluster"),"_",remove = FALSE) %>%
		unite(Ident, c("Cell","Cluster"),sep = "_",remove = FALSE) %>%
		unite(Ident2, c("Cell","Layer","Marker1"),sep = "_",remove = FALSE)

new <- new[order(as.numeric(as.character(new$Cluster))), ]

current.cluster.ids <- new$Cluster
new.cluster.ids <- as.character(new$variable)

SME_Data_Integrated@active.ident <- plyr::mapvalues(x = SME_Data_Integrated@active.ident, 
													from = current.cluster.ids, 
													to = new.cluster.ids)

pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP_Labels.pdf", width = 7, height = 6)
DimPlot(object = SME_Data_Integrated, reduction = "umap", label = TRUE, pt.size = 0.5)+
theme(legend.position="none")
dev.off()

SME_Data_Integrated@meta.data$Cell <- SME_Data_Integrated@active.ident
SME_Data_Integrated@meta.data <- SME_Data_Integrated@meta.data %>%
                                       rownames_to_column("TMP") %>%
                                       mutate(Class = case_when(grepl("Exc", Cell) ~ "Glutamatergic", 
                                                            grepl("Inh", Cell) ~ "Gabaergic",
                                                            grepl("Astro|Olig|OPC", Cell) ~ "NonNeuronal")) %>%
                                       column_to_rownames("TMP")

SME_Data_Integrated@meta.data <- SME_Data_Integrated@meta.data %>%
                                       rownames_to_column("TMP") %>%
                                       mutate(Definition = case_when(grepl("Exc", Cell) ~ "Excitatory", 
                                                            grepl("Inh", Cell) ~ "Inhibitory",
                                                            grepl("Astro", Cell) ~ "Astrocytes", 
                                                            grepl("Olig", Cell) ~ "Oligodendrocytes",
                                                            grepl("OPC", Cell) ~ "OPC")) %>%
                                       column_to_rownames("TMP")

# Reorder levels

SME_Data_Integrated@meta.data$Cell <- factor(SME_Data_Integrated@meta.data$Cell, levels = c(
"Exc_L2-3_LINC00507_FREM3",
"Exc_L3-5_RORB_TWIST2",
"Exc_L3-5_RORB_ESR1",
"Exc_L4-6_RORB_SEMA3E",
"Exc_L4-6_FEZF2_IL26",
"Exc_L5-6_THEMIS_CRABP1",
"Exc_L5-6_THEMIS_C1QL3",
"Exc_L5-6_FEZF2_ABO",
"Inh_L1-2_LAMP5_DBP",
"Inh_L1-2_SST_BAGE2",
"Inh_L1-3_SST_CALB1",
"Inh_L2-4_VIP_SPAG17",
"Inh_L2-4_VIP_CBLN1",
"Inh_L2-5_PVALB_SCUBE3",
"Inh_L2-6_LAMP5_CA1",
"Inh_L4-6_PVALB_SULF1",
"Inh_L5-6_SST_KLHDC8A",
"Astro_L1-6_FGFR3_SLC14A1",
"Oligo_L1-6_OPALIN_NA",
"OPC_L1-6_PDGFRA_NA"))

save(SME_Data_Integrated, file = "integration_scRNAseq/OUTPUTS/SME_Data_Integrated_Labels.RData")

# other stats
prop_cell <- SME_Data_Integrated@meta.data %>% 
  	ggplot(aes(x=Class, fill=Cell)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("Prop Cells")+
theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_prop_cell.pdf", plot = prop_cell, width = 6, height = 4, units = "in", dpi = 150)


prop_cell_persubj <- SME_Data_Integrated@meta.data %>% 
  	ggplot(aes(x=orig.ident, fill=Class)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("Prop Cells per Subject")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_prop_cell_subj.pdf", plot = prop_cell_persubj, width = 6, height = 4, units = "in", dpi = 150)


pie <- table(SME_Data_Integrated@meta.data$Class) %>% 
	   as.data.frame()  %>% 
	   arrange(desc(Freq)) %>% 
	   mutate(percent = scales::percent(Freq/sum(Freq))) %>%
	   ggplot(aes(x = "", y = Freq, fill = fct_inorder(Var1))) +
       geom_bar(width = 1, stat = "identity") +
       coord_polar("y", start = 0) +
       geom_label_repel(aes(label = percent), size=5, show.legend = F, nudge_x = 1) +
       guides(fill = guide_legend(title = "Cell Class"))
ggsave("integration_scRNAseq/PLOTS/SME_Integration_prop_cell_pie.pdf", plot = pie, width = 6, height = 4, units = "in", dpi = 150)

# UMAP Alternative
umap <- as.data.frame(Embeddings(SME_Data_Integrated, reduction = "umap"))
meta <- as.data.frame(SME_Data_Integrated@meta.data)

df <- cbind(umap,meta)%>% 
	group_by(Cell) %>% 
	mutate(N = n()) %>%
  	ungroup() %>% 
  	arrange(Definition) %>%
  	as.data.frame()  


label <- data.frame(Cell=levels(df$Cell),label=levels(df$Cell))

label_2 <- df %>% 
  group_by(Cell) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  left_join(label) %>%
  as.data.frame() %>%
  mutate(Definition = case_when(grepl("Exc", Cell) ~ "Excitatory", 
                       grepl("Inh", Cell) ~ "Inhibitory",
                       grepl("Astro", Cell) ~ "Astrocytes", 
                       grepl("Olig", Cell) ~ "Oligodendrocytes",
                       grepl("OPC", Cell) ~ "OPC"))


colors <- cbind(umap,meta)%>% 
			group_by(Cell) %>%
			summarize(N = n()) %>%
			arrange(desc(N)) %>%
			mutate(Definition = case_when(grepl("Exc", Cell) ~ "Excitatory", 
                       grepl("Inh", Cell) ~ "Inhibitory",
                       grepl("Astro", Cell) ~ "Astrocytes", 
                       grepl("Olig", Cell) ~ "Oligodendrocytes",
                       grepl("OPC", Cell) ~ "OPC")) %>%
            mutate(Class = case_when(grepl("Exc", Cell) ~ "Glutamatergic", 
                        grepl("Inh", Cell) ~ "Gabaergic",
                        grepl("Astro|Olig|OPC", Cell) ~ "NonNeuronal"))	%>%		
  	as.data.frame()  

l <- split(colors,colors$Class)

#l[[1]]$Color <- colorRampPalette(rev(brewer.pal(4,"pal_peach")))(9)
#l[[2]]$Color <- colorRampPalette(rev(brewer.pal(4,"pal_karpfenblau")))(8)
#l[[3]]$Color <- colorRampPalette(rev(brewer.pal(3,"pal_pinky")))(3)

l[[1]]$Color <- rev(usecol(pal_bordeaux, 9))
l[[2]]$Color <- rev(usecol(pal_karpfenblau, 8))
l[[3]]$Color <- rev(usecol(pal_seegruen, 3))


colors <- do.call(rbind,l)
colors <- colors[match(label_2$label, colors$Cell),]


pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP_Alternative.pdf", width = 5, height = 5)
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] 
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Cell),size=0.5) +
theme(legend.position="none") +
ggrepel::geom_text_repel(data = label_2, aes(label = label),
							color = "black",
							#fontface = 'bold',
							segment.colour = "grey60",
						    box.padding = unit(0.25, "lines"),
						    point.padding = unit(0.5, "lines"),
						    nudge_x = .15,
						    nudge_y = 1,
						    size = 2.5) + 
		#scale_color_viridis(discrete=TRUE,option="inferno")
		#scale_colour_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(31))
		scale_colour_manual(values = colors$Color)
dev.off()


# UMAP Alternative by Cell Class
pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP_ByCellClass_Alternative.pdf", width = 6, height = 3)
ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Class)) +
ggrastr::geom_point_rast(size=0.5) +
facet_wrap(.~Class) +
theme_classic() +
theme(legend.position="none")+ 
scale_colour_manual(values = c("blue","red","green"))
dev.off()

# UMAP Alternative by Cell Class
pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP_BySubject_Alternative.pdf", width = 6, height = 4)
ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Batch)) +
ggrastr::geom_point_rast(size=0.5) +
facet_wrap(.~orig.ident,ncol=3) +
theme_classic() +
theme(legend.position="none")+
scale_colour_manual(values = c("#E7B800", "#FC4E07"))
dev.off()

pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP_StackedBar.pdf", width = 8, height = 4)
ggplot(SME_Data_Integrated@meta.data, aes(x=Cell, fill=orig.ident)) + 
geom_bar(position = "fill")+
rotate_x_text(angle = 45)
dev.off()

pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_NCell_Pheatmap.pdf", width = 8, height = 4)
n_cells <- FetchData(SME_Data_Integrated, vars = c("ident", "orig.ident")) %>%
        		count(ident, orig.ident) %>%
        		spread(ident, n) %>%
        		as.data.frame() %>%
        		column_to_rownames("orig.ident") %>%
				pheatmap::pheatmap(border_color = "black",cluster_cols=FALSE)
dev.off()


# Recalculate DGE
Integrated.Markers.Labels <- FindAllMarkers(object = SME_Data_Integrated,
						min.cells.gene = 5,
						only.pos = TRUE, 
						min.pct = 0.25, 
						thresh.use = 0.25)



Integrated.Markers.Labels$cluster <- factor(Integrated.Markers.Labels$cluster, levels = c(
"Exc_L2-3_LINC00507_FREM3",
"Exc_L3-5_RORB_TWIST2",
"Exc_L3-5_RORB_ESR1",
"Exc_L4-6_RORB_SEMA3E",
"Exc_L4-6_FEZF2_IL26",
"Exc_L5-6_THEMIS_CRABP1",
"Exc_L5-6_THEMIS_C1QL3",
"Exc_L5-6_FEZF2_ABO",
"Inh_L1-2_LAMP5_DBP",
"Inh_L1-2_SST_BAGE2",
"Inh_L1-3_SST_CALB1",
"Inh_L2-4_VIP_SPAG17",
"Inh_L2-4_VIP_CBLN1",
"Inh_L2-5_PVALB_SCUBE3",
"Inh_L2-6_LAMP5_CA1",
"Inh_L4-6_PVALB_SULF1",
"Inh_L5-6_SST_KLHDC8A",
"Astro_L1-6_FGFR3_SLC14A1",
"Oligo_L1-6_OPALIN_NA",
"OPC_L1-6_PDGFRA_NA"))

Integrated.Markers.Labels.Sign <- Integrated.Markers.Labels %>%
										filter(p_val_adj < 0.05,pct.1 > 0.5,pct.2 < 0.6)


Integrated.Labels.Top.Markers <- tbl_df(Integrated.Markers.Labels.Sign) %>% 
                  group_by(cluster) %>% 
                  top_n(n = 5, wt = avg_logFC) %>%
                  as.data.frame()

save(Integrated.Markers.Labels,Integrated.Markers.Labels.Sign, Integrated.Labels.Top.Markers,
	file = "integration_scRNAseq/OUTPUTS/SME_Integrated_scRNAseq_Labels_DGEdata.RData")

df <- Integrated.Markers.Labels.Sign %>% select(gene,cluster) %>% rename(Gene = gene)
write.table(df,"integration_scRNAseq/OUTPUTS/SME_Integrated_Labels_DGE.txt",sep="\t",quote=F,row.names=F)

GeneSets <- split(df,df$cluster)
save(GeneSets, file = "integration_scRNAseq/OUTPUTS/SME_Integrated_scRNA_GeneSets.RData")


# Remake teh violin
SME_Data_Integrated@active.ident <- factor(SME_Data_Integrated@active.ident, levels = c(
"Exc_L2-3_LINC00507_FREM3",
"Exc_L3-5_RORB_TWIST2",
"Exc_L3-5_RORB_ESR1",
"Exc_L4-6_RORB_SEMA3E",
"Exc_L4-6_FEZF2_IL26",
"Exc_L5-6_THEMIS_CRABP1",
"Exc_L5-6_THEMIS_C1QL3",
"Exc_L5-6_FEZF2_ABO",
"Inh_L1-2_LAMP5_DBP",
"Inh_L1-2_SST_BAGE2",
"Inh_L1-3_SST_CALB1",
"Inh_L2-4_VIP_SPAG17",
"Inh_L2-4_VIP_CBLN1",
"Inh_L2-5_PVALB_SCUBE3",
"Inh_L2-6_LAMP5_CA1",
"Inh_L4-6_PVALB_SULF1",
"Inh_L5-6_SST_KLHDC8A",
"Astro_L1-6_FGFR3_SLC14A1",
"Oligo_L1-6_OPALIN_NA",
"OPC_L1-6_PDGFRA_NA"))


genes <- c("SLC17A7",
           "CBLN2", 
           "ATP10A",
           "CUX2",
           "RORB",
           "ADAM12",
           "HS3ST2",
           "FOXP2",
           "HTR2C")

pdf("integration_scRNAseq/PLOTS/Markers_Excitatory_Violin.pdf",width=10,height=6)
StackedVlnPlot(SME_Data_Integrated,features = genes) 
dev.off()

pdf("integration_scRNAseq/PLOTS/Markers_Excitatory_Cells_Violin.pdf",width=4,height=5)
ident <- rownames(SME_Data_Integrated@meta.data[SME_Data_Integrated@meta.data$Class == "Glutamatergic",])
exc <- SubsetData(SME_Data_Integrated, cells = ident)
StackedVlnPlot(exc,features = genes) 
dev.off()


genes <- c(
           "GAD1",
           "ADRA1B",
           "RELN",           
           "SST",
           "CCK",
           "LAMP5",
           "PAWR",
           "ROR2",
           "ADAM12",
           "FGD5"
           )

pdf("integration_scRNAseq/PLOTS/Markers_Inhibitory_Violin.pdf",width=10,height=6)
StackedVlnPlot(SME_Data_Integrated,features = genes) 
dev.off()

pdf("integration_scRNAseq/PLOTS/Markers_Cells_Inhibitory_Violin.pdf",width=4,height=6)
ident <- rownames(SME_Data_Integrated@meta.data[SME_Data_Integrated@meta.data$Class == "Gabaergic",])
inh <- SubsetData(SME_Data_Integrated, cells = ident)
StackedVlnPlot(inh,features = genes) 
dev.off()

# Non Neuronal
genes <- c(
           "SLC1A2",
           "FGFR3",
           "ADGRV1",
           "MOBP",
           "MOG",
           "ABCA8",
           "ST18",
           "VCAN",
           "SMOC1")

pdf("integration_scRNAseq/PLOTS/Markers_NonNeuronal_Violin.pdf",width=10,height=6)
StackedVlnPlot(SME_Data_Integrated,features = genes) 
dev.off()

pdf("integration_scRNAseq/PLOTS/Markers_Cells_NonNeuronal_Violin.pdf",width=3,height=6)
ident <- rownames(SME_Data_Integrated@meta.data[SME_Data_Integrated@meta.data$Class == "NonNeuronal",])
nonneu <- SubsetData(SME_Data_Integrated, cells = ident)
StackedVlnPlot(nonneu,features = genes)
dev.off()



genes <- c("SLC17A7",
           "CUX2",           
           "RORB",
           "CBLN2",
           "FOXP2",                       
           "HS3ST2",
           "GAD1",
           "RELN",                      
           "ADRA1B",
           "SST",
           "CCK",
           "LAMP5",
           "SLC14A1",
           "FGFR3",
           "MOBP",
           "MOG",
           "VCAN",
           "SMOC1")

pdf("integration_scRNAseq/PLOTS/Markers_AllMarkers_Violin.pdf",width=5,height=8)
StackedVlnPlot(SME_Data_Integrated,features = genes) 
dev.off()

# Cool Genes
genes <- c(
           "IL1RAPL2",
           "CCDC85A",
           "SMAD3",
           "NEDD4L",
           "TRIO",
           "SYBU")

pdf("integration_scRNAseq/PLOTS/Markers_Cells_CoolGenes_Violin.pdf",width=5,height=5)
StackedVlnPlot(SME_Data_Integrated,features = genes)+ 
scale_fill_viridis(option="plasma")
dev.off()

pdf("integration_scRNAseq/PLOTS/IL1RAPL2_DotPlot.pdf",width=6,height=8)
g <- DotPlot(SME_Data_Integrated, features = c("IL1RAPL2","CCDC85A"))
g$layers[[1]] <- NULL
g <- g + 
    geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
    #scale_color_viridis(option="inferno",direction=-1) + 
    scale_color_gradient(low="white", high="red") +
    guides(color = guide_colorbar(title = 'Average Expression')) + 
    #scale_radius(limits = c(min.pct=0, max.pct=100)) +
    theme(axis.text.x = element_text(angle=45))
print(g)
dev.off()

pdf("integration_scRNAseq/PLOTS/SMAD3_DotPlot.pdf",width=6,height=8)
g <- DotPlot(SME_Data_Integrated, features = c("SMAD3","CCDC85A"))
g$layers[[1]] <- NULL
g <- g + 
    geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
    #scale_color_viridis(option="inferno",direction=-1) + 
    scale_color_gradient(low="white", high="red") +
    guides(color = guide_colorbar(title = 'Average Expression')) + 
    #scale_radius(limits = c(min.pct=0, max.pct=100)) +
    theme(axis.text.x = element_text(angle=45))
print(g)
dev.off()



pdf("integration_scRNAseq/PLOTS/Markers_Cells_IL1RAPL2_Violin.pdf",width=5,height=4)
StackedVlnPlot(SME_Data_Integrated,features = "IL1RAPL2")+ 
scale_fill_viridis(option="plasma")
dev.off()

pdf("integration_scRNAseq/PLOTS/Markers_Cells_SMAD3_Violin.pdf",width=5,height=4)
StackedVlnPlot(SME_Data_Integrated,features = "SMAD3")+ 
scale_fill_viridis(option="plasma")
dev.off()


save(SME_Data_Integrated, file = "integration_scRNAseq/OUTPUTS/SME_Data_Integrated_Labels.RData")



# Create networks
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    wTO  =(cormat)[ut]    
    )
}


Markers_Filt <- Integrated.Markers.Labels.Sign

colnames(SME_Data_Integrated@assays$RNA@data) <- rownames(SME_Data_Integrated@meta.data)

mat <- SME_Data_Integrated@assays$RNA@data %>%
		as.data.frame() %>%
		rownames_to_column("Gene") %>%
		filter(Gene %in% Markers_Filt$gene) %>%
		column_to_rownames("Gene") %>%
		t() %>%
		data.frame() %>%
		rownames_to_column("ID") 
 
pd <- SME_Data_Integrated@meta.data %>%
		rownames_to_column("ID") %>%
		select("ID","Cell")

# Merging data

tmp <- left_join(pd, mat, by = c('ID','ID')) %>%
		select(-ID)

# Create a list of objects by cluster
l <- split(tmp,tmp$Cell)

List_Markers <- split(Markers_Filt,Markers_Filt$cluster)
List_Markers <- List_Markers[match(names(l),names(List_Markers))]

New_Mat <- list()
for (i in 1:length(l))
{
New_Mat[[i]] <- l[[i]] %>%
					select(-Cell) %>% 
					t() %>% 
					as.data.frame() %>%
					rownames_to_column('Gene') %>%
    				filter(Gene %in% List_Markers[[i]]$gene) %>%
    				column_to_rownames('Gene')
}

names(New_Mat) <- names(l)
save(New_Mat,file = "integration_scRNAseq/OUTPUTS/Format_Matrix_ForNet.RData")

# Run correlation g x g in each cluster
List_Cor <- lapply(New_Mat, function(x) cor(t(x),method="pearson"))

# Run the weighted topological overlap
List_wTO <- lapply(List_Cor, function(x) wTO(x,"sign"))
Long_wTO <- lapply(List_wTO, function(x) flattenCorrMatrix(x))
Long_wTO <- lapply(Long_wTO, function(x) x[order(x$wTO, decreasing = TRUE), ])


# Selct top 100 connection to visualized
Top50 <- list()
for (i in 1:length(Long_wTO))
	{
		Top50[[i]] <- Long_wTO[[i]] %>% 
						slice(1:100)
	}
names(Top50) <- names(Long_wTO)

dir.create("integration_scRNAseq/PLOTS/Networks/")
Networks <- list()
layoutFR <- list()
for (i in 1:length(Top50))
	{
		Networks[[i]] <- graph_from_data_frame(Top50[[i]], directed=FALSE)
		Networks[[i]] <- delete.vertices(Networks[[i]],which(degree(Networks[[i]])<1))
		layoutFR[[i]] <- layout_with_fr(Networks[[i]],maxiter = 500)
		V(Networks[[i]])$vertex_degree <-  degree(Networks[[i]])*0.5
		pdf(paste("integration_scRNAseq/PLOTS/Networks/",names(Top50)[[i]], "_CoexpNet.pdf", sep=""),width=10,height=10)
		plot.igraph(Networks[[i]],
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.color="skyblue",
            vertex.size=V(Networks[[i]])$vertex_degree,
            vertex.label.color="black",
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR[[i]],
            edge.color=adjustcolor("grey", alpha.f = .5))
		dev.off()
	}






