library(Seurat)
library(Signac)
library(ggplot2)
library(biomaRt)
library(GenomeInfoDb)
library(GenomicRanges)
library(Matrix)
library(plyr)
library(Matrix.utils)
library(stringi)
library(ggpubr)
library(bedr)
library(reshape2)
library(rio)
set.seed(1234)

### Filter and cluster each sample ###
## Sample 1
counts <- Read10X_h5("~/workdir/STEFANO_SME/cellranger_sample1/cellranger_count_v1.1/outs/raw_peak_bc_matrix.h5")

metadata1 <- read.csv(
  file = "~/workdir/STEFANO_SME/cellranger_sample1/cellranger_count_v1.1/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

sme_atac1 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'sme_atac1',
  min.cells = 1,
  meta.data = metadata1
)

fragpath1 <- "~/workdir/STEFANO_SME/cellranger_sample1/cellranger_count_v1.1/outs/fragments.tsv.gz"

sme_atac1 <- SetFragments(
  object = sme_atac1,
  file = fragpath1
)

# Percentage of reads in peaks
sme_atac1$pct_reads_in_peaks <- sme_atac1$peak_region_fragments / sme_atac1$passed_filters * 100

# Read and filter barcode multiplets
multips1 = read.csv("~/workdir/STEFANO_SME/cellranger_sample1/cellranger_count_v1.1/outs/excluded_barcodes.csv")

# Filter barcode multiplets
sme_atac1[["is_multiplet"]] = 0
sme_atac1[["is_multiplet"]] = ifelse(colnames(sme_atac1) %in% multips1$Excluded.Barcode, 1, 0)


meta = sme_atac1[[]]
meta$Cell = ifelse(meta$passed_filters > 1500 & meta$pct_reads_in_peaks > 15,
		   ifelse(meta$is_multiplet == 0, "Cell", "BC_Multiplet"), "Non-cell")

pdf("InPeaks_InPeaksRatio_sme1.pdf")
ggscatter(data = meta %>% arrange(desc(Cell)),
	 x = 'passed_filters',
	 y = 'pct_reads_in_peaks',
	 xlab = "Total reads in peaks",
	 ylab = "Fraction of reads in peaks",
	 alpha = 1,
	 color = "Cell",
	 palette = c("black", "red", "grey"),
	 size = 1)
dev.off()

# Filter the object
sme_atac1 = subset(sme_atac1, subset = pct_reads_in_peaks > 15 & passed_filters > 1500 & is_multiplet == 0)

# Perform clustering
sme_atac1 <- RunTFIDF(sme_atac1)
sme_atac1 <- FindTopFeatures(sme_atac1, min.cutoff = 'q50')
sme_atac1 <- RunSVD(
  object = sme_atac1,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

sme_atac1 <- RunUMAP(object = sme_atac1, reduction = 'lsi', dims = 1:30)
sme_atac1 <- FindNeighbors(object = sme_atac1, reduction = 'lsi', dims = 1:30)
sme_atac1 <- FindClusters(object = sme_atac1, verbose = FALSE, resolution = 1, algorithm = 3)

pdf("sme_atac1_umap.pdf")
DimPlot(object = sme_atac1, label = TRUE) + NoLegend()
dev.off()

# Save the object
saveRDS(sme_atac1, "~/workdir/STEFANO_SME/sample1_analysis/seur_objects/sme_atac1.RDS")



## Sample 2
counts <- Read10X_h5("~/workdir/STEFANO_SME/cellranger_sample2/outs/raw_peak_bc_matrix.h5")

metadata2 <- read.csv(
  file = "~/workdir/STEFANO_SME/cellranger_sample2/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

fragpath2 <- "~/workdir/STEFANO_SME/cellranger_sample2/outs/fragments.tsv.gz"

sme_atac2 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'sme_atac2',
  min.cells = 1,
  meta.data = metadata2
)

sme_atac2 <- SetFragments(
  object = sme_atac2,
  file = fragpath2
)

# Percentage of reads in peaks
sme_atac2$pct_reads_in_peaks <- sme_atac2$peak_region_fragments / sme_atac2$passed_filters * 100

# Read and filter barcode multiplets
multips2 = read.csv("~/workdir/STEFANO_SME/cellranger_sample2/outs/excluded_barcodes.csv")

# Filter barcode multiplets
sme_atac2[["is_multiplet"]] = 0
sme_atac2[["is_multiplet"]] = ifelse(colnames(sme_atac2) %in% multips2$Excluded.Barcode, 1, 0)


meta = sme_atac2[[]]
meta$Cell = ifelse(meta$passed_filters > 1500 & meta$pct_reads_in_peaks > 15,
		   ifelse(meta$is_multiplet == 0, "Cell", "BC_Multiplet"), "Non-cell")

pdf("InPeaks_InPeaksRatio_sme2.pdf")
ggscatter(data = meta %>% arrange(desc(Cell)),
	 x = 'passed_filters',
	 y = 'pct_reads_in_peaks',
	 xlab = "Total reads in peaks",
	 ylab = "Fraction of reads in peaks",
	 alpha = 1,
	 color = "Cell",
	 palette = c("black", "red", "grey"),
	 size = 1)
dev.off()

# Filter the object
sme_atac2 = subset(sme_atac2, subset = pct_reads_in_peaks > 15 & passed_filters > 1500 & is_multiplet == 0)

# Cluster the cells
sme_atac2 <- RunTFIDF(sme_atac2)
sme_atac2 <- FindTopFeatures(sme_atac2, min.cutoff = 'q50')
sme_atac2 <- RunSVD(
  object = sme_atac2,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

sme_atac2 <- RunUMAP(object = sme_atac2, reduction = 'lsi', dims = 1:30)
sme_atac2 <- FindNeighbors(object = sme_atac2, reduction = 'lsi', dims = 1:30)
sme_atac2 <- FindClusters(object = sme_atac2, verbose = FALSE, resolution = 1, algorithm = 3)

DimPlot(object = sme_atac2, label = TRUE) + NoLegend()
dev.off()
pdf("sme_atac2_umap.pdf", width = 20, height = 10)
CombinePlots(plots = list(p1, p2))
dev.off()

# Save the object
saveRDS(sme_atac2, "~/workdir/STEFANO_SME/allsamples/seurat_objects/sme_atac2.RDS")


## Sample 3
counts <- Read10X_h5("~/workdir/STEFANO_SME/cellranger_sample3/outs/raw_peak_bc_matrix.h5")

metadata3 <- read.csv(
  file = "~/workdir/STEFANO_SME/cellranger_sample3/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

fragpath3 <- "~/workdir/STEFANO_SME/cellranger_sample3/outs/fragments.tsv.gz"

sme_atac3 <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'sme_atac3',
  min.cells = 1,
  meta.data = metadata3
)

sme_atac3 <- SetFragments(
  object = sme_atac3,
  file = fragpath3
)

# Percentage of reads in peaks
sme_atac3$pct_reads_in_peaks <- sme_atac3$peak_region_fragments / sme_atac3$passed_filters * 100

# Read and filter barcode multiplets
multips3 = read.csv("~/workdir/STEFANO_SME/cellranger_sample3/outs/excluded_barcodes.csv")

# Filter barcode multiplets
sme_atac3[["is_multiplet"]] = 0
sme_atac3[["is_multiplet"]] = ifelse(colnames(sme_atac3) %in% multips3$Excluded.Barcode, 1, 0)


meta = sme_atac3[[]]
meta$Cell = ifelse(meta$passed_filters > 1500 & meta$pct_reads_in_peaks > 15,
		   ifelse(meta$is_multiplet == 0, "Cell", "BC_Multiplet"), "Non-cell")

pdf("InPeaks_InPeaksRatio_sme3.pdf")
ggscatter(data = meta %>% arrange(desc(Cell)),
	 x = 'passed_filters',
	 y = 'pct_reads_in_peaks',
	 xlab = "Total reads in peaks",
	 ylab = "Fraction of reads in peaks",
	 alpha = 1,
	 color = "Cell",
	 palette = c("black", "red", "grey"),
	 size = 1)
dev.off()

# Filter the object
sme_atac3 = subset(sme_atac3, subset = pct_reads_in_peaks > 15 & passed_filters > 1500 & is_multiplet == 0)

# Cluster the cells
sme_atac3 <- RunTFIDF(sme_atac3)
sme_atac3 <- FindTopFeatures(sme_atac3, min.cutoff = 'q50')
sme_atac3 <- RunSVD(
  object = sme_atac3,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

sme_atac3 <- RunUMAP(object = sme_atac3, reduction = 'lsi', dims = 1:30)
sme_atac3 <- FindNeighbors(object = sme_atac3, reduction = 'lsi', dims = 1:30)
sme_atac3 <- FindClusters(object = sme_atac3, verbose = FALSE, resolution = 1, algorithm = 3)

pdf("sme_atac3.pdf")
DimPlot(object = sme_atac3, label = TRUE) + NoLegend()
dev.off()

# Save the object
saveRDS(sme_atac3, "~/workdir/STEFANO_SME/allsamples/seurat_objects/sme_atac3.RDS")


# Get union of all peaks
sme_atac1 = readRDS("~/workdir/STEFANO_SME/allsamples/seurat_objects/sme_atac1.RDS")
sme_atac2 = readRDS("~/workdir/STEFANO_SME/allsamples/seurat_objects/sme_atac2.RDS")
sme_atac3 = readRDS("~/workdir/STEFANO_SME/allsamples/seurat_objects/sme_atac3.RDS")


smerows = rbind(rownames(sme_atac1) %>% gsub(":", "-", .) %>% strsplit(., "-") %>% do.call(rbind, .) %>% as.data.frame(),
	    rownames(sme_atac2) %>% gsub(":", "-", .) %>% strsplit(., "-") %>% do.call(rbind, .) %>% as.data.frame(),
	    rownames(sme_atac3) %>% gsub(":", "-", .) %>% strsplit(., "-") %>% do.call(rbind, .) %>% as.data.frame())
smerows$V1 = as.character(smerows$V1)
smerows$V2 = as.numeric(as.character(smerows$V2))
smerows$V3 = as.numeric(as.character(smerows$V3))
colnames(smerows) = c("chr", "start", "end")

smerows = bedr.sort.region(smerows)
smerows_merg = bedr.merge.region(smerows)

comb_peaks = makeGRangesFromDataFrame(smerows_merg)


# Re-make the peak-cell matrices
sme_atac1_cnt <- FeatureMatrix(
  fragments = GetFragments(sme_atac1),
  features = comb_peaks,
  sep = c(":", "-"),
  cells = colnames(sme_atac1)
)

sme_atac2_cnt <- FeatureMatrix(
  fragments = GetFragments(sme_atac2),
  features = comb_peaks,
  sep = c(":", "-"),
  cells = colnames(sme_atac2)
)

sme_atac3_cnt <- FeatureMatrix(
  fragments = GetFragments(sme_atac3),
  features = comb_peaks,
  sep = c(":", "-"),
  cells = colnames(sme_atac3)
)

# Recreate the seurat objects

smeall_atac1 <- CreateSeuratObject(
  counts = sme_atac1_cnt,
  assay = 'peaks',
  project = 'sme_atac1',
  min.cells = 1,
  meta.data = metadata1
)

smeall_atac1 <- SetFragments(
  object = smeall_atac1,
  file = fragpath1
)

smeall_atac2 <- CreateSeuratObject(
  counts = sme_atac2_cnt,
  assay = 'peaks',
  project = 'sme_atac2',
  min.cells = 1,
  meta.data = metadata2
)

smeall_atac2 <- SetFragments(
  object = smeall_atac2,
  file = fragpath2
)


smeall_atac3 <- CreateSeuratObject(
  counts = sme_atac3_cnt,
  assay = 'peaks',
  project = 'sme_atac3',
  min.cells = 1,
  meta.data = metadata3
)

smeall_atac3 <- SetFragments(
  object = smeall_atac3,
  file = fragpath3
)

# First, cluster each sample separately
all_seur = c(smeall_atac1, smeall_atac2, smeall_atac3)
plt = list()
for(i in 1:length(all_seur)){

	print(i)

	all_seur[[i]] <- RunTFIDF(all_seur[[i]])
	all_seur[[i]] <- FindTopFeatures(all_seur[[i]], min.cutoff = 'q50')
	all_seur[[i]] <- RunSVD(
	  object = all_seur[[i]],
	  assay = 'peaks',
	  reduction.key = 'LSI_',
	  reduction.name = 'lsi'
	)

	all_seur[[i]] <- RunUMAP(object = all_seur[[i]], reduction = 'lsi', dims = 1:30)
	all_seur[[i]] <- FindNeighbors(object = all_seur[[i]], reduction = 'lsi', dims = 1:30)
	all_seur[[i]] <- FindClusters(object = all_seur[[i]], verbose = FALSE, resolution = 1)

	plt[[i]] = DimPlot(object = all_seur[[i]], label = TRUE) + NoLegend()

	print(i)
}

saveRDS(all_seur, "~/workdir/STEFANO_SME/allsamples/seurat_objects/all_seur.RDS")

# Give each sample's cells unique names
all_seur[[1]] = RenameCells(object = all_seur[[1]], new.names = paste0("S1_", colnames(x = all_seur[[1]])))
all_seur[[2]] = RenameCells(object = all_seur[[2]], new.names = paste0("S2_", colnames(x = all_seur[[2]])))
all_seur[[3]] = RenameCells(object = all_seur[[3]], new.names = paste0("S3_", colnames(x = all_seur[[3]])))
names(all_seur) = c("s1", "s2", "s3")

smeall_atac1 = all_seur[[1]]
smeall_atac2 = all_seur[[2]]
smeall_atac3 = all_seur[[3]]

# Choose top 10 percent variable features for integration
smebig <- merge(smeall_atac1, merge(smeall_atac2, smeall_atac3))
smebig = FindTopFeatures(smebig, min.cutoff = 'q90')
peaks_use <- VariableFeatures(smebig)
peaks_use = intersect(peaks_use, rownames(smeall_atac1)) # Making sure the peak is not entirely zero for any samples.
peaks_use = intersect(peaks_use, rownames(smeall_atac2))
peaks_use = intersect(peaks_use, rownames(smeall_atac3))

length(peaks_use)

# Find integration anchors between datasets
anchors <- FindIntegrationAnchors(
  object.list = all_seur,
  anchor.features = peaks_use,
  assay = c('peaks', 'peaks', 'peaks'),
  dims = 1:30,
  k.filter = NA
)

gc()

# Integrate data and create a new merged object
integrated <- IntegrateData(
  anchorset = anchors,
  preserve.order = F,
  dims = 1:30
)


# Run LSI again on the adjusted matrix and perform dimentionality reduction and clustering.
integrated <- RunSVD(
  object = integrated,
  n = 30,
  reduction.name = 'integratedLSI'
)

integrated <- RunUMAP(
  object = integrated,
  dims = 1:30,
  reduction = 'integratedLSI'
)

integrated <- FindNeighbors(object = integrated, reduction = 'integratedLSI', dims = 1:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, resolution = 1)

saveRDS(integrated, "~/workdir/STEFANO_SME/allsamples_integrated_66k.RDS")

pdf("Integrated_sme_clusters.pdf")
DimPlot(integrated, pt.size = 0.1, label = T, label.size = 7) +
theme(axis.text.x = element_text(size=20),
	axis.text.y = element_text(size=20),
	axis.title = element_text(size=20),
	legend.position= "none")

dev.off()

pdf("Integrated_sme_samples.pdf")
DimPlot(integrated, pt.size = 0.1, label = T, group.by = 'orig.ident', label.size = 7) +
theme(axis.text.x = element_text(size=20),
	axis.text.y = element_text(size=20),
	axis.title = element_text(size=20),
	legend.position= "none")

dev.off()


# Stacked barplot of each sample per cluster
meta = integrated[[]]
hmm = lapply(names(table(meta$seurat_clusters)), 
		function(x){table(meta[meta$seurat_clusters == x, "orig.ident"])})
hmm = do.call(rbind, hmm) %>% as.data.frame()
hmm$cluster = names(table(meta$seurat_clusters))
hmm = melt(hmm, id.vars='cluster')
colnames(hmm)[2] = "orig.ident"


# Add all cell numbers
hmm = table(meta$orig.ident) %>% as.data.frame() %>%
	cbind(cluster = "Total", .) %>%
	dplyr::rename(orig.ident = Var1, value = Freq) %>% rbind(hmm, .)

pdf("clusters_stack_samples.pdf", width = 15, height = 5)
ggplot(hmm, aes(x = cluster, y = value, fill = orig.ident)) + 
    geom_bar(position = "fill",stat = "identity") +
    xlab("Clusters") + 
    ylab("Percentage") +
    theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=20),
	axis.text.y = element_text(size=20),
	axis.title = element_text(size=20),
	legend.position= "none") +
    scale_y_continuous(labels = scales::percent_format())

dev.off()


### ANNOTATION ##
# Read integrated seurat object and gene activity matrix
integrated = readRDS("~/workdir/STEFANO_SME/allsamples/seurat_objects/allsamples_integrated.RDS")
geneact = readRDS("~/workdir/STEFANO_SME/allsamples/seurat_objects/allsamples_merged_geneact_matrix.RDS")

# Match gene activity matrix cell names and order
int_names = colnames(integrated@assays$peaks@counts)
geneact = geneact[,match(int_names, colnames(geneact))]

# add the gene RNA matrix to the Seurat object as a new assay, and normalize it
integrated[['RNA']] <- CreateAssayObject(counts = geneact)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)

DefaultAssay(integrated) = "RNA"

# Integrate with RNA
load("~/workdir/STEFANO_SME/scrnaseq_from_stefano/SME_Data_Integrated_Labels.RData")
SME_Data_Integrated = FindVariableFeatures(SME_Data_Integrated, nfeatures = 2000)

# Check the common features to be used for label transfer
sum(VariableFeatures(SME_Data_Integrated) %in% rownames(integrated@assays$RNA@counts))

# Find anchors
transfer.anchors <- FindTransferAnchors(reference = SME_Data_Integrated, query = integrated,
					dims = 1:20,features = VariableFeatures(SME_Data_Integrated),
					reference.assay = "RNA", query.assay = "RNA",
					reduction = "cca")

# Make cell type predictions based on anchors
celltype.predictions <- TransferData(anchorset = transfer.anchors,
				     refdata = Idents(SME_Data_Integrated), 
				     weight.reduction = integrated[["integratedLSI"]])

# Add predictions to meta data
integrated <- AddMetaData(integrated, metadata = celltype.predictions)

pdf("predictionscores.pdf")
hist(celltype.predictions$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()

## Plot all cells ##
integrated_filt <- subset(integrated, subset = prediction.score.max > 0)

# Make the colors match
integrated_filt$predicted.id <- factor(integrated_filt$predicted.id, levels = levels(SME_Data_Integrated))

# Plot scATACseq and scRNAseq side by side
p1 <- DimPlot(integrated_filt, group.by = "predicted.id", label = TRUE, repel = TRUE) +
		ggtitle("scATAC-seq cells") + 
		NoLegend() +
		scale_colour_hue(drop = FALSE)

p2 <- DimPlot(SME_Data_Integrated, group.by = "Cell", label = TRUE, repel = TRUE) +
		ggtitle("scRNA-seq cells") +
		NoLegend()

pdf("Integrated_annot.pdf", width = 20, height = 10)
CombinePlots(plots = list(p1, p2))
dev.off()


# Filter clusters. This part is empirical based on lack of clear prediction
annot_conf = integrated
keep_clusters = as.character(c(0,3,4,5,6,7,8,9,10,13,14,16,19,21,24,28,29))

# Update filtered object
Idents(annot_conf) = "seurat_clusters"
annot_conf = subset(annot_conf, idents = keep_clusters)


# Recluster after removing
annot_conf <- RunSVD(
  object = annot_conf,
  n = 30,
  reduction.name = 'integratedLSI'
)

annot_conf <- RunUMAP(
  object = annot_conf,
  dims = 1:30,
  reduction = 'integratedLSI'
)

annot_conf <- FindNeighbors(object = annot_conf, reduction = 'integratedLSI', dims = 1:30)
annot_conf <- FindClusters(object = annot_conf, verbose = FALSE, resolution = 1, graph.name = 'integrated_snn')

pdf("filt_reclustered_sample.pdf")
DimPlot(annot_conf, group.by = "orig.ident", label = TRUE, repel = TRUE)
dev.off()

pdf("filt_reclustered.pdf")
DimPlot(annot_conf, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
dev.off()



meta = annot_conf[[]]
hmm = lapply(names(table(meta$predicted.id)), 
		function(x){table(meta[meta$predicted.id == x, "orig.ident"])})
hmm = do.call(rbind, hmm) %>% as.data.frame()
hmm$cluster = names(table(meta$predicted.id))
hmm = melt(hmm, id.vars='cluster')
colnames(hmm)[2] = "orig.ident"


# Add all cell numbers
hmm = table(meta$orig.ident) %>% as.data.frame() %>%
	cbind(cluster = "Total", .) %>%
	dplyr::rename(orig.ident = Var1, value = Freq) %>% rbind(hmm, .)

pdf("annotated_prediction_stack_sample.pdf", width = 15, height = 5)
ggplot(hmm, aes(x = cluster, y = value, fill = orig.ident)) + 
    geom_bar(position = "fill",stat = "identity") +
    xlab("Clusters") + 
    ylab("Percentage") +
    theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=20),
	axis.text.y = element_text(size=20),
	axis.title = element_text(size=20)) +
    scale_y_continuous(labels = scales::percent_format())

dev.off()




# Replot with the filtered cells
annot_conf_filt <- subset(annot_conf, subset = prediction.score.max > 0)
annot_conf_filt$predicted.id <- factor(annot_conf_filt$predicted.id,
						levels = levels(SME_Data_Integrated))

p1 <- DimPlot(annot_conf_filt, group.by = "predicted.id", label = TRUE, repel = TRUE) +
		ggtitle("scATAC-seq cells") + 
		NoLegend() +
		scale_colour_hue(drop = FALSE)

p2 <- DimPlot(SME_Data_Integrated, group.by = "Cell", label = TRUE, repel = TRUE) +
		ggtitle("scRNA-seq cells") + 
		NoLegend()

pdf("integrated_annot_filt_reclustered.pdf", width = 20, height = 10)
CombinePlots(plots = list(p1, p2))
dev.off()


# Heatmap of prediction by cluster #

# Calculate total number of cells for given predicted-cluster id match
filtlabel = annot_conf[[c("predicted.id", "seurat_clusters")]]
filtlabel$total = 1
tmp = aggregate(. ~ predicted.id + seurat_clusters, filtlabel, sum)

# Add total number of cells with the given predicted id
predid_tab = table(filtlabel$predicted.id)
matchind = match(tmp$predicted.id, names(predid_tab))
tmp$predicted_id_tot = predid_tab[matchind]

# Add total number of cells with the given cluster id
seurid_tab = table(filtlabel$seurat_clusters)
matchind = match(tmp$seurat_clusters, names(seurid_tab))
tmp$seurat_id_tot = seurid_tab[matchind]

# Calculate percentage of label for the given cluster. Keep only the ones that are > 20%
tmp$perc_inclus = tmp$total / tmp$seurat_id_tot
tmp2 = tmp[tmp$perc_inclus > 0.05,]

# Heatmap of annotation
toplot = tmp2[, c('predicted.id', 'seurat_clusters', 'perc_inclus')]
toplot$predicted.id = factor(toplot$predicted.id)
toplot$seurat_clusters = factor(toplot$seurat_clusters)
toplot$perc_inclus = as.numeric(as.character(toplot$perc_inclus))

pdf("predicted_heatmap.pdf", width = 20, height =10)
ggplot(toplot, aes(seurat_clusters, predicted.id, fill = perc_inclus)) +
	geom_tile(color = "white") +	
	scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
	midpoint = 0.05, limit = c(0,1)) +
	theme_classic() +
	theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
		rotate_x_text(90)
dev.off()

# Give each cluster names based on the predicted ID
map = setNames(c("Oligo_L1-6_OPALIN", "OPC_L1-6_PDGFRA", "Exc_L2-3_LINC00507_FREM3",
		 "Oligo_L1-6_OPALIN", "Astro_L1-6_FGFR3", "Exc_L2-3_LINC00507_FREM3",
		 "Exc_L3-5_RORB", "Inh_L1-2_SST-LAMP5", "Inh_L5-6_SST_KLHDC8A",
		 "Inh_L4-6_PVALB_SURF1", "Inh_L2-4_VIP", "Exc_L4-6_RORB_SEMA3A",
		 "Exc_L5-6_THEMIS_C1QL3", "Exc_L5-6_FEZF2_ABO", "Exc_L5-6_THEMIS_CRABP1",
		 "Exc_L5-6_FEZF2_ABO", "OPC_L1-6_PDGFRA", "Exc_L5-6_FEZF2_IL26"), c(0:17))

annot_conf[["annotation"]] <- map[unlist(annot_conf[["seurat_clusters"]])]
Idents(annot_conf) = annot_conf[["annotation"]]

pdf("FINAL_integrated_annot_filt_reclustered.pdf")
DimPlot(annot_conf, group.by = "annotation", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
dev.off()

saveRDS(annot_conf, "~/workdir/STEFANO_SME/annot_conf.RDS")


### Filter module genes and extract module clusters ###
mods = import("~/workdir/STEFANO_SME/sample1_analysis/table_stefano/Table_S3.xlsx")
m4 = mods[mods$ModuleName == "WM4",]
m12 = mods[mods$ModuleName == "WM12",]
m21 = mods[mods$ModuleName == "WM21",]

scATACseq = readRDS("~/workdir/STEFANO_SME/allsamples/seurat_objects/annot_conf.RDS")
DefaultAssay(scATACseq) = "RNA"

## Filter module genes for presence in their own clusters in scATACseq ##

# Module 4 and Module 12
# Filter for presence in atac-seq (>30%)
cl_names = names(table(Idents(scATACseq)))
ids = cl_names[grepl("Exc|Inh", cl_names)]
ExcInh_atac = subset(scATACseq, idents = ids)

m4_ov = m4$Gene[m4$Gene %in% rownames(ExcInh_atac)]
m12_ov = m12$Gene[m12$Gene %in% rownames(ExcInh_atac)]

m4_atac = ExcInh_atac@assays$RNA@counts[m4_ov,]
m12_atac = ExcInh_atac@assays$RNA@counts[m12_ov,]

m4_zero = apply(as.matrix(m4_atac), 1, function(x){sum(x == 0)})
m12_zero = apply(as.matrix(m12_atac), 1, function(x){sum(x == 0)})

m4_zero_filt = m4_zero[((ncol(m4_atac) - m4_zero) / ncol(m4_atac)) > 0.3]
m12_zero_filt = m12_zero[((ncol(m12_atac) - m12_zero) / ncol(m12_atac)) > 0.3]


# Write cells in wanted clusters. This is to call the peaks again.
to_peak_call = WhichCells(scATACseq, idents = names(table(Idents(ExcInh_atac))))
write.table(to_peak_call, "~/workdir/STEFANO_SME/allsamples/ExcInh_cells.txt", sep = "\n", quote = F,
	    col.names = F, row.names = F)

# Module 21
# Filter for presence in atac-seq (>30%)
cl_names = names(table(Idents(scATACseq)))
ids = cl_names[grepl("Oligo|OPC", cl_names)]
OligOPC_atac = subset(scATACseq, idents = ids)
m21_ov = m21$Gene[m21$Gene %in% rownames(OligOPC_atac)]
m21_atac = OligOPC_atac@assays$RNA@counts[m21_ov,]
m21_zero = apply(as.matrix(m21_atac), 1, function(x){sum(x == 0)})
m21_zero_filt = m21_zero[((ncol(m21_atac) - m21_zero) / ncol(m21_atac)) > 0.3]

# Save all surviving module genes
genes = c(names(m4_zero_filt), names(m12_zero_filt), names(m21_zero_filt))
class = c(rep("WM4", length(m4_zero_filt)), rep("WM12", length(m12_zero_filt)), rep("WM21", length(m21_zero_filt)))
towrite_df = data.frame(Gene = genes, Class = class)

export(towrite_df, "~/workdir/STEFANO_SME/allsamples/Mods_ATAC_0.3.xlsx")


# Write cells in wanted clusters. This is to call the peaks again.
to_peak_call = WhichCells(scATACseq, idents = names(table(Idents(OligOPC_atac))))
write.table(to_peak_call, "~/workdir/STEFANO_SME/allsamples/OligOPC_cells.txt", sep = "\n", quote = F,
	    col.names = F, row.names = F)


## Find TFs with motifs and sufficiently accessible in ExcInh and/or OligOPC clusters ##

# Use all TFs from CISBP for now. Later we will filter for only confident TF motifs.
motif_dat = read.table("~/workdir/motif_databases/CISBP_HUMAN_2019/TF_Information_all_motifs_plus.txt",
		header = T, sep = "\t", na.strings = "NA", fill = T)
all_TF_names = unique(motif_dat$TF_Name)

# Filter for presence in atac-seq. M4 and M12
TF_ov = all_TF_names[all_TF_names %in% rownames(ExcInh_atac@assays$RNA@counts)]
TF_atac = ExcInh_atac@assays$RNA@counts[TF_ov,]
zero_TF = apply(as.matrix(TF_atac), 1, function(x){sum(x == 0)})
m4_12_zero_TF_filt = zero_TF[((ncol(TF_atac) - zero_TF) / ncol(TF_atac)) > 0.3]

# Filter for presence in atac-seq. M21
TF_ov = all_TF_names[all_TF_names %in% rownames(OligOPC_atac@assays$RNA@counts)]
TF_atac = OligOPC_atac@assays$RNA@counts[TF_ov,]
zero_TF = apply(as.matrix(TF_atac), 1, function(x){sum(x == 0)}) # Count #of cells for which given TF is not detected
m21_zero_TF_filt = zero_TF[((ncol(TF_atac) - zero_TF) / ncol(TF_atac)) > 0.3] # Keep the ones detected at least 30% of the cells

genes = c(names(m4_12_zero_TF_filt), names(m21_zero_TF_filt))
class = c(rep("Exc_Inh", length(m4_12_zero_TF_filt)), rep("Olig_OPC", length(m21_zero_TF_filt)))
towrite_df = data.frame(Gene = genes, Class = class)

export(towrite_df, "~/workdir/STEFANO_SME/allsamples/TFs_Annot.xlsx")


