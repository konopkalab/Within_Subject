library(Seurat)
library(Signac)
library(ggplot2)
library(biomaRt)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)

####################################################
######  Process Cellranger's peak-cell matrix  #####
####################################################
counts <- Read10X_h5("~/konopkaLab/pr3/signac_tmp/sme_atac/filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "~/konopkaLab/pr3/signac_tmp/sme_atac/singlecell.csv",
  header = TRUE,
  row.names = 1
)

sme_atac <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'sme_atac',
  min.cells = 1,
  meta.data = metadata
)

fragment.path <- "~/konopkaLab/pr3/signac_tmp/sme_atac/fragments.tsv.gz"

sme_atac <- SetFragments(
  object = sme_atac,
  file = fragment.path
)

#To plot band density of read sizes
sme_atac <- NucleosomeSignal(object = sme_atac)
#Percentage of reads in peaks
sme_atac$pct_reads_in_peaks <- sme_atac$peak_region_fragments / sme_atac$total * 100
#Percentage of blacklisted fragments in peaks
sme_atac$blacklist_ratio <- sme_atac$blacklist_region_fragments / sme_atac$peak_region_fragments

#Plot the initial quality controls
plot1 <- VlnPlot(
  object = sme_atac,
  features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1) + NoLegend()

plot2_a <- VlnPlot(
  object = sme_atac,
  features = 'peak_region_fragments',
  pt.size = 0.1, log = TRUE) + NoLegend()

plot2_b <- FeatureScatter(sme_atac,"peak_region_fragments",'nucleosome_signal', pt.size = 0.1) + NoLegend()
plot2_c <- FeatureScatter(sme_atac,"peak_region_fragments",'blacklist_ratio', pt.size = 0.1) + NoLegend()
plot2 <- CombinePlots(plots = list(plot2_a,plot2_b,plot2_c), ncol = 3)

pdf("~/konopkaLab/pr3/signac_tmp/sme_atac/sme_atac_qc.pdf")
CombinePlots(list(plot1,plot2),ncol = 1)
dev.off()

#Calculate mononucleosome / nucleosome-free ratio
sme_atac$nucleosome_group <- ifelse(sme_atac$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

pdf("~/konopkaLab/pr3/signac_tmp/sme_atac/sme_atac_bandingPattern.pdf")
PeriodPlot(object = sme_atac, group.by = 'nucleosome_group')
dev.off()

#Filter for good quality cells
sme_atac <- subset(sme_atac, subset = peak_region_fragments > 1000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10)


# Clustering
sme_atac <- RunTFIDF(sme_atac)
sme_atac <- FindTopFeatures(sme_atac, min.cutoff = 'q0')
sme_atac <- RunSVD(
  object = sme_atac,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

sme_atac <- RunUMAP(object = sme_atac, reduction = 'lsi', dims = 1:30)
sme_atac <- FindNeighbors(object = sme_atac, reduction = 'lsi', dims = 1:30)
sme_atac <- FindClusters(object = sme_atac, verbose = FALSE)

saveRDS(sme_atac, "~/konopkaLab/pr3/signac_tmp/sme_atac/sme_atac_clustered_seurat.RDS")

pdf("~/sme_atac_cellranger_dim30_umap.pdf")
DimPlot(object = sme_atac, label = TRUE) +
	NoLegend() +
	ggtitle("sme_atac-seq sme_atac BA38") + 
	theme(plot.title = element_text(hjust = 0.5))
dev.off()


#Create activity matrix. Do it sequentially to overcome the RAM issue.
sme_atac = readRDS("~/workdir/pr3/sme_atacseq/human/sme_atac_SURGICAL/seurat_objects/sme_atac_clustered_seurat.RDS")
sme_atac_peakmat = sme_atac@assays$peaks@counts
sme_atac_peakmat_1 = sme_atac_peakmat[,1:5000]
sme_atac_peakmat_2 = sme_atac_peakmat[,5001:10000]
sme_atac_peakmat_3 = sme_atac_peakmat[,10001:15000]
sme_atac_peakmat_4 = sme_atac_peakmat[,15001:ncol(sme_atac_peakmat)]


sme_atac_activ_1 = CreateGeneActivityMatrix(peak.matrix = sme_atac_peakmat_1,
					    annotation.file = "/home2/s422159/workdir/reference_genomes/gtf/Homo_sapiens.GRCh38.96.gtf", 
    				 	    seq.levels = c(1:22, "X", "Y"), upstream = 2000,
					    verbose = TRUE)

sme_atac_activ_2 = CreateGeneActivityMatrix(peak.matrix = sme_atac_peakmat_2,
					    annotation.file = "/home2/s422159/workdir/reference_genomes/gtf/Homo_sapiens.GRCh38.96.gtf", 
    				 	    seq.levels = c(1:22, "X", "Y"), upstream = 2000,
					    verbose = TRUE)

sme_atac_activ_3 = CreateGeneActivityMatrix(peak.matrix = sme_atac_peakmat_3,
					    annotation.file = "/home2/s422159/workdir/reference_genomes/gtf/Homo_sapiens.GRCh38.96.gtf", 
    				 	    seq.levels = c(1:22, "X", "Y"), upstream = 2000,
					    verbose = TRUE)


sme_atac_activ_4 = CreateGeneActivityMatrix(peak.matrix = sme_atac_peakmat_4,
					    annotation.file = "/home2/s422159/workdir/reference_genomes/gtf/Homo_sapiens.GRCh38.96.gtf", 
    				 	    seq.levels = c(1:22, "X", "Y"), upstream = 2000,
					    verbose = TRUE)

sme_atac_activ = cbind(sme_atac_activ_1, sme_atac_activ_2, sme_atac_activ_3, sme_atac_activ_4)


sme_atac[["RNA"]] <- CreateAssayObject(counts = sme_atac_activ)

#Normalize the gene activity
sme_atac <- NormalizeData(
  object = sme_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(sme_atac$nCount_RNA)
)

DefaultAssay(sme_atac) <- 'RNA'

pdf("TFs_featurePlot.pdf", width = 10, height = 10)
FeaturePlot(
  object = sme_atac,
  features = c('EZH2', 'SUZ12', 'REST'),
  pt.size = 0.1,
  ncol = 2
)
dev.off()



###########################################
### Integration of sme atac and sme rna ###
###########################################
load("~/workdir/STEFANO_SME/seuObjSME_Norm_PC_Clust_tSNE_UMAP_ReMap.RData")

transfer.anchors <- FindTransferAnchors(reference = seuObjSME_Norm_PC_Clust, query = sme_atac,
					dims = 1:20,features = VariableFeatures(seuObjSME_Norm_PC_Clust),
					reference.assay = "RNA", query.assay = "RNA",
					reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors,
				     refdata = Idents(seuObjSME_Norm_PC_Clust), 
				     weight.reduction = sme_atac[["lsi"]])

sme_atac <- AddMetaData(sme_atac, metadata = celltype.predictions)

seuObjSME_Norm_PC_Clust$cell_type = Idents(seuObjSME_Norm_PC_Clust)

pdf("predictionscores_sme_atac_sme_atac_humanrna_sme_atac_varfeat2k_overlap2k.pdf")
hist(sme_atac$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()

sme_atac_filt <- subset(sme_atac, subset = prediction.score.max > 0.5)
sme_atac_filt$predicted.id <- factor(sme_atac_filt$predicted.id, levels = levels(seuObjSME_Norm_PC_Clust))  # to make the colors match
p1 <- DimPlot(sme_atac_filt, group.by = "predicted.id", label = TRUE, repel = TRUE) +
	ggtitle("scsme_atac-seq cells\nprediction > 0.5") + 
	NoLegend() +
	scale_colour_hue(drop = FALSE) +
	theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(seuObjSME_Norm_PC_Clust, group.by = "cell_type", label = TRUE, repel = TRUE) +
	ggtitle("scRNA-seq cells") + 
	theme(plot.title = element_text(hjust = 0.5)) +
	NoLegend()

pdf("Integr_smeATAC_smeRNA_filt_pred0.5_varfeat2k_overlap2k.pdf", width = 20, height = 10)
CombinePlots(plots = list(p1, p2))
dev.off()


#Annotate cells
map = setNames(c("Oligo_1", "Oligo_2", "Microglia_1",
		 "OPC", "Exc_2", "Exc_4",
		 "Astro", "Exc_3", "Exc_1_1",
		 "Inh_2_1", "Inh_3", "Exc_Unknown",
		 "Inh_1", "Inh_2_2", "Inh_4",
		 "Exc_1_2", "Exc_6", "Microglia_2",
		 "Exc_5_1", "Exc_5_2", "Exc_8",
		 "Exc_7"), c(0:21))

sme_atac[["seurat_clusters"]] <- map[unlist(sme_atac[["seurat_clusters"]])]
Idents(sme_atac) = sme_atac[["seurat_clusters"]]

sme_atac_microglia_removed = subset(sme_atac, idents = c("Microglia_1", "Microglia_2"), invert = T)

saveRDS(sme_atac_microglia_removed,
	file = "~/workdir/STEFANO_SME/sme_atac_microglia_removed_annotated_seur.RDS")





