library(rio)
library(Signac)
library(Seurat)
library(JASPAR2018)
library(TFBSTools)
library("BSgenome.Hsapiens.UCSC.hg38")
library(universalmotif)
library(bedr)
library(Matrix.utils)
library(ggplot2)
library(rio)
library(ggrepel)
library(bedr)

#################################################################
### Filter TF Motifs based on TF's presence in the given data ###
#################################################################

# Get a list of motif position weight matrices from the CISBP database (in jaspar format)
cisbp_hum = readJASPARMatrix("~/workdir/STEFANO_SME/sample1_analysis/EXCINH_bedpe_analysis/motif.jaspar")
cisbp_hum = toPWM(cisbp_hum, type = "prob")

# Take only the TFs in our clusters
TFs = rio::import("~/workdir/STEFANO_SME/allsamples/TFs_Annot.xlsx")
ExcInh_TFs = TFs[TFs$Class == "Exc_Inh", "Gene"]
OligOPC_TFs = TFs[TFs$Class == "Olig_OPC", "Gene"]


# Take only motifs with direct or best prediction TFs.
tmp_list = strsplit(names(cisbp_hum), "_")
theTFs = lapply(tmp_list, function(x){x[2]})
theTFs = unlist(theTFs)

cisbp_hum_ExcInh = cisbp_hum[theTFs %in% ExcInh_TFs]
cisbp_hum_OligOPC = cisbp_hum[theTFs %in% OligOPC_TFs]



#####################################
### WM4 & WM12 - Motif Enrichment ###
#####################################

counts1 <- readMM("~/workdir/STEFANO_SME/allsamples/peaks_ExcInh/sample1_pc_matrix/pc_mat.mtx")
counts1 = t(counts1)
pks1 <- read.table("~/workdir/STEFANO_SME/allsamples/peaks_ExcInh/sample1_pc_matrix/peaks.txt")
brc1 <- read.table("~/workdir/STEFANO_SME/allsamples/peaks_ExcInh/sample1_pc_matrix/barcodes.tsv")

colnames(counts1) = brc1$V1
rownames(counts1) = pks1$V1


# Only keep the upstream peaks per sample

annot_peaks = read.table("~/workdir/STEFANO_SME/allsamples/SME_ExcInh/peak_annotation.tsv", sep = "\t", header = T, stringsAsFactors = F)

tmp = annot_peaks$distance
tmplist = strsplit(tmp, ";")
tmpbool = lapply(tmplist, function(x){sum(as.numeric(x) >= 0) == 0})
tmpbool = unlist(tmpbool)

annot_peaks = annot_peaks[tmpbool,]
tmp = sub("_", ":", annot_peaks$peak)
tmp2 = sub("_", "-", tmp)
annot_peaks$peak = tmp2

# Subset same number of peaks from counts matrix.
# We will switch the rownames of a random count matrix with these peaks. Count matrix is there just so that we can create seurat object.
# Actual motif analysis is not effected by cells at this point as both groups of peaks are already determined (background and module-associated)
counts1 = counts1[1:nrow(annot_peaks),]
rownames(counts1) = annot_peaks$peak

# Create seurat objects
s1_excinh <- CreateSeuratObject(counts = counts1, assay = 'peaks', project = 'ATAC', min.cells = 1)


### Scan the DNA sequence of each peak for the presence of each motif ###

motmat1 <- CreateMotifMatrix(features = StringToGRanges(rownames(s1_excinh), sep = c(":", "-")),
				pwm = cisbp_hum_ExcInh,
				genome = 'hg38',
				sep = c(":", "-"))
saveRDS(motmat1, "~/workdir/STEFANO_SME/allsamples/SME_ExcInh/motmat.RDS")


## Motif Enrichment Analysis ##
tmp = counts1[rownames(motmat1),]
s1_excinh <- CreateSeuratObject(counts = tmp, assay = 'peaks', project = 'ATAC', min.cells = 1)

# Create a new Motif object to store the results
motif1 <- CreateMotifObject(data = motmat1, pwm = cisbp_hum_ExcInh)

# Add the Motif object to the assay
s1_excinh <- AddMotifObject(object = s1_excinh, motif.object = motif1)
s1_excinh <- RegionStats(object = s1_excinh, genome = BSgenome.Hsapiens.UCSC.hg38, sep = c(":", "-"))

# WM4 Motif Enrichment
gene_assoc_peaks = read.table("~/workdir/STEFANO_SME/allsamples/peaks_ExcInh/WM4_specific_peaks_tss_upstr_ExcInh.bed", header = T)
peaks = as.character(gene_assoc_peaks$peak)
peaks = strsplit(peaks, "_")
peaks = lapply(peaks, function(x){paste(paste(x[1], x[2], sep = ":"), x[3], sep = "-")})
peaks = unlist(peaks)

wm4_enriched_motifs_s1 <- FindMotifs(object = s1_excinh, background = rownames(s1_excinh), features = peaks)
wm4_enriched_motifs_s1$p_adjusted = p.adjust(wm4_enriched_motifs_s1$pvalue, method = "BH")
signif = wm4_enriched_motifs_s1[wm4_enriched_motifs_s1$p_adjusted < 0.05,]

rio::export(list(WM4_All_ = wm4_enriched_motifs_s1,
		 WM4_Filtered = signif),
	   file = paste0("~/workdir/STEFANO_SME/allsamples/sample", i, "_peaks_ExcInh/AllPeaksUPS_WM4_ExcInh_EnrichedMotifs.xlsx"))


# WM12 Motif Enrichment
gene_assoc_peaks = read.table("~/workdir/STEFANO_SME/allsamples/peaks_ExcInh/WM12_specific_peaks_tss_upstr_ExcInh.bed", header = T)
peaks = as.character(gene_assoc_peaks$peak)
peaks = strsplit(peaks, "_")
peaks = lapply(peaks, function(x){paste(paste(x[1], x[2], sep = ":"), x[3], sep = "-")})
peaks = unlist(peaks)

WM12_enriched_motifs_s1 <- FindMotifs(object = s1_excinh, background = rownames(s1_excinh), features = peaks)
WM12_enriched_motifs_s1$p_adjusted = p.adjust(WM12_enriched_motifs_s1$pvalue, method = "BH")
signif = WM12_enriched_motifs_s1[WM12_enriched_motifs_s1$p_adjusted < 0.05,]

rio::export(list(WM12_All_ = WM12_enriched_motifs_s1,
		 WM12_Filtered = signif),
	   file = paste0("~/workdir/STEFANO_SME/allsamples/sample", i, "_peaks_ExcInh/AllPeaksUPS_WM12_ExcInh_EnrichedMotifs.xlsx"))


###############################
### WM21 & Motif Enrichment ###
###############################

counts1 <- readMM("~/workdir/STEFANO_SME/allsamples/peaks_OligOPC/sample1_pc_matrix/pc_mat.mtx")
counts1 = t(counts1)
pks1 <- read.table("~/workdir/STEFANO_SME/allsamples/peaks_OligOPC/sample1_pc_matrix/peaks.txt")
brc1 <- read.table("~/workdir/STEFANO_SME/allsamples/peaks_OligOPC/sample1_pc_matrix/barcodes.tsv")

colnames(counts1) = brc1$V1
rownames(counts1) = pks1$V1


# Only keep the upstream peaks per sample

annot_peaks = read.table("~/workdir/STEFANO_SME/allsamples/SME_OligOPC/peak_annotation.tsv", sep = "\t", header = T, stringsAsFactors = F)

tmp = annot_peaks$distance
tmplist = strsplit(tmp, ";")
tmpbool = lapply(tmplist, function(x){sum(as.numeric(x) >= 0) == 0})
tmpbool = unlist(tmpbool)

annot_peaks = annot_peaks[tmpbool,]
tmp = sub("_", ":", annot_peaks$peak)
tmp2 = sub("_", "-", tmp)
annot_peaks$peak = tmp2

# Subset same number of peaks from counts matrix.
# We will switch the rownames of a random count matrix with these peaks. Count matrix is there just so that we can create seurat object.
# Actual motif analysis is not effected by cells at this point as both groups of peaks are already determined (background and module-associated)
counts1 = counts1[1:nrow(annot_peaks),]
rownames(counts1) = annot_peaks$peak

# Create seurat objects
s1_OligOPC <- CreateSeuratObject(counts = counts1, assay = 'peaks', project = 'ATAC', min.cells = 1)


### Scan the DNA sequence of each peak for the presence of each motif ###

motmat1 <- CreateMotifMatrix(features = StringToGRanges(rownames(s1_OligOPC), sep = c(":", "-")),
				pwm = cisbp_hum_OligOPC,
				genome = 'hg38',
				sep = c(":", "-"))
saveRDS(motmat1, "~/workdir/STEFANO_SME/allsamples/SME_OligOPC/motmat1.RDS")


## Motif Enrichment Analysis ##
tmp = counts1[rownames(motmat1),]
s1_OligOPC <- CreateSeuratObject(counts = tmp, assay = 'peaks', project = 'ATAC', min.cells = 1)

# Create a new Motif object to store the results
motif1 <- CreateMotifObject(data = motmat1, pwm = cisbp_hum_OligOPC)

# Add the Motif object to the assay
s1_OligOPC <- AddMotifObject(object = s1_OligOPC, motif.object = motif1)
s1_OligOPC <- RegionStats(object = s1_OligOPC, genome = BSgenome.Hsapiens.UCSC.hg38, sep = c(":", "-"))

# WM21 Motif Enrichment
gene_assoc_peaks = read.table("~/workdir/STEFANO_SME/allsamples/peaks_OligOPC/WM21_specific_peaks_tss_upstr_OligOPC.bed", header = T)
peaks = as.character(gene_assoc_peaks$peak)
peaks = strsplit(peaks, "_")
peaks = lapply(peaks, function(x){paste(paste(x[1], x[2], sep = ":"), x[3], sep = "-")})
peaks = unlist(peaks)

wm21_enriched_motifs_s1 <- FindMotifs(object = s1_OligOPC, background = rownames(s1_OligOPC), features = peaks)
wm21_enriched_motifs_s1$p_adjusted = p.adjust(wm21_enriched_motifs_s1$pvalue, method = "BH")
signif = wm21_enriched_motifs_s1[wm21_enriched_motifs_s1$p_adjusted < 0.05,]

rio::export(list(wm21_All_ = wm21_enriched_motifs_s1,
		 wm21_Filtered = signif),
	   file = paste0("~/workdir/STEFANO_SME/allsamples/sample", i, "_peaks_OligOPC/AllPeaksUPS_wm21_OligOPC_EnrichedMotifs.xlsx"))



## Plot all results ##
data = import_list("~/workdir/STEFANO_SME/allsamples/Combinedpeaks_motif_analysis.xlsx")
data = data[c(1,3,5)]
data[[1]]$module = "WM4"
data[[2]]$module = "WM12"
data[[3]]$module = "WM21"

data = lapply(data, function(x){x$tolabel = 0; x[1:5,]$tolabel = 1; x})

data = lapply(data, function(x){x$tolabel = ifelse(x$p_adjusted < 0.05 & x$tolabel == 1, 1, 0) ; x})


data_df = do.call(rbind, data)
data_df$TFname = do.call(rbind, strsplit(data_df$motif, "_"))[,2]
data_df$log_p_adjusted = -log10(data_df$p_adjusted)
data_df$fold.enrichment = log2(data_df$fold.enrichment + 1)

# Also label SMAD3 and GLIS1
data_df[(data_df$TFname == 'SMAD3' | data_df$TFname == 'GLIS1') &
	(data_df$module == 'WM12' & data_df$p_adjusted < 0.05), 'tolabel'] = 1

data_df$module = factor(data_df$module, levels = c('WM4', 'WM12', 'WM21'))

pdf("motif_enrichment.pdf")
ggplot(data = data_df, aes(,x=log_p_adjusted, y=fold.enrichment, color = module, label = TFname)) +
	geom_point(alpha=0.5) +
	geom_text_repel(data = data_df[data_df$tolabel == 1,], hjust = 1.5, size = 3, fontface = 2, color = 'black') +
	geom_vline(xintercept = 1.3, linetype="dashed", color = "red") +
	theme_classic() +
	xlab(expression(paste(-log[10], FDR))) +
	ylab(expression('Fold Enrichment log'[2])) +
	facet_grid(module ~ .) +
    	theme(axis.text.x = element_text(face = "bold", size = 13),
  	  axis.text.y = element_text(face = "bold", size = 13),
          axis.title=element_text(size=14))
dev.off()


### For IGV ###

# The functions
bedr_fun = function(gene_reg, motif_reg){

	result_int <- bedr(
        input = list(a = gene_reg, b = motif_reg), 
        method = "intersect",
	params = c("-loj", "-header")
        )
	return(result_int)
}

write_bed = function(bedr_out, outname){

	write.table(bedr_out[, c(2,3,4)],
			file = outname,
			sep = "\t",
			quote = F,
			row.names = F,
			col.names = F)

}

# Load the data
motmat = readRDS("~/workdir/STEFANO_SME/allsamples/SME_ExcInh/motmat1.RDS")

### SMAD3 ###
# Filter for the TF
smad3_vec = motmat[, "M09380_SMAD3"]
smad3_vec = smad3_vec[smad3_vec == 1]

# Define 20kb upstream region
#plxna4_20kb_ups = "chr7:132576564-132596564"
cacna_20kb_ups = "chr12:2033563-2053563"
trio_20kb_ups = "chr5:14123702-14143702"
slc7a6_20kb_ups = "chr16:68244516-68264516"
shank2_20kb_ups = "chr11:70661818-70682197"
nedd4l_20kb_ups = "chr18:58024548-58044548"

#inters_pxna4 = bedr_fun(plxna4_20kb_ups, names(smad3_vec))
inters_cacna = bedr_fun(cacna_20kb_ups, names(smad3_vec))
inters_trio = bedr_fun(trio_20kb_ups, names(smad3_vec))
inters_slc7a6 = bedr_fun(slc7a6_20kb_ups, names(smad3_vec))
inters_shank2 = bedr_fun(shank2_20kb_ups, names(smad3_vec))
inters_nedd4l = bedr_fun(nedd4l_20kb_ups, names(smad3_vec))

dir = "~/workdir/STEFANO_SME/allsamples/forIGV/"
write_bed(inters_cacna, paste0(dir, "SMAD3_CACNA_20kb_ups_hits.bed"))
write_bed(inters_trio, paste0(dir, "SMAD3_TRIO_20kb_ups_hits.bed"))
write_bed(inters_slc7a6, paste0(dir, "SMAD3_SLC7A6_20kb_ups_hits.bed"))
write_bed(inters_shank2, paste0(dir, "SMAD3_SHANK2_20kb_ups_hits.bed"))
write_bed(inters_nedd4l, paste0(dir, "SMAD3_NEDD4L_20kb_ups_hits.bed"))


### GLIS1 ###
# Filter for the TF
glis1_vec = motmat[, "M08345_GLIS1"]
glis1_vec = glis1_vec[glis1_vec == 1]

# Define 20kb upstream region
#plxna4_20kb_ups = "chr7:132576564-132596564"
cacna_20kb_ups = "chr12:2033563-2053563"
trio_20kb_ups = "chr5:14123702-14143702"
slc7a6_20kb_ups = "chr16:68244516-68264516"
shank2_20kb_ups = "chr11:70661818-70682197"
nedd4l_20kb_ups = "chr18:58024548-58044548"

#inters_pxna4 = bedr_fun(plxna4_20kb_ups, names(smad3_vec))
inters_cacna = bedr_fun(cacna_20kb_ups, names(glis1_vec))
inters_trio = bedr_fun(trio_20kb_ups, names(glis1_vec))
inters_slc7a6 = bedr_fun(slc7a6_20kb_ups, names(glis1_vec))
inters_shank2 = bedr_fun(shank2_20kb_ups, names(glis1_vec))
inters_nedd4l = bedr_fun(nedd4l_20kb_ups, names(glis1_vec))

dir = "~/workdir/STEFANO_SME/allsamples/forIGV/"
write_bed(inters_cacna, paste0(dir, "GLIS1_CACNA_20kb_ups_hits.bed"))
write_bed(inters_trio, paste0(dir, "GLIS1_TRIO_20kb_ups_hits.bed"))
write_bed(inters_slc7a6, paste0(dir, "GLIS1_SLC7A6_20kb_ups_hits.bed"))
write_bed(inters_shank2, paste0(dir, "GLIS1_SHANK2_20kb_ups_hits.bed"))
write_bed(inters_nedd4l, paste0(dir, "GLIS1_NEDD4L_20kb_ups_hits.bed"))


