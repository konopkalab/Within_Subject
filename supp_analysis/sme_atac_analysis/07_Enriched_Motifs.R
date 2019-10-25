library(rio)
library(Signac)
library(Seurat)
library(JASPAR2018)
library(TFBSTools)
library("BSgenome.Hsapiens.UCSC.hg38")
library(universalmotif)
library(bedr)

#################################################################
### Filter TF Motifs based on TF's presence in the given data ###
#################################################################

# Get a list of motif position weight matrices from the CISBP database (in jaspar format)
cisbp_hum = readJASPARMatrix("~/konopkaLab/pr3/signac_tmp/sme_atac/motif.jaspar")
cisbp_hum = toPWM(cisbp_hum, type = "prob")

#Take only the TFs in our clusters
TFs = rio::import("~/konopkaLab/pr3/signac_tmp/sme_atac/TFs_Annot.xlsx")
ExcInh_TFs = TFs[TFs$Class == "Exc_Inh", "Gene"]
OligOPC_TFs = TFs[TFs$Class == "Olig_OPC", "Gene"]


#Take only motifs with direct or best prediction TFs.
tmp_list = strsplit(names(cisbp_hum), "_")
theTFs = lapply(tmp_list, function(x){x[2]})
theTFs = unlist(theTFs)

cisbp_hum_ExcInh = cisbp_hum[theTFs %in% ExcInh_TFs]
cisbp_hum_OligOPC = cisbp_hum[theTFs %in% OligOPC_TFs]



#####################################
### WM4 & WM12 - Motif Enrichment ###
#####################################

counts <- Read10X_h5("~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/filtered_peak_bc_matrix.h5")

# Only keep the upstream peaks
annot_peaks = read.table("~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/peak_annotation.tsv",
			 sep = "\t", header = T, stringsAsFactors = F)

tmp = annot_peaks$distance
tmplist = strsplit(tmp, ";")
tmpbool = lapply(tmplist, function(x){sum(as.numeric(x) >= 0) == 0})
tmpbool = unlist(tmpbool)

annot_peaks = annot_peaks[tmpbool,]
tmp = sub("_", ":", annot_peaks$peak)
tmp2 = sub("_", "-", tmp)
annot_peaks$peak = tmp2

counts = counts[rownames(counts) %in% annot_peaks$peak,]

metadata <- read.csv(
  file = "~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/singlecell.csv",
  header = TRUE,
  row.names = 1
)

lega_ExcInh <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

fragment.path <- "~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/fragments.tsv.gz"

lega_ExcInh <- SetFragments(
  object = lega_ExcInh,
  file = fragment.path
)


### Scan the DNA sequence of each peak for the presence of each motif ###

motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(lega_ExcInh), sep = c(":", "-")),
  pwm = cisbp_hum_ExcInh,
  genome = 'hg38',
  sep = c(":", "-")
)

saveRDS(motif.matrix, "~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/UPS_PWM_ExcInh_TF_MotifMatrix_Binarized.RDS")

# Create a new Motif object to store the results.
class(cisbp_hum_ExcInh)[1] = "AnyPWM"
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = cisbp_hum_ExcInh
)


# Add the Motif object to the assay
lega_ExcInh[['peaks']] <- AddMotifObject(
  object = lega_ExcInh[['peaks']],
  motif.object = motif
)


lega_ExcInh <- RegionStats(
  object = lega_ExcInh,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)


#WM4 Motif Enrichment
gene_assoc_peaks = read.table("~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/WM4_specific_peaks_upstr_ExcInh.bed", header = T)
peaks = as.character(gene_assoc_peaks$peak)
peaks = strsplit(peaks, "_")
peaks = lapply(peaks, function(x){paste(paste(x[1], x[2], sep = ":"), x[3], sep = "-")})
peaks = unlist(peaks)

wm4_enriched.motifs <- FindMotifs(
  object = lega_ExcInh,
  background = rownames(lega_ExcInh),
  features = peaks
)

wm4_enriched.motifs$p_adjusted = p.adjust(wm4_enriched.motifs$pvalue, method = "BH")
wm4_filter_enriched_motif = wm4_enriched.motifs[wm4_enriched.motifs$p_adjusted < 0.05 &
					wm4_enriched.motifs$fold.enrichment > 1.3,]

rio::export(list(WM4_All_ = wm4_enriched.motifs,
		 WM4_Filtered = wm4_filter_enriched_motif),
	   file = "~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/AllPeaksUPS_WM4_upstr_ExcInh_EnrichedMotifs.xlsx")




#WM12 Motif Enrichment
gene_assoc_peaks = read.table("~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/WM12_specific_peaks_upstr_ExcInh.bed", header = T)
peaks = as.character(gene_assoc_peaks$peak)
peaks = strsplit(peaks, "_")
peaks = lapply(peaks, function(x){paste(paste(x[1], x[2], sep = ":"), x[3], sep = "-")})
peaks = unlist(peaks)

wm12_enriched.motifs <- FindMotifs(
  object = lega_ExcInh,
  background = rownames(lega_ExcInh),
  features = peaks
)


wm12_enriched.motifs$p_adjusted = p.adjust(wm12_enriched.motifs$pvalue, method = "BH")
wm12_filter_enriched_motif = wm12_enriched.motifs[wm12_enriched.motifs$p_adjusted < 0.05 &
					     wm12_enriched.motifs$fold.enrichment > 1.3,]


rio::export(list(WM12_All_ = wm12_enriched.motifs,
		 WM12_Filtered = wm12_filter_enriched_motif),
	   file = "~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/AllPeaksUPS_WM12_upstr_ExcInh_EnrichedMotifs.xlsx")



###############################
### WM21 - Motif Enrichment ###
###############################

### Read and filter data - only Oligo and OPC clusters
counts <- Read10X_h5("~/konopkaLab/pr3/signac_tmp/sme_atac/OligOPC/filtered_peak_bc_matrix.h5")


# Get rid of peaks within any gene
annot_peaks = read.table("~/konopkaLab/pr3/signac_tmp/sme_atac/OligOPC/peak_annotation.tsv",
			 sep = "\t", header = T, stringsAsFactors = F)

tmp = annot_peaks$distance
tmplist = strsplit(tmp, ";")
tmpbool = lapply(tmplist, function(x){sum(as.numeric(x) >= 0) == 0})
tmpbool = unlist(tmpbool)

annot_peaks = annot_peaks[tmpbool,]
tmp = sub("_", ":", annot_peaks$peak)
tmp2 = sub("_", "-", tmp)
annot_peaks$peak = tmp2

counts = counts[rownames(counts) %in% annot_peaks$peak,]


metadata <- read.csv(
  file = "~/konopkaLab/pr3/signac_tmp/sme_atac/OligOPC/singlecell.csv",
  header = TRUE,
  row.names = 1
)

lega_OligOPC <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

fragment.path <- "~/konopkaLab/pr3/signac_tmp/sme_atac/OligOPC/fragments.tsv.gz"

lega_OligOPC <- SetFragments(
  object = lega_OligOPC,
  file = fragment.path
)


# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(lega_OligOPC), sep = c(":", "-")),
  pwm = cisbp_hum_OligOPC,
  genome = 'hg38',
  sep = c(":", "-")
)

saveRDS(motif.matrix, "~/konopkaLab/pr3/signac_tmp/sme_atac/OligOPC/UPS_PWM_OligOPC_TF_MotifMatrix_Binarized.RDS")

# Create a new Mofif object to store the results
class(cisbp_hum_OligOPC)[1] = "AnyPWM"
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = cisbp_hum_OligOPC
)

# Add the Motif object to the assay
lega_OligOPC[['peaks']] <- AddMotifObject(
  object = lega_OligOPC[['peaks']],
  motif.object = motif
)


lega_OligOPC <- RegionStats(
  object = lega_OligOPC,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

#WM21 Motif Enrichment
gene_assoc_peaks = read.table("~/konopkaLab/pr3/signac_tmp/sme_atac/OligOPC/WM21_specific_peaks_OligOPC.bed", header = T)
gene_assoc_peaks = gene_assoc_peaks[gene_assoc_peaks$distance < 0,]
peaks = as.character(gene_assoc_peaks$peak)
peaks = strsplit(peaks, "_")
peaks = lapply(peaks, function(x){paste(paste(x[1], x[2], sep = ":"), x[3], sep = "-")})
peaks = unlist(peaks)

wm21_enriched.motifs <- FindMotifs(
  object = lega_OligOPC,
  background = rownames(lega_OligOPC),
  features = peaks
)


wm21_enriched.motifs$p_adjusted = p.adjust(wm21_enriched.motifs$pvalue, method = "BH")
wm21_filter_enriched_motif = wm21_enriched.motifs[wm21_enriched.motifs$p_adjusted < 0.05 &
					     wm21_enriched.motifs$fold.enrichment > 1.3,]




rio::export(list(WM4_All = wm4_enriched.motifs,
		 WM4_Filtered = wm4_filter_enriched_motif,
		 WM12_All = wm12_enriched.motifs,
		 WM12_Filtered = wm12_filter_enriched_motif,
		 WM21_All = wm21_enriched.motifs,
		 WM21_Filtered = wm21_filter_enriched_motif),
	   file = "~/konopkaLab/pr3/signac_tmp/sme_atac/AllPeaksUPS_WM4-12-21_upstr_EnrichedMotifs.xlsx")




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
motmat = readRDS("~/konopkaLab/pr3/signac_tmp/sme_atac/ExcInh/PWM_ExcInh_TF_MotifMatrix_Binarized.RDS")

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

dir = "~/konopkaLab/pr3/signac_tmp/sme_atac/bed_for_figure/"
write_bed(inters_cacna, paste0(dir, "SMAD3_CACNA_20kb_ups_hits.bed"))
write_bed(inters_trio, paste0(dir, "SMAD3_TRIO_20kb_ups_hits.bed"))
write_bed(inters_slc7a6, paste0(dir, "SMAD3_SLC7A6_20kb_ups_hits.bed"))
write_bed(inters_shank2, paste0(dir, "SMAD3_SHANK2_20kb_ups_hits.bed"))
write_bed(inters_nedd4l, paste(dir, "SMAD3_NEDD4L_20kb_ups_hits.bed"))

### MEF2A ###
# Filter for the TF
mef2a_vec = motmat[, "M08213_MEF2A"]
mef2a_vec = mef2a_vec[mef2a_vec == 1]

# Define 20kb upstream region
il1rapl2_20kb_ups = "chrX:104546315-104566315"
cntnap5_20kb_ups = "chr2:124005287-124025287" #X forw
frmpd4_20kb_ups = "chrX:12118466-12138466"
ccdc85a_20kb_ups = "chr2:56164123-56184123"
pak7_20kb_ups = "chr20:9839041-9859041"
rspo2_20kb_ups ="chr8:108083648-108103648"
lnx1_20kb_ups = "chr4:53591616-53611616" #X rev
kcnq5_20kb_ups = "chr6:72602081-72622081"
dgki_20kb_ups = "chr7:137846536-137866536" #X rev

inters_il1rapl2 = bedr_fun(il1rapl2_20kb_ups, names(mef2a_vec))
inters_cntnap5 = bedr_fun(cntnap5_20kb_ups, names(mef2a_vec))
inters_frmpd4 = bedr_fun(frmpd4_20kb_ups, names(mef2a_vec))
inters_ccdc85a = bedr_fun(ccdc85a_20kb_ups, names(mef2a_vec))
inters_pak7 = bedr_fun(pak7_20kb_ups, names(mef2a_vec))
inters_rspo2 = bedr_fun(rspo2_20kb_ups, names(mef2a_vec))
inters_lnx1 = bedr_fun(lnx1_20kb_ups, names(mef2a_vec))
inters_kcnq5 = bedr_fun(kcnq5_20kb_ups, names(mef2a_vec))
inters_dgki = bedr_fun(dgki_20kb_ups, names(mef2a_vec))

write_bed(inters_cntnap5, "~/konopkaLab/pr3/signac_tmp/sme_atac/Figure/bed_for_figure/MEF2A_CNTNAP5_20kb_ups_hits.bed")
write_bed(inters_lnx1, "~/konopkaLab/pr3/signac_tmp/sme_atac/Figure/bed_for_figure/MEF2A_LNX1_20kb_ups_hits.bed")
write_bed(inters_dgki, "~/konopkaLab/pr3/signac_tmp/sme_atac/Figure/bed_for_figure/MEF2A_DGKI_20kb_ups_hits.bed")

### MEF2C ###
# Filter for the TF
mef2c_vec = motmat[, "M08149_MEF2C"]
mef2c_vec = mef2c_vec[mef2c_vec == 1]

# Define 20kb upstream region
il1rapl2_20kb_ups = "chrX:104546315-104566315"
cntnap5_20kb_ups = "chr2:124005287-124025287"

inters_il1rapl2 = bedr_fun(il1rapl2_20kb_ups, names(mef2c_vec))
inters_cntnap5 = bedr_fun(cntnap5_20kb_ups, names(mef2c_vec))

write_bed(inters_il1rapl2, "~/konopkaLab/pr3/signac_tmp/sme_atac/bed_for_figure/MEF2C_IL1RAPL2_20kb_ups_hits.bed")

write_bed(inters_cntnap5, "~/konopkaLab/pr3/signac_tmp/sme_atac/bed_for_figure/MEF2C_CNTNAP5_20kb_ups_hits.bed")



# Load the data
motmat = readRDS("~/konopkaLab/pr3/signac_tmp/sme_atac/OligOPC/PWM_OligOPC_TF_MotifMatrix_Binarized.RDS")

### SOX5 ###
# Filter for the TF
sox5_vec = motmat[, "M05894_SOX5"]
sox5_vec = sox5_vec[sox5_vec == 1]

# Define 20kb upstream region
shroom4_20kb_ups = "chrX:50591647-50611647"
cntnna3_20kb_ups = "chr10:67696169-67716169"

inters_shroom4 = bedr_fun(shroom4_20kb_ups, names(sox5_vec))
inters_cntnna3 = bedr_fun(cntnna3_20kb_ups, names(sox5_vec))


### SOX6 ###
# Filter for the TF
sox6_vec = motmat[, "M05892_SOX6"]
sox6_vec = sox6_vec[sox6_vec == 1]

# Define 20kb upstream region
shroom4_20kb_ups = "chrX:50591647-50611647"
cntnna3_20kb_ups = "chr10:67696169-67716169"

inters_shroom4 = bedr_fun(shroom4_20kb_ups, names(sox6_vec))
inters_cntnna3 = bedr_fun(cntnna3_20kb_ups, names(sox6_vec))

write_bed(inters_cacna, "SMAD3_CACNA_20kb_ups_hits.bed")
write_bed(inters_shank2, "SMAD3_SHANK2_20kb_ups_hits.bed")





