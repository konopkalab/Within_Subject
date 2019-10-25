library(Seurat)
library(rio)

mods = import("~/workdir/STEFANO_SME/table_stefano/Table_S3.xlsx")
m4 = mods[mods$ModuleName == "WM4",]
m12 = mods[mods$ModuleName == "WM12",]
m21 = mods[mods$ModuleName == "WM21",]

scATACseq = readRDS("~/workdir/STEFANO_SME/seur_objects/sme_atac_microglia_removed_annotated_seur.RDS")
DefaultAssay(scATACseq) = "RNA"

### Filter module genes for presence in their own clusters in scATACseq ###

#Module 4 and Module 12
#Filter for presence in atac-seq (>30%)
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


#Write cells in wanted clusters. This is to call the peaks again.
to_peak_call = WhichCells(scATACseq, idents = names(table(Idents(ExcInh_atac))))
write.table(to_peak_call, "~/workdir/STEFANO_SME/ExcInh_cells.txt", sep = "\n", quote = F,
	    col.names = F, row.names = F)


#Module 21
#Filter for presence in atac-seq (>30%)
cl_names = names(table(Idents(scATACseq)))
ids = cl_names[grepl("Oligo|OPC", cl_names)]
OligOPC_atac = subset(scATACseq, idents = ids)
m21_ov = m21$Gene[m21$Gene %in% rownames(OligOPC_atac)]
m21_atac = OligOPC_atac@assays$RNA@counts[m21_ov,]
m21_zero = apply(as.matrix(m21_atac), 1, function(x){sum(x == 0)})
m21_zero_filt = m21_zero[((ncol(m21_atac) - m21_zero) / ncol(m21_atac)) > 0.3]

genes = c(names(m4_zero_filt), names(m12_zero_filt), names(m21_zero_filt))
class = c(rep("WM4", length(m4_zero_filt)), rep("WM12", length(m12_zero_filt)), rep("WM21", length(m21_zero_filt)))
towrite_df = data.frame(Gene = genes, Class = class)

export(towrite_df, "~/workdir/STEFANO_SME/Mods_ATAC_0.3.xlsx")


#Write cells in wanted clusters. This is to call the peaks again.
to_peak_call = WhichCells(scATACseq, idents = names(table(Idents(OligOPC_atac))))
write.table(to_peak_call, "~/workdir/STEFANO_SME/OligOPC_cells.txt", sep = "\n", quote = F,
	    col.names = F, row.names = F)


### Find TFs with motifs and sufficiently expressed in ExcInh and/or OligOPC clusters ###

#Use all TFs from CISBP for now. Later we will filter for only confident TF motifs.
motif_dat = read.table("~/workdir/motif_databases/CISBP_HUMAN_2019/TF_Information_all_motifs_plus.txt",
		header = T, sep = "\t", na.strings = "NA", fill = T)
all_TF_names = unique(motif_dat$TF_Name)

#Filter for presence in atac-seq. M4 and M12
TF_ov = all_TF_names[all_TF_names %in% rownames(ExcInh_atac@assays$RNA@counts)]
TF_atac = ExcInh_atac@assays$RNA@counts[TF_ov,]
zero_TF = apply(as.matrix(TF_atac), 1, function(x){sum(x == 0)})
m4_12_zero_TF_filt = zero_TF[((ncol(TF_atac) - zero_TF) / ncol(TF_atac)) > 0.3]

#Filter for presence in atac-seq. M21
TF_ov = all_TF_names[all_TF_names %in% rownames(OligOPC_atac@assays$RNA@counts)]
TF_atac = OligOPC_atac@assays$RNA@counts[TF_ov,]
zero_TF = apply(as.matrix(TF_atac), 1, function(x){sum(x == 0)})
m21_zero_TF_filt = zero_TF[((ncol(TF_atac) - zero_TF) / ncol(TF_atac)) > 0.3]

genes = c(names(m4_12_zero_TF_filt), names(m21_zero_TF_filt))
class = c(rep("Exc_Inh", length(m4_12_zero_TF_filt)), rep("Olig_OPC", length(m21_zero_TF_filt)))
towrite_df = data.frame(Gene = genes, Class = class)

export(towrite_df, "~/workdir/STEFANO_SME/TFs_Annot.xlsx")


