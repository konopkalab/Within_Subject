library(rio)

## WM12 and 4 ##
peaks = read.table("~/workdir/STEFANO_SME/SME_EXCINH/outs/peak_annotation.tsv", sep = "\t", header = T,
			stringsAsFactors = F)
mod_genes = import("~/workdir/STEFANO_SME/Mods_ATAC_0.3.xlsx")
wm4_genes = mod_genes[mod_genes$Class == "WM4", "Gene"]
wm12_genes = mod_genes[mod_genes$Class == "WM12", "Gene"]

# Get only upstream of the peaks
peaks_wm4 = peaks[(peaks$gene %in% wm4_genes) & peaks$distance < 0, ]
peaks_wm12 = peaks[(peaks$gene %in% wm12_genes) & peaks$distance < 0, ]

write.table(peaks_wm4, file = "~/workdir/STEFANO_SME/WM4_specific_peaks_tss_upstr_ExcInh.bed",
		quote = F, row.names = F, sep = "\t")

write.table(peaks_wm12, file = "~/workdir/STEFANO_SME/WM12_specific_peaks_tss_upstr_ExcInh.bed",
		quote = F, row.names = F, sep = "\t")



## WM21 ##
peaks = read.table("~/workdir/STEFANO_SME/SME_OligOPC/outs/peak_annotation.tsv", sep = "\t", header = T,
			stringsAsFactors = F)

wm21_genes = mod_genes[mod_genes$Class == "WM21", "Gene"]

# Get only upstream of the peaks
peaks_wm21 = peaks[(peaks$gene %in% wm21_genes) & peaks$distance < 0, ]

write.table(peaks_wm21, file = "~/workdir/STEFANO_SME/WM21_specific_peaks_tss_upstr_OligOPC.bed",
		quote = F, row.names = F, sep = "\t")
