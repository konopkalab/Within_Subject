library(rio)

## WM12 and 4 ##

# Loop through each sample
for(i in 1:3){
	peaks = read.table(paste0("~/workdir/STEFANO_SME/allsamples/cellrangerout_ExcInh_sample", i, "/peak_annotation.tsv"),
				sep = "\t", header = T, stringsAsFactors = F)
	mod_genes = import("~/workdir/STEFANO_SME/allsamples/Mods_ATAC_0.3.xlsx")
	wm4_genes = mod_genes[mod_genes$Class == "WM4", "Gene"]
	wm12_genes = mod_genes[mod_genes$Class == "WM12", "Gene"]

	# Get only upstream of the peaks
	peaks_wm4 = peaks[(peaks$gene %in% wm4_genes) & peaks$distance < 0, ]
	peaks_wm12 = peaks[(peaks$gene %in% wm12_genes) & peaks$distance < 0, ]

	write.table(peaks_wm4, paste0("~/workdir/STEFANO_SME/allsamples/sample", i, "_peaks_ExcInh/WM4_specific_peaks_tss_upstr_ExcInh.bed"),
			quote = F, row.names = F, sep = "\t")

	write.table(peaks_wm12, paste0("~/workdir/STEFANO_SME/allsamples/sample", i, "_peaks_ExcInh/WM12_specific_peaks_tss_upstr_ExcInh.bed"),
			quote = F, row.names = F, sep = "\t")
}


## WM21 ##
# Loop through each sample
for(i in 1:3){
	peaks = read.table(paste0("~/workdir/STEFANO_SME/allsamples/cellrangerout_OligOPC_sample", i, "/peak_annotation.tsv"),
				sep = "\t", header = T, stringsAsFactors = F)
	mod_genes = import("~/workdir/STEFANO_SME/allsamples/Mods_ATAC_0.3.xlsx")
	wm21_genes = mod_genes[mod_genes$Class == "WM21", "Gene"]

	# Get only upstream of the peaks
	peaks_wm21 = peaks[(peaks$gene %in% wm21_genes) & peaks$distance < 0, ]

	write.table(peaks_wm21, paste0("~/workdir/STEFANO_SME/allsamples/sample", i, "_peaks_OligOPC/WM21_specific_peaks_tss_upstr_OligOPC.bed"),
			quote = F, row.names = F, sep = "\t")

}


