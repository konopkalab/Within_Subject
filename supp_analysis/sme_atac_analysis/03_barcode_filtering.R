library(gdata)


# Load all fragments
dir = "~/workdir/pr3/atacseq/human/LEGA_SURGICAL/"
allbarcs = readRDS(paste0(dir, "fragments.RDS"))

# Filter for neuronal fragments (WM4 and WM12 enriched)
sub_barcs = read.table("~/workdir/STEFANO_SME/ExcInh_cells.txt", header = F, stringsAsFactors = F)
sub_fragm = allbarcs[allbarcs$V4 %in% sub_barcs$V1,]

saveRDS(sub_fragm, file = "~/workdir/STEFANO_SME/ExcInh_cells_fragments.RDS")
write.table(sub_fragm, file = "~/workdir/STEFANO_SME/ExcInh_cells_fragments.txt", row.names = F, quote = F, col.names = F, sep = "\t")


# Filter for oligodendrocyte and OPC fragments (WM21 enriched)
sub_barcs = read.table("~/workdir/STEFANO_SME/OligOPC_cells_cells.txt", header = F, stringsAsFactors = F)
sub_fragm = allbarcs[allbarcs$V4 %in% sub_barcs$V1,]

saveRDS(sub_fragm, file = "~/workdir/STEFANO_SME/OligOPC_cells_fragments.RDS")
write.table(sub_fragm, file = "~/workdir/STEFANO_SME/OligOPC_cells_fragments.txt", row.names = F, quote = F, col.names = F, sep = "\t")


# Adjust neuronal fragments
fragm = readRDS("~/workdir/STEFANO_SME/ExcInh_cells_fragments.RDS")

fragm1 = fragm
fragm2 = fragm

posmat1 = as.matrix(fragm1[, 2:3])
for(i in 1:nrow(posmat1)){

	posmat1[i,2] = posmat1[i,1] + 100
	posmat1[i,1] = posmat1[i,1] - 100
}

fragm1[, 2:3] = posmat1


posmat2 = as.matrix(fragm2[, 2:3])
for(i in 1:nrow(posmat2)){

	posmat2[i,1] = posmat2[i,2] - 100
	posmat2[i,2] = posmat2[i,2] + 100
}

fragm2[, 2:3] = posmat2

smooth_sum_200 = interleave(fragm1, fragm2)
options(scipen = 999)

write.table(smooth_sum_200, file = "~/workdir/STEFANO_SME/ExcInh_cells_adj_fragments.txt",
		sep = "\t", col.names = F, row.names = F, quote = F)


# Adjust oligodendrocyte and OPC fragments
fragm = readRDS("~/workdir/STEFANO_SME/OligOPC_cells_fragments.RDS")

fragm1 = fragm
fragm2 = fragm

posmat1 = as.matrix(fragm1[, 2:3])
for(i in 1:nrow(posmat1)){

	posmat1[i,2] = posmat1[i,1] + 100
	posmat1[i,1] = posmat1[i,1] - 100
}

fragm1[, 2:3] = posmat1


posmat2 = as.matrix(fragm2[, 2:3])
for(i in 1:nrow(posmat2)){

	posmat2[i,1] = posmat2[i,2] - 100
	posmat2[i,2] = posmat2[i,2] + 100
}

fragm2[, 2:3] = posmat2

smooth_sum_200 = interleave(fragm1, fragm2)
options(scipen = 999)

write.table(smooth_sum_200, file = "~/workdir/STEFANO_SME/OligOPC_cells_adj_fragments.txt",
		sep = "\t", col.names = F, row.names = F, quote = F)

