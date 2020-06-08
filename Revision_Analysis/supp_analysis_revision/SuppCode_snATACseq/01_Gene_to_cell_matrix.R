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
set.seed(1234)


fragpath1 <- "~/workdir/STEFANO_SME/cellranger_sample1/cellranger_count_v1.1/outs/fragments.tsv.gz"
fragpath2 <- "~/workdir/STEFANO_SME/cellranger_sample2/outs/fragments.tsv.gz"
fragpath3 <- "~/workdir/STEFANO_SME/cellranger_sample3/outs/fragments.tsv.gz"


### GENE TO CELL MATRIX ###
hum_gtf_file = "~/workdir/reference_genomes/cellranger_reference/Hsa_GRCh38/HomSap_GRCh38/genes/genes.gtf"

# Keep only protein coding genes
hum_gtf = read.table(hum_gtf_file, sep = "\t", skip = 5, stringsAsFactors = F)
hum_genes = hum_gtf[(hum_gtf$V3 == "gene") & grepl("protein_coding", hum_gtf$V9),]

# Reshape to expand metadata
tmp = hum_genes$V9 %>% gsub(";", "", .) %>% strsplit(., " ") %>% do.call(rbind, .)
hum_genes$V9 = NULL
hum_genes = cbind(hum_genes, tmp)

# Reshape gtf derived df
hum_genes = hum_genes[, c(1,4,5,7,10,14,16,18)]
colnames(hum_genes) = c("chr", "start", "end", "strand",
			"gene_id", "gene_name", "gene_source", "gene_biotype")

# Turn into granges object
gtf_gr = makeGRangesFromDataFrame(hum_genes, keep.extra.columns = T)
gtf_gr_coords <- Extend(x = gtf_gr, upstream = 2000, downstream = 0)

# Split into list
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
gtf_coords = chunk(gtf_gr_coords, 10)

origcells = rownames(meta[meta$orig.ident == 'sme_atac1',]) %>% gsub("S1_", "", .)

# Gene to cell matrix for s1
gene_act_list = mclapply(gtf_coords, mc.cores = 10, function(x){gene_activities <- FeatureMatrix(
										  fragments = fragpath1,
										  features = x,
										  cells = origcells,
										  chunk = 5
										  )})

# Bind into one matrix
gene_activities = rbind.fill.matrix(gene_act_list)
gene_activities[is.na(gene_activities)] = 0

# Assign rownames
gene_coord_only = lapply(gene_act_list, function(x){rownames(x)})
gene_coord_only = unlist(gene_coord_only)
rownames(gene_activities) = gene_coord_only

# convert rownames from chromsomal coordinates into gene names
gene_key <- gtf_gr_coords$gene_name
names(gene_key) <- GRangesToString(grange = gtf_gr_coords)
rownames(gene_activities) <- gene_key[rownames(gene_activities)]

s1_geneact = gene_activities


# Gene to cell matrix for s2 #
origcells = rownames(meta[meta$orig.ident == 'sme_atac2',]) %>% gsub("S2_", "", .)

gene_act_list = mclapply(gtf_coords, mc.cores = 10, function(x){gene_activities <- FeatureMatrix(
										  fragments = fragpath2,
										  features = x,
										  cells = origcells,
										  chunk = 5
										  )})

# Bind into one matrix
gene_activities = rbind.fill.matrix(gene_act_list)
gene_activities[is.na(gene_activities)] = 0

# Assign rownames
gene_coord_only = lapply(gene_act_list, function(x){rownames(x)})
gene_coord_only = unlist(gene_coord_only)
rownames(gene_activities) = gene_coord_only

# convert rownames from chromsomal coordinates into gene names
gene_key <- gtf_gr_coords$gene_name
names(gene_key) <- GRangesToString(grange = gtf_gr_coords)
rownames(gene_activities) <- gene_key[rownames(gene_activities)]

s2_geneact = gene_activities


# Gene to cell matrix for s3 #
origcells = rownames(meta[meta$orig.ident == 'sme_atac3',]) %>% gsub("S3_", "", .)

gene_act_list = mclapply(gtf_coords, mc.cores = 10, function(x){gene_activities <- FeatureMatrix(
										  fragments = fragpath3,
										  features = x,
										  cells = origcells,
										  chunk = 5
										  )})

# Bind into one matrix
gene_activities = rbind.fill.matrix(gene_act_list)
gene_activities[is.na(gene_activities)] = 0

# Assign rownames
gene_coord_only = lapply(gene_act_list, function(x){rownames(x)})
gene_coord_only = unlist(gene_coord_only)
rownames(gene_activities) = gene_coord_only

# convert rownames from chromsomal coordinates into gene names
gene_key <- gtf_gr_coords$gene_name
names(gene_key) <- GRangesToString(grange = gtf_gr_coords)
rownames(gene_activities) <- gene_key[rownames(gene_activities)]

s3_geneact = gene_activities


# Save all in a list
all = list(s1_geneact, s2_geneact, s3_geneact)
saveRDS(all, "~/workdir/STEFANO_SME/all_geneact_matrices.RDS")

# Convert cell names
colnames(s1_geneact) = gsub("^", "S1_", colnames(s1_geneact))
colnames(s2_geneact) = gsub("^", "S2_", colnames(s2_geneact))
colnames(s3_geneact) = gsub("^", "S3_", colnames(s3_geneact))

# Combine them in one matrix
commongenes = Reduce(intersect, list(rownames(s1_geneact),
				     rownames(s2_geneact),
				     rownames(s3_geneact)))

mergedmat = cbind(s1_geneact[commongenes, ], s2_geneact[commongenes, ], s3_geneact[commongenes, ])

saveRDS(mergedmat, "~/workdir/STEFANO_SME/allsamples_merged_geneact_matrix.RDS")

