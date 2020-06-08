#!/bin/bash

## This script calls peaks for a subset of cells.
## Peak calling will be done using all cutsites instead of fragments.
## We will use BED format to account for both cutsites.


# Load modules
module load macs/2.1.2
module load bedops
module load bedtools/2.29.0

cd ~/workdir/STEFANO_SME/sample1_analysis/peak_bam

# Separate header and the rest of sam file
samtools view -@ 23 possorted_bam.bam > nohead.sam
samtools view -H possorted_bam.bam > onlyhead.sam

# Divide the sam file to parallalize across 23 cpus.
parts=$(echo `wc -l < nohead.sam` '/ 23' | bc)
split -l $parts -d nohead.sam temp
rm nohead.sam


# Function to convert sam to bed while keeping the cell barcode.
sam2bed_par(){
cat onlyhead.sam $1 > $1".head"
sam2bed --do-not-sort < $1".head" > $1".init.bed"
cut -f1-3,6 $1".init.bed" > $1".tmp1"
sed 's/.*CR:Z:\(.\{16\}\).*/\1/' $1".init.bed" > $1".tmp2"
paste $1".tmp1" $1".tmp2" > $1".final.bed"
rm $1".head" $1".init.bed" $1".tmp1" $1".tmp2"
}


# Run the function in parallel
for k in temp*
do
sam2bed_par $k &
done

wait

# Generate final bed file with barcode and strand information
cat *.final.bed > allmapped.bed
LC_ALL=C sort -k1,1 -k2,2n allmapped.bed > allmapped_sorted.bed
rm temp*
rm allmapped.bed


######### In R #######################
## Filter all reads for ExcInh/OligOPC barcodes in R per sample
library(data.table)
library(tidyr)
smp = 'sample2'
allbed = fread(paste0(smp, "_allmapped_sorted.bed"))
allbed$V5 = gsub("^", "S2_", allbed$V5) %>% gsub("$", "-1", .)

# ExcInh cells
cells = fread("ExcInhcells.txt", header = F)
ov = allbed[allbed$V5 %in% cells$V1,]
write.table(ov, paste0(smp, "_ExcInh.bed"), row.names = F, col.names = F, quote = F, sep = '\t')

# OligOPC cells
cells = fread("OligOPCcells.txt", header = F)
ov = allbed[allbed$V5 %in% cells$V1,]
write.table(ov, paste0(smp, "_OligOPC.bed"), row.names = F, col.names = F, quote = F, sep = '\t')

######################################


# Merge all samples per cell group
cat *_ExcInh.bed > all_ExcInh.bed
cat *_OligOPC.bed > all_OligOPC.bed

# Call peaks using BED file
macs2 callpeak  -t all_ExcInh.bed \
		-f BED \
                -g hs \
                --outdir peaks_ExcInh \
                -n hsa \
                --nomodel \
                --keep-dup all \
                --extsize 200 \
		--shift -100

macs2 callpeak  -t all_OligOPC.bed \
		-f BED \
                -g hs \
                --outdir peaks_OligOPC \
                -n hsa \
                --nomodel \
                --keep-dup all \
                --extsize 200 \
		--shift -100


# Create sorted peaks file to run cellranger reanalyze
cd ~/workdir/STEFANO_SME/allsamples/peaks_ExcInh/peaks/
cut -f1-3 sme_ExcInh_peaks.narrowPeaks > simple_peaks.bed
grep -v 'random\|Un\|M' simple_peaks.bed | sort -k1,1V -k2,2n -k3,3n > simple_peaks_sorted.bed

cd ~/workdir/STEFANO_SME/allsamples/peaks_OligOPC/peaks/
cut -f1-3 sme_OligOPC_peaks.narrowPeaks > simple_peaks.bed
grep -v 'random\|Un\|M' simple_peaks.bed | sort -k1,1V -k2,2n -k3,3n > simple_peaks_sorted.bed

# Run cellranger reanalyze to create peak annotation file.
cd /home2/s422159/workdir/STEFANO_SME/allsamples
~/programs/cellranger-atac-1.1.0/./cellranger-atac reanalyze --id=SME_ExcInh \
							     --peaks=peaks_ExcInh/peaks/simple_peaks_sorted.bed\
							     --fragments=~/workdir/STEFANO_SME/sample1_analysis/fragments/test_sorted.bed.gz \
							     --reference=~/workdir/reference_genomes/cellranger_reference/refdata-cellranger-atac-GRCh38-1.1.0/


cd /home2/s422159/workdir/STEFANO_SME/allsamples
~/programs/cellranger-atac-1.1.0/./cellranger-atac reanalyze --id=SME_OligOPC \
							     --peaks=peaks_OligOPC/peaks/simple_peaks_sorted.bed\
							     --fragments=~/workdir/STEFANO_SME/sample1_analysis/fragments/test_sorted.bed.gz \
							     --reference=~/workdir/reference_genomes/cellranger_reference/refdata-cellranger-atac-GRCh38-1.1.0/

# END OF SCRIPT

