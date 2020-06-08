#!/bin/bash

module load UCSC_userApps/v317
module load bedtools
#module load bedops
module load samtools/gcc/1.8

cd /home2/s422159/workdir/STEFANO_SME/allsamples/forIGV/
ln -s ~/workdir/STEFANO_SME/allsamples/all_ExcInh.bed .
ln -s ~/workdir/STEFANO_SME/allsamples/all_OligOPC.bed .

# Create bw file for peaks
fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' ../peaks_ExcInh/simple_peaks_sorted2.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > peaks.bedGraph
bedGraphToBigWig peaks.bedGraph hg38.chrom.sizes peaks.bw


### SHANK2 ###

# First get all reads from the 20kb upstream region, convert to bigwig files.
int="chr11\t70661818\t70682197"
echo -e "$int" > tmp.bed

grep "chr11" all_ExcInh.bed > tmp2.bed
bedtools sort -i tmp2.bed > tmp3.bed
bedtools intersect -sorted -a tmp.bed -b tmp3.bed > up_20kb.bed

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes shank2_up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_SHANK2_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_SHANK2_20kb_ups_hits.bedGraph
bedGraphToBigWig SMAD3_SHANK2_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_SHANK2_20kb_ups_hits.bw

awk '{print}' GLIS1_SHANK2_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > GLIS1_SHANK2_20kb_ups_hits.bedGraph
bedGraphToBigWig GLIS1_SHANK2_20kb_ups_hits.bedGraph hg38.chrom.sizes GLIS1_SHANK2_20kb_ups_hits.bw

### CACNA1C ###

# First get all reads from the 20kb upstream region, convert to bigwig files.
int="chr12\t2033563\t2053563"
echo -e "$int" > tmp.bed

grep "chr12" all_ExcInh.bed > tmp2.bed
bedtools sort -i tmp2.bed > tmp3.bed
bedtools intersect -sorted -a tmp.bed -b tmp3.bed > up_20kb.bed

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes cacna1c_up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_CACNA_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_CACNA1C_20kb_ups_hits.bedGraph
bedGraphToBigWig SMAD3_CACNA1C_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_CACNA1C_20kb_ups_hits.bw

awk '{print}' GLIS1_CACNA_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > GLIS1_CACNA1C_20kb_ups_hits.bedGraph
bedGraphToBigWig GLIS1_CACNA1C_20kb_ups_hits.bedGraph hg38.chrom.sizes GLIS1_CACNA1C_20kb_ups_hits.bw




### TRIO ###

# First get all reads from the 20kb upstream region, convert to bigwig files.
int="chr5\t14123702\t14143702"
echo -e "$int" > tmp.bed

grep "chr5" all_ExcInh.bed > tmp2.bed
bedtools sort -i tmp2.bed > tmp3.bed
bedtools intersect -sorted -a tmp.bed -b tmp3.bed > up_20kb.bed

awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes TRIO_up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_TRIO_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_TRIO_20kb_ups_hits.bedGraph
bedGraphToBigWig SMAD3_TRIO_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_TRIO_20kb_ups_hits.bw

awk '{print}' GLIS1_TRIO_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > GLIS1_TRIO_20kb_ups_hits.bedGraph
bedGraphToBigWig GLIS1_TRIO_20kb_ups_hits.bedGraph hg38.chrom.sizes GLIS1_TRIO_20kb_ups_hits.bw




### NEDD4L ###
# First get all reads from the 20kb upstream region, convert to bigwig files.
int="chr18\t58024548\t58044548"
echo -e "$int" > tmp.bed

grep "chr18" all_ExcInh.bed > tmp2.bed
bedtools sort -i tmp2.bed > tmp3.bed
bedtools intersect -sorted -a tmp.bed -b tmp3.bed > up_20kb.bed

awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes NEDD4L_up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_NEDD4L_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_NEDD4L_20kb_ups_hits.bedGraph
bedGraphToBigWig SMAD3_NEDD4L_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_NEDD4L_20kb_ups_hits.bw

awk '{print}' GLIS1_NEDD4L_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > GLIS1_NEDD4L_20kb_ups_hits.bedGraph
bedGraphToBigWig GLIS1_NEDD4L_20kb_ups_hits.bedGraph hg38.chrom.sizes GLIS1_NEDD4L_20kb_ups_hits.bw



