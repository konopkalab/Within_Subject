#!/bin/bash
#
# CREATED USING THE BIOHPC PORTAL on Fri Jun 21 2019 22:56:40 GMT-0500 (Central Daylight Time)
#
# This file is batch script used to run commands on the BioHPC cluster.
# The script is submitted to the cluster using the SLURM `sbatch` command.
# Lines starting with # are comments, and will not be run.
# Lines starting with #SBATCH specify options for the scheduler.
# Lines that do not start with # or #SBATCH are commands that will run.

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name 09_

# Name of the SLURM partition that this job should run on.
#SBATCH -p super       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Memory (RAM) requirement/limit in MB.
#SBATCH --mem 254976      # Memory Requirement (MB)

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 2-0:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o 09_%j.out
#SBATCH -e 09_%j.err

module load UCSC_userApps/v317
module load bedtools
module load samtools/gcc/1.8

### SHANK2 ###

# First get all reads from the 20kb upstream region, convert to bigwig files.
cd /home2/s422159/workdir/STEFANO_SME/STEFANO_SME_WM4/outs/tmp
int="chr11\t70661818\t70682197"
echo -e "$int" > tmp.bed

grep "chr11" fragments.tsv > tmp_shank2.tsv
bedtools intersect -a tmp.bed -b tmp_shank2.tsv > ~/workdir/STEFANO_SME/IGV_related/shank2/up_20kb.bed


cd ~/workdir/STEFANO_SME/IGV_related/shank2

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_SHANK2_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_SHANK2_20kb_ups_hits.bedGraph

bedGraphToBigWig SMAD3_SHANK2_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_SHANK2_20kb_ups_hits.bw


### CACNA1C ###


module load UCSC_userApps/v317
module load bedtools
module load samtools/gcc/1.8


# First get all reads from the 20kb upstream region, convert to bigwig files.
cd /home2/s422159/workdir/STEFANO_SME/STEFANO_SME_WM4/outs/tmp
int="chr12\t2033563\t2053563"
echo -e "$int" > tmp.bed

grep "chr12" fragments.tsv > tmp_cacna1c.tsv
bedtools intersect -a tmp.bed -b tmp_cacna1c.tsv > ~/workdir/STEFANO_SME/IGV_related/cacna1c/cacna1c_up_20kb.bed


cd ~/workdir/STEFANO_SME/IGV_related/cacna1c

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' cacna1c_up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > cacna1c_up_20kb.bedGraph

bedGraphToBigWig cacna1c_up_20kb.bedGraph hg38.chrom.sizes cacna1c_up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_CACNA_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_CACNA_20kb_ups_hits.bedGraph

bedGraphToBigWig SMAD3_CACNA_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_CACNA_20kb_ups_hits.bw


### SLC7A6 ###


module load UCSC_userApps/v317
module load bedtools
module load samtools/gcc/1.8

# First get all reads from the 20kb upstream region, convert to bigwig files.
cd /home2/s422159/workdir/STEFANO_SME/STEFANO_SME_WM4/outs/tmp
int="chr16\t68244516\t68264516"
echo -e "$int" > tmp2.bed

grep "chr16" fragments.tsv > tmp_slc7a6.tsv
bedtools intersect -a tmp2.bed -b tmp_slc7a6.tsv > ~/workdir/STEFANO_SME/IGV_related/slc7a6/up_20kb.bed


cd ~/workdir/STEFANO_SME/IGV_related/slc7a6

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_SLC7A6_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_SLC7A6_20kb_ups_hits.bedGraph

bedGraphToBigWig SMAD3_SLC7A6_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_SLC7A6_20kb_ups_hits.bw



### TRIO ###

module load UCSC_userApps/v317
module load bedtools
module load samtools/gcc/1.8

# First get all reads from the 20kb upstream region, convert to bigwig files.
cd /home2/s422159/workdir/STEFANO_SME/STEFANO_SME_WM4/outs/tmp
int="chr5\t14123702\t14143702"
echo -e "$int" > tmp3.bed

grep "chr5" fragments.tsv > tmp_trio.tsv
bedtools intersect -a tmp3.bed -b tmp_trio.tsv > ~/workdir/STEFANO_SME/IGV_related/trio/up_20kb.bed


cd ~/workdir/STEFANO_SME/IGV_related/trio

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_TRIO_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_TRIO_20kb_ups_hits.bedGraph

bedGraphToBigWig SMAD3_TRIO_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_TRIO_20kb_ups_hits.bw



### NEDD4L ###

module load UCSC_userApps/v317
module load bedtools
module load samtools/gcc/1.8

# First get all reads from the 20kb upstream region, convert to bigwig files.
cd /home2/s422159/workdir/STEFANO_SME/STEFANO_SME_WM4/outs/tmp
int="chr18\t58024548\t58044548"
echo -e "$int" > tmp4.bed

grep "chr18" fragments.tsv > tmp_nedd4l.tsv
bedtools intersect -a tmp4.bed -b tmp_nedd4l.tsv > ~/workdir/STEFANO_SME/IGV_related/nedd4l/up_20kb.bed


cd ~/workdir/STEFANO_SME/IGV_related/nedd4l

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' SMAD3_NEDD4L_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > SMAD3_NEDD4L_20kb_ups_hits.bedGraph

bedGraphToBigWig SMAD3_NEDD4L_20kb_ups_hits.bedGraph hg38.chrom.sizes SMAD3_NEDD4L_20kb_ups_hits.bw


##### MEF2A HITS #####

### CNTNAP5 ###

# First get all reads from the 20kb upstream region, convert to bigwig files.
cd /home2/s422159/workdir/STEFANO_SME/STEFANO_SME_WM4/outs/tmp
int="chr2\t124005287\t124025287"
echo -e "$int" > tmp.bed

grep "chr2" fragments.tsv > tmp_cntnap5.tsv
mkdir ~/workdir/STEFANO_SME/IGV_related/cntnap5
bedtools intersect -a tmp.bed -b tmp_cntnap5.tsv > ~/workdir/STEFANO_SME/IGV_related/cntnap5/up_20kb.bed


cd ~/workdir/STEFANO_SME/IGV_related/cntnap5

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' MEF2A_CNTNAP5_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > MEF2A_CNTNAP5_20kb_ups_hits.bedGraph

bedGraphToBigWig MEF2A_CNTNAP5_20kb_ups_hits.bedGraph hg38.chrom.sizes MEF2A_CNTNAP5_20kb_ups_hits.bw



### LNX1 ###

# First get all reads from the 20kb upstream region, convert to bigwig files.
cd /home2/s422159/workdir/STEFANO_SME/STEFANO_SME_WM4/outs/tmp
int="chr4\t53591616\t53611616"
echo -e "$int" > tmp.bed

grep "chr4" fragments.tsv > tmp_lnx1.tsv
mkdir ~/workdir/STEFANO_SME/IGV_related/lnx1
bedtools intersect -a tmp.bed -b tmp_lnx1.tsv > ~/workdir/STEFANO_SME/IGV_related/lnx1/up_20kb.bed


cd ~/workdir/STEFANO_SME/IGV_related/lnx1

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' MEF2A_LNX1_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > MEF2A_LNX1_20kb_ups_hits.bedGraph

bedGraphToBigWig MEF2A_LNX1_20kb_ups_hits.bedGraph hg38.chrom.sizes MEF2A_LNX1_20kb_ups_hits.bw



### DGKI ###

# First get all reads from the 20kb upstream region, convert to bigwig files.
cd /home2/s422159/workdir/STEFANO_SME/STEFANO_SME_WM4/outs/tmp
int="chr7\t137846536\t137866536"
echo -e "$int" > tmp.bed

grep "chr7" fragments.tsv > tmp_dgki.tsv
mkdir ~/workdir/STEFANO_SME/IGV_related/dgki
bedtools intersect -a tmp.bed -b tmp_dgki.tsv > ~/workdir/STEFANO_SME/IGV_related/dgki/up_20kb.bed


cd ~/workdir/STEFANO_SME/IGV_related/dgki

fetchChromSizes hg38 > hg38.chrom.sizes
awk '{print}' up_20kb.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > up_20kb.bedGraph

bedGraphToBigWig up_20kb.bedGraph hg38.chrom.sizes up_20kb.bw

# Then get motif occurrence - peak intersections (20kb upstream), convert to bigwig files.

awk '{print}' MEF2A_DGKI_20kb_ups_hits.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > MEF2A_DGKI_20kb_ups_hits.bedGraph

bedGraphToBigWig MEF2A_DGKI_20kb_ups_hits.bedGraph hg38.chrom.sizes MEF2A_DGKI_20kb_ups_hits.bw

