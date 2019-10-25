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
#SBATCH --job-name pr3_sc4

# Name of the SLURM partition that this job should run on.
#SBATCH -p super       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Memory (RAM) requirement/limit in MB.
#SBATCH --mem 254976      # Memory Requirement (MB)

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 1-0:0:00
# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o pr3_sc4_%j.out
#SBATCH -e pr3_sc4_%j.err

# Create cellranger style fragments file
cd /home2/s422159/workdir/STEFANO_SME/fragments
bedtools sort -i ExcInh_adj_fragments.txt > ExcInh_adj_fragments_sorted.bed
bgzip ExcInh_adj_fragments_sorted.bed
tabix -p bed ExcInh_adj_fragments_sorted.bed.gz

bedtools sort -i OligOPC_adj_fragments.txt > OligOPC_adj_fragments.bed
bgzip OligOPC_adj_fragments_sorted.bed
tabix -p bed OligOPC_adj_fragments_sorted.bed.gz


# Sort the peaks
cd ../peak_calling/ExcInh_adj_fragments_peaks
sort -k1,1V -k2,2n -k3,3n simple_bed_peaks.bed > simple_bed_peaks_cellranger_sort.bed

cd ../OligOPC_adj_fragments_peaks
sort -k1,1V -k2,2n -k3,3n simple_bed_peaks.bed > simple_bed_peaks_cellranger_sort.bed


# Cellranger reanalyze to create peak annotation file
cd /home2/s422159/workdir/STEFANO_SME
~/programs/cellranger-atac-1.1.0/./cellranger-atac reanalyze --id=SME_EXCINH \
							     --peaks=peak_calling/ExcInh_adj_fragments_peaks/simple_bed_peaks_cellranger_sort.bed\
							     --fragments=fragments/ExcInh_adj_fragments_sorted.bed.gz \
							     --reference=~/workdir/reference_genomes/cellranger_reference/refdata-cellranger-atac-GRCh38-1.1.0/


~/programs/cellranger-atac-1.1.0/./cellranger-atac reanalyze --id=SME_OligOPC \
							     --peaks=peak_calling/OligOPC_adj_fragments_peaks/simple_bed_peaks_cellranger_sort.bed\
							     --fragments=fragments/OligOPC_adj_fragments_sorted.bed.gz \
							     --reference=~/workdir/reference_genomes/cellranger_reference/refdata-cellranger-atac-GRCh38-1.1.0/


# END OF SCRIPT

