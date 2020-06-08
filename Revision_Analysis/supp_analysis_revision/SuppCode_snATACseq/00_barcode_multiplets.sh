#!/bin/bash
dir=~/onlybiohpc/STEFANO_SME/BAM_based_analysis

source ~/programs/cellranger-atac-1.1.0/sourceme.bash

# Human #
# First sample
cd  ~/workdir/STEFANO_SME/cellranger_sample1/cellranger_count_v1.1/outs/
python $dir/clean_barcode_multiplets_1.1.py ~/workdir/STEFANO_SME/cellranger_sample1/cellranger_count_v1.1/outs/

# Second sample
cd ~/workdir/STEFANO_SME/cellranger_sample2/outs/
python $dir/clean_barcode_multiplets_1.1.py ~/workdir/STEFANO_SME/cellranger_sample2/outs/

# Third sample
cd ~/workdir/STEFANO_SME/cellranger_sample3/outs/
python $dir/clean_barcode_multiplets_1.1.py ~/workdir/STEFANO_SME/cellranger_sample3/outs/



