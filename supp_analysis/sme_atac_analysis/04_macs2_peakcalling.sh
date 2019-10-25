#!/bin/bash

module load macs/2.1.2

cd ~/workdir/STEFANO_SME/fragments

macs2 callpeak  -t OligOPC_adj_fragments.txt \
                -f BEDPE \
                -g hs \
                --outdir OligOPC_adj_fragments_peaks \
                -n ba38_OligOPC_fragments \
                --nomodel \
                --keep-dup all \
                --extsize 1


macs2 callpeak  -t ExcInh_adj_fragments.txt \
                -f BEDPE \
                -g hs \
                --outdir ExcInh_adj_fragments_peaks \
                -n ba38_ExcInh_fragments \
                --nomodel \
                --keep-dup all \
                --extsize 1

