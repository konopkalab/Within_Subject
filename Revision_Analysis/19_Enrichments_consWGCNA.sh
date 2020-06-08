# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776, brain expressed = 15585, WGCNA list = 6029)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir enrichments_wgcna/

cp rawdata/geneset/*.RData enrichments_wgcna/
cp UTILS/Enrichment.r enrichments_wgcna/
cp processing_wgcna/ModuleOutput_WITHIN.txt enrichments_wgcna/

cd enrichments_wgcna/

mkdir STATS/

# scRNA healty
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l Allen_MultiReg_CellMarkers_GeneSet.RData -p -b 6029 -o STATS/Allen_Markers -W 7 -H 4
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l GeneSets_BBLake_2018.RData -p -b 6029 -o STATS/BBLake_Markers -W 7 -H 4
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l SME_Integrated_scRNA_GeneSets.RData -p -b 6029 -o STATS/BA38_Integrated_Markers -W 7 -H 4

# scRNA Disorders
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l ALZ_SingleCell_DEGs.RData -p -b 15585 -o STATS/ALZ_SingleCell -W 7 -H 4
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l ASD_SingleCell_DEGs.RData -p -b 15585 -o STATS/ASD_SingleCell -W 7 -H 4

# Neuropsy
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l PsychENCODE_DEGs.RData -p -b 15585 -o STATS/PSY_DEGS -W 12 -H 2
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l PsychEncode_Modules.RData -p -b 15585 -o STATS/PSY_MODS -W 7 -H 5
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l GeneSets_RareVar.RData -p -b 15585 -o STATS/RARE_VAR -W 7 -H 3
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l GeneSets_UltraRare_CaseSpecific.RData -p -b 15585 -o STATS/SCZ_UltraRare -W 7 -H 2
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l ASD_SFARI.RData -p -b 15585 -o STATS/ASD_Sfari -W 12 -H 2
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l GeneSets_UltraRare_SCZspecific.RData -p -b 15585 -o STATS/SCZspecific_UltraRare -W 7 -H 2
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l GeneSets_DeNovo_Var.RData -p -b 15585 -o STATS/DeNovo_Var -W 7 -H 3

# SME
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l New_Genes_SMEws.RData -p -b 15585 -o STATS/WithinSubj_SME -W 7 -H 2.5
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l GeneSets_CerCort.RData -p -b 15585 -o STATS/CerCort_SME -W 7 -H 2
Rscript Enrichment.r -g ModuleOutput_WITHIN.txt -l OtherGenes_SME.RData -p -b 15585 -o STATS/Behavioral_Genes -W 7 -H 2

rm *.RData

