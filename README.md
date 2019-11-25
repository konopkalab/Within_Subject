Within Subject data analysis
==========================

This repository contains analysis code for the Within Subject project carried out by researchers at the [Konopka Lab, UTSW](http://konopkalab.org/).

## Citation

If you use anything in this repository please cite the following publication:

Pre-print URL: https://www.biorxiv.org/content/early/2019/11/25/853531.full.pdf

## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`processing_qc`](processing_qc/) | Output data from initial processing and quality check. | 01_Data_processing_QC.R |
| [`processing_memory`](processing_memory/) | Output data from memory (SME) analysis. | 02_SME_Analysis.R |
| [`processing_math`](processing_math/) | Output data from math task analysis. | 03_MATH_Analysis.R |
| [`processing_mri`](processing_mri/) | Output data from MRI (thickness) analysis. | 04_MRI_Analysis.R |
| [`processing_behavior`](processing_behavior/) | Output data from behavioral analysis. | 05_BEHAVIORAL_Analysis.R |
| [`final_visualizations`](final_visualizations/) | Some visualization and data integration. | 06_Visualizations.R |
| [`enrichments_SME`](enrichments_SME/) | Enrichment analysis for SME genes. | 07_Cross-Enrich_SME.R |
| [`processing_wgcna`](processing_wgcna/) | Output data from the consensus WGCNA analysis. | 08_consWGCNA.R |
| [`enrichments_wgcna`](enrichments_wgcna/) | Enrichment analysis for co-expression modules. | 09_Enrichments_consWGCNA.R |
| [`magma_wgcna`](magma_wgcna/) | GWAS enrichment for the co-expression modules. | 10_consWGCNA_Magma.sh |
| [`processing_scRNAseq`](processing_scRNAseq/) | Output data from single-nuclei RNA-seq analysis. | 11_SingleCell_Analysis.R |
| [`supp_tables`](supp_tables/) | Databases and supplementary tables. | 12_Database.R |
| [`networking`](networking/) | Output data from PPI netowrk analysis. | 13_PPI_Networks.R |
| [`supp_analysis`](supp_analysis/) | Supplementary data and analysis. |
| [`Shiny_App`](Shiny_App/) | Shiny app for SME - Gene Scatterplot. |

![](ScreenShot.png)]
