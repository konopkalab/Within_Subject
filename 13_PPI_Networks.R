# Covariate plot
rm(list=ls())
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pls))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(effsize))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rio))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(STRINGdb))
suppressPackageStartupMessages(library(igraph))

# Create Output directories
folder_names <- c("networking")
sapply(folder_names, dir.create)

#Load the final gene list of differentially correlated genes
finalgenes = read.table("processing_wgcna/ModuleOutput_WITHIN.txt", header = T)

names <- c("WM4","WM11","WM12","WM19","WM21","WM22")

# Create string db
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=400, input_directory="" )

# Loop across all the cool modules
input <- list()
mapped <- list()
hits <- list()
enrich <- list()
pdfname <- list()
for(i in 1:length(names))
	{
		pdfname[[i]]<-paste0("networking/",names[[i]],"_PPI_network.pdf") 
		pdf(file=pdfname[[i]],,width = 5,height=5)
		input[[i]] <- finalgenes %>% filter(ModuleName==names[[i]])
		mapped[[i]] <- string_db$map(input[[i]], "Gene", removeUnmappedRows = TRUE )
		enrich[[i]] <- string_db$get_ppi_enrichment(mapped[[i]]$STRING_id[1:nrow(input[[i]])])
 		hits[[i]] <- mapped[[i]]$STRING_id[1:nrow(input[[i]])]
		string_db$plot_network(hits[[i]])
		dev.off()
	}



