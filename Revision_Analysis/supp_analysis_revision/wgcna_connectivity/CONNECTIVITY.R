### Merging multiple dataframe by rownames
## list of tables with the gene counts

library(ggpubr)
library(ggplot2)
library(ggsignif)
library(tidyverse)

mod <- read.table("ModuleOutput_WITHIN.txt",header=T)
load("SME_QuantReg.RData")
sme <- SME_Cor_QuantReg_Sign

names <- c("WM4","WM11","WM12","WM19","WM21","WM22")
input <- mod %>% 
			filter(ModuleName%in%names) %>% 
			mutate(Class=case_when(Gene %in% sme$Gene ~ "SME", !(Gene %in% sme$Gene) ~ "BKG")) %>%
			mutate(LOG = log2(kWithin))


input$ModuleName <- factor(input$ModuleName,levels=c("WM4","WM11","WM12","WM19","WM21","WM22"))

pdf("Connecivity_SME_genes.pdf",width=3,height=5)
my_comparisons <- list( c("SME", "BKG"))
ggboxplot(input, "Class", "LOG",
fill = "Class", 
palette = c("grey", "red"),facet.by="ModuleName",notch = TRUE)+ 
stat_compare_means(comparisons = my_comparisons,label = "p.signif",method="t.test",label.y=9)+ # Add pairwise comparisons p-value
xlab("")+
ylab("log2(Connectivity)")+
theme_classic()+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
theme(legend.position="none")
dev.off()
