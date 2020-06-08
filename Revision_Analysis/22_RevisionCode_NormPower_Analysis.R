# Load libraries
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(sqldf))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(made4))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(VennDiagram))
source("UTILS/Utils.R")

# Create Output directories
dir.create("normpower_results")

# Load inputs
load(here("rawdata","expAdj.RData"))
Lateral_resected_NormPower <- read.table("rawdata/NormPower_LR_Electrode.txt",header=T,sep="\t")

# Average
Lateral_resected_Average  <- Lateral_resected_NormPower %>%
                                  select(-Region,-Electrode) %>% 
                                  group_by(ID) %>%
                                  summarise_all(funs(mean)) %>% 
                                  as.data.frame() %>%
                                  column_to_rownames("ID") %>%
                                  t()

# LR
res <- list()
for(i in 1:nrow(Lateral_resected_Average))
{
res[[i]] <- FastCor(within_data,Lateral_resected_Average[i,],method="spearman",alternative="two.sided",cores=12,override=TRUE)%>%
                      as.data.frame() %>%
                      rownames_to_column('Gene') %>%
                      mutate(FDR = p.adjust(Pval, method="BH")) 
}

SME_cor_LR_NormPower <- data.frame()
for(i in 1:length(res)){
	res[[i]]$Waves <- rep(paste(rownames(Lateral_resected_Average)[i]),nrow(res[[i]]))
}

SME_cor_LR_NormPower <- do.call(rbind,res)
SME_cor_LR_NormPower <- SME_cor_LR_NormPower[c(1,5,2,3,4)]



colnames(SME_cor_LR_NormPower) <- c("Gene","Waves","Rho_NormPower","Pval_NormPower","FDR_NormPower")
save(SME_cor_LR_NormPower, file = "normpower_results/SME_cor_LR_NormPower.RData" )

# Load original data
load("processing_memory_QuantReg/SME_QuantReg.RData")

Comb <- Reduce(dplyr::full_join, list(SME_Cor_QuantReg, SME_cor_LR_NormPower))

# Compare with original
pdf(file="normpower_results/GeneCorrelations_ScatterCor.pdf",width=8,height=8)
ggscatter(Comb, x = "Rho", y = "Rho_NormPower",
   color = "lightgray", size = 0.5, 
   add = "reg.line",  
   add.params = list(color = "red", fill = "blue"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.sep = "\n"),
   ellipse = TRUE)+ 
xlab("Observed rho") + ylab("NormPower rho") +
geom_hline(yintercept = 0, linetype="dotted", color = "black", size=1)+ 
geom_vline(xintercept = 0, linetype="dotted", color = "black", size=1)+
facet_wrap(~Waves)
dev.off()

# 

sign <- Reduce(dplyr::left_join, list(SME_Cor_QuantReg_Sign, SME_cor_LR_NormPower)) %>%
            mutate(log10_Obs = -log10(Pval), log10_NormPower = -log10(Pval_NormPower)) 

pdf(file="normpower_results/GeneCorrelations_SMEgenes_ScatterCor_AllElectrodes.pdf",width=8,height=6)
ggscatter(sign, x = "Rho", y = "Rho_NormPower",
   color = "lightgray", size = 0.5, 
   add = "reg.line",  
   add.params = list(color = "red", fill = "blue"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.sep = "\n"),
   ellipse = TRUE)+ 
xlab("Observed rho") + ylab("NormPower rho") +
geom_hline(yintercept = 0, linetype="dotted", color = "black", size=1)+ 
geom_vline(xintercept = 0, linetype="dotted", color = "black", size=1)+
facet_wrap(~Waves)
dev.off()

pdf(file="normpower_results/Log10_Pvalue_SMEgenes_ScatterCor_AllElectrodes.pdf",width=8,height=6)
ggscatter(sign, x = "log10_Obs", y = "log10_NormPower",
   color = "lightgray", size = 0.5, 
   add = "reg.line",  
   add.params = list(color = "steelblue", fill = "gray"), 
   conf.int = TRUE, 
   cor.coef = TRUE, 
   cor.coeff.args = list(method = "pearson", label.sep = "\n"),
   ellipse = FALSE)+ 
xlab("-log10(observed P-value)") + ylab("-log10(NormPower P-value)") +
geom_hline(yintercept = 1.3, linetype="dotted", color = "black", size=1)+ 
geom_vline(xintercept = 1.3, linetype="dotted", color = "black", size=1)+
facet_wrap(~Waves)
dev.off()

sign <- sign %>%
        mutate(NormPower_Sign = case_when(Pval_NormPower < 0.05 ~ "NormPower_Sign", Pval_NormPower > 0.05  ~ "SME_Sign"))

pdf("normpower_results/Barplot_SME_NormPower_Saved.pdf",width=4,height=3,useDingbats=FALSE)
df <- reshape2::melt(table(sign$Waves,sign$NormPower_Sign))
df$Var1 <- factor(df$Var1,levels = c("Delta","Theta","Alpha","Beta","Gamma","High.Gamma"))
ggbarplot(df, 
  "Var1", 
  "value",
     fill = "Var2", 
     color = "Var2", 
     palette = c("red","black")) + 
  xlab("") + 
  ylab("# of Observed SME genes")+
  theme_classic()+
  theme(legend.position="right")+
  rotate_x_text(angle = 45)
dev.off()






# Venns
a <- Comb %>% filter(Waves == "Delta", Pval < 0.05)
b <- Comb %>% filter(Waves == "Delta", Pval_NormPower < 0.05)

ven <- venn.diagram(x= list(
Obs = a$Gene,
NormPow = b$Gene),
filename = NULL,
fill = c("salmon","olivedrab"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2,
lty = "dashed",
main="Delta")

pdf(file="normpower_results/Venn_Delta.pdf",width=4,height=4)
grid.draw(ven)
dev.off()


# Venns
a <- Comb %>% filter(Waves == "Theta", Pval < 0.05)
b <- Comb %>% filter(Waves == "Theta", Pval_NormPower < 0.05)

ven <- venn.diagram(x= list(
Obs = a$Gene,
NormPow = b$Gene),
filename = NULL,
fill = c("salmon","olivedrab"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2,
lty = "dashed",
main="Theta")

pdf(file="normpower_results/Venn_Theta.pdf",width=4,height=4)
grid.draw(ven)
dev.off()


# Venns
a <- Comb %>% filter(Waves == "Alpha", Pval < 0.05)
b <- Comb %>% filter(Waves == "Alpha", Pval_NormPower < 0.05)

ven <- venn.diagram(x= list(
Obs = a$Gene,
NormPow = b$Gene),
filename = NULL,
fill = c("salmon","olivedrab"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2,
lty = "dashed",
main="Alpha")

pdf(file="normpower_results/Venn_Alpha.pdf",width=4,height=4)
grid.draw(ven)
dev.off()

# Venns
a <- Comb %>% filter(Waves == "Beta", Pval < 0.05)
b <- Comb %>% filter(Waves == "Beta", Pval_NormPower < 0.05)

ven <- venn.diagram(x= list(
Obs = a$Gene,
NormPow = b$Gene),
filename = NULL,
fill = c("salmon","olivedrab"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2,
lty = "dashed",
main="Beta")

pdf(file="normpower_results/Venn_Beta.pdf",width=4,height=4)
grid.draw(ven)
dev.off()

# Venns
a <- Comb %>% filter(Waves == "Gamma", Pval < 0.05)
b <- Comb %>% filter(Waves == "Gamma", Pval_NormPower < 0.05)

ven <- venn.diagram(x= list(
Obs = a$Gene,
NormPow = b$Gene),
filename = NULL,
fill = c("salmon","olivedrab"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2,
lty = "dashed",
main="Gamma")

pdf(file="normpower_results/Venn_Gamma.pdf",width=4,height=4)
grid.draw(ven)
dev.off()

# Venns
a <- Comb %>% filter(Waves == "High.Gamma", Pval < 0.05)
b <- Comb %>% filter(Waves == "High.Gamma", Pval_NormPower < 0.05)

ven <- venn.diagram(x= list(
Obs = a$Gene,
NormPow = b$Gene),
filename = NULL,
fill = c("salmon","olivedrab"),
alpha = 0.50, 
cex = 1, 
fontfamily = "serif", 
cat.cex = 1.5,
cat.fontfamily = "serif",
cat.dist = 0.2,
margin = 0.2,
lty = "dashed",
main="High.Gamma")

pdf(file="normpower_results/Venn_HighGamma.pdf",width=4,height=4)
grid.draw(ven)
dev.off()