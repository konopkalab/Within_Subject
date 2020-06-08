suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(ggthemes))
source("UTILS/Utils.R")

# 
load("processing_memory_QuantReg/SME_QuantReg.RData")

df <- SME_Cor_QuantReg_Sign


# PIE CHART LR
df <- as.data.frame(table(SME_Cor_QuantReg_Sign$Waves))
colnames(df) <- c("Waves","Counts")

# Mutate the table in proportion
df <- df %>%
  arrange(desc(Waves)) %>%
  mutate(prop = round(Counts*100/sum(Counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)

pdf("processing_memory_QuantReg/PIE_QR.pdf",width=4,height=4,useDingbats=FALSE)
ggplot(df, aes(x = "", y = prop, fill = Waves)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  coord_polar("y", start = 0)+
  theme_void()+
scale_fill_manual(values=c("blue", "red","orange","magenta","green","grey60"))+
ggtitle("QR Genes per Waves: P < 0.05")
dev.off()

# N genes LR
df <- as.data.frame(table(SME_Cor_QuantReg_Sign$Waves))
df$Var1 <- factor(df$Var1,levels = c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta"))

pdf("processing_memory_QuantReg/Barplot_QR_Ngenes_005.pdf",width=4,height=4,useDingbats=FALSE)
ggbarplot(df, x = "Var1", y = "Freq",
          fill = "Var1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c("green","magenta", "red","blue","grey60","orange"),            
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45,           # Rotate vertically x axis texts
          label = TRUE, label.pos = "out",
          orientation = "horiz")+
xlab("Waves")+ 
ylab("# of Genes")+
theme(legend.position="none")+
ylim(0,1000)
dev.off()

# Stacked
df <- as.data.frame(table(SME_Cor_QuantReg_Sign$Waves))
df$Var1 <- factor(df$Var1,levels = c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta"))
df$Class <- "FDR"

tmp1 <- SME_Cor_QuantReg %>%
			filter(Pval < 0.05,BootP < 0.05,PermP < 0.05, BootP_All < 0.05)

df1 <- as.data.frame(table(tmp1$Waves))
df1$Var1 <- factor(df1$Var1,levels = c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta"))
df1$Freq <- df1$Freq - df$Freq
df1$Class <- "Pvalue"

df <- rbind(df,df1)

pdf("processing_memory_QuantReg/Barplot_QR_Ngenes_Stacked.pdf",width=5,height=4,useDingbats=FALSE)
ggbarplot(df, x = "Var1", y = "Freq",
          fill = "Class",               # change fill color by cyl
          color = "Var1",            # Set bar border colors to white
          palette = c("green","magenta", "red","blue","grey60","orange"),            
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45,           # Rotate vertically x axis texts
          label = TRUE, label.pos = "out",
          orientation = "horiz")+
xlab("Waves")+ 
ylab("# of Genes")+
theme(legend.position="right")+
ylim(0,1000)
dev.off()

# Intersection
pdf("processing_memory_QuantReg/WavesIntersection_BothCorrleation_LR.pdf",width=6,height=4,useDingbats=FALSE)
l <- split(as.character(SME_Cor_QuantReg_Sign$Gene),SME_Cor_QuantReg_Sign$Waves)
Waves <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Waves, ToTGene))
names(metadata) <- c("Waves", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),keep.order = F,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), list(type = "matrix_rows", 
    column = "sets", colors = c(Alpha = "blue", Beta = "red", Theta = "grey",Delta = "orange",Gamma = "purple",High.Gamma = "green"), 
    alpha = 0.5))))
dev.off()

# R squared plot
tmp <- SME_Cor_QuantReg_Sign
tmp$Rsq <- tmp$Rho^2
tmp$Waves <- factor(tmp$Waves, levels = c("Delta","Theta","Alpha","Beta","Gamma","High.Gamma"))

pdf("processing_memory_QuantReg/Rsq_Plot.pdf",width=5,height=4,useDingbats=FALSE)
ggline(tmp, 
    x = "Waves", 
    y = "Rsq", 
        add = c("mean_se"),
        color = "Waves", palette = c("orange","grey60","blue", "red","magenta","green"))+
    labs(title="",x="", y = "rho^2")+
    theme_classic()+
    theme(legend.position="none")+
    rotate_x_text(angle = 45)+
      stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.",label.y = 0.4)
dev.off()


pdf("processing_memory_QuantReg/Rsq_Plot_Violin.pdf",width=3,height=3,useDingbats=FALSE)
ggviolin(tmp, 
    x = "Waves", 
    y = "Rsq", 
        add = c("mean_sd"),
        color = "Waves", palette = c("orange","grey60","blue", "red","magenta","green"))+
    labs(title="",x="", y = "rho^2")+
    theme_classic()+
    theme(legend.position="none")+
    rotate_x_text(angle = 45)+
      stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.",label.y = 0.8)
dev.off()