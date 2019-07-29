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
source("Utils.R")

# Covariates
dir.create("final_visualizations")
load(here("rawdata","expData.RData"))
load(here("rawdata","expAdj.RData"))
load(here("rawdata","SME_Values.RData"))
memory.genes <- read.table("processing_memory/SME_Significant_Genes.txt")
pheno.ws$Batch <- as.factor(pheno.ws$Batch)

PCA <- prcomp(t(within_data))
PCi<-data.frame(PCA$x,Labels = pheno.ws$UT.code)
eig <- (PCA$sdev)^2
variance <- eig*100/sum(eig)

pcaExp_WS <- ggscatter(PCi, x = "PC1", y = "PC2",color = "steelblue", size = 3,repel = TRUE,label="Labels")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()+
xlim(-50,50) + 
ylim(-40,40) 

PCA2<- prcomp(t(boot_data_all))
PCi2<-data.frame(PCA2$x,Labels = rownames(pheno_all),Class=pheno_all$Class)
eig <- (PCA2$sdev)^2
variance <- eig*100/sum(eig)

pcaExp_All <- ggscatter(PCi2, 
              x = "PC1", 
              y = "PC2",
              color = "Class",
              size = 3,
              shape="Class",
              palette = c("#00AFBB", "steelblue", "#FC4E07"))+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()+
xlim(-100,100) + 
ylim(-50,50) +
theme(legend.title = element_text(size=12, color = "black"),
           legend.justification=c(0,1), 
           legend.position=c(0.05, 0.95),
           legend.background = element_blank(),
           legend.key = element_blank())

plot2by2 <- cowplot::plot_grid(pcaExp_WS, pcaExp_All,labels=c("A", "B"), ncol = 2)
cowplot::save_plot("final_visualizations/PCAs.pdf", plot2by2, ncol = 2,base_height=4,base_width=4)

# Demo analysis
pdf("final_visualizations/Demo_Pairs.pdf",width=12,height=10)
p <- ggpairs(pheno_all, aes(color = Class),columns = names(pheno_all)[2:10]) + 
theme_classic()
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c("#00AFBB", "steelblue", "#FC4E07")) +
        scale_color_manual(values=c("#00AFBB", "steelblue", "#FC4E07"))  
  }
}
print(p)
dev.off()

# Var Exp


exp <- exp[,match(rownames(pheno_all),colnames(exp))]
demo <- pheno_all[c(-11)] %>% droplevels()
demo[is.na(demo)] <- 0
var <- VarExp(exp,demo,10,FALSE)

plot.dat <- data.frame(covar=names(var), prop=var)

pdf("final_visualizations//Variance_Explained_All.pdf",width=3,height=4)
ggbarplot(plot.dat, "covar", "prop", orientation = "horiz",color = "steelblue",label = TRUE, lab.vjust=-0.5,lab.hjust = -0.2,lab.nb.digits = 2)+
scale_x_discrete(limits=names(var))+
ylim(0,0.5)+
xlab("")+
ylab("")+
theme_hc()
dev.off()



# SME distribution 
lr <- as.data.frame(Lateral_resected)
lr$SME <- rownames(lr)
tmp.m <- melt(lr)
tmp.m$SME <- factor(tmp.m$SME, levels = c("Delta","Theta","Alpha","Beta","Gamma","High.Gamma"))

pdf("final_visualizations/plotSME_Distribution.pdf",width=3,height=4,useDingbats=FALSE)
ggplot(tmp.m, aes(x = value, y = SME)) +
geom_density_ridges(aes(fill = SME),scale = 2, alpha = 0.7) +
scale_fill_manual(values = c("orange","grey60","blue", "red","magenta","green"))+
labs(title="",x="", y = "")+
xlim(0,4)+
theme_classic()+
theme(legend.position="none")
dev.off()

# Corr mat
lr <- Lateral_resected
lr <- lr[match(c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta"),rownames(lr)),]

pdf("final_visualizations/SME_CorMat.pdf",width=4,height=4,useDingbats=FALSE)
lr %>% 
    t() %>%
    cor(method="pearson") %>%
    corrplot::corrplot(type="lower",method="pie")
dev.off()


#### BOXPLOT SME Subject
pdf("final_visualizations/plotSME_Boxplot.pdf",width=4,height=4,useDingbats=FALSE)
ggboxplot(tmp.m, x = "SME",
          y = "value",
          combine = TRUE,
          color = "SME", palette = c("green","magenta", "red","blue","grey60","orange"),
          ylab = "Expression",
          add = "jitter",                               
          add.params = list(size = 0.1, jitter = 0.2),  
          )+
ggtitle("SME Distribution")+
xlab("Waves")+ 
ylab("SME")+
theme_classic()+
theme(legend.position="none")+
rotate_x_text(angle = 45)
dev.off()

# PIE CHART LR
df <- as.data.frame(table(memory.genes$Waves))
colnames(df) <- c("Waves","Counts")

# Mutate the table in proportion
df <- df %>%
  arrange(desc(Waves)) %>%
  mutate(prop = round(Counts*100/sum(Counts), 1),
         lab.ypos = cumsum(prop) - 0.5*prop)

pdf("final_visualizations/PIE_LR.pdf",width=4,height=4,useDingbats=FALSE)
ggplot(df, aes(x = "", y = prop, fill = Waves)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  coord_polar("y", start = 0)+
  theme_void()+
scale_fill_manual(values=c("blue", "red","orange","magenta","green","grey60"))+
ggtitle("LR Genes per Waves: P < 0.05")
dev.off()

# N genes LR
df <- as.data.frame(table(memory.genes$Waves))
df$Var1 <- factor(df$Var1,levels = c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta"))

pdf("final_visualizations/Barplot_LR_Ngenes_005.pdf",width=4,height=4,useDingbats=FALSE)
ggbarplot(df, x = "Var1", y = "Freq",
          fill = "Var1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c("green","magenta", "red","blue","grey60","orange"),            
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45,           # Rotate vertically x axis texts
          label = TRUE, label.pos = "out")+
ggtitle("Genes Distribution P 0.05")+
xlab("Waves")+ 
ylab("# of Genes")+
theme(legend.position="none")
dev.off()

# N genes MR
tab <- memory.genes
pos <- tab[tab$Rho > 0,]
neg <- tab[tab$Rho < 0,]
df1 <- as.data.frame(table(pos$Waves))
df2 <- as.data.frame(table(neg$Waves))
df1$Class <- rep("POS",nrow(df1))
df2$Class <- rep("NEG",nrow(df2))

df <- rbind(df1,df2)
df$Var1 <- factor(df$Var1,levels = c("High.Gamma","Gamma","Beta","Alpha","Theta","Delta"))

pdf("final_visualizations/Barplot_LR_Ngenes_005_BothCor.pdf",width=4,height=4,useDingbats=FALSE)
ggbarplot(df, x = "Var1", y = "Freq",
          fill = "Var1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c("green","magenta", "red","blue","grey60","orange"),            
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45,           # Rotate vertically x axis texts
          label = TRUE, label.pos = "out")+
ggtitle("Genes Distribution P 0.05")+
xlab("Waves")+ 
ylab("# of Genes")+
theme(legend.position="none")+
ylim(0,600)+
facet_wrap(~Class)
dev.off()

# Scatter matrix MR
lr <- as.data.frame(t(Lateral_resected))

pdf("final_visualizations/ScatterMatrix_LR.pdf",width=5,height=5,useDingbats=FALSE)
pairs.panels(lr, 
             method = "spearman", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             smoother=FALSE)
dev.off()

# Intersection
pdf("final_visualizations/WavesIntersection_BothCorrleation_LR.pdf",width=6,height=4,useDingbats=FALSE)
l <- split(as.character(memory.genes$Gene),memory.genes$Waves)
Waves <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Waves, ToTGene))
names(metadata) <- c("Waves", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), list(type = "matrix_rows", 
    column = "sets", colors = c(Alpha = "blue", Beta = "red", Theta = "grey",Delta = "orange",Gamma = "purple",High.Gamma = "green"), 
    alpha = 0.5))))
dev.off()

# R squared plot
tmp <- memory.genes
tmp$Rsq <- tmp$Rho^2
tmp$Waves <- factor(tmp$Waves, levels = c("Delta","Theta","Alpha","Beta","Gamma","High.Gamma"))

pdf("final_visualizations/Rsq_Plot.pdf",width=5,height=4,useDingbats=FALSE)
ggline(tmp, 
    x = "Waves", 
    y = "Rsq", 
        add = c("mean_se"),
        color = "Waves", palette = c("orange","grey60","blue", "red","magenta","green"))+
    labs(title="",x="", y = "rho^2")+
    theme_classic()+
    theme(legend.position="none")+
    rotate_x_text(angle = 45)+
      stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y = 0.4)
dev.off()


# Scatter plot candidates
lr <- as.data.frame(Lateral_resected)
lr$SME <- rownames(lr)
tmp.m <- melt(lr)
tmp.m$SME <- factor(tmp.m$SME, levels = c("Delta","Theta","Alpha","Beta","Gamma","High.Gamma"))

IL1RAPL2 <- as.numeric(within_data["IL1RAPL2",])
SMAD3 <- as.numeric(within_data["SMAD3",])

tmp.2 <- tmp.m %>%
          split(tmp.m$SME) %>%
          map(~mutate(., IL1RAPL2 = IL1RAPL2, SMAD3 = SMAD3)) %>%
          bind_rows()

pdf("final_visualizations/IL1RAPL2.pdf",width=5,height=5,useDingbats=FALSE)
ggscatter(tmp.2, 
          x = "IL1RAPL2", 
          y = "value",
          add = "reg.line",                         # Add regression line
          color = "SME", palette = c("orange","grey60","blue", "red","magenta","green"),           # Color by groups "cyl"
          shape = "SME",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = TRUE,
          facet.by="SME"                                # Add marginal rug
          )+
  stat_cor(aes(color = SME),method = "spearman",label.y=3.5)+
    labs(title="",x="Adjusted Expression", y = "SME")+
    theme_classic()+
      theme(legend.position="none")
dev.off()

pdf("final_visualizations/SMAD3.pdf",width=5,height=5,useDingbats=FALSE)
ggscatter(tmp.2, 
          x = "SMAD3", 
          y = "value",
          add = "reg.line",                         # Add regression line
          color = "SME", palette = c("orange","grey60","blue", "red","magenta","green"),           # Color by groups "cyl"
          shape = "SME",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = TRUE,
          facet.by="SME"                                # Add marginal rug
          )+
  stat_cor(aes(color = SME),method = "spearman",label.y=3.5)+
    labs(title="",x="Adjusted Expression", y = "SME")+
    theme_classic()+
      theme(legend.position="none")
dev.off()

