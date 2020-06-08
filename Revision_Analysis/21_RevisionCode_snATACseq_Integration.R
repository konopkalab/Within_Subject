# Single Cell Analysis
rm(list = ls())
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(wTO))
suppressPackageStartupMessages(library(unikn))
source("UTILS/Utils.R")

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 5000 * 1024^2)

# UMAP Alternative
df <- readRDS("https://www.dropbox.com/s/cgkqacfekdce2z0/sme_atac_meta_arranged.RDS?dl=1")
df$pct_reads_in_peaks <- df$peak_region_fragments / df$passed_filters * 100
df$blacklist_ratio <- df$blacklist_region_fragments / df$peak_region_fragments

df <- df %>%
      mutate(Class = case_when(grepl("Exc", annotation) ~ "Glutamatergic", 
                        grepl("Inh", annotation) ~ "Gabaergic",
                        grepl("Astro|Olig|OPC", annotation) ~ "NonNeuronal"))


label <- data.frame(annotation=levels(as.factor(df$annotation)),label=levels(as.factor(df$annotation)))

label_2 <- df %>% 
  group_by(annotation) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  left_join(label) %>%
  as.data.frame() %>%
  mutate(Definition = case_when(grepl("Exc", annotation) ~ "Excitatory", 
                       grepl("Inh", annotation) ~ "Inhibitory",
                       grepl("Astro", annotation) ~ "Astrocytes", 
                       grepl("Olig", annotation) ~ "Oligodendrocytes",
                       grepl("OPC", annotation) ~ "OPC"))


colors <- df %>% 
			group_by(annotation) %>%
			summarize(N = n()) %>%
			arrange(desc(N)) %>%
			mutate(Definition = case_when(grepl("Exc", annotation) ~ "Excitatory", 
                       grepl("Inh", annotation) ~ "Inhibitory",
                       grepl("Astro", annotation) ~ "Astrocytes", 
                       grepl("Olig", annotation) ~ "Oligodendrocytes",
                       grepl("OPC", annotation) ~ "OPC")) %>%
            mutate(Class = case_when(grepl("Exc", annotation) ~ "Glutamatergic", 
                        grepl("Inh", annotation) ~ "Gabaergic",
                        grepl("Astro|Olig|OPC", annotation) ~ "NonNeuronal"))	%>%		
  	as.data.frame()  

l <- split(colors,colors$Class)

#l[[1]]$Color <- colorRampPalette(rev(brewer.pal(4,"pal_peach")))(9)
#l[[2]]$Color <- colorRampPalette(rev(brewer.pal(4,"pal_karpfenblau")))(8)
#l[[3]]$Color <- colorRampPalette(rev(brewer.pal(3,"pal_pinky")))(3)

l[[1]]$Color <- rev(usecol(pal_bordeaux, 4))
l[[2]]$Color <- rev(usecol(pal_karpfenblau, 7))
l[[3]]$Color <- rev(usecol(pal_seegruen, 3))


colors <- do.call(rbind,l)
colors <- colors[match(label_2$label, colors$annotation),]


pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP_Alternative_scATAC.pdf", width = 5, height = 5)
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] 
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = annotation),size=0.5) +
ggrepel::geom_text_repel(data = label_2, aes(label = label),
							color = "black",
							#fontface = 'bold',
							segment.colour = "grey60",
						    box.padding = unit(0.25, "lines"),
						    point.padding = unit(0.5, "lines"),
						    nudge_x = .15,
						    nudge_y = 1,
						    size = 2.5) + 
    theme_classic()+
    theme(legend.position="none") +
		#scale_color_viridis(discrete=TRUE,option="inferno")
		#scale_colour_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(31))
		scale_colour_manual(values = colors$Color)
dev.off()


pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP_BySubject_Alternative_scATAC.pdf", width = 5, height = 3)
ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=orig.ident)) +
ggrastr::geom_point_rast(size=0.5) +
facet_wrap(.~orig.ident,ncol=3) +
theme_classic() +
theme(legend.position="none")
#scale_colour_manual(values = c("#E7B800", "#FC4E07"))
dev.off()

# UMAP Alternative by Cell Class
pdf("integration_scRNAseq/PLOTS/SME_Data_Integrated_UMAP_ByCellClass_Alternative_scATAC.pdf", width = 6, height = 3)
ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Class)) +
ggrastr::geom_point_rast(size=0.5) +
facet_wrap(.~Class) +
theme_classic() +
theme(legend.position="none")+ 
scale_colour_manual(values = c("blue","red","green"))
dev.off()

pie <- table(df$Class) %>% 
     as.data.frame()  %>% 
     arrange(desc(Freq)) %>% 
     mutate(percent = scales::percent(Freq/sum(Freq))) %>%
     ggplot(aes(x = "", y = Freq, fill = fct_inorder(Var1))) +
       geom_bar(width = 1, stat = "identity") +
       coord_polar("y", start = 0) +
       geom_label_repel(aes(label = percent), size=5, show.legend = F, nudge_x = 1) +
       guides(fill = guide_legend(title = "Cell Class"))+
       theme_classic()
ggsave("integration_scRNAseq/PLOTS/SME_Integration_prop_cell_pie_scATAC.pdf", plot = pie, width = 6, height = 4, units = "in", dpi = 150)

ngenes <- df %>% 
    ggboxplot( 
      x = "orig.ident", 
      y = "nFeature_peaks",
      color = "orig.ident",
      palette ="Paired") + 
      theme_classic() +
    rotate_x_text(angle = 45) + 
    xlab("") + 
    ylab("# of Genes")+
    theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_nGENES_scATAC.pdf", plot = ngenes, width = 4, height = 4, units = "in", dpi = 150)


numis <- df %>% 
    ggboxplot( 
      x = "orig.ident", 
      y = "nCount_peaks",
      color = "orig.ident",
      palette ="Paired") + 
      theme_classic() +
    rotate_x_text(angle = 45) + 
    xlab("") + 
    ylab("# of UMI")+
    theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_nUMI_scATAC.pdf", plot = numis, width = 4, height = 4, units = "in", dpi = 150)


pct_reads <- df %>% 
    ggboxplot( 
      x = "orig.ident", 
      y = "pct_reads_in_peaks",
      color = "orig.ident",
      palette ="Paired") + 
      theme_classic() +
    rotate_x_text(angle = 45) + 
    xlab("") + 
    ylab("% Reads in Peaks")+
    theme(legend.position="none")+
    ylim(0,100)
ggsave("integration_scRNAseq/PLOTS/SME_Integration_PctReadsPeaks_scATAC.pdf", plot = pct_reads, width = 3, height = 4, units = "in", dpi = 150)


peak_region <- df %>% 
    ggboxplot( 
      x = "orig.ident", 
      y = "peak_region_fragments",
      color = "orig.ident",
      palette ="Paired") + 
      theme_classic() +
    rotate_x_text(angle = 45) + 
    xlab("") + 
    ylab("Total reads in Peaks")+
    theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_TotalReadsPeaks_scATAC.pdf", plot = peak_region, width = 3, height = 4, units = "in", dpi = 150)


nCount_peaks <- df %>% 
    ggboxplot( 
      x = "orig.ident", 
      y = "nCount_peaks",
      color = "orig.ident",
      palette ="Paired") + 
      theme_classic() +
    rotate_x_text(angle = 45) + 
    xlab("") + 
    ylab("# of Peaks")+
    theme(legend.position="none")
ggsave("integration_scRNAseq/PLOTS/SME_Integration_NumberOfPeaks_scATAC.pdf", plot = nCount_peaks, width = 3, height = 4, units = "in", dpi = 150)

pdf("integration_scRNAseq/PLOTS/SME_Integrated_PCTvsPeaks_scATAC.pdf", width=5,height=4)
ggplot(df, aes(x=nCount_peaks, y=pct_reads_in_peaks, color=orig.ident)) +
ggrastr::geom_point_rast(size=0.5) + 
theme_classic()+
ylim(0,100) + 
scale_colour_brewer("Colors in Paired", palette="Paired")
dev.off()





