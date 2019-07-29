library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
enableWGCNAThreads()

#cool <- read.table("coolmod.txt")
p <- read.table("modTraitP_WITHIN.txt")
rownames(p) <- gsub("ME","",rownames(p))
#p <- as.matrix(p[rownames(p) %in% cool$V1,])
cor <- read.table("modTraitCor_WITHIN.txt")
rownames(cor) <- gsub("ME","",rownames(cor))
#cor <- as.matrix(cor[rownames(cor) %in% cool$V1,])

rownames(p) <- paste("W","M",1:nrow(p),sep="")
rownames(cor) <- paste("W","M",1:nrow(cor),sep="")
p$Module <- rownames(p)
cor$Module <- rownames(cor)

p <- melt(p)
cor <- melt(cor)
p$cor <- cor$value
p$log <- -log10(p$value)
p$Direction <- ifelse(p$cor > 0,"Pos","Neg")
p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))
p$variable <- gsub("_LR","",p$variable)

p$log[p$log < 1.3] <- NA
p$abs <- abs(p$cor)
p$abs[p$abs < 0.5] <- NA

pdf("Module_Traits_Bubble.pdf",width=8,height=1.8)
ggscatter(p, 
  			x = "variable",
  			y = "Module",
   			size="abs",
   			color="Direction",
   			palette=c("blue","red"),
   			alpha = 0.8,
   			xlab = "",
            ylab = "",) +
   			theme_minimal() +
   			rotate_x_text(angle = 45)+
        coord_flip()
dev.off()

