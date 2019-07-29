library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(WGCNA)

files=list.files(pattern="STATISTICS")
tmp=as.data.frame(lapply(files,read.table,sep="\t",header=T)[[1]])

tmp <- tmp[c(1,7,8)]
df <- melt(tmp)
df$log <- -log10(df$value)
df$VARIABLE <- factor(df$VARIABLE,levels=paste("WM",1:26,sep=""))
df$Sample <- factor(df$Sample,levels=rev(levels(df$Sample)))

df$log[df$log < 1.3] <- NA


pdf("MAGMA_BUBBLE.pdf",width=8,height=3.5)
ggscatter(df, 
  			x = "Sample",
  			y = "VARIABLE",
   			size="log",
   			color="log",
   			alpha = 0.8,
   			xlab = "",
            ylab = "",) +
   			theme_minimal() + 
   			gradient_color(c("red", "darkred"))+
   			rotate_x_text(angle = 45)+
        coord_flip()
dev.off()


