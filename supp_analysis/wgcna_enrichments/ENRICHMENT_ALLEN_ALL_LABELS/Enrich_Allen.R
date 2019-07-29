library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
enableWGCNAThreads()

# Brain Exp
rm(list=ls())
wgcna = list.files(pattern = 'ModuleOutput')
tab=read.table(wgcna,sep="\t",header=T)
tab <- tab[c(1,4)]
colnames(tab)=c("Gene","DEFINITION")
tab$DEFINITION <- factor(tab$DEFINITION,levels=paste("WM",1:length(unique(tab$DEFINITION)),sep=""))

Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
load("Allen_New_Clusters.RData")
#for(i in 1:length(GeneSets))
#{
#	GeneSets[[i]] <- GeneSets[[i]][GeneSets[[i]]$Gene %in% tab$Gene,]
#}
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- 14373-Genes$Freq-nrow(GeneSets[[i]]) #Protein Coding bkg = 19901 Ortho bkg = 14373 Brain Expressed bkg = 15585
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file="AllenFull_Hsap.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,7]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='bonferroni') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

p <- as.data.frame(FisherAdj)
cor <- as.data.frame(FisherOR)


p$Module <- rownames(p)
cor$Module <- rownames(cor)

p <- melt(p)
cor <- melt(cor)
p$cor <- cor$value
p$log <- -log10(p$value)
p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))

p$log[p$log < 1.3] <- NA
p$abs <- abs(p$cor)
p$OR<- ifelse(is.na(p$log), p$log, p$abs)


#p$variable <- factor(p$variable,levels=c("Endo_L2-6","Micro_L1-3","Astro_L1-6","Astro_L1-2","Oligo_L1-6","OPC_L1-6","Exc_L6","Exc_L5-6","Exc_L4-6","Exc_L4-5","Exc_L3-5","Exc_L3-4","Exc_L2-4","Exc_L2-3","Exc_L2",
#"Inh_L5-6","Inh_L4-6","Inh_L4-5","Inh_L3-6","Inh_L3-5","Inh_L2-6","Inh_L2-5","Inh_L2-4","Inh_L2-3","Inh_L1-4","Inh_L1-3","Inh_L1-2","Inh_L1"))

pdf("Modules_Allen_Bubble.pdf",width=8,height=10)
ggscatter(p, 
  			x = "variable",
  			y = "Module",
   			size="OR",
   			color="log",
   			alpha = 0.8,
   			xlab = "",
            ylab = "",) +
   			theme_minimal() +
   			rotate_x_text(angle = 45)+
        coord_flip()+
  		scale_size(range = c(2, 10))+ 
  		gradient_color(c("orange","darkorange"))
dev.off()





