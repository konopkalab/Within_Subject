library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
library(igraph)
library(tidyverse)

exp <- read.table("logCPM_ALL_QuantNorm_LMreg.txt")

mod <- read.table("ModuleOutput_WITHIN.txt",header=T)
#selection <- c("WM4","WM11","WM12","WM19","WM21","WM22")
selection <- "WM4"

mod <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- mod %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("MODULES/WM4_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()

# WM11
selection <- "WM11"
mod <- read.table("ModuleOutput_WITHIN.txt",header=T)

mod <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- mod %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("MODULES/WM11_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()


# WM11
selection <- "WM12"
mod <- read.table("ModuleOutput_WITHIN.txt",header=T)

mod <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- mod %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("MODULES/WM12_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()

# WM19
selection <- "WM19"
mod <- read.table("ModuleOutput_WITHIN.txt",header=T)

mod <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- mod %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("MODULES/WM19_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()

# WM21
selection <- "WM21"
mod <- read.table("ModuleOutput_WITHIN.txt",header=T)

mod <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- mod %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("MODULES/WM21_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()

# WM22
selection <- "WM22"
mod <- read.table("ModuleOutput_WITHIN.txt",header=T)

mod <- mod[mod$ModuleName %in% selection,] %>%
            droplevels()

hubs <- mod %>% 
            group_by(ModuleName) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

df <- exp[rownames(exp) %in% hubs$Gene,]
adjMat <- bicor(t(df))
adjMat[abs(adjMat)<0.5] <- 0
adjMat <- abs(adjMat)
adjMat <- adjMat[match(hubs$Gene,rownames(adjMat)), match(hubs$Gene,colnames(adjMat))]
modColors <- as.matrix(t(hubs$ModuleColor))


valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE)
g1=delete.vertices(g1,which(degree(g1)<1))
#layoutFR <- layout.fruchterman.reingold(g1,dim=2)
layoutFR <- layout_with_lgl(g1,maxiter = 500)



pdf("MODULES/WM22_IGRAPH.pdf",width=3.5,height=3.5)
plot.igraph(g1,
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hubs$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()
