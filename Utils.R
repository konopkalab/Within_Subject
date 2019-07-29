##########################################################
## Functions correlation and brain waves         	##
## FastCor = multicore rapid correlation		##
## VarExp = variance explained by PCAs 			##
## plotVarExp = Plot variance explained  		##
##########################################################


#### Rapid correlation
# exp = gene x sample matrix
# sme = sample x oscillation (one by one as loop)
FastCor <- function(exp,sme,method="pearson",alternative="two.sided",cores=1,override=FALSE,...){
	suppressPackageStartupMessages(library(parallel))
	mat = as.matrix(exp)
	nr = nrow(mat)
	if(nr > 1000){
		if(override){
			show("Wait for it!")
		} else {
			stop("Set override to TRUE")
		}
	}
	if(is.null(rownames(mat))){
		nm = as.character(1:nr)
	} else { 
		nm = rownames(mat)
	}
	corrAndP = mclapply(1:nr,function(i){
			res = cor.test(mat[i,],as.numeric(sme),method=method,alternative=alternative,exact=FALSE)
			cbind(res$estimate,res$p.value)
		},mc.cores=cores,...)
	Rho = unlist(sapply(corrAndP,function(i){i[,1]}))
	Pval  = unlist(sapply(corrAndP,function(i){i[,2]}))
	results = cbind(Rho,Pval)
	rownames(results) <- rownames(exp)
	return(results)
}

#### Variance Explained
# counts = gene x sample matrix
# meta = sample x covariate matrix
# threshold = number of PCA to consider (e.g. 5)
# inter = interaction between covariates (e.g. Age:Sex)

VarExp <- function(counts, meta, threshold, inter){
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(optimx))
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)

  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc <- 3}

  pred.list <- colnames(meta)
  meta <- droplevels(meta)

  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}

  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit,control = lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}

# Bar plot for data visualization
plotVarExp <- function(pvca.res, title){
  suppressPackageStartupMessages(library(ggplot2))
  plot.dat <- data.frame(eff=names(pvca.res), prop=pvca.res)
  p <- ggplot2::ggplot(plot.dat, aes(x=eff, y=prop))
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::geom_bar(stat="identity", fill="steelblue", colour="steelblue")
  p <- p + ggplot2::geom_text(aes(label=round(prop,3), y=prop+0.04), size=4)
  p <- p + ggplot2::scale_x_discrete(limits=names(pvca.res))
  p <- p + ggplot2::scale_y_continuous(limits = c(0,1))
  p <- p + ggplot2::labs(x= "Effects", y= "WAPV")
  p <- p + ggplot2::theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
}

# Fisher Exact multi-factor
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

### For Single Cell

##-------------------------------------------------------
## FUNCTION TO SEURATIFY DATA
##-------------------------------------------------------
seuratify <- function(countData, sampleName)
 {
 seuratObject <- CreateSeuratObject(counts = countData, project = sampleName)
 # save(seuratObject, countData, file = paste(sampleName, "_DATA.RData", sep = ""))
 return(seuratObject)
 }

##-------------------------------------------------------
## FUNCTION TO CREATE BIN HISTOGRAM
##-------------------------------------------------------
plotbinhist <- function(seuratObject)
 {
 brks <- unique(c(seq(0, 1000, by = 100), seq(1000, 25000, by = 1000)))
 # brks <- c(0, 200, 400, 600, 800, 1000, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 100000)
 umiBins <- cut(seuratObject@meta.data$nCount_RNA, brks, include.lowest = T, right = FALSE)

 binData <- cbind(summary(umiBins))
 colnames(binData) <- c("nUMI")
 binData <- as.data.frame(binData)
 binData$bins <- row.names(binData)
 binData2 <- melt(binData)
 binData2$bins <- factor(binData2$bins, levels = row.names(binData))

 plotBins <- ggplot(binData2, aes(bins, value, fill = variable)) +
             geom_bar(stat = "identity", position = "dodge") +
             # facet_grid(rows = vars(variable)) +
             facet_grid(. ~ variable) +
             # coord_flip() +
             theme_bw() +
             theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
             geom_text_repel(aes(label = value), vjust = -1, angle = 90) +
             ylim(0, 50000) +
             labs(title = seuratObject@project.name, x = "Bins", y = "Number of Cells") +
             NULL
 # plotBins
 ggsave(paste(seuratObject@project.name, "_UMI_HIST.pdf", sep = ""), plot = plotBins, width = 16, height = 8, units = "in", dpi = 300)  
 }

##-------------------------------------------------------
## MITO GENES AND QC PLOTS
##-------------------------------------------------------
plotqcmito <- function(seuratObject)
 {
 ## Identify mitochondrial genes and plot percent mito genes along with nGenes and nUMI
 mito.genes <- grep(pattern = "^MT-", x = rownames(x = seuratObject@assays$RNA@data), value = TRUE)
 percent.mito <- Matrix::colSums(seuratObject@assays$RNA@data[mito.genes, ])/Matrix::colSums(seuratObject@assays$RNA@data)
 seuratObject <- AddMetaData(object = seuratObject, metadata = percent.mito, col.name = "pMito")
 
 ## QC Plots before removing mito genes
 plot_nU <- VlnPlot(object = seuratObject, features = "nCount_RNA", pt.size = 0)
 plot_nG <- VlnPlot(object = seuratObject, features = "nFeature_RNA", pt.size = 0)
 plot_pM <- VlnPlot(object = seuratObject, features = "pMito", pt.size = 0)
 plotQC1 <- grid.arrange(plot_nU, plot_nG, plot_pM, ncol = 3)
 ggsave(paste(seuratObject@project.name, "_QC_1.pdf", sep = ""), plot = plotQC1, width = 12, height = 4, units = "in", dpi = 300)

 plot_nUpM <- FeatureScatter(object = seuratObject, feature1 = "nCount_RNA", feature2 = "pMito", cex.use = 1)
 plot_nUnG <- FeatureScatter(object = seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cex.use = 1)
 plotQC2 <- grid.arrange(plot_nUpM, plot_nUnG, ncol = 2)
 ggsave(paste(seuratObject@project.name, "_QC_2.pdf", sep = ""), plot = plotQC2, width = 10, height = 4, units = "in", dpi = 300)

 return(seuratObject)
 }

##-------------------------------------------------------
## FUNCTION TO CREATE HISTOGRAM
##-------------------------------------------------------
plotnumihist <- function(all.data, meta.data, prefix)
 {
 ## Generate histograms for UMI data
 dataHist <- as.data.frame(colSums(all.data))
 colnames(dataHist) <- "UMI_Count"

 hist.data <- merge(dataHist, meta.data, by = "row.names")
 row.names(hist.data) <- hist.data$Row.names
 hist.data$Row.names <- NULL

 ## Plot histograms UMI Distribution (auto-scale)
 plotHist <- ggplot(hist.data, aes(x = log10(UMI_Count), fill = LIBRARY)) + geom_histogram(bins = 100) + facet_grid(. ~ LIBRARY) + labs(title = prefix, x = "Log10(Number of UMIs)", y = "Number of Cells") + theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 90))
 ggsave(paste(prefix, "_HIST.pdf", sep = ""), plot = plotHist, width = 12, height = 4, units = "in", dpi = 150)
 }


