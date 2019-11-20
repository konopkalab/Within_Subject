
#' Load Libraries
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(DMwR))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(d3heatmap))
suppressPackageStartupMessages(library(xlsx))


dir.create("Dge_Results")
# Data Exp
exp <- read.table("exp_data.txt")
exp <- exp[,grep("Lega|Tamm",names(exp))]
#exp$EpBB04<-NULL
filter=apply(exp, 1, function(x) (all(x[grep("Lega",names(x))] > 0) | all(x[grep("Tamm",names(x))] > 0)))
exp <- exp[filter,]

pd = read.table("Demographic.txt",header=T,sep="\t",row.names=1)
pd <- droplevels(pd[rownames(pd)%in%colnames(exp),])


#' Mutate covariate into factorial
pd$Diagnosis = as.factor(pd$Diagnosis)
pd$Sex = as.factor(pd$Sex)
pd$Batch = as.factor(pd$Batch)
pd$Hemisphere = as.factor(pd$Hemisphere)


# 
Data=t(exp)
outputSV <- data.frame(matrix(nrow=ncol(Data), ncol=3, dimnames = list(colnames(Data), c("logFC", "Pval", "Warning"))))
outputSV[,] <- NA

#" This is the formula I am using
print(formula(paste("Data[,i] ~ ", paste(colnames(pd),collapse = " + "))))

#' Looping on the Data
for (i in 1:ncol(Data)) {	
			Model=tryCatch(lm(as.formula(paste("Data[,i] ~ ", paste(colnames(pd),collapse = " + "))), data = pd),warning =  function(w) w)
				if (i %% 1000 == 0) {cat(paste("Done on gene ",i,"\n",sep=""))}
				if(typeof(Model) == "list"){
    				coefs = data.frame(coef(summary(Model)))
    				t_value = coefs["Diagnosis", "t.value"]
    				outputSV[i,"Pval"] = 2 * (1 - pnorm(abs(t_value)))
    				outputSV[i,"logFC"]= coefs["Diagnosis", "Estimate"]
    				} else {
    				outputSV[i,"Warning"] = as.character(Model)
    				outputSV[i, "logFC"] = 0
    				outputSV[i,"Pval"] = 1
  			}
		}
outputSV$FDR <- p.adjust(outputSV$Pval, method="BH")
outputSV$Warning <- NULL
signSV <- outputSV[outputSV$FDR< 0.05,]
write.table(outputSV, "Dge_Results/Dge_WS_FrozBrainHealthy.txt",sep="\t",quote=F)
openxlsx::write.xlsx(outputSV, file = "Dge_Results/Dge_WS_FrozBrainHealthy.xlsx", colNames = TRUE, borders = "columns",sheetName="WS_Tamm",row.names=T)


# For Shiny

p <- exp
pd_sva=pd[c(-1)]
pd_sva=droplevels(pd_sva)
avebeta.lm<-lapply(1:nrow(p), function(x){
  lm(unlist(p[x,])~.,data=pd_sva)
})
residualsSV<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residualsSV<-do.call(rbind, residualsSV)
adj.residualsSV<-residualsSV+matrix(apply(p, 1, mean), nrow=nrow(residualsSV), ncol=ncol(residualsSV))
rownames(adj.residualsSV)<-rownames(p)
rownames(residualsSV) <- rownames(p)
write.table(adj.residualsSV, "Shiny/expAdj_WS_FrozBrainHealthy.txt",sep="\t",quote=F)

pd$Type <- ifelse(pd$Diagnosis == "Lega","Type1","Type2")
write.table(pd, "Shiny/WS_FrozBrainHealthy_pheno.txt",sep="\t",quote=F)





