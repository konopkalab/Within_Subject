#!/usr/bin/env Rscript
#Check if the optparse library is installed, required to parse arguments
`%notin%` <- function(x,y) !(x %in% y) 

if(suppressPackageStartupMessages(require("optparse"))){
	print("optparse is loaded correctly");
} else {
	print ("trying to install optparse")
	install.packages("optparse",repos="http://cran.rstudio.com/");
	if(suppressPackageStartupMessages(require("optparse"))){
		print ("optparse installed and loaded");
	} else{
		stop("could not install optparse");
	}
}

if(suppressPackageStartupMessages(require("tidyverse"))){
	print("tidyverse is loaded correctly");
} else {
	print ("trying to install tidyverse")
	install.packages("tidyverse",repos="http://cran.rstudio.com/");
	if(suppressPackageStartupMessages(require("tidyverse"))){
		print ("tidyverse installed and loaded");
	} else{
		stop("could not install tidyverse");
	}
}

if(suppressPackageStartupMessages(require("reshape2"))){
	print("reshape2 is loaded correctly");
} else {
	print ("trying to install reshape2")
	install.packages("reshape2",repos="http://cran.rstudio.com/");
	if(suppressPackageStartupMessages(require("reshape2"))){
		print ("reshape2 installed and loaded");
	} else{
		stop("could not install reshape2");
	}
}

#Perform the argument checking before spending time in loading DESeq2
option_list <- list(
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
	make_option(c("-q", "--quietly"), action="store_false",	dest="verbose", help="Print little output"),
	make_option(c("-g", "--genes"), type="character",  help="The input. List of objects. Each object should be of the format <Gene> <Association>", dest="inputFile"),
	make_option(c("-l", "--list"), type="character",  help="The input for the overlap. Should be of the format <Gene> <Association>  ", dest="inputDisease"),
	make_option(c("-p", "--plot"), action="store_true", default=FALSE, help="Generate the bubble chart for enrichment"),
	make_option(c("-b", "--bkg"), type="integer", default=15585, help="Number of genes as background", metavar="number"),
	make_option(c("-o", "--out"), type="character", help="Output prefix for the Rdata of enrichment stats", metavar="output"),
	make_option(c("-d", "--dir"), type="character", help="Directory where to save", metavar="directory"),
	make_option(c("-W", "--width"), type="integer", default=4, help="Width for pdf", metavar="number"),
	make_option(c("-H", "--height"), type="integer", default=5, help="Height for pdf", metavar="number")

)
 
opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$inputFile)){
	stop("Require input!");
}

if(is.null(opt$inputDisease)){
	stop("Require Input for the enrichment!");
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

if(file.access(opt$inputFile) ==-1){
	stop(sprintf("Input does not exist", opt$inputFile));
} else{
	inputTable = read.table(opt$inputFile,sep="\t",header=T);
	inputTable <- inputTable %>%
			        select(Gene,ModuleName) %>%
			        rename(DEFINITION = ModuleName)
	inputTable$DEFINITION <- factor(inputTable$DEFINITION,levels=paste("WM",1:length(unique(inputTable$DEFINITION)),sep=""))
	Genes=as.data.frame(table(inputTable$DEFINITION))
}
print("Input obtained")

if(file.access(opt$inputDisease) ==-1){
	stop(sprintf("Input Disease does not exist", opt$inputDisease));
} else{
	inputDisease = load(opt$inputDisease) %>%
            get()
	GeneSets <- inputDisease
}
print("Input obtained, GeneSet generate")

#print(head(inputTable));
#print(sampleInfo);

ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=inputTable[inputTable$Gene %in% GeneSets[[i]]$Gene,]
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
INT[[i]]$d <- opt$bkg-Genes$Freq-nrow(GeneSets[[i]]) #Our = 15193 #BrainSpan = 15585
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
save(FisherMat,file=paste(opt$out, "_Stats.RData", sep=""))

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

#library(reshape2)
p <- FisherP %>% 
		as.data.frame() %>%
		rownames_to_column("Module") %>%
		melt() %>%
		group_by(Module) %>%
		mutate(FDR = p.adjust(value,method="BH")) %>%
		mutate(log = -log10(FDR)) %>%
		rename(Pval = value) %>%
		as.data.frame()
		
OR <- FisherOR %>%
		as.data.frame() %>%
		rownames_to_column("Module") %>%
		melt()%>%
		rename(OR = value)

p <- Reduce(dplyr::left_join, list(p, OR))
p$OR[!is.finite(p$OR)] <- max(p$OR[is.finite(p$OR)])
p$log[p$log < 1.3] <- NA
p$OR<- ifelse(is.na(p$log), p$log, p$OR)
p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))

#Graph plotting related codes
if(opt$plot){
	#Check the required packages
	if(suppressPackageStartupMessages(require(ggpubr))){
		print("ggpubr is loaded correctly");
	} else {
		print ("trying to install ggpubr")
		install.packages("ggpubr",repos="http://cran.rstudio.com/");
		if(suppressPackageStartupMessages(require("ggpubr"))){
			print ("ggpubr installed and loaded");
		} else{
			print("could not install ggpubr, will not plot the sample clustering plot");
			opt$plot = FALSE;
		}
	}
	if(suppressPackageStartupMessages(require(RColorBrewer))){
		print("RColorBrewer is loaded correctly");
	} else {
		print ("trying to install RColorBrewer")
		install.packages("RColorBrewer",repos="http://cran.rstudio.com/");
		if(suppressPackageStartupMessages(require("RColorBrewer"))){
			print ("RColorBrewer installed and loaded");
		} else if(opt$plot){
			print("could not install RColorBrewer, will not plot the sample clustering plot");
			opt$plot=FALSE;
		} else{
			print("could not install RColorBrewer, it is also required to plot the sample clustering plot");
		}
	}
	if(opt$plot){
		suppressPackageStartupMessages(require(ggpubr))
		#pdf(paste(opt$out, "_Enrichment.pdf", sep=""))
		ggscatter(p, 
        	x = "variable",
        	y = "Module",
        	size="OR",
        	color="log",
        	alpha = 0.8,
        	xlab = "",
            ylab = "") +
        		theme_minimal()+
        		rotate_x_text(angle = 45)+
        		coord_flip()+
      			scale_size(range = c(2, 10))+ 
      			gradient_color(c("lightblue","darkblue"))+
        #geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=3.5,ymax=4.5),colour="cyan",fill=NA,alpha=0.5) +
        #geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=10.5,ymax=11.5),colour="greenyellow",fill=NA,alpha=0.5) +
        #geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=11.5,ymax=12.5),colour="grey60",fill=NA,alpha=0.5) +
        #geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=18.5,ymax=19.5),colour="pink",fill=NA,alpha=0.5) +
        #geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=20.5,ymax=21.5),colour="red",fill=NA,alpha=0.5) +
        #geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=21.5,ymax=22.5),colour="steelblue",fill=NA,alpha=0.5)
        annotate("rect", xmin=-Inf,xmax=Inf,ymin=3.5,ymax=4.5,fill="cyan",alpha=0.2) +
        annotate("rect", xmin=-Inf,xmax=Inf,ymin=10.5,ymax=11.5,fill="greenyellow",alpha=0.2) +
        annotate("rect", xmin=-Inf,xmax=Inf,ymin=11.5,ymax=12.5,fill="grey60",alpha=0.2) +
        annotate("rect", xmin=-Inf,xmax=Inf,ymin=18.5,ymax=19.5,fill="pink",alpha=0.2) +
        annotate("rect", xmin=-Inf,xmax=Inf,ymin=20.5,ymax=21.5,fill="red",alpha=0.2) +
        annotate("rect", xmin=-Inf,xmax=Inf,ymin=21.5,ymax=22.5,fill="steelblue",alpha=0.2)

		ggsave(paste(opt$out, "_Enrichment.pdf", sep=""),width=opt$width,height=opt$height)		
	}	
}


