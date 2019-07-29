library(ggpubr)
library(tidyverse)

files = list.files(pattern = '*.txt')
names <- gsub( "MF_|.txt", "", files )
GeneSets = lapply(files, read.table,header=T,sep="\t")
names(GeneSets) <- names
GeneSets <- GeneSets[c("WM4","WM11","WM12","WM19","WM21","WM22")]


filt <- vector("list", length = length(GeneSets))
names(filt) <- names(GeneSets)
class <- names(GeneSets)
for (i in 1:length(GeneSets))
    {
      filt[[i]] <- GeneSets[[i]] %>% 
      				filter(Count > 5 & 
       				Pvalue < 0.05) %>%
        			select(Term,Pvalue,OddsRatio) %>% 
        			mutate(log = -log10(Pvalue)) %>%
        			as.data.frame()
	}

for (i in 1:length(GeneSets))
    {
	filt[[i]]$Class <- factor(rep(class[i],nrow(filt[[i]])))
	}

df <- do.call(rbind,filt)


top_labelled <- tbl_df(df) %>% 
                  group_by(Class) %>% 
                  top_n(n = 3, wt = abs(log))


pdf("GO_enrichment.pdf",width=10,height=2,useDingbats=FALSE)
ggscatter(top_labelled, x = "log", y = "OddsRatio",
   color = "Class", palette = c("cyan", "greenyellow", "grey60","pink","red","royalblue"),size = 2,
   label = "Term", repel = TRUE,font.label = c(8, "plain"))+
xlab("-log10(p-value)")+ 
ylab("Odds Ratio") +
facet_wrap(~Class,ncol=6,nrow=2,scales="free")+
theme_classic()+
theme(legend.position="none")

dev.off()
