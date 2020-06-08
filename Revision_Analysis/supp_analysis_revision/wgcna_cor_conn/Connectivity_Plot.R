library(tidyverse)
library(ggpubr)

load("SME_QuantReg.RData")
sme <- SME_Cor_QuantReg_Sign
mod <- read.table("ModuleOutput_WITHIN.txt",header=T)

wm4 <- mod %>% 
		filter(ModuleName == "WM4")

df <- merge(wm4,sme,by="Gene",all=F) %>%
		filter(Waves == "Delta",Rho > 0)


pdf("WM4_Connectivity.pdf",width=6,height=5)
ggscatter(df, 
		   	x = "Rho", 
		   	y = "kWithin",
   			color = "Rho", 
   			size = "kWithin",
   			label = "Gene", 
   			repel = TRUE, 
   			ylab = "Connectivity",
          	xlab = "Spearman's rho") +
          	theme(legend.position="right")
dev.off()

# 
wm12 <- mod %>% 
		filter(ModuleName == "WM12")

df <- merge(wm12,sme,by="Gene",all=F) %>%
		filter(Waves == "Delta",Rho > 0)

pdf("WM12_Connectivity.pdf",width=6,height=5)
ggscatter(df, 
		   	x = "Rho", 
		   	y = "kWithin",
   			color = "Rho", 
   			size = "kWithin",
   			label = "Gene", 
   			repel = TRUE, 
   			ylab = "Connectivity",
          	xlab = "Spearman's rho") +
          	theme(legend.position="right")
dev.off()


wm21 <- mod %>% 
		filter(ModuleName == "WM21")

df <- merge(wm21,sme,by="Gene",all=F) %>%
		filter(Waves == "Delta",Rho < 0, kWithin > 70) %>%
		mutate(Abs = abs(Rho))

pdf("WM21_Connectivity.pdf",width=6,height=5)
ggscatter(df, 
		   	x = "Abs", 
		   	y = "kWithin",
   			color = "Abs", 
   			size = "kWithin",
   			label = "Gene", 
   			repel = TRUE, 
   			ylab = "Connectivity",
          	xlab = "Spearman's rho") +
          	theme(legend.position="right")
dev.off()



