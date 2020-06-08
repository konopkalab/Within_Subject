### Stefano Berto, 10/2018 - 
### Code for SME - gene expression correlation analysis
### Perform expression filtering, regression, correlation and cross-validations]

# Load libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(fastpos))
suppressPackageStartupMessages(library(tidyverse))
source("UTILS/Utils.R")

# Create Output directories
dir.create("simulation")

# Load inputs
load(here("rawdata","expAdj.RData"))
load(here("rawdata","SME_Values.RData"))

mc_sim <- function(correlation, sample, replication, confidence) {
      r <- abs(correlation)
      n <- sample
      conf <- confidence
      rep <- replication

      z <- .5*log((1+r)/(1-r))
      se <- 1/sqrt(n-2)
      zvec <- rnorm(rep)*se+z
      ivec <- (exp(2*zvec)-1)/(exp(2*zvec)+1)
      low <- (1-conf/100)/2
      upp <- ((1-conf/100)/2)+(conf/100)
      LowQuant <- quantile(ivec,low)
      UpperQuant <- quantile(ivec,upp)
      pLQ <- 2*pt(-abs(LowQuant/(sqrt(1-LowQuant^2)/(n-2))),df=n-2)
      pUQ <- 2*pt(-abs(UpperQuant/(sqrt(1-UpperQuant^2)/(n-2))),df=n-2)
      out = list(Stata = data.frame(LQ = LowQuant, UQ = UpperQuant, Pval_LQ = pLQ, Pval_UQ=pUQ), histo = ivec)

    return(out)
}

df12 <- mc_sim(0.59,12,1000,95)
df16 <- mc_sim(0.59,16,1000,95)
df30 <- mc_sim(0.59,30,1000,95)
df50 <- mc_sim(0.59,50,1000,95)
df100 <- mc_sim(0.59,100,1000,95)
df500 <- mc_sim(0.59,500,1000,95)
df1000 <- mc_sim(0.59,1000,1000,95)

data <- data.frame(Sample = as.factor(c(rep("Size12",1000),
                              rep("Size16",1000),
                              rep("Size30",1000),
                              rep("Size50",1000),
                              rep("Size100",1000),
                              rep("Size500",1000),
                              rep("Size1000",1000))),
                   Simulation = c(df12$histo,
                                  df16$histo,
                                  df30$histo,
                                  df50$histo,
                                  df100$histo,
                                  df500$histo,
                                  df1000$histo))

data$Sample <- factor(data$Sample,levels = c("Size12","Size16","Size30","Size50","Size100","Size500","Size1000"))


pdf("simulation/histogram_MonteCarlo.pdf",width=5,height=3)
a <- gghistogram(data, 
              x = "Simulation",
              add = "mean", 
              rug = FALSE,
              color = "Sample",
              #fill = "Sample",
              bins=50, 
              palette = "Paired")+
            theme(legend.position="right")+
            rotate_x_text(angle = 45)+
            geom_vline(xintercept = 0.59, linetype="dotted", color = "black", size=1)
print(a)
dev.off()

# Pvals
log <- c(-log10(df12$Stata$Pval_LQ),
  -log10(df16$Stata$Pval_LQ),
  -log10(df30$Stata$Pval_LQ),
  -log10(df50$Stata$Pval_LQ))

p <- c(df12$Stata$Pval_LQ,
        df16$Stata$Pval_LQ,
        df30$Stata$Pval_LQ,
        df50$Stata$Pval_LQ)

n <- c(12,16,30,50)

data <- data.frame(Size = n, Pval=signif(p,3),log=log)
pdf("simulation/Lines_MonteCarlo_pval.pdf",width=3,height=3)
b <- ggscatter(data, 
          x = "Size", 
          y = "log",
          shape = 21, 
          size = 3,
          label = "Pval", repel = TRUE)+
          theme(legend.position="none")+
          rotate_x_text(angle = 45)+
          geom_hline(yintercept = 1.3, linetype="dotted", color = "red", size=1)+
          geom_vline(xintercept = 16, linetype="dotted", color = "steelblue", size=1) +
          xlab("Sample Size") + 
          ylab("-log10(P-value)")
print(b)
dev.off()

plot2by2 <- cowplot::plot_grid(a,b,labels=c("A","B"), ncol = 2)
cowplot::save_plot("simulation/Simulations_stats.pdf", plot2by2, ncol = 2,base_height=3,base_width=4)

# Calculate Point Of Stability for each of the gene
load(here("processing_memory","SME_Genes.RData"))

df <- find_critical_pos(rho = SME_cor_LR$Rho,
          precision = 0.2, 
          sample_size_min = 12,
          sample_size_max = 200,
          confidence_levels = 0.95,
          n_studies = 1000)

colnames(df)[2] <- "POS"

write.table(df,"simulation/SME_PointOfStability.txt",sep="\t",quote=F)

openxlsx::write.xlsx(df, 
                     file = "simulation/SME_PointOfStability.xlsx", 
                     colNames = TRUE, 
                     borders = "columns",
                     sheetName="POS")

pdf("simulation/SME_PointOfStability.pdf",width=5,height=3.5)
gghistogram(df, 
            x = "POS",
            add = "median", 
            rug = FALSE,
            color = "navy",
            fill="steelblue",
            bins = 50,
            add_density = TRUE)
dev.off()

df$ABS <- abs(df$rho_pop)

pdf("simulation/SME_Genes_PointOfStability.pdf",width=5,height=3.5)
ggscatter(df, 
            x = "POS",
            add = "median", 
            rug = FALSE,
            color = "navy",
            fill="steelblue",
            bins = 50,
            add_density = TRUE)
dev.off()

# Combine data
load("processing_memory/SME_Genes.RData")
pos <- read.table("simulation/SME_PointOfStability.txt")

SME_cor_LR_PoS <- cbind(SME_cor_LR,pos)

write.table(SME_cor_LR_PoS,"simulation/SME_cor_LR_PoS.txt",sep="\t",quote=F)
save(SME_cor_LR_PoS, file ="simulation/SME_Genes_PointOfStability.RData")

sign <- SME_cor_LR_PoS %>% 
          filter(Pval < 0.05, PermP < 0.05, BootP_All < 0.05,BootP < 0.05)

openxlsx::write.xlsx(sign, 
                     file = "simulation/SME_Significant_PointOfStability.xlsx", 
                     colNames = TRUE, 
                     borders = "columns",
                     sheetName="POS")

