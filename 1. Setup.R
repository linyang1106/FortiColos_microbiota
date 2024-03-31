#setup
library(tidyverse)
library(vegan)
library(BiocManager)
library(phyloseq)
library(devtools)
library(readxl)
library(dplyr)
library(tibble)
library(data.table)
library(tableone)
library(ggplot2)
library(broom)
library(readxl)
library(dplyr)
library(tibble)
library(data.table)
# DESeq2
library(DESeq2)
library(RColorBrewer)
library(ComplexHeatmap)
# correlation
library(Hmisc)
library(ggpubr)
#figure
library(cowplot)

# color
mycolors_for <- c("yellowgreen", "#ffcc1d")
mycolors_BM <- c("#3282b8","#FA7F6F" )
mycolors_hos <- c("#46ACC6", "#1A86A3", "#0E5B76", "#023047","#FFCA5F", "#FEA809", "#FB8402",  "#CD5C08")
mycolors_FT <- c("#D5F0C1","#AAD9BB","#99BC85","gray50")

# theme
# without legend
mytheme_n <- theme_bw()+theme(plot.title = element_text(hjust=0,size=10),
                           legend.position ="none",
                           axis.line = element_line(colour="black"),
                           axis.text = element_text(size = 8, color="black"),
                           axis.title = element_text(size = 10, color = "black"),
                           panel.grid = element_blank(),
                           panel.border = element_blank())
# with legend
mytheme_r <- theme_bw()+theme(plot.title = element_text(hjust=0,size=10),
                              legend.position ="right",
                              axis.line = element_line(colour="black"),
                              axis.text = element_text(size = 8, color="black"),
                              axis.title = element_text(size = 10, color = "black"),
                              panel.grid = element_blank(),
                              panel.border = element_blank())
# abundance plot 
mytheme_abundance <- theme(plot.title = element_text(hjust=0.5, size=12),
                           legend.key.size = unit(0.15,"inches"),
                           legend.position ="right",
                           legend.text = element_text( size = 8,),
                           legend.text.align = 0,
                           legend.background = element_blank(),
                           strip.background = element_blank(),
                           strip.placement = "outside",
                           strip.text = element_text( size = 10),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           axis.text = element_text(size = 10, color="black"),
                           axis.line = element_line(color = "black", size = 0.5),
                           panel.grid = element_blank())
guides(fill=guide_legend(ncol = 1),legend.key.size = unit(0.1, "inches"))
#S1 distribution
themeS1 <- theme(plot.title = element_text(hjust=0.5,  size=10),
      legend.position ="right",
      legend.key.size = unit(0.2,"inches"),
      legend.text = element_text(size = 8),
      legend.background = element_blank(),
      legend.title = element_text(size=8),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.title = element_text(size = 10, color = "black"),
      axis.text = element_text(color = "black", size = 8),
      strip.background = element_blank(),
      strip.text = element_text(size=10))
