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
# 
library(DESeq2)
library(RColorBrewer)
library(ComplexHeatmap)

# color and theme
mycolors_for <- c("yellowgreen", "#ffcc1d")
mycolors_BM <- c()
mycolors_hos <- c()
mycolors_pro <- c()
mytheme<- theme_bw()+theme(plot.title = element_text(hjust=0,size=10),
                           legend.position ="none",
                           axis.line = element_line(colour="black"),
                           axis.text = element_text(size = 8, color="black"),
                           axis.title = element_text(size = 10, color = "black"),
                           panel.grid = element_blank(),
                           panel.border = element_blank())
