setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
#################################################################################
###################################################### Alpha diversity, data for supplemantary figure 3a, 3b ##############################################
rich <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon"))
mapping <- data.frame(sample_data(ps_rarefied))
index <- match(rownames(rich), rownames(mapping))
rich_full <- cbind.data.frame(rich, mapping[index, ])
rich_0 <- subset(rich_full, FT == "0")
rich_1 <- subset(rich_full, FT == "1")
rich_2 <- subset(rich_full, FT == "2")

# Shannon index
lm_shannon0 <- lm(Shannon~ fortification + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, rich_FT0)
shannon0 <- as.data.frame(tidy(lm_shannon0))
lm_shannon1 <- lm(Shannon~ fortification + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, rich_FT1)
shannon1 <- as.data.frame(tidy(lm_shannon1))
lm_shannon2 <- lm(Shannon~fortification + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, rich_FT2)
shannon2 <- as.data.frame(tidy(lm_shannon2))
BM <- data.frame(group1 = c("C","C","C"), group2=c("V", "V", "V"))
shannon_BM <- as.data.frame(rbind(shannon0[2,],shannon1[2,],shannon2[2,]))
shannon_BM <- shannon_BM[,-1]
row.names(shannon_BM) <- c("FT0","FT1","FT2")
shannon_BM <- cbind(BM,shannon_BM)
print(shannon_BM, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

# Observed zOTUs
lm_observed0 <- lm(Observed~fortification + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, rich_0)
observed0 <- as.data.frame(tidy(lm_observed0))
lm_observed1 <- lm(Observed~fortification  + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, rich_1)
observed1 <- as.data.frame(tidy(lm_observed1))
lm_observed2 <- lm(Observed~ fortification  + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, rich_2)
observed2 <- as.data.frame(tidy(lm_observed2))
observed_BM <- as.data.frame(rbind(observed0[2,],observed1[2,],observed2[2,]))
observed_BM <- observed_BM[,-1]
row.names(observed_BM) <- c("FT0","FT1","FT2")
observed_BM <- cbind(BM,observed_BM)
print(observed_BM, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

############################################################## Beta diversity, data for figure 2a-c #########################################################
# adonis 
adonis <- function(ps, method_c) {
  set.seed(123)
  dist <- phyloseq::distance(ps, method = method_c)
  metadata <- data.frame(sample_data(ps))
  re_adonis <- adonis2(dist ~ fortification + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, data = metadata, permutations = 999, by="margin")
}
adonis0w <-adonis(ps_FT0, "Weighted Unifrac")
adonis0u <-adonis(ps_FT0, "Unweighted Unifrac")
adonis1w <-adonis(ps_FT1, "Weighted Unifrac")
adonis1u <-adonis(ps_FT1, "Unweighted Unifrac")
adonis2w <-adonis(ps_FT2, "Weighted Unifrac")
adonis2u <-adonis(ps_FT2, "Unweighted Unifrac")
# extract results of BM
weighted_pBM <- as.data.frame(rbind(adonis0w$`Pr(>F)`[2],adonis1w$`Pr(>F)`[2],adonis2w$`Pr(>F)`[2]))
weighted_rBM <- as.data.frame(rbind(adonis0w$R2[2],adonis1w$R2[2],adonis2w$R2[2]))
unweighted_pBM <- as.data.frame(rbind(adonis0u$`Pr(>F)`[2],adonis1u$`Pr(>F)`[2],adonis2u$`Pr(>F)`[2]))
unweighted_rBM <- as.data.frame(rbind(adonis0u$R2[2],adonis1u$R2[2],adonis2u$R2[2]))
adonis_BM <- cbind(weighted_pBM,weighted_rBM,unweighted_pBM,unweighted_rBM)
row.names(adonis_BM) <- c("FT0","FT1","FT2")
colnames(adonis_BM) <-c("p_w","r2_w","p_uw","r2_uw")

######################################################### dominated phylum type, data for figure 2d ##########################################################
# calculate the relative abundance of phylum
ps.rel <- transform_sample_counts(ps_rarefied, function(x) x/sum(x)*100)
tax <- data.table::as.data.table(tax_table(ps.rel),keep.rownames = TRUE)
# change to charater for easy adjust level
glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
# get the relative abundance
ps.melt_phylum <- ps.melt %>%
  dplyr::group_by(FullID,Phylum) %>%
  dplyr::mutate(mean=mean(Abundance))
ps.melt_phylum <- ps.melt_phylum[order(ps.melt_phylum$mean, decreasing = TRUE),]
# keep the highest abundance phylum
ps.melt_phylum <- ps.melt_phylum[!duplicated(ps.melt_phylum$FullID),]
# change those with the highest abundance value < 50 as mixed type
ps.melt_phylum$Phylum[ps.melt_phylum$mean < 50] <- "mixed"
table(ps.melt_phylum$FT,ps.melt_phylum$Phylum,ps.melt_phylum$delivery_mode)
# add new variable "type"
for (i in 1:dim(ps.melt_phylum)[1]) {
  if(ps.melt_phylum$Phylum[i] == "Proteobacteria") {ps.melt_phylum$type[i] <- "Proteobacteria"
  } else if (trimws(ps.melt_phylum$Phylum[i]) == "Firmicutes") {ps.melt_phylum$type[i] <- "Firmicutes"
  } else {ps.melt_phylum$type[i] <- "Other"
  }
}  
metadata_type <- ps.melt_phylum

#### type and delivery mode
library(nnet)
fitt0 <- multinom(type~delivery_mode, data = subset(metadata_type, FT=="0"))
summary(fitt0)
z0 <- summary(fitt0)$ coefficients/summary(fitt0)$standard.errors
p0 <- (1-pnorm(abs(z0),0,1))*2
fitt1 <- multinom(type~delivery_mode, data = subset(metadata_type, FT=="1"))
summary(fitt1)
z1 <- summary(fitt1)$ coefficients/summary(fitt1)$standard.errors
p1 <- (1-pnorm(abs(z1),0,1))*2
fitt2 <- multinom(type~delivery_mode, data = subset(metadata_type, FT=="2"))
summary(fitt2)
z2 <- summary(fitt2)$ coefficients/summary(fitt2)$standard.errors
p2 <- (1-pnorm(abs(z2),0,1))*2

# another method
# add a new variable "typep": whether is Proteobacteria dominated
for (i in 1:dim(metadata_type)[1]) {
  if(metadata_type$Phylum[i] == "Proteobacteria") {metadata_type$typep[i] <- "Yes"
  } else {metadata_type$typep[i] <- "No"
  }
}  
# add a new variable "typef": whether is Firmicutes dominated
for (i in 1:dim(metadata_type)[1]) {
  if(metadata_type$Phylum[i] == "Firmicutes") {metadata_type$typef[i] <- "Yes"
  } else {metadata_type$typef[i] <- "No"
  }
} 
# data preparation 
library(scales)
metadata_type$group2 <- factor(metadata_type$group2, levels = c("C_0","V_0","C_7","V_7","C_14","V_14"))
metadata_type$typep <- as.factor(metadata_type$typep)
metadata_type$typef <- as.factor(metadata_type$typef)
type0 <- subset(metadata_type,FT=="0")
type1 <- subset(metadata_type,FT=="1")
type2 <- subset(metadata_type,FT=="2")
typeCS<- subset(metadata_type, delivery_mode == "C")
typeVD<- subset(metadata_type, delivery_mode == "V")
library(Publish)
library(nnet)
#Proteobacteria
fitp0_p<-glm(typep~delivery_mode + hospital + fortification + GA + DOL + SGA + antibiotics,data = type0, family = "binomial")
publish(fitp0_p)
proteobacteria0 <- as.data.frame(tidy(fitp0_p))
fitp1_p<-glm(typep~delivery_mode + hospital + fortification + GA + DOL + SGA + antibiotics,data = type1, family = "binomial")
proteobacteria1 <- as.data.frame(tidy(fitp1_p))
fitp2_p<-glm(typep~delivery_mode + hospital + fortification + GA + DOL + SGA + antibiotics,data = type2, family = "binomial")
proteobacteria2 <- as.data.frame(tidy(fitp2_p))
deliverygroup <- data.frame(group1 = c("C","C","C"), group2=c("V", "V", "V"))
proteobacteria_p <- as.data.frame(rbind(proteobacteria0[2,],proteobacteria1[2,],proteobacteria2[2,]))
proteobacteria_p <- proteobacteria_p[,-1]
row.names(proteobacteria_p) <- c("FT0","FT1","FT2")
proteobacteria_p <- cbind(deliverygroup,proteobacteria_p)
print(proteobacteria_p, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)
#Firmicutes
fitf0_f<-glm(typef~delivery_mode + hospital + fortification + GA + DOL + SGA + antibiotics,data = type0, family = "binomial")
firmicutes0 <- as.data.frame(tidy(fitf0_f))
fitf1_f<-glm(typef~delivery_mode + hospital + fortification + GA + DOL + SGA + antibiotics,data = type1, family = "binomial")
firmicutes1 <- as.data.frame(tidy(fitf1_f))
fitf2_f<-glm(typef~delivery_mode + hospital + fortification + GA + DOL + SGA + antibiotics,data = type2, family = "binomial")
firmicutes2 <- as.data.frame(tidy(fitf2_f))
deliverygroup <- data.frame(group1 = c("C","C","C"), group2=c("V", "V", "V"))
firmicutes_p <- as.data.frame(rbind(firmicutes0[2,],firmicutes1[2,],firmicutes2[2,]))
firmicutes_p <- firmicutes_p[,-1]
row.names(firmicutes_p) <- c("FT0","FT1","FT2")
firmicutes_p <- cbind(deliverygroup,firmicutes_p)
print(firmicutes_p, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

################################################################### DESeq2, data for figure 2e #################################################################
library(stringr)
library(DESeq2)
ps.full <- phyloseq(OTU, TAX, sample_data(metadata_rarefied), TREE)
ps_FT0f <- subset_samples(ps.full, FT == "0")
ps_FT1f <- subset_samples(ps.full, FT == "1")
ps_FT2f <- subset_samples(ps.full, FT == "2")
### FT0
ps.sample <- ps_FT0f
ps.sample <- tax_glom(ps.sample, taxrank = "Genus", NArm = FALSE)
path_table <- "/desep2/"
dir.create(path_table)
ps <- ps.sample
ps.ds <- phyloseq_to_deseq2(ps, ~ hospital_catch + GA + DOL + SGA + antibiotics + fortification + delivery_mode)
ps.ds <-  DESeq(ps.ds, test="Wald", fitType="parametric", sfType="poscounts")
# result
res = results(ps.ds, cooksCutoff = FALSE, pAdjustMethod = "none")
sigtab = res
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)
# select appearance > 0.5
ps.sample_C <- prune_samples(sample_data(ps.sample)$delivery_mode %in% "C", ps.sample)
ps.sample_V<- prune_samples(sample_data(ps.sample)$delivery_mode %in% "V", ps.sample)
ps.sample_C <- prune_taxa(rowSums(otu_table(ps.sample_C) == 0) < ncol(otu_table(ps.sample_C)) * 0.5, ps.sample)
ps.sample_V <- prune_taxa(rowSums(otu_table(ps.sample_V) == 0) < ncol(otu_table(ps.sample_V)) * 0.5, ps.sample)
otu.f_C <- data.frame(otu_table(ps.sample_C))
otu.f_F <- data.frame(otu_table(ps.sample_V))
otu.f <- rbind(otu.f_C,otu.f_F)
otu.f_rname <- row.names(otu.f)
sigtabf <- sigtab[rownames(sigtab) %in% otu.f_rname,]
# with relative abundance > 1%
ps_rel <- transform_sample_counts(ps.full, function(x) x / sum(x) * 100)
# agglomerate taxa
glom <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
ps_melt <- psmelt(glom)
# change to character for easy-adjusted level
ps_melt$Genus <- as.character(ps_melt$Genus)
ps_melt <- ps_melt %>%
  dplyr::group_by(day, Genus) %>%
  dplyr::mutate(mean=mean(Abundance))
# select grouped mean > 1
keep <- unique(ps_melt$Genus[ps_melt$mean > 1])
sigtabf <- sigtab[sigtab$Genus %in% keep,]
sigtabf <- sigtabf[,-6]
#FDR for multiple test (genera)
library("multtest")
padj <- mt.rawp2adjp(sigtabf$pvalue, proc = "TSBH", alpha = 0.2)
padj <- data.frame(padj$adjp)
colnames(padj)[1] <- 'pvalue'
sigtabf <- merge(sigtabf, padj, by="pvalue")
colnames(sigtabf)[13] <- 'padj'
sigtabf
### FT1 or FT2: change the ps.sample <- ps_FT1f or ps_FT2f
