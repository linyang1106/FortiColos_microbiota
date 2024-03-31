setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
#################################################################################

############################################################# rb-rda, data for Supplementary Table 3  ################################################################
otu_FT0 <- t(as.data.frame(otu_table(ps_FT0)))
env_FT0 <- subset(metadata_FT0,select=c(fortification, GA, SGA, birth_mode, DOL, hospital, antibiotics))
rda_FT0 <- capscale(otu_FT0~., env_FT0, dist="bray", add=TRUE)
dbanov_FT0 <- anova.cca(rda_FT0, by="term", permutations = 999)

otu_FT1 <- t(as.data.frame(otu_table(ps_FT1)))
env_FT1 <- subset(metadata_FT1,select=c(fortification, GA, SGA, birth_mode, DOL, hospital, antibiotics))
rda_FT1 <- capscale(otu_FT1~., env_FT1, dist="bray", add=TRUE)
dbanov_FT1 <- anova.cca(rda_FT1, by="term", permutations = 999)

otu_FT2 <- t(as.data.frame(otu_table(ps_FT2)))
env_FT2 <- subset(metadata_FT2,select=c(fortification, GA, SGA, birth_mode, DOL, hospital, antibiotics))
rda_FT2 <- capscale(otu_FT2~., env_FT2, dist="bray", add=TRUE)
dbanov_FT2 <- anova.cca(rda_FT2, by="term", permutations = 999)

############################################################### Alpha diversity, data for figure 3d-f ##########################################################
rich <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon"))
mapping <- data.frame(sample_data(ps_rarefied))
index <- match(rownames(rich), rownames(mapping))
rich_full <- cbind.data.frame(rich, mapping[index, ])
rich_0 <- subset(rich_full, FT == "0")
rich_1 <- subset(rich_full, FT == "1")
rich_2 <- subset(rich_full, FT == "2")

# Shannon index
lm_shannon0 <- lm(Shannon~ fortification + birth_mode + hospital + GA + DOL + SGA + antibiotics, rich_0)
shannon0 <- as.data.frame(tidy(lm_shannon0))
lm_shannon1 <- lm(Shannon~ fortification + birth_mode + hospital + GA + DOL + SGA + antibiotics, rich_1)
shannon1 <- as.data.frame(tidy(lm_shannon1))
lm_shannon2 <- lm(Shannon~fortification + birth_mode + hospital + GA + DOL + SGA + antibiotics, rich_2)
shannon2 <- as.data.frame(tidy(lm_shannon2))
fortification <- data.frame(group1 = c("BC","BC","BC"), group2=c("CF", "CF", "CF"))
shannon_p <- as.data.frame(rbind(shannon0[2,],shannon1[2,],shannon2[2,]))
shannon_p <- shannon_p[,-1]
row.names(shannon_p) <- c("FT0","FT1","FT2")
shannon_p <- cbind(fortification,shannon_p)
print(shannon_p, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

# Observed zOTUs
lm_observed0 <- lm(Observed~fortification + birth_mode + hospital + GA + DOL + SGA + antibiotics, rich_0)
observed0 <- as.data.frame(tidy(lm_observed0))
lm_observed1 <- lm(Observed~fortification  + birth_mode + hospital + GA + DOL + SGA + antibiotics, rich_1)
observed1 <- as.data.frame(tidy(lm_observed1))
lm_observed2 <- lm(Observed~ fortification  + birth_mode + hospital + GA + DOL + SGA + antibiotics, rich_2)
observed2 <- as.data.frame(tidy(lm_observed2))
observed_p <- as.data.frame(rbind(observed0[2,],observed1[2,],observed2[2,]))
observed_p <- observed_p[,-1]
row.names(observed_p) <- c("FT0","FT1","FT2")
observed_p <- cbind(fortification,observed_p)
print(observed_p, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

############################################################### Beta diversity, data for figure 3a-c ##########################################################
################## adonis 
adonis <- function(ps, method_c) {
  set.seed(123)
  dist <- phyloseq::distance(ps, method = method_c)
  metadata <- data.frame(sample_data(ps))
  re_adonis <- adonis2(dist ~ fortification + birth_mode + hospital + GA + DOL + SGA + antibiotics, data = metadata, permutations = 999, by="margin")
}
adonis0w <-adonis(ps_FT0, "Weighted Unifrac")
adonis0u <-adonis(ps_FT0, "Unweighted Unifrac")
adonis1w <-adonis(ps_FT1, "Weighted Unifrac")
adonis1u <-adonis(ps_FT1, "Unweighted Unifrac")
adonis2w <-adonis(ps_FT2, "Weighted Unifrac")
adonis2u <-adonis(ps_FT2, "Unweighted Unifrac")
adonis0b <-adonis(ps_FT0, "bray")
adonis0j <-adonis(ps_FT0, "jaccard")
adonis1b <-adonis(ps_FT1, "bray")
adonis1j <-adonis(ps_FT1, "jaccard")
adonis2b <-adonis(ps_FT2, "bray")
adonis2j <-adonis(ps_FT2, "jaccard")
# extract results of fortifier
weighted_pfor <- as.data.frame(rbind(adonis0w$`Pr(>F)`[1],adonis1w$`Pr(>F)`[1],adonis2w$`Pr(>F)`[1]))
weighted_rfor <- as.data.frame(rbind(adonis0w$R2[1],adonis1w$R2[1],adonis2w$R2[1]))
unweighted_pfor <- as.data.frame(rbind(adonis0u$`Pr(>F)`[1],adonis1u$`Pr(>F)`[1],adonis2u$`Pr(>F)`[1]))
unweighted_rfor <- as.data.frame(rbind(adonis0u$R2[1],adonis1u$R2[1],adonis2u$R2[1]))
adonis_for <- cbind(weighted_pfor,weighted_rfor,unweighted_pfor,unweighted_rfor)
row.names(adonis_for) <- c("FT0","FT1","FT2")
colnames(adonis_for) <-c("p_w","r2_w","p_uw","r2_uw")

################## sub-analysis by birth mode (CS or BV), data for figure supplementary 5
# subset ps
ps_CS <- subset_samples(ps_rarefied, birth_mode == "C")
ps_CS_FT0 <- subset_samples(ps_CS, FT == "0")
ps_CS_FT1 <- subset_samples(ps_CS, FT == "1")
ps_CS_FT2 <- subset_samples(ps_CS, FT == "2")

ps_VB <- subset_samples(ps_rarefied, birth_mode == "V")
ps_VB_FT0 <- subset_samples(ps_VB, FT == "0")
ps_VB_FT1 <- subset_samples(ps_VB, FT == "1")
ps_VB_FT2 <- subset_samples(ps_VB, FT == "2")

# adonis
adonis <- function(ps, method_c) {
  set.seed(123)
  dist <- phyloseq::distance(ps, method = method_c)
  metadata <- data.frame(sample_data(ps))
  re_adonis <- adonis2(dist ~ fortification+ hospital + GA + DOL + SGA + antibiotics, data = metadata, permutations = 999, by="margin")
}
adonisCS0w <-adonis(ps_CS_FT0, "Weighted Unifrac")
adonisCS0u <-adonis(ps_CS_FT0, "Unweighted Unifrac")
adonisCS1w <-adonis(ps_CS_FT1, "Weighted Unifrac")
adonisCS1u <-adonis(ps_CS_FT1, "Unweighted Unifrac")
adonisCS2w <-adonis(ps_CS_FT2, "Weighted Unifrac")
adonisCS2u <-adonis(ps_CS_FT2, "Unweighted Unifrac")

adonisVB0w <-adonis(ps_VB_FT0, "Weighted Unifrac")
adonisVB0u <-adonis(ps_VB_FT0, "Unweighted Unifrac")
adonisVB1w <-adonis(ps_VB_FT1, "Weighted Unifrac")
adonisVB1u <-adonis(ps_VB_FT1, "Unweighted Unifrac")
adonisVB2w <-adonis(ps_VB_FT2, "Weighted Unifrac")
adonisVB2u <-adonis(ps_VB_FT2, "Unweighted Unifrac")

##################################################################### DESeq2, data for figure 3f ##########################################################
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
ps.ds <- phyloseq_to_deseq2(ps, ~ hospital + GA + DOL + SGA + antibiotics + birth_mode + fortification)
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
ps.sample_BC <- prune_samples(sample_data(ps.sample)$fortification %in% "BC", ps.sample)
ps.sample_CF<- prune_samples(sample_data(ps.sample)$fortification %in% "FM85", ps.sample)
ps.sample_BC <- prune_taxa(rowSums(otu_table(ps.sample_BC) == 0) < ncol(otu_table(ps.sample_BC)) * 0.5, ps.sample)
ps.sample_CF <- prune_taxa(rowSums(otu_table(ps.sample_CF) == 0) < ncol(otu_table(ps.sample_CF)) * 0.5, ps.sample)
otu.f_BC <- data.frame(otu_table(ps.sample_BC))
otu.f_CF <- data.frame(otu_table(ps.sample_CF))
otu.f <- rbind(otu.f_BC,otu.f_CF)
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
