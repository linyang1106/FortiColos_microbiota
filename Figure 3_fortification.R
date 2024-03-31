setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
##################################################################################
# alpha diversity
rich <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon"))
mapping <- data.frame(sample_data(ps_rarefied))
index <- match(rownames(rich), rownames(mapping))
rich_full <- cbind.data.frame(rich, mapping[index, ])
rich_0 <- subset(rich_full, FT == "0")
rich_1 <- subset(rich_full, FT == "1")
rich_2 <- subset(rich_full, FT == "2")
#figure 3d
p_shannon <- ggplot(rich_full, aes(x=FT, y=Shannon, fill=fortification)) + 
  stat_boxplot(geom = "errorbar", position = position_dodge(), width=0.6, linetype = 1)+
  geom_boxplot(position = position_dodge(), width=0.6)+
  scale_y_continuous(breaks = seq(0,3,1),labels = c(" 0.00","1.00","2.00","3.00"))+
  labs(x="FT", y = "Shannon index", title = "")+
  scale_fill_manual("fortification", values = mycolors_for) + 
  mytheme_n
#figure 3e
p_observed <- ggplot(rich_full, aes(x=FT, y=Observed, fill=fortification)) + 
  stat_boxplot(geom = "errorbar", position = position_dodge(), width=0.6, linetype = 1)+
  geom_boxplot(position = position_dodge(), width=0.6)+
  labs(x="FT", y = "Observed zOTUs", title = "")+
  scale_fill_manual("fortification", values = mycolors_for) + 
  mytheme_n

#PcoA
plot_sub_ordi <- function(ps_rarefied_f, method_c, color_c) {
  metadata <- data.frame(sample_data(ps_rarefied_f))
  # distance and pca
  dist <- phyloseq::distance(ps_rarefied_f, method = method_c)
  pcoa <- cmdscale(dist, k = 3, eig = T)
  pcs <- as.data.frame(pcoa$points)
  colnames(pcs) <- c("x", "y", "z")
  eig <- pcoa$eig
  points <- cbind.data.frame(pcs, metadata)
  var_color <- eval(substitute(get(color_c)), points)
  # adonis
  set.seed(123)
  adonis_sub <- adonis2(dist ~ fortification + birth_mode + hospital + GA + DOL + SGA + antibiotics, data = points, permutations = 999, by="margin")
  adonis_R2 <- adonis_sub$R2
  adonis_p <- adonis_sub$`Pr(>F)`
  # centroid
  centroids <- aggregate(cbind(x,y) ~ var_color,data = points, mean)
  colnames(centroids) <- c(color_c, "x", "y")
  points <- merge(points, centroids, by = color_c,
                  suffixes = c("", ".centroid"))
  # re-evaluate
  var_color<-eval(substitute(get(color_c)), points)
  
  p_ordi <- ggplot(points, aes(x = x, y = y, fill = var_color, color = var_color)) +
    scale_fill_manual("fortification", values = mycolors_for) + 
    scale_color_manual("fortification", values = mycolors_for) +
    geom_point(alpha = 0.8, size = 2, aes(color = var_color)) + 
    stat_ellipse(aes(fill=var_color),level = 0.8, geom = "polygon", alpha = 1/12) +
    labs(x = paste0("PCoA 1 (", round(100 * eig[1] / sum(eig), 1), " %)", sep=""),
         y = paste0("PCoA 2 (", round(100 * eig[2] / sum(eig), 1), " %)", sep="")) +
    ggtitle(paste(value = "FT",format(metadata$FT[1]),
                  ", R2 =", format(adonis_R2[1],digits = 1),", P =",format.pval(adonis_p[1],1,0.01,nsmall=1))) +
    mytheme_n
}
#figure 3a
p0w <- plot_sub_ordi(ps_FT0, "Weighted Unifrac", "fortification")
#figure 3b
p1w <- plot_sub_ordi(ps_FT1, "Weighted Unifrac", "fortification")
#figure 3c
p2w <- plot_sub_ordi(ps_FT2, "Weighted Unifrac", "fortification")

#figure 3f
ps_keep <- prune_samples(colSums(otu_table(ps_rarefied)) >= lib_size, ps_rarefied)
ps_keep <- prune_taxa(rowSums(otu_table(ps_keep)) > 0, ps_keep)
ps_rel <- transform_sample_counts(ps_keep, function(x) x / sum(x) * 100)
glom_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
ps_melt_genus <- psmelt(glom_genus)
ps_melt_genus <- as.data.frame(ps_melt_genus)
stapabundance <-ps_melt_genus %>% filter(Genus == "Staphylococcus")
stapabundance$day <- factor(stapabundance$day, level= c("0","7","14"))
pbox_sta <- ggplot(stapabundance, aes(x = FT, y = Abundance, fill=fortification)) +
  stat_boxplot(geom = "errorbar", position = position_dodge(), width=0.6, linetype = 1)+
  geom_boxplot(position = position_dodge(),width=0.6,
               aes(fill = fortification), alpha = 1) +
  labs(x = "FT", y = "Relative abundance(%)",title = substitute(paste(italic('Staphylococcus')))) +
  scale_fill_manual("fortification", values = mycolors_for) + 
  scale_y_continuous(breaks = seq(0,100,25),labels = c(" 0.0","25.0","50.0","75.0","100.0"))+
  mytheme_n

#figure 3g
source("Figure rabuplots_FortiColos.R")
p_relative_for

#label
p0w2 <- p0w + labs(tag = "a") + theme(plot.tag.position = "topleft")
p1w2 <- p1w + labs(tag = "b") + theme(plot.tag.position = "topleft")
p2w2 <- p2w + labs(tag = "c") + theme(plot.tag.position = "topleft")
p_shannon2 <- p_shannon + labs(tag = "d") + theme(plot.tag.position = "topleft")
p_observed2 <- p_observed + labs(tag = "e") + theme(plot.tag.position = "topleft")
pbox_sta2 <- pbox_sta + labs(tag = "f") + theme(plot.tag.position = "topleft")
p_relative_for2 <- p_relative_for + labs(tag = "g") + theme(plot.tag.position = "topleft")

library("cowplot")
pdf("Figrue3_fortification.pdf",height = 15, width=12,bg="white")
ggdraw()+
  draw_plot(p0w2,x=0,y=0.8, width = 0.32,height = 0.2)+
  draw_plot(p1w2,x=0.34,y=0.8, width = 0.32,height = 0.2)+
  draw_plot(p2w2,x=0.68,y=0.8, width = 0.32,height = 0.2)+
  draw_plot(p_shannon2,x=0,y=0.58,width = 0.32,height = 0.2)+
  draw_plot(p_observed2,x=0.34,y=0.58,width = 0.32,height = 0.2)+
  draw_plot(pbox_sta,x=0.68,y=0.58,width = 0.32,height = 0.2)+
  draw_plot(p_relative_for2,x=0,y=0,width = 0.99,height = 0.56)
dev.off()
