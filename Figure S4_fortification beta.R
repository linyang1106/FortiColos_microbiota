setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
##################################################################################
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
    ggtitle(paste(value = method_c ,paste(value = ", FT",format(metadata$FT[1])),
                  ", R2 =", format(adonis_R2[1],digits = 1),", P =",format(adonis_p[1], digits = 1))) +
    mytheme
}
p0u <- plot_sub_ordi(ps_FT0, "Unweighted Unifrac", "fortification")
p1u <- plot_sub_ordi(ps_FT1, "Unweighted Unifrac", "fortification")
p2u <- plot_sub_ordi(ps_FT2, "Unweighted Unifrac", "fortification")
p0b <- plot_sub_ordi(ps_FT0, "bray", "fortification")
p1b <- plot_sub_ordi(ps_FT1, "bray", "fortification")
p2b <- plot_sub_ordi(ps_FT2, "bray", "fortification")
p0j <- plot_sub_ordi(ps_FT0, "jaccard", "fortification")
p1j <- plot_sub_ordi(ps_FT1, "jaccard", "fortification")
p2j <- plot_sub_ordi(ps_FT2, "jaccard", "fortification")

#label
p0u2 <- p0u + labs(tag = "a") + theme(plot.tag.position = "topleft")
p1u2 <- p1u + labs(tag = "b") + theme(plot.tag.position = "topleft")
p2u2 <- p2u + labs(tag = "c") + theme(plot.tag.position = "topleft")
p0b2 <- p0b + labs(tag = "d") + theme(plot.tag.position = "topleft")
p1b2 <- p1b + labs(tag = "e") + theme(plot.tag.position = "topleft")
p2b2 <- p2b + labs(tag = "f") + theme(plot.tag.position = "topleft")
p0j2 <- p0j + labs(tag = "g") + theme(plot.tag.position = "topleft")
p1j2 <- p1j + labs(tag = "h") + theme(plot.tag.position = "topleft")
p2j2 <- p2j + labs(tag = "i") + theme(plot.tag.position = "topleft")

library("cowplot")
pdf("FigrueS3_fortification_beta.pdf",height = 10, width=12,bg="white")
ggdraw()+
  draw_plot(p0u2,x=0,y=0.68, width = 0.32,height = 0.32)+
  draw_plot(p1u2,x=0.34,y=0.68, width = 0.32,height = 0.32)+
  draw_plot(p2u2,x=0.68,y=0.68, width = 0.32,height = 0.32)+
  draw_plot(p0b2,x=0,y=0.33,width = 0.32,height = 0.32)+
  draw_plot(p1b2,x=0.34,y=0.33,width = 0.32,height = 0.32)+
  draw_plot(p2b2,x=0.68,y=0.33,width = 0.32,height = 0.32)+
  draw_plot(p0j2,x=0,y=0,width = 0.32,height = 0.32)+
  draw_plot(p1j2,x=0.34,y=0,width = 0.32,height = 0.32)+
  draw_plot(p2j2,x=0.68,y=0,width = 0.32,height = 0.32)
dev.off()
