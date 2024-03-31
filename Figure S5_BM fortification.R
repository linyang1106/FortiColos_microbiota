setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4. Microbiota_fortifier.R")
##################################################################################
#PcoA CS S3a-c
# subset ps
ps_CS <- subset_samples(ps_rarefied, delivery_mode == "C")
ps_CS_FT0 <- subset_samples(ps_CS, FT == "0")
ps_CS_FT1 <- subset_samples(ps_CS, FT == "1")
ps_CS_FT2 <- subset_samples(ps_CS, FT == "2")

ps_VB <- subset_samples(ps_rarefied, delivery_mode == "V")
ps_VB_FT0 <- subset_samples(ps_VB, FT == "0")
ps_VB_FT1 <- subset_samples(ps_VB, FT == "1")
ps_VB_FT2 <- subset_samples(ps_VB, FT == "2")

#PCoA
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
  adonis_sub <- adonis2(dist ~ fortification + hospital_catch + GA + DOL + SGA + antibiotics, data = points, permutations = 999, by="margin")
  adonis_R2 <- adonis_sub$R2
  adonis_p <- adonis_sub$`Pr(>F)`
  # centroid
  centroids <- aggregate(cbind(x,y) ~ var_color,data = points, mean)
  colnames(centroids) <- c(color_c, "x", "y")
  points <- merge(points, centroids, by = color_c,
                  suffixes = c("", ".centroid"))
  # re-evaluate
  var_color<-eval(substitute(get(color_c)), points)
  
  p_ordi <- ggplot(points, aes(x = x, y = y,  fill = var_color, color = var_color, shape="fortification")) +
    scale_fill_manual("fortification", values = mycolors2) + 
    scale_color_manual("fortification", values = mycolors2) +
    geom_point(alpha = 0.8, size = 2, aes(color = var_color)) + 
    scale_shape_manual(values=c(17,17))+
    stat_ellipse(aes(fill=var_color),level = 0.8, geom = "polygon", alpha = 1/12) +
    labs(x = paste0("PCoA 1 (", round(100 * eig[1] / sum(eig), 1), " %)", sep=""),
         y = paste0("PCoA 2 (", round(100 * eig[2] / sum(eig), 1), " %)", sep="")) +
    ggtitle(paste("CS, ",value = "FT",format(metadata$FT[1]),
                  ", R2 =", format(adonis_R2[1],digits = 1),", P =",format.pval(adonis_p[1],1,0.01,nsmall=1))) +
    mytheme_n
}

pCS0w <- plot_sub_ordi(ps_CS_FT0, "Weighted Unifrac", "fortification")
pCS1w <- plot_sub_ordi(ps_CS_FT1, "Weighted Unifrac", "fortification")
pCS2w <- plot_sub_ordi(ps_CS_FT2, "Weighted Unifrac", "fortification")

# PcoA VB S3d-f
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
  adonis_sub <- adonis2(dist ~ fortification + hospital_catch + GA + DOL + SGA + antibiotics, data = points, permutations = 999, by="margin")
  adonis_R2 <- adonis_sub$R2
  adonis_p <- adonis_sub$`Pr(>F)`
  # centroid
  centroids <- aggregate(cbind(x,y) ~ var_color,data = points, mean)
  colnames(centroids) <- c(color_c, "x", "y")
  points <- merge(points, centroids, by = color_c,
                  suffixes = c("", ".centroid"))
  # re-evaluate
  var_color<-eval(substitute(get(color_c)), points)
  
  p_ordi <- ggplot(points, aes(x = x, y = y, fill = var_color, color = var_color, shape="fortification")) +
    scale_fill_manual("fortification", values = mycolors2) + 
    scale_color_manual("fortification", values = mycolors2) +
    geom_point(alpha = 0.8, size = 2, aes(color = var_color)) + 
    scale_shape_manual(values=c(15,15))+
    stat_ellipse(aes(fill=var_color),level = 0.8, geom = "polygon", alpha = 1/12) +
    labs(x = paste0("PCoA 1 (", round(100 * eig[1] / sum(eig), 1), " %)", sep=""),
         y = paste0("PCoA 2 (", round(100 * eig[2] / sum(eig), 1), " %)", sep="")) +
    ggtitle(paste("VB, ",value = "FT",format(metadata$FT[1]),
                  ", R2 =", format(adonis_R2[1],digits = 1),", P =",format.pval(adonis_p[1],1,0.01,nsmall=1))) +
    mytheme_n
}

pVB0w <- plot_sub_ordi(ps_VB_FT0, "Weighted Unifrac", "fortification")
pVB1w <- plot_sub_ordi(ps_VB_FT1, "Weighted Unifrac", "fortification")
pVB2w <- plot_sub_ordi(ps_VB_FT2, "Weighted Unifrac", "fortification")

#label
pCS0w2 <- pCS0w + labs(tag = "a") + theme(plot.tag.position = "topleft")
pCS1w2 <- pCS1w + labs(tag = "b") + theme(plot.tag.position = "topleft")
pCS2w2 <- pCS2w + labs(tag = "c") + theme(plot.tag.position = "topleft")
pVB0w2 <- pVB0w + labs(tag = "d") + theme(plot.tag.position = "topleft")
pVB1w2 <- pVB1w + labs(tag = "e") + theme(plot.tag.position = "topleft")
pVB2w2 <- pVB2w + labs(tag = "f") + theme(plot.tag.position = "topleft")
library("cowplot")
pdf("FigrueS4_BMfortification.pdf",height = 6.5, width=12,bg="white")
ggdraw()+
  draw_plot(pCS0w2,x=0,y=0.5,width = 0.32,height = 0.5)+
  draw_plot(pCS1w2,x=0.34,y=0.5,width = 0.32,height = 0.5)+
  draw_plot(pCS2w2,x=0.68,y=0.5,width = 0.32,height = 0.5)+
  draw_plot(pVB0w2,x=0,y=0,width = 0.32,height = 0.5)+
  draw_plot(pVB1w2,x=0.34,y=0,width = 0.32,height = 0.5)+
  draw_plot(pVB2w2,x=0.68,y=0,width = 0.32,height = 0.5)
dev.off()
