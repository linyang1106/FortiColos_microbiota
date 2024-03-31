setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("5. Microbiota_birth mode.R")
##################################################################################
# PCoA
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
  adonis_sub <- adonis2(dist ~ birth_mode + fortification + hospital + GA + DOL + SGA + antibiotics, data = points, permutations = 999, by="margin")
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
    scale_fill_manual("birth_mode", values = mycolors_BM) + 
    scale_color_manual("birth_mode", values = mycolors_BM) +
    geom_point(alpha = 0.8, size = 2, aes(color = var_color)) + 
    stat_ellipse(aes(fill=var_color),level = 0.8, geom = "polygon", alpha = 1/12) +
    labs(x = paste0("PCoA 1 (", round(100 * eig[1] / sum(eig), 1), " %)", sep=""),
         y = paste0("PCoA 2 (", round(100 * eig[2] / sum(eig), 1), " %)", sep=""),
         title = paste0(color_c, "R2=", round(adonis_R2[1],1),
                        ", P=", format(adonis_p[1], digits = 1))) +
    ggtitle(paste(value = "FT",format(metadata$FT[1]),
                  ", R2=", format(adonis_R2[1],digits = 1),", P=",format.pval(adonis_p[1],1,0.01,nsmall=2))) +
    mytheme_n
}
# figure 2a
p0wBM <- plot_sub_ordi(ps_FT0, "Weighted Unifrac", "birth_mode")
# figure 2b
p1wBM <- plot_sub_ordi(ps_FT0, "Weighted Unifrac", "birth_mode")
# figure 2c
p2wBM <- plot_sub_ordi(ps_FT0, "Weighted Unifrac", "birth_mode")

# figure 2d
newlabel <-c("C"="CS",
             "V"="VB")
metadata_rarefied$type <- factor(metadata_rarefied$type,level=c("Other","Firmicutes","Proteobacteria"))
p_type <- ggplot(data=metadata_rarefied, aes(x=FT,fill=type))+
  geom_bar(position = "fill",color="black",width=0.8)+
  labs(x="FT", y = "Percentage(%)", title = NULL)+
  facet_wrap(~birth_mode,labeller=as_labeller(newlabel))+
  scale_y_continuous(labels = c("0 ","25 ","50 ","75 ","100 "))+
  scale_x_discrete(labels = c("0","1","2","0","1","2"))+
  scale_fill_manual(values = c("#d57149","#c3aed6","#8675a9"),
                    breaks=c("Other","Firmicutes","Proteobacteria"),
                    labels=c("Other",
                             substitute(paste(italic('Firmicutes'),"-dominating")),
                             substitute(paste(italic('Proteobacteria'),"-dominating"))))+
  mytheme_abundance

# figure 2e
ps_keep <- prune_samples(colSums(otu_table(ps_rarefied)) >= lib_size, ps_rarefied)
ps_keep <- prune_taxa(rowSums(otu_table(ps_keep)) > 0, ps_keep)
ps_rel <- transform_sample_counts(ps_keep, function(x) x / sum(x) * 100)
glom_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
ps_melt_genus <- psmelt(glom_genus)
ps_melt_genus <- as.data.frame(ps_melt_genus)
stapabundance <-ps_melt_genus %>% filter(Genus == "Staphylococcus")
stapabundance$day <- factor(stapabundance$FT, level= c("0","1","2"))

pbox_sta <- ggplot(stapabundance, aes(x = FT, y = Abundance, fill=birth_mode)) +
  stat_boxplot(geom = "errorbar", position = position_dodge(), width=0.6, linetype = 1)+
  geom_boxplot(position = position_dodge(),width=0.6,
               aes(fill = birth_mode), alpha = 1) +
  labs(x = "FT", y = "Relative abundance(%)", title = substitute(paste(italic('Staphylococcus')))) +
  scale_fill_manual("birth_mode", values = mycolors_BM) + 
  scale_y_continuous(breaks = seq(0,100,25),labels = c(" 0.0","25.0","50.0","75.0","100.0"))+
  mytheme_n

# figure 2f
source("Figure rabuplots_FortiColos.R")
p_relative_BM

#label
p0wBM2 <- p0wBM + labs(tag = "a") + theme(plot.tag.position = "topleft")
p1wBM2 <- p1wBM + labs(tag = "b") + theme(plot.tag.position = "topleft")
p2wBM2 <- p2wBM + labs(tag = "c") + theme(plot.tag.position = "topleft")
p_type2 <- p_type + labs(tag = "d") + theme(plot.tag.position = "topleft")
pbox_sta2 <- pbox_sta + labs(tag = "e") + theme(plot.tag.position = "topleft")
p_relative_BM2 <- p_relative_BM + labs(tag = "f") + theme(plot.tag.position = "topleft")

library("cowplot")
pdf("p_Figrue2_birthmode.pdf",height = 15, width=12,bg="white")
ggdraw()+
  draw_plot(p0wBM2,x=0,y=0.8, width = 0.32,height = 0.2)+
  draw_plot(p1wBM2,x=0.34,y=0.8, width = 0.32,height = 0.2)+
  draw_plot(p2wBM2,x=0.68,y=0.8, width = 0.32,height = 0.2)+
  draw_plot(p_type2,x=0,y=0.58,width = 0.53,height = 0.2)+
  draw_plot(pbox_sta2,x=0.68,y=0.58,width = 0.32,height = 0.2)+
  draw_plot(p_relative_BM2,x=0,y=0,width = 0.99,height = 0.56)
dev.off()
