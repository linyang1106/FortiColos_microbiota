setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
##################################################################################
mycolors_hos <- c("#46ACC6", "#1A86A3", "#0E5B76", "#023047","#FFCA5F", "#FEA809", "#FB8402",  "#CD5C08")

# figure S2a
rich <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon"))
mapping <- data.frame(sample_data(ps_rarefied))
index <- match(rownames(rich), rownames(mapping))
rich_full <- cbind.data.frame(rich, mapping[index, ])

p_shannon_h <- ggplot(rich_full, aes(x=FT, y=Shannon, fill=hospital)) + 
  stat_boxplot(geom = "errorbar", position = position_dodge(), width=0.6, linetype = 1)+
  geom_boxplot(position = position_dodge(), width=0.6)+
  labs(x="FT", y = "Shannon index", title = "")+
  scale_fill_manual("hospital", values = mycolors_hos) + 
  mytheme_r

# figure S2c-e
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
  adonis_sub <- adonis2(dist ~ hospital +birth_mode + fortification +  GA + DOL + SGA + antibiotics, data = points, permutations = 999, by="margin")
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
    scale_fill_manual("hospital", values = mycolors_hos) + 
    scale_color_manual("hospital", values = mycolors_hos) +
    geom_point(alpha = 0.8, size = 2, aes(color = var_color)) + 
    stat_ellipse(aes(fill=var_color),level = 0.8, geom = "polygon", alpha = 1/12) +
    labs(x = paste0("PCoA 1 (", round(100 * eig[1] / sum(eig), 1), " %)", sep=""),
         y = paste0("PCoA 2 (", round(100 * eig[2] / sum(eig), 1), " %)", sep="")) +
    mytheme_n
}
# figure S2c
p0w_h <- plot_sub_ordi(ps_FT0, "Weighted Unifrac", "hospital")
# figure S2d
p1w_h <- plot_sub_ordi(ps_FT1, "Weighted Unifrac", "hospital")
# figure S2e
p2w_h <- plot_sub_ordi(ps_FT2, "Weighted Unifrac", "hospital")

# figure S2b
# include large palette
ps_rel <- transform_sample_counts(ps_rarefied, function(x) x / sum(x) * 100)
# agglomerate taxa
glom <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
ps_melt <- psmelt(glom)
# change to character for easy-adjusted level
ps_melt$Genus <- as.character(ps_melt$Genus)
ps_melt <- ps_melt %>%
  dplyr::group_by(day, Genus) %>%
  dplyr::mutate(mean=mean(Abundance))
# select grouped median > 1
keep <- unique(ps_melt$Genus[ps_melt$mean > 1])
ps_melt$Genus[!(ps_melt$Genus %in% keep)] <- "<1% mean abund."
# to get the same rows together
ps_melt_sum <- ps_melt %>%
  dplyr::group_by(FullID, hospital, Genus,FT) %>%
  dplyr::summarise(Abundance = sum(Abundance))
ps_melt_sum$FT <- factor(ps_melt_sum$FT, levels = c("0", "1","2"))
ps_melt_sum <- ps_melt_sum%>%  group_by(Genus)%>%  mutate(mean=mean(Abundance))%>%  ungroup()%>%  arrange(desc(mean))%>%  mutate(Genus=factor(Genus, levels = unique(Genus)))
# build color palette
library(paletteer)
newlabel_FT <-c("0"="FT0", "1"="FT1","2"="FT2")
d_palettes<- palettes_d_names
paletteer_d("ggthemes::Tableau_20",n=20)
p_bar_genus_h <- ggplot(
  ps_melt_sum,
  aes(x = hospital, y = Abundance, fill = Genus)) +
  geom_bar(stat = "summary",position = "stack", width = 0.8, aes(fill = Genus), fun = "mean") +
  labs(x = "Hospital", y = "Relative Abundance (%)") +
  facet_wrap(~FT,labeller=as_labeller(newlabel_FT)) +
  scale_fill_manual(
    values = c( "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", "#F1CE63FF", "#499894FF", "#86BCB6FF", "#E15759FF",
                "#FF9D9AFF", "#79706EFF", "#BAB0ACFF", "#D37295FF", "#FABFD2FF", "#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF"), 
    labels=c(substitute(paste(italic('Staphylococcus'))),
             substitute(paste("Unclassified",italic(' Enterobacteriaceae'))),
             substitute(paste(italic('Enterococcus'))),
             substitute(paste(italic('Bifidobacterium'))),
             substitute(paste("<1% mean abund.")),
             substitute(paste(italic('Clostridium'))),
             substitute(paste('lactobacilli')),
             substitute(paste(italic('Veillonella'))),
             substitute(paste(italic('Streptococcus'))),
             substitute(paste(italic('Klebsiella'))),
             substitute(paste('Unclassified',italic(' Clostridiaceae'))),
             substitute(paste(italic('Corynebacterium'))),
             substitute(paste(italic('Haemophilus'))),
             substitute(paste(italic('Serratia '))))) +
  mytheme_abundance

# label
p_shannon_h2 <- p_shannon_h + labs(tag = "a") + theme(plot.tag.position = "topleft")
p_bar_genus_h2 <- p_bar_genus_h + labs(tag = "b") + theme(plot.tag.position = "topleft")
p0w_h2 <- p0w_h + labs(tag = "c") + theme(plot.tag.position = "topleft")
p1w_h2 <- p1w_h + labs(tag = "d") + theme(plot.tag.position = "topleft")
p2w_h2 <- p2w_h + labs(tag = "e") + theme(plot.tag.position = "topleft")
# merge
library("cowplot")
pdf("FigrueS2_hospital.pdf",height = 7, width=12,bg="white")
ggdraw()+
  draw_plot(p_shannon_h2,x=0,y=0.50, width = 0.42,height = 0.5)+
  draw_plot(p_bar_genus_h2,x=0.42,y=0.50,width = 0.58,height = 0.5)+
  draw_plot(p0w_h2,x=0,y=0, width = 0.33,height = 0.5)+
  draw_plot(p1w_h2,x=0.33,y=0, width = 0.33,height = 0.5)+
  draw_plot(p2w_h2,x=0.66,y=0, width = 0.33,height = 0.5)
dev.off()

# further labelled by Adobe illustrator