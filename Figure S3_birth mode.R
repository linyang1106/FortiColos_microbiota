setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
##################################################################################
# figure S3A
# calculate shannon index
rich <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon"))
mapping <- data.frame(sample_data(ps_rarefied))
index <- match(rownames(rich), rownames(mapping))
rich_full <- cbind.data.frame(rich, mapping[index, ])
rich_full$BMlabel <- factor(rich_full$BMlabel, levels=c("CS","VB"))
rich_full$birth_mode <- factor(rich_full$birth_mode, levels = c("C", "V"))
# plot
#figure S4a
p_shannon_BM <- ggplot(rich_full, aes(x=FT, y=Shannon, fill=birth_mode)) + 
  stat_boxplot(geom = "errorbar", position = position_dodge(), width=0.6, linetype = 1)+
  geom_boxplot(position = position_dodge(), width=0.6)+
  scale_y_continuous(breaks = seq(0,3,1),labels = c(" 0.00","1.00","2.00","3.00"))+
  labs(x="FT", y = "Shannon", title = "")+
  scale_fill_manual("BMlabel", values = mycolors_BM) + 
  mytheme_n
#figure S4b
p_observed_BM <- ggplot(rich_full, aes(x=FT, y=Observed, fill=BMlabel)) + 
  stat_boxplot(geom = "errorbar", position = position_dodge(), width=0.6, linetype = 1)+
  geom_boxplot(position = position_dodge(), width=0.6)+
  labs(x="FT", y = "Observed zOTUs", title = "")+
  scale_fill_manual("Birth mode", values = mycolors_BM) + 
  mytheme_r

#figure S4c
#Phylum
# include large palette
ps_rel <- transform_sample_counts(ps_rarefied, function(x) x / sum(x) * 100)
# agglomerate taxa
glom <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
ps_melt <- psmelt(glom)
# change to character for easy-adjusted level
ps_melt$Phylum <- as.character(ps_melt$Phylum)
ps_melt <- ps_melt %>%
  dplyr::group_by(day, Phylum) %>%
  dplyr::mutate(mean=mean(Abundance))
# select grouped median > 1
keep <- unique(ps_melt$Phylum[ps_melt$mean > 1])
ps_melt$Phylum[!(ps_melt$Phylum %in% keep)] <- "<1% mean abund."
# to get the same rows together
ps_melt_sum <- ps_melt %>%
  dplyr::group_by(FullID, birth_mode, Phylum,FT) %>%
  dplyr::summarise(Abundance = sum(Abundance))
ps_melt_sum$FT <- factor(ps_melt_sum$FT, levels = c("0", "1","2"))
# order
ps_melt_sum <- ps_melt_sum%>%  group_by(Phylum)%>%  mutate(mean=mean(Abundance))%>%  ungroup()%>%  arrange(desc(mean))%>%  mutate(Phylum=factor(Phylum, levels = unique(Phylum)))
# plot
newlabel_BM <-c("C"="CS", "V"="VB")
p_bar_phylum_BM <- ggplot(
  ps_melt_sum,
  aes(x = FT, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "summary",position = "stack", width = 0.8, aes(fill = Phylum), fun = "mean") +
  labs(x = "FT", y = "Relative Abundance (%)") +
  facet_wrap(~birth_mode,labeller=as_labeller(newlabel_BM))+
  scale_fill_manual(breaks=c("<1% mean abund.","Firmicutes","Proteobacteria", "Actinobacteria", "Bacteroidetes"),
                    values = c("#8CD17DFF","#aa96da","#424874","#bbded6","#71c9ce"),
                    labels=c(substitute(paste("<1% mean abund.")),
                             substitute(paste(italic('Actinobacteria'))),
                             substitute(paste(italic('Bacterioidetes'))),
                             substitute(paste(italic('Firmicutes'))),
                             substitute(paste(italic('Proteobacteria')))))+
  mytheme_abundance

# label
p_shannon_BM2 <- p_shannon_BM + labs(tag = "a") + theme(plot.tag.position = "topleft")
p_observed_BM2 <- p_observed_BM + labs(tag = "b") + theme(plot.tag.position = "topleft")
p_bar_phylum_BM2 <- p_bar_phylum_BM + labs(tag = "c") + theme(plot.tag.position = "topleft")

# merge
pdf("FigrueS3_birthmode.pdf",height = 8, width=10,bg="white")
ggdraw()+
  draw_plot(p_shannon_BM2,x=0,y=0.50, width = 0.43,height = 0.5)+
  draw_plot(p_observed_BM2,x=0.46,y=0.50,width = 0.51,height = 0.5)+
  draw_plot(p_bar_phylum_BM2,x=0,y=0, width = 1,height = 0.5)
dev.off()

