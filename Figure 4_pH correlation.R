setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("7. Growth.R")
##################################################################################
#figure 4a
genera_growth_for
genera_growth_sta <- na.omit(genera_growth_for[, c("fortification","Rsta02", "Dw02")])

p_wsta02 <- ggscatter(genera_growth_sta,
                       x="Rsta02",y="Dw02",
                       add = "reg.line",
                       conf.int = TRUE,
                       size=1.5,
                       conf.int.level = 0.95,
                       xlab = "Ratio of log10-transformed relative Abundance of Staphylococcus",
                       ylab = "Weight delta Z-score",
                       add.params = list(color = "black",
                                         size = 1,
                                         linetype = 1,
                                         fill = "gray30"),
                       color = "fortification",
                       palette = c("yellowgreen", "#ffcc1d"))+
  labs(x="Ratio of log10-transformed relative Abundance of Staphylococcus",
       ylab = "Weight delta Z-score",
       title = substitute(paste("n = 131, R = -0.28, q < 0.05")))+
  mytheme_r

source("6. pH.R")
# figure 4b
p_ph <- ggplot(metadata_rarefied, aes(x=FT, y=pH, fill=fortification)) + 
  stat_boxplot(geom = "errorbar", position = position_dodge(), width=0.6, linetype = 1)+
  geom_boxplot(position = position_dodge(), width=0.6)+
  labs(x="FT", y = "pH", title = "")+
  scale_fill_manual(values=mycolors_for,
                    name="Fortification",
                    breaks=c("BC","FM85"),
                    labels=c("BC","CF"))+ mytheme_n
# figure 4c-d
my_scaled_dataframe <- cbind(meta_tab$fortification,data.frame(my_scaled_data))# data after CLR 
# Figure 4c
p_ph_staphylococcus <- ggscatter(my_scaled_dataframe,
                                 x="Staphylococcus",y="pH",
                                 conf.int = TRUE,
                                 size=1.5,
                                 conf.int.level = 0.95,
                                 color= "meta_tab$fortification", 
                                 palette = c("yellowgreen", "#ffcc1d"),
                                 add.params = list(color = "black",
                                                   size = 1,
                                                   linetype = 1,
                                                   fill = "gray30"),
                                 add = "reg.line",)+
  labs(x="Relative abundance (CLR transformed)",
       y="pH (CLR transformed)",
       title = substitute(paste(italic('Staphylococcus,'),
                                " n = 538, R = 0.16, q < 0.01")))+
  mytheme_n
#Figure 4d
p_ph_corynebacterium <- ggscatter(my_scaled_dataframe,
                                  x="Corynebacterium",y="pH",
                                  conf.int = TRUE,
                                  size=1.5,
                                  conf.int.level = 0.95,
                                  color= "meta_tab$fortification", 
                                  palette = c("yellowgreen", "#ffcc1d"),
                                  add.params = list(color = "black",
                                                    size = 1,
                                                    linetype = 1,
                                                    fill = "gray30"),
                                  add = "reg.line",)+
  labs(x="Relative abundance (CLR transformed)", 
       y="pH (CLR transformed)", 
       title = substitute(paste(italic('Corynebacterium,'), 
                                " n = 394, R = 0.15, q = 0.02")))+
  mytheme_n
#Figure 4e
p_ph_bifidobacterium <- ggscatter(my_scaled_dataframe,
                                  x="Bifidobacterium",y="pH",
                                  conf.int = TRUE,
                                  size=1.5,
                                  conf.int.level = 0.95,
                                  color= "meta_tab$fortification", 
                                  palette = c("yellowgreen", "#ffcc1d"),
                                  add.params = list(color = "black",
                                                    size = 1,
                                                    linetype = 1,
                                                    fill = "gray30"),
                                  add = "reg.line",)+
  labs(x="Relative abundance (CLR transformed)", 
       y="pH (CLR transformed)",
       title = substitute(paste(italic('Bifidobaterium,'),
                                " n = 248, R = -0.23, q < 0.01")))+
  mytheme_n

#label
p_wsta022 <- p_wsta02 + labs(tag = "a") + theme(plot.tag.position = "topleft")
p_ph2 <- p_ph + labs(tag = "b") + theme(plot.tag.position = "topleft")
p_ph_staphylococcus2 <- p_ph_staphylococcus + labs(tag = "c") + theme(plot.tag.position = "topleft")
p_ph_bifidobacterium2 <- p_ph_bifidobacterium + labs(tag = "e") + theme(plot.tag.position = "topleft")
p_ph_corynebacterium2 <- p_ph_corynebacterium + labs(tag = "d") + theme(plot.tag.position = "topleft")

#merge
pdf("Figrue4_correlation.pdf",height = 6, width=12,bg="white")
ggdraw()+
  draw_plot(p_wsta0142,x=0,y=0.5, width = 0.4,height = 0.5)+
  draw_plot(p_ph2,x=0.5,y=0.5, width = 0.4,height = 0.5)+
  draw_plot(p_ph_staphylococcus2,x=0.33,y=0, width = 0.32,height = 0.5)+
  draw_plot(p_ph_bifidobacterium2,x=0,y=0, width = 0.32,height = 0.5)+
  draw_plot(p_ph_corynebacterium2,x=0.66,y=0,width = 0.33,height = 0.5)
dev.off()
