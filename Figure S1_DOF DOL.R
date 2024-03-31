setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
##################################################################################
# distribution of DOF of each samples
date_sample <- read.table(file = "data/Sample_date.tab", sep = "\t", header = T, row.names = 1, comment.char = "")
# DOF: day of fortification (sample_date - fortification_date); DOL: day of life (sample_date - birth_date)
for (i in 1:dim(date_sample)[1]) {
  if(date_sample$DOF[i] <= 0) {date_sample$FT[i] <- "0"
  } else if (date_sample$DOF[i] > 3 & date_sample$DOF[i] <=10) {date_sample$FT[i] <- "1"
  } else if (date_sample$DOF[i] > 10 &date_sample$DOF[i] <=17) {date_sample$FT[i] <- "2"
  } else {date_sample$FT[i] <- "Excluded"
  }
} 
date_sample$FT <- factor(date_sample$FT, levels= c("0","1","2","Excluded"))

# Figure S1a DOF of all samples (DOF: day of fortification)
mycolors_FT <- c("#D5F0C1","#AAD9BB","#99BC85","gray50")
p_DOF <- ggplot(date_sample, aes(x=DOF,y=frequency(DOF),fill=FT))+
  geom_bar(stat = "identity")+
  scale_fill_manual(name="FT",breaks=c("0","1","2","Excluded"), values = mycolors_FT)+
  labs(x = "DOF(days)", y = "Frequency") +
  scale_x_continuous(expand = c(0,0),breaks = seq(-10,35,5),limits = c(-10,35))+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,80,10),limits = c(0,80))+
  ggtitle("Distribution of samples by DOF (days)")+
  themeS1

# Figure S2a distribution of DOL of each DOF group (DOL: day of life)
date_sample_clean <- subset(date_sample,date_sample$FT != "Excluded")
p_DOL_DOF <-ggplot(date_sample_clean, aes(x=DOL, y= frequency(DOL),fill=FT)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(name="FT",breaks=c("0","1","2"), values = mycolors_FT)+
  labs(x = "DOL(days)", y = "Frequency") +
  scale_x_continuous(expand = c(0,0),breaks = seq(0,45,5),limits = c(0,45))+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,40,5),limits = c(0,40))+
  ggtitle("Distribution of samples by DOL (days)")+
  themeS1

#label
p_DOF2 <- p_DOF + labs(tag = "a") + theme(plot.tag.position = "topleft")
p_DOL_DOF2 <- p_DOL_DOF + labs(tag = "b") + theme(plot.tag.position = "topleft")
#merge
library("cowplot")
pdf("FigrueS1_inclusion.pdf",height = 8, width=12,bg="white")
ggdraw()+
  draw_plot(p_DOF2,x=0,y=0.5, width = 1,height = 0.5)+
  draw_plot(p_DOL_DOF2,x=0,y=0, width = 1,height = 0.5)
dev.off()
