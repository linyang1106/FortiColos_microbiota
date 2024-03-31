setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
##################################################################################
# Rhea correlation, ref. https://github.com/Lagkouvardos/Rhea
################### data preparation ##########################
### top 3 genera and alpha diversity
## genera
ps.rarefied.rel <- transform_sample_counts(ps_rarefied, function(x) x/sum(x)*100)
psdat.gen <- tax_glom(ps.rarefied.rel, taxrank = "Genus")
psdat.gen
ps.melt <- psmelt(psdat.gen)
# change to charater for easy adjust level
ps.melt$Genus <- as.character(ps.melt$Genus)
head(ps.melt) 
# long format to wide format
library(reshape2)
m_Tab <- dcast(ps.melt, Genus~Sample, value.var= "Abundance",fun.aggregate = sum)
rownames(m_Tab) <- m_Tab$Genus
m_Tab <- m_Tab[,-1]
colSums(m_Tab)
m_Tabs <- m_Tab[which(rownames(m_Tab) %in% c("Staphylococcus","Unclassified  Enterobacteriaceae","Enterococcus")),]
m_Tablog <- log10(m_Tabs)
m_Tablog <- t(m_Tablog)
## alpha diversity
rich <- estimate_richness(ps_rarefied, measures = c("Observed", "Shannon"))
index_gs <- match(rownames(m_Tablog),rownames(rich))
### microbial file
m_Tabgs <- cbind(m_Tablog, rich[index_gs,])

# metadata
mapping <- data.frame(sample_data(ps_rarefied))
correlationRaw <- as.data.frame(cbind(rownames(metadata_rarefied),metadata_rarefied$record_id,metadata_rarefied$FT,metadata_rarefied$fortification))
rownames(correlationRaw) <- correlationRaw[,1]
correlationRaw <- subset(correlationRaw, select=-V1)
colnames(correlationRaw)[1]<-"record_id"
colnames(correlationRaw)[2]<-"FT"
colnames(correlationRaw)[3]<-"fortification"
original_levels <- levels(metadata_rarefied$FT)
correlationRaw$FT <- factor(metadata_rarefied$FT, levels = original_levels)
# only use the sample with correlation data
mapping <- mapping[rownames(mapping) %in% rownames(correlationRaw),]
# match the mapping table and the correlation data
index <- match(rownames(mapping), rownames(correlationRaw))
mapping_merge <- correlationRaw[index,]
rownames(mapping_merge) <- rownames(mapping)

# match the microbial table
index <- match(rownames(mapping_merge),rownames(m_Tabgs))
meta_tab <- cbind(mapping_merge, m_Tabgs[index,])
head(meta_tab)
# seperate FT
meta_tab_0 <- subset(meta_tab, FT==0)
meta_tab_1 <- subset(meta_tab, FT==1)
meta_tab_2 <- subset(meta_tab, FT==2)
meta_tab_merge <- merge(meta_tab_0,meta_tab_1,by = "record_id",all=TRUE)
meta_tab_merge <- merge(meta_tab_merge,meta_tab_2,by = "record_id",all=TRUE)
rownames(meta_tab_merge) <- meta_tab_merge$record_id
meta_tab_merge <- meta_tab_merge[,-c(2,9,10,16,17)]
colnames(meta_tab_merge) <- c("record_id","fortification","ent0","sta0","unent0","Obs0","Sha0","ent1","sta1","unent1","Obs1","Sha1","ent2","sta2","unent2","Obs2","Sha2")
# delta alpha diversity
meta_tab_merge$DObs01 <- (meta_tab_merge$Obs1)-(meta_tab_merge$Obs0)
meta_tab_merge$DObs02 <- (meta_tab_merge$Obs2)-(meta_tab_merge$Obs0)
meta_tab_merge$DObs12 <- (meta_tab_merge$Obs2)-(meta_tab_merge$Obs1)
meta_tab_merge$DSha01 <- (meta_tab_merge$Sha1)-(meta_tab_merge$Sha0) 
meta_tab_merge$DSha02 <- (meta_tab_merge$Sha2)-(meta_tab_merge$Sha0)
meta_tab_merge$DSha12 <- (meta_tab_merge$Sha2)-(meta_tab_merge$Sha1)
# ratio abundance
meta_tab_merge$Rent01 <- meta_tab_merge$ent1/meta_tab_merge$ent0 
meta_tab_merge$Rent02 <- meta_tab_merge$ent2/meta_tab_merge$ent0 
meta_tab_merge$Rent12 <- meta_tab_merge$ent2/meta_tab_merge$ent1
meta_tab_merge$Rsta01 <- meta_tab_merge$sta1/meta_tab_merge$sta0 
meta_tab_merge$Rsta02 <- meta_tab_merge$sta2/meta_tab_merge$sta0 
meta_tab_merge$Rsta12 <- meta_tab_merge$sta2/meta_tab_merge$sta1
meta_tab_merge$Runent01 <- meta_tab_merge$unent1/meta_tab_merge$unent0 
meta_tab_merge$Runent02 <- meta_tab_merge$unent2/meta_tab_merge$unent0 
meta_tab_merge$Runent12 <- meta_tab_merge$unent2/meta_tab_merge$unent1
# keep calculated data
my_dataall <- meta_tab_merge[,-c(3:17)]
## remove outliner by 1.5 IQR method
my_dataall_cleanIQR<-my_dataall
variable_names <- c("DObs01","DObs02","DObs12",
                    "DSha01","DSha02","DSha12",
                    "Rent01","Rent02","Rent12",
                    "Rsta01","Rsta02","Rsta12",
                    "Runent01","Runent02","Runent12")
for (variable_name in variable_names) {
  # replace outliner as NA
  Q1 <- quantile(my_dataall_cleanIQR[, variable_name], 0.25, na.rm = TRUE)
  Q3 <- quantile(my_dataall_cleanIQR[, variable_name], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_limit <- Q1 - 1.5 * IQR
  upper_limit <- Q3 + 1.5 * IQR
  # save as new datafrome
  my_dataall_cleanIQR[, variable_name] <- ifelse(my_dataall_cleanIQR[, variable_name] < lower_limit | my_dataall_cleanIQR[, variable_name] > upper_limit, NA, my_dataall_cleanIQR[, variable_name])
}

# input delta growth data: delta weight, length and head circumference of each two time points (Z-score)
growth <- read_excel("data/delta weight length head.xlsx",sheet="growth012")
growth$Dw01 <- (growth$weightz1)-(growth$weightz0)
growth$Dw02 <- (growth$weightz2)-(growth$weightz0)
growth$Dw12 <- (growth$weightz2)-(growth$weightz1)
growth$Dl01 <- (growth$lengthz1)-(growth$lengthz0)
growth$Dl02 <- (growth$lengthz2)-(growth$lengthz0)
growth$Dl12 <- (growth$lengthz2)-(growth$lengthz1)
growth$Dh01 <- (growth$headz1)-(growth$headz0)
growth$Dh02 <- (growth$headz2)-(growth$headz0)
growth$Dh12 <- (growth$headz2)-(growth$headz1)
genera_growth <- merge(my_dataall_cleanIQR,growth,by="record_id",all=TRUE)
rownames(genera_growth) <- genera_growth$record_id
genera_growth_for <- genera_growth[,-c(1,18:26)]
genera_growth <- genera_growth[,-c(1:2,18:26)]
### Spearman correlation
my_rcorr <- rcorr(as.matrix(genera_growth), type = "spearman")
# result of FT0-1
rcorrr01 <- data.frame(my_rcorr$r)[c(1,4,7,10,13),c(16,19,22)]
rcorrn01 <- data.frame(my_rcorr$n)[c(1,4,7,10,13),c(16,19,22)]
rcorrp01 <- data.frame(my_rcorr$P)[c(1,4,7,10,13),c(16,19,22)]
# result of FT0-2
rcorrr02 <- data.frame(my_rcorr$r)[c(2,5,8,11,14),c(17,20,23)]
rcorrn02 <- data.frame(my_rcorr$n)[c(2,5,8,11,14),c(17,20,23)]
rcorrp02 <- data.frame(my_rcorr$P)[c(2,5,8,11,14),c(17,20,23)]
# result of FT1-2
rcorrr12 <- data.frame(my_rcorr$r)[c(3,6,9,12,15),c(18,21,24)]
rcorrn12 <- data.frame(my_rcorr$n)[c(3,6,9,12,15),c(18,21,24)]
rcorrp12 <- data.frame(my_rcorr$P)[c(3,6,9,12,15),c(18,21,24)]
# summary result
rcorrr_merge <- cbind(t(rcorrr01),t(rcorrr02))
rcorrr_merge <- as.data.frame(t(cbind(rcorrr_merge,t(rcorrr12))))
r <- unlist(rcorrr_merge)
rcorrn_merge <- cbind(t(rcorrn01),t(rcorrn02))
rcorrn_merge <- as.data.frame(t(cbind(rcorrn_merge,t(rcorrn12))))
n <- unlist(rcorrn_merge)
rcorrp_merge <- cbind(t(rcorrp01),t(rcorrp02))
rcorrp_merge <- as.data.frame(t(cbind(rcorrp_merge,t(rcorrp12))))
p <- unlist(rcorrp_merge)
cor_merge <- as.data.frame(cbind(r,n,p))
rownames(cor_merge) <- c("DObs_w01","DSha_w01","Rent_w01","Rsta_w01","Runent_w01","DObs_w02","DSha_w02","Rent_w02","Rsta_w02","Runent_w02","DObs_w12","DSha_w12","Rent_w12","Rsta_w12","Runent_w12",
                         "DObs_l01","DSha_l01","Rent_l01","Rsta_l01","Runent_l01","DObs_l02","DSha_l02","Rent_l02","Rsta_l02","Runent_l02","DObs_l12","DSha_l12","Rent_l12","Rsta_l12","Runent_l12",
                         "DObs_h01","DSha_h01","Rent_h01","Rsta_h01","Runent_h01","DObs_h02","DSha_h02","Rent_h02","Rsta_h02","Runent_h02","DObs_h12","DSha_h12","Rent_h12","Rsta_h12","Runent_h12")
# FDR adjustment
library(multtest)
p_adj <- mt.rawp2adjp(as.numeric(cor_merge[,3]), proc = "BH", alpha = 0.2) 
p_adj$adjp <- round(p_adj$adjp[order(p_adj$index),], 3)
p_value <-as.data.frame(p_adj$adjp)  
cor_merge$adjp <- p_value$BH
# subset adjp <0.05
p_plot <- subset(cor_merge,cor_merge$adjp < 0.05)
print(p_plot)
