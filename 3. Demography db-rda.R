setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2. Build phyloseq.R")
#########################################################################demographic table
# Rarefy and remove low reads sample
depth_r <- colSums(otu_table(ps))
lib_size <- 10000
ps_rarefied <- rarefy_even_depth(ps, rngseed = 1,
                                 sample.size = lib_size, replace = F)
metadata_rarefied <- data.frame(sample_data(ps_rarefied))
#infants metadata
metadata_infant <- metadata_rarefied[!duplicated(metadata_rarefied$record_id),]
#subset group
ps_FT0 <- subset_samples(ps_rarefied, FT == "0")
ps_FT1 <- subset_samples(ps_rarefied, FT == "1")
ps_FT2 <- subset_samples(ps_rarefied, FT == "2")
metadata_FT0 <- data.frame(sample_data(ps_FT0))
metadata_FT1 <- data.frame(sample_data(ps_FT1))
metadata_FT2 <- data.frame(sample_data(ps_FT2))

#baseline
################################### Table 1 fortifier
myvar_all <- c("GA", "BW", "SGA","PMA", "delivery_mode", "multiplebirth_merge", "sex", "hospital_catch","region_catch", "probiotics_before")
myvar_FT <- c("antibiotics","MMpercent3d")
catvar <- c( "SGA", "delivery_mode", "multiplebirth_merge", "sex", "hospital_catch","region_catch", "antibiotics", "probiotics_before", "fortification")
nonvar <- c("hospital_catch")
#basic characteristics
tab_all <- CreateTableOne(vars = myvar_all, strata = "fortification", data = metadata_infant, factorVars = catvar)
sepline_all <- print(tab_all, showAllLevels = TRUE, exact = nonvar,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)
#interventions
tab_FT0 <- CreateTableOne(vars = myvar_FT, strata = "fortification", data = metadata_FT0, factorVars = catvar)
sepline_FT0 <- print(tab_FT0, showAllLevels = TRUE, exact = nonvar,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

tab_FT1 <- CreateTableOne(vars = myvar_FT, strata = "fortification", data = metadata_FT1, factorVars = catvar)
sepline_FT1 <- print(tab_FT1, showAllLevels = TRUE, exact = nonvar,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

tab_FT2 <- CreateTableOne(vars = myvar_FT, strata = "fortification", data = metadata_FT2, factorVars = catvar)
sepline_FT2 <- print(tab_FT2, showAllLevels = TRUE, exact = nonvar,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

#################################### Supplementary Table 1 birth mode
myvar_all_BM <- c("GA", "BW", "SGA","PMA", "fortification", "multiplebirth_merge", "sex", "hospital_catch","region_catch", "probiotics_before")
tab_all_BM <- CreateTableOne(vars = myvar_all_BM, strata = "delivery_mode", data = metadata_infant, factorVars = catvar)
sepline_all_BM <- print(tab_all_BM, showAllLevels = TRUE, exact = nonvar,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)
#interventions
tab_FT0_BM <- CreateTableOne(vars = myvar_FT, strata = "delivery_mode", data = metadata_FT0, factorVars = catvar)
sepline_FT0_BM <- print(tab_FT0_BM, showAllLevels = TRUE, exact = nonvar,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

tab_FT1_BM <- CreateTableOne(vars = myvar_FT, strata = "delivery_mode", data = metadata_FT1, factorVars = catvar)
sepline_FT1 <- print(tab_FT1_BM, showAllLevels = TRUE, exact = nonvar,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

tab_FT2_BM <- CreateTableOne(vars = myvar_FT, strata = "delivery_mode", data = metadata_FT2, factorVars = catvar)
sepline_FT2_BM <- print(tab_FT2_BM, showAllLevels = TRUE, exact = nonvar,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

##################################### Supplementary Table 2 hospital and probiotics
probiotics_infant <- select(metadata_rarefied, record_id, hospital_catch, probiotics_before)
probiotics_infant <- probiotics_infant[!duplicated(probiotics_infant$record_id),]
table(probiotics_infant$hospital_catch,probiotics_infant$probiotics_before)

#############################################################################db-rda
##################################### Supplementary Table 3 rb-rda
otu_FT0 <- t(as.data.frame(otu_table(ps_FT0)))
env_FT0 <- subset(metadata_FT0,select=c(fortification, GA, SGA, delivery_mode, DOL, hospital_catch, antibiotics))
rda_FT0 <- capscale(otu_FT0~., env_FT0, dist="bray", add=TRUE)
dbanov_FT0 <- anova.cca(rda_FT0, by="term", permutations = 999)

otu_FT1 <- t(as.data.frame(otu_table(ps_FT1)))
env_FT1 <- subset(metadata_FT1,select=c(fortification, GA, SGA, delivery_mode, DOL, hospital_catch, antibiotics))
rda_FT1 <- capscale(otu_FT1~., env_FT1, dist="bray", add=TRUE)
dbanov_FT1 <- anova.cca(rda_FT1, by="term", permutations = 999)

otu_FT2 <- t(as.data.frame(otu_table(ps_FT2)))
env_FT2 <- subset(metadata_FT2,select=c(fortification, GA, SGA, delivery_mode, DOL, hospital_catch, antibiotics))
rda_FT2 <- capscale(otu_FT2~., env_FT2, dist="bray", add=TRUE)
dbanov_FT2 <- anova.cca(rda_FT2, by="term", permutations = 999)
