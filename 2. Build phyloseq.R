setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("1. Setup.R")
########################################################################## build phyloseq

otu <- read.table(file = "data/zOTU_table_GG.txt", sep = "\t", header = T, row.names = 1, comment.char = "")
metadata <-mapping_sample <- read.table("data/mapping_sample.txt", sep = "\t", header = T, row.names = 1, comment.char = "")

# remove unclassified
otu <- otu[otu$taxonomy != "Unassigned", ]

taxonomy <- otu[, "taxonomy", drop = FALSE]
otu <- subset(otu, select = -taxonomy)

# clean the taxonomy, rm [, ], [a-z]__
taxonomy$taxonomy <- gsub("[a-z]__", "", taxonomy$taxonomy)
taxonomy$taxonomy <- gsub("\\[|\\]", "", taxonomy$taxonomy)
taxonomy$taxonomy <- gsub("\\ ", "", taxonomy$taxonomy)
tax <- taxonomy %>%
  dplyr::select(taxonomy) %>%
  separate(taxonomy,
           c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           ";") %>%
  mutate_all(na_if, "") # non-value string to NA

# complete species name
tax$Species <- mapply(
  function(x, y) ifelse(! is.na(y), paste(x, y, sep = " "), y),
  tax$Genus, tax$Species)

for (i in 1:nrow(tax)) {
  if (is.na(tax[i, 2])) {
    kingdom <- paste("Unclassified ", tax[i, 1], sep = " ")
    tax[i, 2:7] <- kingdom
  } else if (is.na(tax[i, 3])) {
    phylum <- paste("Unclassified ", tax[i, 2], sep = " ")
    tax[i, 3:7] <- phylum
  } else if (is.na(tax[i, 4])) {
    class <- paste("Unclassified ", tax[i, 3], sep = " ")
    tax[i, 4:7] <- class
  } else if (is.na(tax[i, 5])) {
    order <- paste("Unclassified ", tax[i, 4], sep = " ")
    tax[i, 5:7] <- order
  } else if (is.na(tax[i, 6])) {
    family <- paste("Unclassified ", tax[i, 5], sep = " ")
    tax[i, 6:7] <- family
  } else if (is.na(tax[i, 7])) {
    tax$Species[i] <- paste("Unclassified ", tax$Genus[i], sep = " ")
  }
}

# build phyloseq object
OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
SAMPLE <- sample_data(metadata)
TREE <- read_tree_greengenes("./data/zOTU.tree")
## merge the data
ps <- phyloseq(OTU, TAX, SAMPLE, TREE)
