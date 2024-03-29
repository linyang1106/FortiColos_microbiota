setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("3. Demography.R")
######################################################################## pH, data for figure 4b ##########################################################
lm_ph0 <- lm(pH~fortification + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, metadata_FT0)
ph0 <- as.data.frame(tidy(lm_ph0))
lm_ph1 <- lm(pH~fortification + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, metadata_FT1)
ph1 <- as.data.frame(tidy(lm_ph1))
lm_ph2 <- lm(pH~fortification + delivery_mode + hospital_catch + GA + DOL + SGA + antibiotics, metadata_FT2)
ph2 <- as.data.frame(tidy(lm_ph2))
fortification <- data.frame(group1 = c("BC","BC","BC"), group2=c("CF", "CF", "CF"))
ph_p <- as.data.frame(rbind(ph0[2,],ph1[2,],ph2[2,]))
ph_p <- ph_p[,-1]
row.names(ph_p) <- c("FT0","FT1","FT2")
ph_p <- cbind(fortification,ph_p)
print(ph_p, quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

################################################################ pH correlation, data for figure 4c-e ####################################################
# Rhea correlation, ref. https://github.com/Lagkouvardos/Rhea
################### data preparation ##########################
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
#combine the file for correlation use
mapping <- data.frame(sample_data(ps_rarefied))
correlationRaw <- as.data.frame(cbind(rownames(metadata_rarefied),metadata_rarefied$pH,metadata_rarefied$fortification))
rownames(correlationRaw) <- correlationRaw[,1]
correlationRaw <- subset(correlationRaw, select=-V1)
colnames(correlationRaw)[1]<-"pH"
colnames(correlationRaw)[2]<-"fortification"
head(mapping)
head(m_Tab)
head(correlationRaw)
# only use the sample with correlation data
mapping <- mapping[rownames(mapping) %in% rownames(correlationRaw),]
# match the mapping table and the correlation data
index <- match(rownames(mapping), rownames(correlationRaw))
mapping_merge <- correlationRaw[index,]
rownames(mapping_merge) <- rownames(mapping)
#match the microbial table
m_Tab <- t(m_Tab)
index <- match(rownames(mapping_merge),rownames(m_Tab))
meta_tab <- cbind(mapping_merge, m_Tab[index,])
head(meta_tab)
write.table(meta_tab, "data/correlation_meta_Genusgroup.tab", sep = "\t", col.names = NA, quote = F)
################### Pearson correlation ##########################
#ref. https://github.com/Lagkouvardos/Rhea
setwd("data")                     #<--- CHANGE ACCORDINGLY !!!
#' Please give the file name of the table containing the variables for analysis
input_file <- "correlation_meta_Genusgroup.tab"              #<--- CHANGE ACCORDINGLY !!!
#' Please give the position where the taxonomic variables (OTUs or taxonomic groups) start!!
otu_variables_start <- 2                                        #<--- CHANGE ACCORDINGLY !!!
# Set the cutoff for significance
signf_cutoff <- 0.05
# Calculate correlation among taxonomic variables
# 1 = calculate correlations within OTUs or taxa;
# 0 = NO test within taxonomic variables
includeTax <- 0
# Calculate correlation among meta-variables
# 1 = calculate correlations within meta-variables
# 0 = NO test within meta-variables
includeMeta <- 0
# Handling of missing values for meta-variables
# 1 = missing values are filled with the mean for the corresponding variable
# 0 = NO imputation (replacing missing data with substituted values)
fill_NA <- 0
# Treat zeros in taxonomic variables as missing values
# 1 = Consider taxonomic zeros as missing values
# 0 = Keep zeros for the calculation of correlations
replace_zeros <- 1
# Set a cutoff for the minimum number of values (prevalence) for a given taxonomic variable to be considered for calculation
prevalence_exclusion <- 0.3
# Set a cutoff for the minimal number of pairs observations required for calculation of correlations
min_pair_support <- 200
# Set a significance cutoff for graphical output
plot_pval_cutoff <- 0.05
# Set a correlation coefficient cutoff for graphical output
plot_corr_cutoff <- 0.5

###### Main Script ######
# Check if required packages are already installed, and install if missing
packages <-c("Hmisc","corrplot") 
# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack)
  } 
}
# Applying the installation on the list of packages
lapply(packages, InsPack)
# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)
# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))

###################Read input table####################
# Load the tab-delimited file containing the values to be checked (rownames in the first column)
my_data <-
  read.table (
    file = input_file,
    check.names = FALSE,
    header = TRUE,
    dec = ".",
    sep = "\t",
    row.names = 1,
    comment.char = ""
  )
# Clean table from empty lines
my_data <- my_data[!apply(is.na(my_data) | my_data=="",1,all),]

#################### Functions #####################
# Function for filling missing values with the mean of the column
fill_NA.mean <- function(vec)
{
  # Calculate mean value of each column (excluding missing values)
  m <- mean(vec, na.rm = TRUE)
  # Replace missing values with mean
  vec[is.na(vec)] <- m
  # Return the new input data frame
  return(vec)
}
# Function to logarithmically normalized OTU values
log_ratio <- function(data)
{
  # Compute the logarithmus
  log_data <- log(data)
  # Calculate exponential function of column-wise mean values for finite log transformed data
  gm <- exp(mean(log_data[is.finite(log_data)]))
  # Compute the logarithmus
  log_gm <- log(gm)
  # Take the difference of both log-transformed datasets
  data <- log_data - log_gm
  # Return the new OTU table
  return(data)
}

###################### MAIN PROGRAM #####################
my_data <- as.data.frame(apply(my_data,2,as.numeric))
first_OTU <- colnames(my_data)[otu_variables_start]
# Split the meta and taxonomic parts of the table
# Choose the continuous scaled variables
my_meta_data <- my_data[1:otu_variables_start - 1]
# Choose the taxonomic variables
my_otu_data <- my_data[otu_variables_start:dim(my_data)[2]]
# Process the meta measurements according to selection
if (fill_NA == 0) {
  # Do not do anything, just rename the file
  my_meta_fixed =  my_meta_data
}
if (fill_NA == 1) {
  # Fill non-zero missing meta-values with the mean of the column (optional)
  # Apply the previously implemented function 'fill_Na.mean' to the meta-data subset
  my_meta_fixed = apply(my_meta_data, 2, fill_NA.mean)
}
# The maximal number of absence
prevalence_cutoff <- dim(my_otu_data)[1] - (prevalence_exclusion*dim(my_otu_data)[1])
# Count how many missing values are found for each OTU 
na_count <-sapply(my_otu_data, function(y) sum(length(which(is.na(y)))))
# Count how many zeros are found for each OTU 
zero_count <-sapply(my_otu_data, function(y) sum(length(which(y==0))))
prevalence_count <- na_count + zero_count
# A new OTU-table is generated, where the number of missing values is below the set cutoff
my_otu_data <- my_otu_data[, prevalence_count <= prevalence_cutoff ]
# If the parameter is set, zeros are replaced with missing values
if (replace_zeros == 1) {
  my_otu_data[my_otu_data==0] <- NA
}
# Replace zeros with 0.0001 to avoid infinite number when calculating logarithmus (log(0)=-INF)
my_otu_data[my_otu_data==0] <- 0.0001
# Transform compositional data by log ratio transformation
my_otu_fixed = apply(my_otu_data, 2, log_ratio)
# Merge the meta- and OTU data in one table
transformed_data <- cbind(my_meta_fixed, my_otu_fixed)
# Centre and scale the values
my_scaled_data <- scale(transformed_data, center = TRUE, scale = TRUE)
# Calculate all pairwise correlations using Pearson correlation method
my_rcorr <- rcorr(as.matrix(my_scaled_data), type = "pearson")
# Generate vector with variable names
var_names <- row.names(my_rcorr$r)
# Depending on which parameters were set at the beginning, one query type is selected
# In each query type, three matrices are generated: p-value matrix, correlation matrix, support matrix
# All possibles pairs are saved in a vector (pairs)
if(includeTax==1 & includeMeta==0){
  # Correlation among OTUs and NO correlation among meta-variables
  row_names <- var_names[c(otu_variables_start:dim(my_rcorr$r)[1])]
  col_names <- var_names
  pairs <-expand.grid(row_names, col_names)
  my_cor_matrix <- my_rcorr$r[c(otu_variables_start:dim(my_rcorr$r)[1]),]
  my_pvl_matrix <-my_rcorr$P[c(otu_variables_start:dim(my_rcorr$P)[1]),]
  my_num_matrix <- my_rcorr$n[c(otu_variables_start:dim(my_rcorr$n)[1]),]
  # Set variable for plotting
  diagonale=0
} else if(includeTax==1 & includeMeta==1){
  # Correlation among OTUs and correlation among meta-variables
  row_names <-var_names
  col_names <- var_names
  pairs <-expand.grid(row_names, col_names)
  my_cor_matrix <- my_rcorr$r
  my_pvl_matrix <-my_rcorr$P
  my_num_matrix <- my_rcorr$n
  # Set variable for plotting
  diagonale=0
} else if (includeTax==0 & includeMeta==1) {
  # NO correlation among OTUs and correlation among meta-variables
  row_names <-var_names[c(1:(otu_variables_start - 1))]
  col_names <- var_names
  pairs <-expand.grid(row_names, col_names)
  my_cor_matrix <- my_rcorr$r[c(1:(otu_variables_start - 1)),]
  my_pvl_matrix <-my_rcorr$P[c(1:(otu_variables_start - 1)),]
  my_num_matrix <- my_rcorr$n[c(1:(otu_variables_start - 1)),]
  # Set variable for plotting
  diagonale=1
} else {
  # NO correlation among OTUs and NO correlation among meta-variables
  row_names <- var_names[c(1:(otu_variables_start - 1))]
  col_names <- var_names[otu_variables_start:dim(my_rcorr$r)[1]]
  pairs <-expand.grid(row_names, col_names)
  my_cor_matrix <- my_rcorr$r[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$r)[1])]
  my_pvl_matrix <-my_rcorr$P[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$P)[1])]
  my_num_matrix <- my_rcorr$n[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$n)[1])]
  # Set variable for plotting
  diagonale=1
}
# Select the corresponding p-value for each pair 
p_vector <- as.vector(my_pvl_matrix)
# Select the corresponding correlation coefficient for each pair 
c_vector <- as.vector(my_cor_matrix)
# Select the corresponding number of observations for each pair 
n_vector <- as.vector(my_num_matrix)
# Generate matrix with the pairwise comparisons
my_pairs <-
  matrix(ncol = 5,
         c(
           as.character(pairs[, 2]),
           as.character(pairs[, 1]),
           c_vector,
           p_vector,
           n_vector
         ))
# Delete all pairs with insufficient number of pairs
my_pairs <- subset(my_pairs, as.numeric(my_pairs[,5]) > min_pair_support)
# Adjust p-value for multiple testing using the Benjamin-Hochberg method
pVal_BH <- round(p.adjust(my_pairs[,4], method = "BH"), 4)
# Add the corrected p-value in the table 
my_pairs <- cbind(my_pairs,as.numeric(pVal_BH))
# Remove similar pairs (values along the diagonal)
my_pairs <- my_pairs[!as.character(my_pairs[, 1]) == as.character(my_pairs[, 2]),]
# Remove duplicate pairs
my_pairs <- my_pairs[!duplicated(my_pairs[, 3]), ]
# Created matrix columns represent correlation coefficients, p-values, number of observations, and corrected p-values
# Rows represent the pairs
matrix_names <- list(c(rep("",times=dim(my_pairs)[1])),
                     c(
                       "variable1",
                       "variable2",
                       "correlation",
                       "pValue",
                       "support",
                       "Corrected"
                     ))
dimnames(my_pairs) <- matrix_names
# Create subset of pairs with significant p-values
my_pairs_cutoff <- my_pairs[as.numeric(my_pairs[, 4]) <= signf_cutoff, ]
# Convert to matrix
my_pairs_cutoff <- matrix(my_pairs_cutoff,ncol=6,dimnames = list(c(rep("",times=dim(my_pairs_cutoff)[1])),c("variable1","variable2","correlation","pvalue","support","corrected pvalue")))
# Create subset of significant pairs with strong correlation (above 0.5)
my_pairs_cutoff_corr <- my_pairs_cutoff[abs(as.numeric(my_pairs_cutoff[, 3])) >= 0, ]
# Remove columns containing no information
my_cor_matrix <- my_cor_matrix[, colSums(is.na(my_cor_matrix)) != nrow(my_cor_matrix)]
# Missing values in the correlation matrix are set to zero
my_cor_matrix[is.na(my_cor_matrix)] <- 0
########################## result ###########################
my_pairs_cutoff
