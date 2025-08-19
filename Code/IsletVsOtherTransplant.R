# File Info ----
# Author: Jyotirmoy Roy
# Title: Analysis for Islet Transplant Rejection
# Date Created: October 2024
# Info: 
# Ref: 
# Load libraries ----
# install.packages("pacman") # Install pacman if you haven't before
# Install the following packages if you haven't, then load them as follows:
pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)


# MSIGDBR Pathways ----
# Needs msigdbr package: https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
msigdbr_collections() # Take a look at all the pathway groups in the msigdbr database
sets_hallmark <- msigdbr(species="Mus musculus", category="H") # Large df w/ categories
pwl_hallmark <- split(sets_hallmark$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_hallmark$gs_name) # Pathway names
sets_reactome <- msigdbr(species="Mus musculus", subcategory="CP:REACTOME") # Large df w/ categories
pwl_reactome <- split(sets_reactome$gene_symbol, # Genes to split into pathways, by ensembl
                      sets_reactome$gs_name) # Pathway names
kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names
biocarta_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:BIOCARTA") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names
pwl_msigdbr <- c(pwl_hallmark, pwl_kegg) # Compile them all
length(pwl_msigdbr)









# USE CASE ONE ----------------------------------------------------------------- 
# Organize Data
setwd("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch1")
getwd()

#Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch1/Metadata_Batch1.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch2/Metadata_Batch2.csv", sep=",", header=T) # Metadata file
meta_STx_HTx_batch <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Russ Transplant Data/Sequencing Analyses Script/STx_HTx_metadata_Updated.csv", sep=",", header=T) # Metadata file

meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_STx_HTx_batch<-as.data.frame(meta_STx_HTx_batch)
# Merge metadata by columns (i.e., add samples from Batch 2, Batch 1 and STx/HTx)
meta_combined <- rbind(meta_batch1, meta_batch2,meta_STx_HTx_batch)

# Preview the combined metadata
head(meta_combined)


#Counts Data Importing
counts_batch1 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch1/IsTx_gene_expected_count_annot_batch1.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch1 <- na.omit(counts_batch1)

counts_batch2 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch2/IsTx_gene_expected_count_annot_batch2.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch2 <- na.omit(counts_batch2)

counts_STx_HTx <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Russ Transplant Data/Sequencing Analyses Script/Raw_counts_STx_HTx.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_STx_HTx <- na.omit(counts_STx_HTx)
# Install biomaRt if it's not installed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}

# Load biomaRt library
library(biomaRt)

# Connect to Ensembl Biomart for mouse genes
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Extract Ensembl gene IDs from the first column of your counts data
ensembl_gene_ids <- counts_STx_HTx[, 1]

# Retrieve gene symbols corresponding to the Ensembl gene IDs
gene_info <- getBM(
  filters = "ensembl_gene_id", 
  attributes = c("ensembl_gene_id", "mgi_symbol"), 
  values = ensembl_gene_ids, 
  mart = mart
)

# Merge the gene symbols with the original counts data
counts_STx_HTx <- merge(gene_info, counts_STx_HTx, by.x = "ensembl_gene_id", by.y = colnames(counts_STx_HTx)[1])

# Optionally, remove the Ensembl gene IDs and keep only gene symbols
counts_STx_HTx <- counts_STx_HTx[, -1]
colnames(counts_STx_HTx)[1] <- "GeneSymbol"
# Remove any rows with NA values (if necessary)
counts_STx_HTx <- na.omit(counts_STx_HTx)
# View the updated dataframe
head(counts_STx_HTx)


#Remove duplicate names
counts_batch1 <- counts_batch1[!duplicated(counts_batch1[, 1]), ]
genes <- counts_batch1[, 1]
rownames(counts_batch1) <- genes
counts_batch1 <- counts_batch1[, -1]

counts_batch2 <- counts_batch2[!duplicated(counts_batch2[, 1]), ]
genes <- counts_batch2[, 1]
rownames(counts_batch2) <- genes
counts_batch2 <- counts_batch2[, -1]

counts_STx_HTx <- counts_STx_HTx[!duplicated(counts_STx_HTx[, 1]), ]
genes <- counts_STx_HTx[, 1]
rownames(counts_STx_HTx) <- genes
counts_STx_HTx <- counts_STx_HTx[, -1]
#Combine data
# Step 2: Identify common genes between the three datasets
common_genes <- Reduce(intersect, list(rownames(counts_batch1), rownames(counts_batch2), rownames(counts_STx_HTx)))

# Step 3: Subset each dataset to include only the common genes
counts_batch1_common <- counts_batch1[common_genes, ]
counts_batch2_common <- counts_batch2[common_genes, ]
counts_STx_HTx_common <- counts_STx_HTx[common_genes, ]

# Step 4: Combine the datasets by columns
# Assuming you want to combine them as different samples/conditions
combined_counts <- cbind(counts_batch1_common, counts_batch2_common, counts_STx_HTx_common)

# Check the combined result
head(combined_counts)

#### Start Analysis ####

combined_counts <- combined_counts[, meta_combined$Samples] # Put counts & metadata in same order
case1_f1 <- flexiDEG.function1(combined_counts, meta_combined, # Run Function 1
                               convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                               batches = T, quality = T, variance = T,use_pseudobulk = F) # Select filters: 2, 0, 15

# Saving case1_f1 dataframe as a CSV file
#write.csv(case1_f1, file = "Normalized_Counts_ITx_HTx_STx.csv", row.names = TRUE)

# Saving meta_combined dataframe as a CSV file
#write.csv(meta_combined, file = "Combined_Metadata_ITx_HTx_STx.csv", row.names = FALSE)
library(dplyr)

meta_combined <- meta_combined %>%
  mutate(Batch = ifelse(Batch %in% c("Islet Transplant 1", "Islet Transplant 2"), "Islet Transplant", Batch),
         Group = paste(Group, Batch, sep = "_"))


# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_combined)
unique_colSide <- unique(colSide)

# Double Volcano Plot ----
library(dplyr)
# Create separate design matrices for Islet and Other Transplants
meta_combined <- meta_combined %>%
  mutate(Batch = ifelse(Batch %in% c("Islet Transplant 1", "Islet Transplant 2"), "Islet Transplant", Batch))

# A) Biology Agnostic ---- 
getwd() 
sample_set <- "Case1ag" 
currentDate <- Sys.Date()
save_name <- paste(currentDate, sample_set)
case1_f2 <- flexiDEG.function2(case1_f1, meta_combined) # Run Function 2
case1_f3 <- flexiDEG.function3(case1_f2, meta_combined, fdr_cutoff = 1, logfc_cutoff = 2.5) # Run Function 3       ++++ Doesn't seem to be working correctly
# Gene Clustering
# Double Volcano
case1_f4 <- flexiDEG.function4(case1_f2, meta_combined,validation_option = 1) # Run Function 4
ENplots <- flexiDEG.ENplots(case1_f1, case1_f4, colSide, unique_colSide) # Generate PCA plots
ggarrange(plotlist = ENplots, ncol=5, nrow=4) # Plots in 5 cols & 4 rows
# Collect EN results
EN1 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[1]])]), ]) 
EN.95 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[2]])]), ]) 
EN.9 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[3]])]), ]) 
EN.85 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[4]])]), ]) 
EN.8 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[5]])]), ]) 
EN.75 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[6]])]), ]) 
EN.7 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[7]])]), ]) 
EN.65 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[8]])]), ]) 
EN.6 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[9]])]), ]) 
EN.55 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[10]])]), ]) 
EN.5 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[11]])]), ]) 
EN.45 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[12]])]), ]) 
EN.4 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[13]])]), ]) 
EN.35 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[14]])]), ]) 
EN.3 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[15]])]), ]) 
EN.25 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[16]])]), ]) 
EN.2 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[17]])]), ]) 
EN.15 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[18]])]), ]) 
EN.1 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[19]])]), ]) 
EN.05 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[20]])]), ]) 
heatmap.2(as.matrix(EN.9), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,10), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN1), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=colSide)

# B) Pre-hoc Biology ---- 
getwd() 
sample_set <- "Case1ph" 
currentDate <- Sys.Date()
save_name <- paste(currentDate, sample_set)
case1_f1<-as.matrix(case1_f1)# Change dataframe to matrix
is.matrix(case1_f1) # GSVA needs data as matrix w/ genes as rownames; Must be TRUE to proceed

#sets_celltype <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
#wl_celltype <- split(sets_celltype$gene_symbol, # Genes to split into pathways, by ensembl
#                      sets_celltype$gs_name) # Pathway names
pwl_msigdbr <- c(pwl_hallmark, pwl_kegg) # Compile them all

case1_f1<-as.matrix(case1_f1)# Change dataframe to matrix
is.matrix(case1_f1) # GSVA needs data as matrix w/ genes as rownames; Must be TRUE to proceed
gsvapar<-gsvaParam(case1_f1, c(pwl_hallmark), maxDiff=TRUE,minSize=3)
gsva_case1 <- gsva(gsvapar)

gsva_case1<-as.data.frame(gsva_case1) # Convert to dataframe
gsvaf2 <- flexiDEG.function2(gsva_case1, meta_combined) # Run Function 2
dev.off()
heatmap.2(as.matrix(gsva_case1), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,70), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both",
          cexRow = 1.5,  # Increase font size for row labels
          cexCol = 1.5)   # Increase font size for column labels)

#gsvaf3 <- flexiDEG.function3(gsva, meta_case1) # Run Function 3
#genes_cons <- flexiDEG.sharedgenes(pwl_msigdbr, gsvaf3, cutoff = 1) # Conserved genes
# Pathway clustering
# Double volcano
# Differential networks (DNApath)
gsvaf4 <- flexiDEG.function4(gsvaf2, meta_case1,validation_option = 1) # Run Function 4
ENplots <- flexiDEG.ENplots(gsva_case1, gsvaf4, colSide, unique_colSide) # Generate PCA plots
ggarrange(plotlist = ENplots, ncol=5, nrow=4) # Plots in 5 cols & 4 rows
# Collect EN results
EN1 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[1]])]), ])
EN.95 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[2]])]), ])
EN.9 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[3]])]), ])
EN.85 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[4]])]), ])
EN.8 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[5]])]), ])
EN.75 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[6]])]), ])
EN.7 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[7]])]), ])
EN.65 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[8]])]), ])
EN.6 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[9]])]), ])
EN.55 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[10]])]), ])
EN.5 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[11]])]), ])
EN.45 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[12]])]), ])
EN.4 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[13]])]), ])
EN.35 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[14]])]), ])
EN.3 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[15]])]), ])
EN.25 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[16]])]), ])
EN.2 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[17]])]), ])
EN.15 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[18]])]), ])
EN.1 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[19]])]), ])
EN.05 <- na.omit(gsva_case1[unique(rownames(gsva_case1)[as_vector(gsvaf4[[20]])]), ])
heatmap.2(as.matrix(EN.9), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN.1), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=colSide)














# USE CASE TWO ----------------------------------------------------------------- 
# Organize Data
meta_case2 <- read.table(".\\Data\\case2_meta.csv", sep=",", header=T) # Metadata file
meta_case2 <- as.data.frame(meta_case2)
rownames(meta_case2) <- meta_case2$Samples
counts_case2 <- as.data.frame(read.table(".\\Data\\case2_raw.csv", sep=",", header=T)) # Raw counts file
genes <- counts_case2[, 1]
counts_case2 <- counts_case2[, -1]
rownames(counts_case2) <- genes
counts_case2 <- counts_case2[, meta_case2$Samples] # Put counts & metadata in same order
f1 <- flexiDEG.function1(counts_case2, meta_case2, # Run Function 1
                         convert_genes = F, exclude_riken = F, exclude_pseudo = F,
                         batches = T, quality = T, variance = T)
# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_case2)
unique_colSide <- unique(colSide)

# A) Biology Agnostic ---- 
getwd() 
sample_set <- "Case2ag" 
currentDate <- Sys.Date()
save_name <- paste(currentDate, sample_set)
f2 <- flexiDEG.function2(f1, meta_case2) # Run Function 2
f3 <- flexiDEG.function3(f2, meta_case2) # Run Function 3
# Gene Clustering
# Double Volcano (Function 6)
f4 <- flexiDEG.function4(f2, meta_case2, validation_option = 2) # Run Function 4
ENplots <- flexiDEG.ENplots(f1, f4, colSide, unique_colSide) # Generate PCA plots
ggarrange(plotlist = ENplots, ncol=5, nrow=4) # Plots in 5 cols & 4 rows
# Collect EN results
EN1 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[1]])]), ]) 
EN.95 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[2]])]), ]) 
EN.9 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[3]])]), ]) 
EN.85 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[4]])]), ]) 
EN.8 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[5]])]), ]) 
EN.75 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[6]])]), ]) 
EN.7 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[7]])]), ]) 
EN.65 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[8]])]), ]) 
EN.6 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[9]])]), ]) 
EN.55 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[10]])]), ]) 
EN.5 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[11]])]), ]) 
EN.45 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[12]])]), ]) 
EN.4 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[13]])]), ]) 
EN.35 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[14]])]), ]) 
EN.3 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[15]])]), ]) 
EN.25 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[16]])]), ]) 
EN.2 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[17]])]), ]) 
EN.15 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[18]])]), ]) 
EN.1 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[19]])]), ]) 
EN.05 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[20]])]), ]) 
heatmap.2(as.matrix(bc_EN.8), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(bc_EN.8), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=colSide)

# B) Pre-hoc Biology ---- 
getwd() 
sample_set <- "Case2ph" 
currentDate <- Sys.Date()
save_name <- paste(currentDate, sample_set)
is.matrix(f1) # GSVA needs data as matrix w/ genes as rownames; Must be TRUE to proceed
gsva_case2 <- gsva(f1, c(pwl_msigdbr), method = "gsva",
                   kcdf = "Gaussian", # "Gaussian" for continuous counts; "Poisson" integers
                   min.sz = 15, max.sz = 500, # Min & max gene set size
                   mx.diff = TRUE, # Compute Gaussian-distributed scores
                   verbose = TRUE) # Progress bar
gsvaf2 <- flexiDEG.function2(gsva_case2, meta_case2) # Run Function 2
gsvaf3 <- flexiDEG.function3(gsvaf2, meta_case2) # Run Function 3
genes_cons <- flexiDEG.sharedgenes(pwl_msigdbr, gsvaf3, cutoff = 1) # Conserved genes
# Pathway clustering
# Double volcano
# Differential networks (DNApath)
gsvaf4 <- flexiDEG.function4(gsvaf3, meta_case2) # Run Function 4
ENplots <- flexiDEG.ENplots(gsva_case2, gsvaf4, colSide, unique_colSide) # Generate PCA plots
ggarrange(plotlist = ENplots, ncol=5, nrow=4) # Plots in 5 cols & 4 rows
# Collect EN results
EN1 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[1]])]), ])
EN.95 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[2]])]), ])
EN.9 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[3]])]), ])
EN.85 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[4]])]), ])
EN.8 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[5]])]), ])
EN.75 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[6]])]), ])
EN.7 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[7]])]), ])
EN.65 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[8]])]), ])
EN.6 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[9]])]), ])
EN.55 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[10]])]), ])
EN.5 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[11]])]), ])
EN.45 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[12]])]), ])
EN.4 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[13]])]), ])
EN.35 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[14]])]), ])
EN.3 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[15]])]), ])
EN.25 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[16]])]), ])
EN.2 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[17]])]), ])
EN.15 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[18]])]), ])
EN.1 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[19]])]), ])
EN.05 <- na.omit(gsva_case2[unique(rownames(gsva_case2)[as_vector(gsvaf4[[20]])]), ])
heatmap.2(as.matrix(EN.8), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN.8), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=colSide)





# USE CASE THREE ----------------------------------------------------------------- 
# Organize Data
setwd("/Users/lonnieshea/Desktop/Jyotirmoy/ElasticNetPackage/Data/")
meta_case3 <- read.table("case3_meta.csv", sep=",", header=T) # Metadata file
meta_case3 <- as.data.frame(meta_case3)
rownames(meta_case3) <- meta_case3$Samples
counts_case3 <- as.data.frame(read.table("case3_raw.csv", sep=",", header=T)) # Raw counts file
genes <- counts_case3[, 1]
counts_case3 <- counts_case3[, -1]
rownames(counts_case3) <- genes
counts_case3 <- counts_case3[, meta_case3$Samples] # Put counts & metadata in same order
f1 <- flexiDEG.function1(counts_case3, meta_case3, # Run Function 1
                         convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                         batches = T, quality = T, variance = T)
# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_case3)
unique_colSide <- unique(colSide)

# A) Biology Agnostic ---- 
getwd() 
sample_set <- "Case3ag" 
currentDate <- Sys.Date()
save_name <- paste(currentDate, sample_set)
f2 <- flexiDEG.function2(f1, meta_case3) # Run Function 2
f3 <- flexiDEG.function3(f2, meta_case3) # Run Function 3
# Gene Clustering
# Double Volcano (Function 6)
f4 <- flexiDEG.function4(f2, meta_case3,validation_option = 2) # Run Function 4

EN1 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[1]])]), ])
EN2 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[2]])]), ])
EN3 <- na.omit(f1[unique(rownames(f1)[as_vector(f4[[3]])]), ])

png(paste("Heatmap_BiologyAgonistic", names(f4)[1], ".png", sep=""), width = 800, height = 800)
heatmap(as.matrix(EN1), scale="row", col=coul, key= T, xlab="", ylab="", ColSideColors=colSide, trace="none",
        key.title=NA, 
        key.ylab=NA, keysize=0.3, dendrogram="both",
        main=paste("Alpha:", names(f4)[1]))
dev.off()

png(paste("Heatmap_BiologyAgonistic", names(f4)[2], ".png", sep=""), width = 800, height = 800)
heatmap(as.matrix(EN2), scale="row", col=coul, key= T, xlab="", ylab="", 
        margins=c(7,7), ColSideColors=colSide, trace="none", key.title=NA, 
        key.ylab=NA, keysize=0.8, dendrogram="both",
        main=paste("Alpha:", names(f4)[2]))
dev.off()

png(paste("Heatmap_BiologyAgonistic", names(f4)[3], ".png", sep=""), width = 800, height = 800)
heatmap(as.matrix(EN3), scale="row", col=coul, key= T, xlab="", ylab="", 
        margins=c(7,7), ColSideColors=colSide, trace="none", key.title=NA, 
        key.ylab=NA, keysize=0.8, dendrogram="both",
        main=paste("Alpha:", names(f4)[3]))
dev.off()
# Create ggbiplot plot with automated color scheme
install_github("vqv/ggbiplot")
library(ggbiplot)

plot<- ggbiplot(prcomp(t(EN1), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, var.scale=1, circle=T) + 
  theme_classic() + geom_point(size=3, color=colSide) + 
  scale_color_manual(name="Group", values=unique_colSide)+
  labs(title = paste("Alpha:", names(f4)[1]))

# Convert ggbiplot plot to a gtable object
plot_gtable <- ggplotGrob(plot)
# Create a new plot with the PDF device
grid.newpage()
grid.draw(plot_gtable)
getwd()
ggsave(paste("PCA_BiologyAgonistic_", names(f4)[1], ".png", sep=""), plot, width = 6, height = 6, units = "in", dpi = 300)



plot<- ggbiplot(prcomp(t(EN2), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, var.scale=1, circle=T) + 
  theme_classic() + geom_point(size=3, color=colSide) + 
  scale_color_manual(name="Group", values=unique_colSide)+
  labs(title = paste("Alpha:", names(f4)[2]))

# Convert ggbiplot plot to a gtable object
plot_gtable <- ggplotGrob(plot)
# Create a new plot with the PDF device
grid.newpage()
grid.draw(plot_gtable)
ggsave(paste("PCA_BiologyAgonistic_", names(f4)[2], ".png", sep=""), plot, width = 6, height = 6, units = "in", dpi = 300)

plot<- ggbiplot(prcomp(t(EN3), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, var.scale=1, circle=T) + 
  theme_classic() + geom_point(size=3, color=colSide) + 
  scale_color_manual(name="Group", values=unique_colSide)+
  labs(title = paste("Alpha:", names(f4)[3]))

# Convert ggbiplot plot to a gtable object
plot_gtable <- ggplotGrob(plot)
# Create a new plot with the PDF device
grid.newpage()
grid.draw(plot_gtable)
ggsave(paste("PCA_BiologyAgonistic_", names(f4)[3], ".png", sep=""), plot, width = 6, height = 6, units = "in", dpi = 300)


# B) Pre-hoc Biology ---- 
getwd() 
sample_set <- "Case3ph" 
currentDate <- Sys.Date()
save_name <- paste(currentDate, sample_set)
library(org.Hs.eg.db)
Mm <- org.Hs.eg.db
my.symbols <- rownames(f1)


# Subset the Mm data frame based on the selected columns
#Table <- Mm[, c("ENTREZID", "SYMBOL")]

Table<-AnnotationDbi::select(Mm, 
                             keys = my.symbols,
                             columns = c("ENSEMBL", "SYMBOL"),
                             keytype = "SYMBOL")

Table <- Table %>%
  distinct(SYMBOL, .keep_all = TRUE)
# Filter columns based on non-NA entrezgene values
filtered_rows <- my.symbols[!is.na(Table$ENSEMBL)]

# Convert character column names to numeric indicesd
filtered_indices <- match(filtered_rows, rownames(f1))

# Create a new data frame with only selected columns
filtered_f1 <- f1[c(filtered_indices),]

# Rename columns to entrezgene IDs
rownames(filtered_f1) <- Table$ENSEMBL[!is.na(Table$ENSEMBL)]



f1_m<-as.matrix(filtered_f1)
is.matrix(f1_m) # GSVA needs data as matrix w/ genes as rownames; Must be TRUE to proceed


gsva_case3 <- gsva(f1_m, c(pwl_msigdbr), method = "gsva",
                   kcdf = "Gaussian", # "Gaussian" for continuous counts; "Poisson" integers
                   min.sz = 15, max.sz = 500, # Min & max gene set size
                   mx.diff = TRUE, # Compute Gaussian-distributed scores
                   verbose = TRUE) #3Progress bar
gsva_case3<- as.data.frame(gsva_case3)
gsvaf2 <- flexiDEG.function2(gsva_case3, meta_case3) # Run Function 2
gsvaf3 <- flexiDEG.function3(gsvaf2, meta_case3) # Run Function 3
genes_cons <- flexiDEG.sharedgenes(pwl_msigdbr, gsvaf2, cutoff = 1) # Conserved genes
# Pathway clustering
# Double volcano
# Differential networks (DNApath)
gsvaf4 <- flexiDEG.function4(gsvaf2, meta_case3) # Run Function 4
# Create ggbiplot plot with automated color scheme
install_github("vqv/ggbiplot")
library(ggbiplot)

for (i in 1:length(gsvaf4)) {
  # Generate ENi
  EN <- na.omit(gsva_case3[unique(rownames(gsva_case3)[as_vector(gsvaf4[[i]])]), ])
  # Plot heatmap
  png(paste("Heatmap_PreHocBiology", names(gsvaf4)[i], ".png", sep=""), width = 800, height = 800)
  heatmap(as.matrix(EN), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7, 40), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both",
          main=paste("Alpha:", names(gsvaf4)[i]))
  dev.off()
  # Plot PCA
  plot <- ggbiplot(prcomp(t(EN), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, var.scale=1, circle=T) + 
    theme_classic() + geom_point(size=3, color=colSide) + 
    scale_color_manual(name="Group", values=unique_colSide) +
    labs(title = paste("Alpha:", names(gsvaf4)[i]))
  # Convert ggbiplot plot to a gtable object
  plot_gtable <- ggplotGrob(plot)
  # Save PCA plot
  ggsave(paste("PCA_PreHocBiology_", names(gsvaf4)[i], ".png", sep=""), plot, width = 6, height = 6, units = "in", dpi = 300)
}

# Use Genes from PreHoc Biology to Run Elastic Net For Alpha=1
EN_1 <- na.omit(gsva_case3[unique(rownames(gsva_case3)[as_vector(gsvaf4[[3]])]), ])
genes_cons_prehoc <- flexiDEG.sharedgenes(pwl_msigdbr, EN_1, cutoff = 5)

f1_m<-as.data.frame(f1_m)
f2<-na.omit(f1_m[(genes_cons_prehoc), ])
f6 <- flexiDEG.function4(f2, meta_case3) # Run Function 4
# Create ggbiplot plot with automated color scheme
