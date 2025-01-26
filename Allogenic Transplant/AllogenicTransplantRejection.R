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

library(dplyr)
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
pwl_msigdbr <- c(pwl_hallmark, pwl_reactome, pwl_kegg, pwl_biocarta) # Compile them all
length(pwl_msigdbr)








# 1. Load the Data ----
# Organize Data
setwd("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch1")
getwd()

#Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch1/Metadata_Batch1.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch2/Metadata_Batch2.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch3/Metadata_Batch3.csv", sep=",", header=T) # Metadata file

meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)

# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3)

# Preview the combined metadata
head(meta_combined)
unique(meta_combined$Group)

#Counts Data Importing
counts_batch1 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch1/IsTx_gene_expected_count_annot_batch1.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch1 <- na.omit(counts_batch1)

counts_batch2 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch2/IsTx_gene_expected_count_annot_batch2.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch2 <- na.omit(counts_batch2)

counts_batch3 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch3/IsTx_gene_expected_count_annot_batch3.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch3 <- na.omit(counts_batch3)

#Remove duplicate names
counts_batch1 <- counts_batch1[!duplicated(counts_batch1[, 1]), ]
genes <- counts_batch1[, 1]
rownames(counts_batch1) <- genes
counts_batch1 <- counts_batch1[, -1]

counts_batch2 <- counts_batch2[!duplicated(counts_batch2[, 1]), ]
genes <- counts_batch2[, 1]
rownames(counts_batch2) <- genes
counts_batch2 <- counts_batch2[, -1]

counts_batch3 <- counts_batch3[!duplicated(counts_batch3[, 1]), ]
genes <- counts_batch3[, 1]
rownames(counts_batch3) <- genes
counts_batch3 <- counts_batch3[, -1]

#Combine data
# First merge counts_batch1 and counts_batch2
combined_counts <- merge(counts_batch1, counts_batch2, by = "row.names", all = TRUE)
# Rename the Row.names column back
rownames(combined_counts) <- combined_counts$Row.names
combined_counts$Row.names <- NULL
# Then merge the result with counts_batch3
combined_counts <- merge(combined_counts, counts_batch3, by = "row.names", all = TRUE)
# Rename the Row.names column back
rownames(combined_counts) <- combined_counts$Row.names
combined_counts$Row.names <- NULL

# Identify the samples to keep (not "Technical Rejection")
samples_to_keep <- meta_combined$Sample[meta_combined$Group != "Technical Rejection"]

# Filter the meta_combined data
meta_combined <- meta_combined[meta_combined$Group != "Technical Rejection", ]

# Filter the combined_counts data to keep only columns corresponding to samples_to_keep
combined_counts <- combined_counts[, colnames(combined_counts) %in% samples_to_keep]


# Preview the combined dataset
head(combined_counts)


# 2.Preprocessing- DESeq ----

# Select columns in combined_counts that match the remaining sample names in meta_combined
combined_counts <- combined_counts[, meta_combined$Samples]  # Ensure Sample_IDs match column names in combined_counts

case1_f1 <- flexiDEG.function1(combined_counts, meta_combined, # Run Function 1
                               convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                               batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 2, 0, 15
rows_to_remove <- grep("^Gm[0-9]", rownames(case1_f1))
# Remove those rows from case1_f1
case1_f1 <- case1_f1[-rows_to_remove, ]
# meta_combined <- meta_combined %>%
#   mutate(Group = case_when(
#     Group %in% c("Healthy_7", "Healthy_14") ~ "Early",
#     Group == "Healthy_28" ~ "Intermediate",
#     Group %in% c("Healthy_42", "Healthy_56") ~ "Late",
#     TRUE ~ Group  # Keep other values unchanged
#   ))
# Saving case1_f1 dataframe as a CSV file
write.csv(case1_f1, file = "/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Combined_Raw_Counts_ITx_Filtered.csv", row.names = TRUE)

# Saving meta_combined dataframe as a CSV file
write.csv(meta_combined, file = "/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Combined_Metadata_ITx.csv", row.names = FALSE)

# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_combined)
unique_colSide <- unique(colSide)

#3.DESEQ Analysis-Rejection with Anti-CD40L vs Without Anti-CD40L ----
unique(meta_combined$Group)
# Identify the samples to keep (not "Technical Rejection")
samples_to_keep <- meta_combined$Sample[meta_combined$Group != "Tolerance"]

# Filter the meta_combined data
meta_combined_rejected <- meta_combined[meta_combined$Group != "Tolerance", ]

# Filter the combined_counts data to keep only columns corresponding to samples_to_keep
rejected_counts <- case1_f1[, colnames(case1_f1) %in% samples_to_keep]

# DESeq2 analysis for Rejection Anti-CD40L with vs withour--> Cannot Do Batch
dds_rejected <- DESeqDataSetFromMatrix(countData = rejected_counts, colData = meta_combined_rejected, design = ~ Group)
dds_rejected  <- DESeq(dds_rejected)
results_rejected <- as.data.frame(results(dds_rejected, contrast = c("Group", "Rejection(No Anti-CD40L)", "Rejection(Anti-CD40L)")))

# Load required libraries
library(ggplot2)
library(ggrepel)

# Replace NA in pval and padj with 1
results_rejected$pvalue[is.na(results_rejected$pvalue)] <- 1
results_rejected$padj[is.na(results_rejected$padj)] <- 1

# Add -log10(p-value)
results_rejected$logP <- -log10(results_rejected$pvalue)

# Assign significance
results_rejected$Significance <- ifelse(
  results_rejected$padj < 0.05 & results_rejected$log2FoldChange > 1, "Upregulated",
  ifelse(
    results_rejected$padj < 0.05 & results_rejected$log2FoldChange < -1, "Downregulated",
    "Not Significant"
  )
)
#unique(results_rejected$Significance)
# Set custom colors for the plot
# Saving meta_combined dataframe as a CSV file
write.csv(results_rejected, file = "/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/DEG_RejectedControlVsAntiCD40L_ITx.csv", row.names = FALSE)

# Subset data for labeling (most significant genes)
label_data <- subset(results_rejected, padj < 0.05 & abs(log2FoldChange) > 1)

# Set custom colors for the plot
custom_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")

# Create the volcano plot
volcano_plot <- ggplot(results_rejected, aes(x = log2FoldChange, y = logP)) +
  # Add points with different colors for significance
  geom_point(aes(color = Significance), size = 3, alpha = 0.8) +
  # Add labels for the most significant points
  geom_text_repel(
    data = label_data,
    aes(label = rownames(label_data)),
    size = 3,
    max.overlaps = 20
  ) +
  # Customize colors
  scale_color_manual(values = custom_colors) +
  # Add vertical and horizontal lines for thresholds
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # Customize axis labels and title
  labs(
    title = "Rejection Control vs Anti-CD40L",
    x = expression(Log[2] ~ "Fold Change"),  # Subscript for Log2
    y = expression(-Log[10] ~ "(p-value)"),  # Subscript for Log10
    color = "Significance"
  ) +
  # Customize theme for publication quality
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  # Set axis limits (optional)
  coord_cartesian(xlim = c(-6, 6), ylim = c(0, max(results_rejected$logP, na.rm = TRUE)))

# Display the plot
print(volcano_plot)

#4.DESEQ Analysis-Rejection with Anti-CD40L vs Tolerance Anti-CD40L ----
unique(meta_combined$Group)
# Identify the samples to keep (not "Technical Rejection")
samples_to_keep <- meta_combined$Sample[meta_combined$Group != "Rejection(No Anti-CD40L)"]

# Filter the meta_combined data
meta_combined_tol <- meta_combined[meta_combined$Group != "Rejection(No Anti-CD40L)", ]

# Filter the combined_counts data to keep only columns corresponding to samples_to_keep
tol_counts <- case1_f1[, colnames(case1_f1) %in% samples_to_keep]

# DESeq2 analysis for Rejection Anti-CD40L vs Tolerance
dds_tol <- DESeqDataSetFromMatrix(countData = tol_counts, colData = meta_combined_tol, design = ~ Batch+Group)
dds_tol  <- DESeq(dds_tol)
results_tol <- as.data.frame(results(dds_tol, contrast = c("Group", "Rejection(Anti-CD40L)", "Tolerance")))

# Load required libraries
library(ggplot2)
library(ggrepel)

# Replace NA in pval and padj with 1
results_tol$pvalue[is.na(results_tol$pvalue)] <- 1
results_tol$padj[is.na(results_tol$padj)] <- 1

# Add -log10(p-value)
results_tol$logP <- -log10(results_tol$pvalue)

# Assign significance
results_tol$Significance <- ifelse(
  results_tol$padj < 0.05 & results_tol$log2FoldChange > 1, "Upregulated",
  ifelse(
    results_tol$padj < 0.05 & results_tol$log2FoldChange < -1, "Downregulated",
    "Not Significant"
  )
)
unique(results_tol$Significance)
# Set custom colors for the plot
write.csv(results_tol, file = "/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/DEG_RejectedVsTolerance_WithAntiCD40L_ITx.csv", row.names = FALSE)

# Subset data for labeling (most significant genes)
label_data <- subset(results_tol, padj < 0.05 & abs(log2FoldChange) > 1)

# Set custom colors for the plot
custom_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")

# Create the volcano plot
volcano_plot <- ggplot(results_tol, aes(x = log2FoldChange, y = logP)) +
  # Add points with different colors for significance
  geom_point(aes(color = Significance), size = 3, alpha = 0.8) +
  # Add labels for the most significant points
  geom_text_repel(
    data = label_data,
    aes(label = rownames(label_data)),
    size = 3,
    max.overlaps = 30
  ) +
  # Customize colors
  scale_color_manual(values = custom_colors) +
  # Add vertical and horizontal lines for thresholds
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # Customize axis labels and title
  labs(
    title = "Rejection vs Tolerance",
    x = expression(Log[2] ~ "Fold Change"),  # Subscript for Log2
    y = expression(-Log[10] ~ "(p-value)"),  # Subscript for Log10
    color = "Significance"
  ) +
  # Customize theme for publication quality
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  # Set axis limits (optional)
  coord_cartesian(xlim = c(-8, 4), ylim = c(0, max(results_tol$logP, na.rm = TRUE)))

# Display the plot
print(volcano_plot)
unique(meta_combined$Day)

#5. Preprocesing for Elastic Net---- 
case1_EN <- flexiDEG.function1(combined_counts, meta_combined, # Run Function 1
                               convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                               batches = T, quality = T, variance = F,use_pseudobulk = F) # Select filters: 2, 0, 15

rows_to_remove <- grep("^Gm[0-9]", rownames(case1_EN))
# Remove those rows from case1_f1
case1_EN <- case1_EN[-rows_to_remove, ]
# meta_combined <- meta_combined %>%
#   mutate(Group = case_when(
#     Group %in% c("Healthy_7", "Healthy_14") ~ "Early",
#     Group == "Healthy_28" ~ "Intermediate",
#    
case2_EN <- flexiDEG.function2(case1_EN, meta_combined) # Run Function 2

dev.off()
heatmap.2(as.matrix(case2_EN), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both",  cexRow = 1.2,  # Increase the size of y-axis text (row labels)
          cexCol = 1.2   # (Optional) Adjust the size of x-axis text (column labels)
)
ggbiplot(prcomp(t(case2_EN), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=colSide)

#case1_f3 <- flexiDEG.function3(case1_f2, meta_combined, fdr_cutoff = 1, logfc_cutoff = 2.5) # Run Function 3       ++++ Doesn't seem to be working correctly
# Gene Clustering
# Double Volcano
# Identify rows that start with "Gm" followed by any digit (0-9)


# case1_f4 <- flexiDEG.function4(case1_f2, meta_combined,validation_option = 1) # Run Function 4
# ENplots <- flexiDEG.ENplots(case1_f1, case1_f4, colSide, unique_colSide) # Generate PCA plots
# ggarrange(plotlist = ENplots, ncol=5, nrow=4) # Plots in 5 cols & 4 rows
# # Collect EN results
# EN1 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[1]])]), ]) 
# EN.95 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[2]])]), ]) 
# EN.9 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[3]])]), ]) 
# EN.85 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[4]])]), ])
# EN.8 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[5]])]), ]) 
# EN.75 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[6]])]), ]) 
# EN.7 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[7]])]), ]) 
# EN.65 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[8]])]), ]) 
# EN.6 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[9]])]), ])
# EN.55 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[10]])]), ]) 
# EN.5 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[11]])]), ]) 
# EN.45 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[12]])]), ]) 
# EN.4 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[13]])]), ]) 
# EN.35 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[14]])]), ]) 
# EN.3 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[15]])]), ]) 
# EN.25 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[16]])]), ]) 
# EN.2 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[17]])]), ]) 
# EN.15 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[18]])]), ]) 
# EN.1 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[19]])]), ]) 
# EN.05 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[20]])]), ]) 
# EN.0 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[21]])]), ]) 

#6. Pre-hoc Biology-All Groups ---- 
getwd() 
sample_set <- "Case1ph" 
currentDate <- Sys.Date()
save_name <- paste(currentDate, sample_set)

#sets_celltype <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
#wl_celltype <- split(sets_celltype$gene_symbol, # Genes to split into pathways, by ensembl
#                      sets_celltype$gs_name) # Pathway names
pwl_msigdbr <- c(pwl_hallmark, pwl_kegg) # Compile them all

case1_EN<-as.matrix(case1_EN)# Change dataframe to matrix
is.matrix(case1_EN) # GSVA needs data as matrix w/ genes as rownames; Must be TRUE to proceed
gsvapar<-gsvaParam(case1_EN, c(pwl_msigdbr), maxDiff=TRUE,minSize=3)
gsva_case1 <- gsva(gsvapar)

gsva_case1<-as.data.frame(gsva_case1) # Convert to dataframe
# gsvaf2 <- flexiDEG.function2(gsva_case1, meta_combined) # Run Function 2
# dev.off()  # Close any open graphics devices

library(dplyr)

# meta_combined <- meta_combined %>%
#   mutate(Group = case_when(
#     Group %in% c("Healthy_7", "Healthy_14") ~ "Early",
#     Group == "Healthy_28" ~ "Intermediate",
#     Group %in% c("Healthy_42", "Healthy_56") ~ "Late",
#     TRUE ~ Group  # Keep other values unchanged
#   ))
# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_combined)
unique_colSide <- unique(colSide)
dev.off()
heatmap.2(as.matrix(gsva_case1), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,35), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both",
          cexRow = 1.2,  # Increase font size for row labels
          cexCol = 0.8)   # Increase font size for column labels)



#7.CellTypes GSVA----

#BiocManager::install("clusterProfiler")
library(clusterProfiler)
biocarta_gene_sets <- msigdbr(species="Mus musculus", category="C8") # Large df w/ categories
mm_celltype_sets <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                          biocarta_gene_sets$gs_name) # Pathway names

case1_EN<-as.matrix(case1_EN)# Change dataframe to matrix
is.matrix(case1_EN) # GSVA needs data as matrix w/ genes as rownames; Must be TRUE to proceed
gsvapar_celltype<-gsvaParam(case1_EN, c(mm_celltype_sets), maxDiff=TRUE,minSize=3)
gsva_celltype <- gsva(gsvapar_celltype)

gsva_celltype<-as.data.frame(gsva_celltype) # Convert to dataframe
gsvaf_celltype <- flexiDEG.function2(gsva_celltype, meta_combined) # Run Function 2
# dev.off()  # Close any open graphics devices
library(dplyr)

# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_combined)
unique_colSide <- unique(colSide)
dev.off()
heatmap.2(as.matrix(gsvaf_celltype), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,40), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both",
          cexRow = 1,  # Increase font size for row labels
          cexCol = 0.8)   # Increase font size for column labels)


#8. Opertional Tolerance-Elastic Net ----

meta_combined$TimewiseGroup <- ifelse(
  meta_combined$Group == "Tolerance",
  paste(meta_combined$Group, meta_combined$Day, sep = "_"),
  meta_combined$Group
)
unique(meta_combined$TimewiseGroup)

# Filter out samples that are not in the unwanted groups
samples_to_keep <- !meta_combined$Group %in% c("Rejection(No Anti-CD40L)", "Rejection(Anti-CD40L)")

# Update meta_combined to exclude the unwanted samples
meta_combined_filtered <- meta_combined[samples_to_keep, ]

# Update case1_EN to exclude the corresponding columns
case1_EN_filtered <- case1_EN[, samples_to_keep]
case1_EN_filtered <-as.data.frame(case1_EN_filtered)


# Step 1: Identify Day 7 samples
day7_samples <- meta_combined_filtered$Day == "7"

# Step 2: Compute normalization factors
# Assuming `case1_EN_filtered` contains numeric data with samples as columns
# Calculate the mean expression values for Day 7 samples
day7_means <- rowMeans(case1_EN_filtered[, day7_samples], na.rm = TRUE)

# Step 3: Normalize other samples to Day 7
# Divide each column by the Day 7 mean (row-wise operation)
normalized_case1_EN_filtered <- case1_EN_filtered
for (sample in 1:ncol(normalized_case1_EN_filtered)) {
  normalized_case1_EN_filtered[, sample] <- normalized_case1_EN_filtered[, sample] / day7_means
}



# Normalized data is ready in `normalized_case1_EN_filtered`

meta_combined_filtered$Stage <- with(meta_combined_filtered, ifelse(
  TimewiseGroup %in% c("Tolerance_14"), "Early",
  ifelse(TimewiseGroup == "Tolerance_28", "Intermediate",
         ifelse(TimewiseGroup %in% c("Tolerance_42", "Tolerance_56"), "Late",
                ifelse(TimewiseGroup == "Tolerance_70", "End Timepoint", NA)
         )
  )
))


unique(meta_combined_filtered$Stage)
meta_combined_filtered$Group<-meta_combined_filtered$Stage

# Filter out samples that are not in the unwanted groups
samples_to_keep <- !meta_combined_filtered$TimewiseGroup %in% c("Tolerance_7","Tolerance_70")

# Update meta_combined to exclude the unwanted samples
meta_combined_filtered <- meta_combined_filtered[samples_to_keep, ]
meta_combined_filtered 
# Update case1_EN to exclude the corresponding columns
normalized_case1_EN_filtered <- normalized_case1_EN_filtered[, samples_to_keep]
normalized_case1_EN_filtered <-as.data.frame(normalized_case1_EN_filtered)
meta_combined_filtered$Group<-meta_combined_filtered$Day
case2_EN_filtered <- flexiDEG.function2(normalized_case1_EN_filtered, meta_combined_filtered) # Run Function 2
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_combined_filtered)
unique_colSide <- unique(colSide)
unique(meta_combined_filtered$Group)

case4_EN_filtered <- flexiDEG.function4(case2_EN_filtered, meta_combined_filtered,validation_option = 1) # Run Function 4-,validation_option = 1

EN1 <- na.omit(normalized_case1_EN_filtered[unique(rownames(normalized_case1_EN_filtered)[as_vector(case4_EN_filtered[[1]])]), ]) 
EN2 <- na.omit(normalized_case1_EN_filtered[unique(rownames(normalized_case1_EN_filtered)[as_vector(case4_EN_filtered[[2]])]), ]) 
EN3 <- na.omit(normalized_case1_EN_filtered[unique(rownames(normalized_case1_EN_filtered)[as_vector(case4_EN_filtered[[3]])]), ]) 
EN4 <- na.omit(normalized_case1_EN_filtered[unique(rownames(normalized_case1_EN_filtered)[as_vector(case4_EN_filtered[[4]])]), ]) 
EN7 <- na.omit(normalized_case1_EN_filtered[unique(rownames(normalized_case1_EN_filtered)[as_vector(case4_EN_filtered[[7]])]), ]) 

dev.off()
heatmap.2(as.matrix(EN3), scale="row", col=coul_gsva, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colSide, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both",  cexRow = 1.2,  # Increase the size of y-axis text (row labels)
          cexCol = 1.2   # (Optional) Adjust the size of x-axis text (column labels)
)
ggbiplot(prcomp(t(EN7), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=colSide)

#case1_f3 <- flexiDEG.function3(case1_f2, meta_combined, fdr_cutoff = 1, logfc_cutoff = 2.5) # Run Function 3       ++++ Doesn't seem to be working correctly
# Gene Clustering
# Double Volcano
# Identify rows that start with "Gm" followed by any digit (0-9)


# case1_f4 <- flexiDEG.function4(case1_f2, meta_combined,validation_option = 1) # Run Function 4
# ENplots <- flexiDEG.ENplots(case1_f1, case1_f4, colSide, unique_colSide) # Generate PCA plots
# ggarrange(plotlist = ENplots, ncol=5, nrow=4) # Plots in 5 cols & 4 rows
# # Collect EN results
# EN1 <- na.omit(case1_f1[unique(rownames(case1_f1)[as_vector(case1_f4[[1]])]), ]) 
