# File Info ----
# Author: Jyotirmoy Roy
# Title: Analysis for Islet Transplant Rejection
# Date Created: October 2025
# Info: 
# Ref: 
# Load libraries ----
# install.packages("pacman") 
# Install the following packages
pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)

library(patchwork)
library(tibble)
library(Seurat)
library(EnhancedVolcano)
library(stringr)
library(msigdbr)
library(dplyr)
library(clusterProfiler)
library(gghalves)
library(circlize)
library(ComplexHeatmap)
library(tidyr)



# 1. Load the Data ----
# Organize Data
setwd("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/")
getwd()

# Allo Transplant Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch1/Metadata_Batch1.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch2/Metadata_Batch2.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch3/Metadata_Batch3.csv", sep=",", header=T) # Metadata file
meta_batch4 <- read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch4/Metadata_Batch4.csv", sep=",", header=T) # Metadata file

meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)
meta_batch4 <- as.data.frame(meta_batch4)
# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3,meta_batch4)

# Preview the combined metadata
head(meta_combined)
unique(meta_combined$Group)

#All Transplant Counts Data Importing
counts_batch1 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch1/IsTx_gene_expected_count_annot_batch1.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch1 <- na.omit(counts_batch1)

counts_batch2 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch2/IsTx_gene_expected_count_annot_batch2.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch2 <- na.omit(counts_batch2)

counts_batch3 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch3/IsTx_gene_expected_count_annot_batch3.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch3 <- na.omit(counts_batch3)

counts_batch4 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch4/IsTx_gene_expected_count_annot_batch4.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch4 <- na.omit(counts_batch4)


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

counts_batch4 <- counts_batch4[!duplicated(counts_batch4[, 1]), ]
genes <- counts_batch4[, 1]
rownames(counts_batch4) <- genes
counts_batch4 <- counts_batch4[, -1]

#Combine counts data
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
# Then merge the result with counts_batch4
combined_counts <- merge(combined_counts, counts_batch4, by = "row.names", all = TRUE)
# Rename the Row.names column back
rownames(combined_counts) <- combined_counts$Row.names
combined_counts$Row.names <- NULL

# Identify the samples to keep (not "Technical Rejection")
samples_to_keep <- meta_combined$Samples[meta_combined$Group != "Technical Rejection"]

# Filter the meta_combined data
meta_combined <- meta_combined[meta_combined$Group != "Technical Rejection", ]

# Filter the combined_counts data to keep only columns corresponding to samples_to_keep
combined_counts <- combined_counts[, colnames(combined_counts) %in% samples_to_keep]


# Preview the combined dataset
head(combined_counts)


# 2.Preprocessing and Cleaning ----

# Remove zero and low expressed genes
combined_counts <- combined_counts[, meta_combined$Samples]  # Ensure Sample_IDs match column names in combined_counts

IsletTransplantCounts <- flexiDEG.function1(combined_counts, meta_combined, # Run Function 1
                                            convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                            batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0, 0

#Remove undefined and pseudogenes
remove_pattern <- "^Gm[0-9]|^AC[0-9]|^AL[0-9]|^AI[0-9]|^AW[0-9]|^AF[0-9]|^BB[0-9]|^BC[0-9]|^CT[0-9]|^CAAA|^BX[0-9]|^CN[0-9]|^CR[0-9]|^C[0-9]{4,}|^Olfr"
rows_to_remove <- grep(remove_pattern, rownames(IsletTransplantCounts)) #Remove Gm genes
IsletTransplantCounts <- IsletTransplantCounts[-rows_to_remove, ]
# connect to Ensembl mouse database

options(timeout = 120)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "www"
)
# get gene biotypes
gene_info <- getBM(
  attributes = c("mgi_symbol", "gene_biotype"),
  filters = "mgi_symbol",
  values = rownames(IsletTransplantCounts),
  mart = mart
)

pseudo_genes <- gene_info$mgi_symbol[grep("pseudogene", gene_info$gene_biotype)]
# remove them
IsletTransplantCounts <- IsletTransplantCounts[
  !(rownames(IsletTransplantCounts) %in% pseudo_genes), ]


write.csv(IsletTransplantCounts, file = "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Combined/Combined_Raw_Counts_ITx_Filtered.csv", row.names = TRUE)

# Saving meta_combined dataframe as a CSV file
write.csv(meta_combined, file = "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Combined/Combined_Metadata_ITx.csv", row.names = FALSE)

# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_combined)
unique_colSide <- unique(colSide)

# 3.Create Univeral DESqEQ Object ----

IsletTransplantCounts <- as.matrix(IsletTransplantCounts)
storage.mode(IsletTransplantCounts) <- "integer"

dds_IsletTransplant <- DESeqDataSetFromMatrix(IsletTransplantCounts, meta_combined,
                                              design = ~ 1)   # dummy design for now
colData(dds_IsletTransplant)$logIEQ <- log(colData(dds_IsletTransplant)$IEQ)

saveRDS(dds_IsletTransplant, file = "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Robjects/dds_IsletTransplant_master_v2.rds")
unique(colData(dds_IsletTransplant)$Group)

# Check IEQ Distribution----

meta_df <- meta_combined %>%
  as.data.frame() %>%
  as_tibble()

meta_mouse <- meta_df %>%
  mutate(
    MouseID = paste0(as.character(Animal), "_", as.character(Cohort))
  ) %>%
  distinct(MouseID, .keep_all = TRUE)

ggplot(meta_mouse, aes(x = Group, y = IEQ, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_minimal(base_size = 16) +
  labs(y = "IEQ")



# 4. Allogeneic Vs Syngeneic Signature----
dds_IsletTransplant<-readRDS("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Robjects/dds_IsletTransplant_master_v2.rds")
# subset samples of Allo and Syn (without treatment) from Day 7 and Day 14
sel <- colData(dds_IsletTransplant)$Group %in% c("Control Allogeneic","Control Syngeneic") & colData(dds_IsletTransplant)$Day %in% c(7,14)
dds_IsletTransplant_AlloVsSyn <- dds_IsletTransplant[, sel]
dds_IsletTransplant_AlloVsSyn$RecGender <- factor(dds_IsletTransplant_AlloVsSyn$RecGender)
dds_IsletTransplant_AlloVsSyn$Batch       <- factor(dds_IsletTransplant_AlloVsSyn$Batch)
dds_IsletTransplant_AlloVsSyn$LibraryPrep <- factor(dds_IsletTransplant_AlloVsSyn$LibraryPrep)
dds_IsletTransplant_AlloVsSyn$Group <- factor(
  dds_IsletTransplant_AlloVsSyn$Group,
  levels = c("Control Syngeneic", "Control Allogeneic")  # order sets baseline
)

dds_IsletTransplant_AlloVsSyn$Day <- factor(dds_IsletTransplant_AlloVsSyn$Day,
                                            levels = c(7, 14))  # baseline = 7

unique(dds_IsletTransplant_AlloVsSyn$RecGender)

cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
# Check sample distribution
lapply(cd[, c("Batch","LibraryPrep","logIEQ","Day","Group")], function(x) table(x, useNA="ifany"))
table(cd$Group,cd$Day)

#                   7 14
#Control Syngeneic  5  4
#Control Allogeneic 5  4

## DESEQ- Combined Timepoints ----
design(dds_IsletTransplant_AlloVsSyn) <- ~ Batch +  Day + Group 
keep <- rowSums(counts(dds_IsletTransplant_AlloVsSyn) >= 10) >= 4  # Smallest number of samples in a group
dds_IsletTransplant_AlloVsSyn <- dds_IsletTransplant_AlloVsSyn[keep, ]
dds_IsletTransplant_AlloVsSyn <- DESeq(dds_IsletTransplant_AlloVsSyn)
resultsNames(dds_IsletTransplant_AlloVsSyn)

res_ALLOvSYN <- results(dds_IsletTransplant_AlloVsSyn,
                           name = "Group_Control.Allogeneic_vs_Control.Syngeneic")

summary(res_ALLOvSYN)


# Thresholds
q_cut  <- 0.10
fc_cut <- 1

## Prep results as data.frames with a 'gene' column
ALLOvSYN  <- as.data.frame(res_ALLOvSYN);  ALLOvSYN$gene  <- rownames(res_ALLOvSYN)


write.csv(ALLOvSYN, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/DESEQResults_CombinedDays_Allo_vs_Syn.csv", row.names = TRUE)

## Elastic Net Feature Selection----


# DO analysis for timepoint combined
genes_ALLOvSYN  <- as.data.frame(res_ALLOvSYN);  genes_ALLOvSYN$gene  <- rownames(res_ALLOvSYN)

sig_genes_ALLOvSYN <- genes_ALLOvSYN$gene[
  !is.na(genes_ALLOvSYN$pvalue) &
    genes_ALLOvSYN$pvalue <= 0.05 &
    abs(genes_ALLOvSYN$log2FoldChange) >= 0.5
]

length(sig_genes_ALLOvSYN)
# 393 significant genes

vsd <- vst(dds_IsletTransplant_AlloVsSyn, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
design_enet <- model.matrix(~ Group, data = cd)
mat_enet <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  covariates = 
    model.matrix(~ Day, data = cd)[, -1, drop = FALSE]
    ,
  design = design_enet
)

# expression matrix for glmnet: samples x genes
x_AlloVsSyn <- t(mat_enet[sig_genes_ALLOvSYN, , drop = FALSE])

# binary outcome
y_AlloVsSyn <- ifelse(cd$Group == "Control Allogeneic", 1, 0)

# Find best alpha with LOOCV
set.seed(123)
alpha_grid <- seq(0, 1, by = 0.1)
cv_summary <- data.frame()

for (a in alpha_grid) {
  cvfit <- cv.glmnet(
    x = x_AlloVsSyn,
    y = y_AlloVsSyn,
    family = "binomial",
    alpha = a,
    foldid = 1:length(y_AlloVsSyn),   # LOOCV
    type.measure = "deviance",
    standardize = TRUE
  )
  
  cv_summary <- rbind(
    cv_summary,
    data.frame(
      alpha = a,
      cv_error = min(cvfit$cvm),
      lambda_min = cvfit$lambda.min,
      lambda_1se = cvfit$lambda.1se
    )
  )
}

cv_summary[order(cv_summary$cv_error), ]
cv_summary

best_alpha <- cv_summary$alpha[which.min(cv_summary$cv_error)]
best_alpha
#0.5
# Stability selection
set.seed(123)
n_iter <- 1000
n_samples <- nrow(x_AlloVsSyn)


class1_idx <- which(y_AlloVsSyn == unique(y_AlloVsSyn)[1])
class2_idx <- which(y_AlloVsSyn == unique(y_AlloVsSyn)[2])

selected_list <- vector("list", n_iter)

for (i in 1:n_iter) {
  
  # Subsample ~80% of samples each time

  idx1 <- sample(class1_idx, size = round(0.8 * length(class1_idx)))
  idx2 <- sample(class2_idx, size = round(0.8 * length(class2_idx)))
  idx  <- c(idx1, idx2)
  
  x_sub <- x_AlloVsSyn[idx, ]
  y_sub <- y_AlloVsSyn[idx]
  
  cvfit <- cv.glmnet(
    x_sub,
    y_sub,
    family = "binomial",
    alpha = best_alpha,
    foldid = 1:length(y_sub),
    type.measure = "deviance",
    standardize = TRUE
  )
  
  coef_mat <- coef(cvfit, s = "lambda.1se")
  genes <- rownames(coef_mat)[coef_mat[,1] != 0]
  genes <- setdiff(genes, "(Intercept)")
  
  selected_list[[i]] <- genes
}

#Select Genes which appear atleast 70% of times
all_genes <- colnames(x_AlloVsSyn)

freq <- sapply(all_genes, function(g) {
  mean(sapply(selected_list, function(s) g %in% s))
})

freq_table <- data.frame(
  gene = names(freq),
  selection_frequency = freq
)
freq_table <- freq_table[order(freq_table$selection_frequency, decreasing = TRUE), ]
stable_genes <- subset(freq_table, selection_frequency >= 0.7) #Select genes appearing atleast 70 percent of time
stable_genes
# 30 genes -3 sex related genes removed

#Plot Heatmp and PCA using Selected Genes 
stable_gene_names <- stable_genes$gene
# Sex related genes to be removed
remove_genes <- c("Xist", "Eif2s3x", "Kdm6a")

# remove selected genes
stable_gene_names <- setdiff(stable_gene_names, remove_genes)


mat_stable <- mat_enet[stable_gene_names, , drop = FALSE]
mat_scaled <- t(scale(t(mat_stable)))


# Column annotation
annotation_col <- data.frame(
  Day = factor(cd$Day, levels = c("7", "14")),
  Group = factor(
    cd$Group,
    levels = c("Control Syngeneic", "Control Allogeneic")
  )
)

rownames(annotation_col) <- colnames(mat_scaled)

stopifnot(
  identical(rownames(annotation_col), colnames(mat_scaled))
)

# Top annotation
ha <- HeatmapAnnotation(
  df = annotation_col,
  col = list(
    Group = c(
      "Control Allogeneic" = "#FF2400",
      "Control Syngeneic" = "#2E6F40"
    ),
    Day = c(
      "7" = "#F28E6B",
      "14" = "#6FA287"
    )
  ),
  
  border = TRUE,                # <-- adds borders around annotation cells
  
  gp = grid::gpar(
    col = "grey80",              # border color
    lwd = 0.5                   # border width
  ),
  
  annotation_name_gp = grid::gpar(fontsize = 12)
)

Heatmap(
  mat_scaled,
  name = "Z-score",
  
  top_annotation = ha,
  
  # Separate Day 7 and Day 14
  column_split = annotation_col$Day,
  
  # Cluster independently within each day
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  
  column_title = c("Day 7", "Day 14"),
  
  # Cluster genes
  cluster_rows = TRUE,
  
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  show_column_names = FALSE,
  show_row_names = TRUE,
  
  row_names_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  col =  colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  
  # Borders around individual heatmap cells
  rect_gp = grid::gpar(
    col = "grey80",
    lwd = 0.5
  ),
  
  # Gap between Day 7 and Day 14
  column_gap = grid::unit(4, "mm"),
  
  heatmap_legend_param = list(
    title = "Expression\nZ-score"
  )
)


#PCA Analysis
mat_pca <- t(mat_stable)
pca <- prcomp(mat_pca, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = cd$Group
)

# Define colors & shapes
group_colors <- c("Control Allogeneic" = "#E60000", "Control Syngeneic" = "#2E6F40" )
group_shapes <- c("Control Allogeneic" = 25, "Control Syngeneic" = 24)

ggplot(pca_df, aes(PC1, PC2, color = Group,, shape = Group)) +
  geom_point(aes(fill = Group), size = 5, stroke = 1.2) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA of Stable Elastic Net Genes",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1],1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2],1), "%)")
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

### Gene Expressions Over Time----


vsd <- vst(dds_IsletTransplant_AlloVsSyn, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
design_temporal <- model.matrix(
  ~ Group * Day,
  data = cd
)

mat_temporal <- limma::removeBatchEffect(
  mat,
  batch = cd$Batch,
  design = design_temporal
)

mat_plot <- mat_temporal[
  stable_gene_names,
  ,
  drop = FALSE
]

# Check that metadata and expression samples align
stopifnot(
  identical(colnames(mat_plot), rownames(cd))
)


#  Convert expression matrix to long format
plot_df <- as.data.frame(mat_plot) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  left_join(
    cd %>%
      rownames_to_column("Sample") %>%
      select(Sample, Day, Group),
    by = "Sample"
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = c("7", "14")
    ),
    Group = factor(
      Group,
      levels = c(
        "Control Syngeneic",
        "Control Allogeneic"
      )
    )
  )

# Check for missing metadata
stopifnot(
  !any(is.na(plot_df$Day)),
  !any(is.na(plot_df$Group))
)


gene_trajectory_plot <- ggplot(
  plot_df,
  aes(
    x = Day,
    y = Expression,
    color = Group,
    group = Group
  )
) +
  
  # Individual samples
  geom_point(
    aes(
      shape = Group,
      fill = Group
    ),
    position = position_jitter(
      width = 0.06,
      height = 0
    ),
    size = 2,
    alpha = 0.55,
    stroke = 0.5
  ) +
  
  # Group mean connecting Day 7 and Day 14
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.9
  ) +
  
  # Group mean symbols
  stat_summary(
    aes(
      shape = Group,
      fill = Group
    ),
    fun = mean,
    geom = "point",
    size = 3.3,
    stroke = 0.7
  ) +
  
  # Mean ± SEM
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.10,
    linewidth = 0.6
  ) +
  
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_color_manual(
    values = c(
      "Control Syngeneic" = "#2E6F40",
      "Control Allogeneic" = "#FF2400"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Control Syngeneic" = "#2E6F40",
      "Control Allogeneic" = "#FF2400"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Control Syngeneic" = 24,  # upward triangle
      "Control Allogeneic" = 25  # downward triangle
    )
  ) +
  
  labs(
    title = "Temporal Expression- Allogeneic vs Syngeneic Signature",
    x = "Day post-transplant",
    y = "Variance-stabilized (VST) expression",
    color = "Group",
    fill = "Group",
    shape = "Group"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    strip.text = element_text(
      face = "bold.italic",
      colour = "black",
      size = 12
    ),
    
    strip.background = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 0.6
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    panel.spacing = grid::unit(0.8, "lines")
  )

gene_trajectory_plot





## GSEA Analysis----
# Paths to your saved results
day_combined<-"/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/DESEQResults_CombinedDays_Allo_vs_Syn.csv"

# Import
daycombined_ALLOvSYN<- read.csv(day_combined, row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_AlloSyn <- daycombined_ALLOvSYN$stat; names(lfc_vector_AlloSyn) <- rownames(daycombined_ALLOvSYN)

# Drop NAs
lfc_vector_AlloSyn<- lfc_vector_AlloSyn[!is.na(lfc_vector_AlloSyn)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_AlloSyn<- sort(lfc_vector_AlloSyn, decreasing = TRUE)


# --- Collect each set and convert into 2-column (gs_name, gene_symbol) ---

# C8
CellTypeMSigDB_gene_sets <- msigdbr(species="Mus musculus", category="C8")
mm_c8_sets <- split(CellTypeMSigDB_gene_sets$gene_symbol, CellTypeMSigDB_gene_sets$gs_name)
mm_c8_df <- data.frame(
  gs_name = rep(names(mm_c8_sets), sapply(mm_c8_sets, length)),
  gene_symbol = unlist(mm_c8_sets)
)


# Hallmark
hallmark <- msigdbr(species = "Mus musculus", category  = "H")
mm_h_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
mm_h_df <- data.frame(
  gs_name = rep(names(mm_h_sets), sapply(mm_h_sets, length)),
  gene_symbol = unlist(mm_h_sets)
)


# KEGG
kegg_all <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_LEGACY")
mm_kegg_sets <- split(kegg_all$gene_symbol, kegg_all$gs_name)
mm_kegg_df <- data.frame(
  gs_name = rep(names(mm_kegg_sets), sapply(mm_kegg_sets, length)),
  gene_symbol = unlist(mm_kegg_sets)
)
mm_all_df <- rbind(mm_c8_df, mm_h_df, mm_kegg_df)

# Combined

gsea_results_daycombined_allosyn <- GSEA(
  geneList      = lfc_vector_AlloSyn,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_daycombined_allosyn_df <- as.data.frame(gsea_results_daycombined_allosyn)


# Full results Day Combined
write.csv(gsea_results_daycombined_allosyn_df,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/GSEAResults_AlloVsSyn_DayCombined.csv",
          row.names = FALSE)


# Combined Day Plot

library(forcats)
pathways_of_interest <- unique(c(
  # Up in Allo
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "DESCARTES_FETAL_LIVER_LYMPHOID_CELLS",
  "DESCARTES_FETAL_PANCREAS_LYMPHOID_CELLS",
  "KEGG_RETINOL_METABOLISM",
  "KEGG_TYROSINE_METABOLISM",
  "KEGG_TRYPTOPHAN_METABOLISM",
  "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION",
  "KEGG_BUTANOATE_METABOLISM",
  "KEGG_STEROID_BIOSYNTHESIS",
  #Down in Allo
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "AIZARANI_LIVER_C23_KUPFFER_CELLS_3",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_E2F_TARGETS",
  "HAY_BONE_MARROW_NEUTROPHIL",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY"
))

# Subset and prepare plotting data
plot_df <- gsea_results_daycombined_allosyn_df %>%
  filter(ID %in% pathways_of_interest) %>%
  mutate(
    neglog10_padj = -log10(p.adjust),
    neglog10_padj = ifelse(is.infinite(neglog10_padj), NA, neglog10_padj),
    Direction = ifelse(NES >= 0, "Up in Allo", "Down in Allo")
  ) %>%
  arrange(NES) %>%
  mutate(
    ID = factor(ID, levels = ID)
  )

# Plot
ggplot(plot_df, aes(x = NES, y = ID, size = neglog10_padj, color = NES)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](adjusted~italic(p))), range = c(3, 10)) +
  scale_color_gradient2(
    low = "forestgreen",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "NES"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 16) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Allogeneic vs Syngeneic"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

## Islet Single Cell Mapping----

IsletScRNA = readRDS("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Islet ScRNASeq/islet_graft_seurat_v8.rds")
DimPlot(IsletScRNA,label = TRUE, label.box = T,label.size = 8,repel = T,pt.size = 0.9)+
  NoAxes() +NoLegend()
DimPlot(
  IsletScRNA,
  label = TRUE,
  label.box = TRUE,
  label.size = 8,
  repel = TRUE,
  group.by = "condition",
  pt.size = 0.9,
  cols = c(
    "Allogeneic" = "red",
    "Syngeneic" = "forestgreen"
  )
) +
  NoAxes() +
  NoLegend()
count_df <- IsletScRNA@meta.data %>%
  dplyr::select(celltype, condition) %>%
  dplyr::filter(condition %in% c("Allogeneic", "Syngeneic")) %>%
  dplyr::group_by(celltype, condition) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")


count_df$celltype <- factor(count_df$celltype)
count_df$condition <- factor(count_df$condition, levels = c("Syngeneic", "Allogeneic"))

ggplot(count_df, aes(x = celltype, y = n, fill = condition)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(
    values = c(
      "Syngeneic" = "forestgreen",
      "Allogeneic" = "red"
    )
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "Cell Types",
    y = "Number of cells",
    title = "Islet Graft Cell Counts"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    legend.title = element_blank()
  )

### IN Allo vs Syn EN Score ----
#Downregulated and Upregulated wrt Allo vs Syn


downregulated_Signature <- c("Myo18b","Cd59b","Rep15","Shisa9","Rragb","Gdf3",
                             "Lncpint","Ifitm5","Malat1","Myo3b")
upregulated_Signature <- c("Hipk4","S1pr5","Trpm6","Tctn2","Lhfpl4","Ido2","Fsip1","Emx2os","Hlf","Gfra1","Bmp5","Col6a5","Myh3","Cbs","Vtn","Dnah8","Sybu")

downregulated_Signature_use <- intersect(downregulated_Signature, rownames(IsletScRNA))
upregulated_Signature_use <- intersect(upregulated_Signature, rownames(IsletScRNA))

combined_signature_use<-c(upregulated_Signature_use,downregulated_Signature_use)
DotPlot(IsletScRNA, features = combined_signature_use) +
  scale_color_gradient(low = "grey", high = "red") +
  scale_size(range = c(0, 8), limits = c(0, 100)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

# IN based Allogeneic  Score
DefaultAssay(IsletScRNA)
#SCT Assay used
expr <- GetAssayData(IsletScRNA,layer = "data")
up_score <- colMeans(expr[upregulated_Signature_use, ])
down_score <- colMeans(expr[downregulated_Signature_use, ])
IsletScRNA$ElasticNetSignature <- up_score - down_score



plot_df <- IsletScRNA@meta.data %>%
  dplyr::select(celltype, condition, ElasticNetSignature) %>%
  filter(!is.na(celltype), !is.na(condition), !is.na(ElasticNetSignature)) %>%
  filter(condition %in% c("Allogeneic", "Syngeneic"))

plot_df$condition <- factor(plot_df$condition, levels = c("Syngeneic", "Allogeneic"))

stat_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(ElasticNetSignature ~ condition)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = 0.15,#max(ElasticNetSignature, na.rm = TRUE) + 0.15 * diff(range(ElasticNetSignature, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = paste0("FDR=", signif(p_adj, 2))
  )
stat_df <- stat_df %>%
  mutate(label = paste0("FDR=", signif(p_adj, 2)))


ggplot(plot_df, aes(x = condition, y = ElasticNetSignature, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.15, size = 0.1, alpha = 0.0) +
  facet_wrap(~celltype, nrow = 2) +
  geom_text(
    data = stat_df,
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Niche-Based Signature Score",
    title = "Niche-Based Signature Score Across Islet Graft Cell Types",
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )+
  scale_fill_manual(
    values = c(
      "Syngeneic" = "forestgreen",
      "Allogeneic" = "red"
    )
  ) +coord_cartesian(ylim = c(-1.2, 0.2))


### ssGSEA Pathway Score ----
colnames(mm_all_df)
pathways_of_interest_updated <- unique(c(
  # Up in Allo
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "KEGG_RETINOL_METABOLISM",
  "KEGG_TYROSINE_METABOLISM",
  "KEGG_TRYPTOPHAN_METABOLISM",
  "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION",
  "KEGG_BUTANOATE_METABOLISM",
  "KEGG_STEROID_BIOSYNTHESIS",
  #Down in Allo
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_E2F_TARGETS",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY"
))


gene_sets_list <- mm_all_df %>%
  filter(gs_name %in% pathways_of_interest_updated) %>%
  distinct(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe()

expr <- GetAssayData(IsletScRNA, assay = "SCT", layer = "data")
gene_sets_list <- lapply(gene_sets_list, function(gs) intersect(gs, rownames(expr)))
gene_sets_list <- gene_sets_list[sapply(gene_sets_list, length) >= 5]
sapply(gene_sets_list, length)

ssgsea_par <- ssgseaParam(
  exprData = expr,
  geneSets = gene_sets_list,
  assay = NA_character_,
  annotation = NULL,
  minSize = 5,
  maxSize = 500,
  normalize = TRUE,
  checkNA = "auto",
  use = "everything",
  verbose = TRUE
)
ssgsea_scores <- gsva(ssgsea_par)

ssgsea_df <- as.data.frame(t(ssgsea_scores)) %>%
  rownames_to_column("cell")

meta_df <- IsletScRNA@meta.data %>%
  rownames_to_column("cell")

meta_ssgsea <- meta_df %>%
  left_join(ssgsea_df, by = "cell")


summary_long <- meta_ssgsea %>%
  filter(condition %in% c("Allogeneic", "Syngeneic")) %>%
  select(celltype, condition, all_of(rownames(ssgsea_scores))) %>%
  pivot_longer(
    cols = all_of(rownames(ssgsea_scores)),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  group_by(celltype, condition, pathway) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    .groups = "drop"
  )


heatmap_df <- summary_long %>%
  select(celltype, condition, pathway, mean_score) %>%
  pivot_wider(names_from = condition, values_from = mean_score) %>%
  mutate(delta_Allo_minus_Syn = Allogeneic - Syngeneic)

heatmap_mat <- heatmap_df %>%
  select(celltype, pathway, delta_Allo_minus_Syn) %>%
  pivot_wider(names_from = pathway, values_from = delta_Allo_minus_Syn) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

stats_long <- meta_ssgsea %>%
  filter(condition %in% c("Allogeneic", "Syngeneic")) %>%
  select(celltype, condition, all_of(rownames(ssgsea_scores))) %>%
  pivot_longer(
    cols = all_of(rownames(ssgsea_scores)),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  group_by(celltype, pathway) %>%
  summarise(
    p_value = tryCatch(wilcox.test(score ~ condition)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

sig_mat <- stats_long %>%
  select(celltype, pathway, sig) %>%
  pivot_wider(names_from = pathway, values_from = sig) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

col_fun <- colorRamp2(
  c(min(heatmap_mat, na.rm = TRUE), 0, max(heatmap_mat, na.rm = TRUE)),
  c("forestgreen", "white", "red")
)

max_val <- max(abs(heatmap_mat), na.rm = TRUE)

# Visulaize all the pathways
ht<-Heatmap(
  heatmap_mat,
  name = "Allo - Syn",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  row_names_side = "left",
  column_names_rot = 45,
  heatmap_legend_param = list(title = "ssGSEA score diff"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(sig_mat[i, j], x, y, gp = grid::gpar(fontsize = 12))
  }
)

draw(
  ht,
  padding = unit(c(30, 40, 10, 10), "mm")  # top, right, bottom, left
)


## Auto-Allo Genes Signature Cross Preservation----

# subset samples of D7 and D14
sel <- colData(dds_IsletTransplant)$Group %in% c("Control Allogeneic","Control Syngeneic") & colData(dds_IsletTransplant)$Day %in% c(7,14)
dds_IsletTransplant_AlloVsSyn <- dds_IsletTransplant[, sel]
dds_IsletTransplant_AlloVsSyn$RecGender <- factor(dds_IsletTransplant_AlloVsSyn$RecGender)
dds_IsletTransplant_AlloVsSyn$Batch       <- factor(dds_IsletTransplant_AlloVsSyn$Batch)
dds_IsletTransplant_AlloVsSyn$LibraryPrep <- factor(dds_IsletTransplant_AlloVsSyn$LibraryPrep)
dds_IsletTransplant_AlloVsSyn$Group <- factor(
  dds_IsletTransplant_AlloVsSyn$Group,
  levels = c("Control Syngeneic", "Control Allogeneic")  # order sets baseline
)

dds_IsletTransplant_AlloVsSyn$Day <- factor(dds_IsletTransplant_AlloVsSyn$Day,
                                            levels = c(7, 14))  # baseline = 7


design(dds_IsletTransplant_AlloVsSyn) <- ~ Batch +  Day + Group 
dds_IsletTransplant_AlloVsSyn <- DESeq(dds_IsletTransplant_AlloVsSyn)
resultsNames(dds_IsletTransplant_AlloVsSyn)


vsd <- vst(dds_IsletTransplant_AlloVsSyn, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
design_enet <- model.matrix(~ Group, data = cd)
mat_enet <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  covariates = 
    model.matrix(~ Day, data = cd)[, -1, drop = FALSE]
  ,
  design = design_enet
)



#Plot Heatmp and PCA using Selected Genes 
Auto_allo_signature <- c("Rsph10b","S1pr5","Erdr1","Pla2g4b","Clec2g","Ctsg","Syce1","Fndc7",
                         "Epgn","Tmem132e","Stkld1","Osbpl6","Sycp3")
setdiff(Auto_allo_signature, rownames(mat_enet))
mat_stable <- mat_enet[Auto_allo_signature, , drop = FALSE]
mat_scaled <- t(scale(t(mat_stable)))
annotation_col <- data.frame(
  Day = cd$Day,
  Group = cd$Group
)
rownames(annotation_col) <- colnames(mat_scaled)


# Column annotation
annotation_col <- data.frame(
  Day = factor(cd$Day, levels = c("7", "14")),
  Group = factor(
    cd$Group,
    levels = c("Control Syngeneic", "Control Allogeneic")
  )
)

rownames(annotation_col) <- colnames(mat_scaled)

stopifnot(
  identical(rownames(annotation_col), colnames(mat_scaled))
)

# Top annotation
ha <- HeatmapAnnotation(
  df = annotation_col,
  col = list(
    Group = c(
      "Control Syngeneic" = "forestgreen",
      "Control Allogeneic" = "red"
    ),
    Day = c(
      "7" = "#F28E6B",
      "14" = "#6FA287"
    )
  ),
  
  border = TRUE,                # <-- adds borders around annotation cells
  
  gp = grid::gpar(
    col = "grey80",              # border color
    lwd = 0.5                   # border width
  ),
  
  annotation_name_gp = grid::gpar(fontsize = 12)
)


Heatmap(
  mat_scaled,
  name = "Z-score",
  top_annotation = ha,
  # Separate Day 7 and Day 14
  column_split = annotation_col$Day,
  # Cluster independently within each day
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_title = c("Day 7", "Day 14"),
  # Cluster genes
  cluster_rows = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  col = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  # Borders around individual heatmap cells
  rect_gp = grid::gpar(
    col = "grey80",
    lwd = 0.5
  ),
  # Gap between Day 7 and Day 14
  column_gap = grid::unit(4, "mm"),
  heatmap_legend_param = list(
    title = "Expression\nZ-score"
  )
)

#PCA Analysis
mat_pca <- t(mat_stable)
pca <- prcomp(mat_pca, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = cd$Group
)

# Define colors & shapes
group_colors <- c( "Control Syngeneic" = "forestgreen",
                   "Control Allogeneic" = "red" )
group_shapes <- c("Control Syngeneic" = 24,  # upward triangle
                  "Control Allogeneic" = 25  # downward triangle
                  )

ggplot(pca_df, aes(PC1, PC2, color = Group,, shape = Group)) +
  geom_point(aes(fill = Group), size = 5, stroke = 1.2) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "Auto-Allo Genes from NOD: Syn vs Allo B6",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1],1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2],1), "%)")
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )


### Gene Expressions Over Time

vsd <- vst(dds_IsletTransplant_AlloVsSyn, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
design_temporal <- model.matrix(
  ~ Group * Day,
  data = cd
)

mat_temporal <- limma::removeBatchEffect(
  mat,
  batch = cd$Batch,
  design = design_temporal
)

mat_plot <- mat_temporal[
  Auto_allo_signature,
  ,
  drop = FALSE
]

# Check that metadata and expression samples align
stopifnot(
  identical(colnames(mat_plot), rownames(cd))
)


#  Convert expression matrix to long format
plot_df <- as.data.frame(mat_plot) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  left_join(
    cd %>%
      rownames_to_column("Sample") %>%
      select(Sample, Day, Group),
    by = "Sample"
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = c("7", "14")
    ),
    Group = factor(
      Group,
      levels = c(
        "Control Syngeneic",
        "Control Allogeneic"
      )
    )
  )

# Check for missing metadata
stopifnot(
  !any(is.na(plot_df$Day)),
  !any(is.na(plot_df$Group))
)


gene_trajectory_plot <- ggplot(
  plot_df,
  aes(
    x = Day,
    y = Expression,
    color = Group,
    group = Group
  )
) +
  
  # Individual samples
  geom_point(
    aes(
      shape = Group,
      fill = Group
    ),
    position = position_jitter(
      width = 0.06,
      height = 0
    ),
    size = 2,
    alpha = 0.55,
    stroke = 0.5
  ) +
  
  # Group mean connecting Day 7 and Day 14
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.9
  ) +
  
  # Group mean symbols
  stat_summary(
    aes(
      shape = Group,
      fill = Group
    ),
    fun = mean,
    geom = "point",
    size = 3.3,
    stroke = 0.7
  ) +
  
  # Mean ± SEM
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.10,
    linewidth = 0.6
  ) +
  
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_color_manual(
    values = c(
      "Control Syngeneic" = "forestgreen", 
      "Control Allogeneic" = "red" 
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Control Syngeneic" = "forestgreen", 
      "Control Allogeneic" = "red" 
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Control Syngeneic" = 24,  # upward triangle
      "Control Allogeneic" = 25  # downward triangle
    )
  ) +
  
  labs(
    title = "Temporal Expression- Auto-Allo NOD Signature in Allo vs Syn B6",
    x = "Day post-transplant",
    y = "Variance-stabilized (VST) expression",
    color = "Group",
    fill = "Group",
    shape = "Group"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    strip.text = element_text(
      face = "bold.italic",
      colour = "black",
      size = 12
    ),
    
    strip.background = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 0.6
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    panel.spacing = grid::unit(0.8, "lines")
  )

gene_trajectory_plot


# 5. Allo+Anti-CD40L-Rejection vs Acceptance-----

# subset samples of D7 and D14
sel <- colData(dds_IsletTransplant)$Group %in% c("Acceptance","Rejection") & colData(dds_IsletTransplant)$Day %in% c(7,14)
dds_IsletTransplant_RejVsAccep <- dds_IsletTransplant[, sel]

dds_IsletTransplant_RejVsAccep$Batch       <- factor(dds_IsletTransplant_RejVsAccep$Batch)
dds_IsletTransplant_RejVsAccep$LibraryPrep <- factor(dds_IsletTransplant_RejVsAccep$LibraryPrep)
dds_IsletTransplant_RejVsAccep$Group <- factor(
  dds_IsletTransplant_RejVsAccep$Group,
  levels = c("Acceptance", "Rejection")  # order sets baseline
)
levels(dds_IsletTransplant_RejVsAccep$Group)

dds_IsletTransplant_RejVsAccep$Day <- factor(dds_IsletTransplant_RejVsAccep$Day,
                                            levels = c(7, 14))  # baseline = 7
levels(dds_IsletTransplant_RejVsAccep$Day)

cd <- as.data.frame(colData(dds_IsletTransplant_RejVsAccep))
table(cd$Group,cd$Day)
#            7 14
#Acceptance 14 13
#Rejection   2  3
#Design Formula
keep <- rowSums(counts(dds_IsletTransplant_RejVsAccep) >= 10) >= 5  # Total Samples in Rejection group used to account for group imbalance
dds_IsletTransplant_RejVsAccep <- dds_IsletTransplant_RejVsAccep[keep, ]
design(dds_IsletTransplant_RejVsAccep) <- ~ Batch + Day + Group 
dds_IsletTransplant_RejVsAccep <- DESeq(dds_IsletTransplant_RejVsAccep)
design(dds_IsletTransplant_RejVsAccep)
resultsNames(dds_IsletTransplant_RejVsAccep)

# Combined Timepoints
res_REJvACCEP_ALL <- results(dds_IsletTransplant_RejVsAccep,
                          name = "Group_Rejection_vs_Acceptance")

summary(res_REJvACCEP_ALL)

ALL_REJVSACCEP  <- as.data.frame(res_REJvACCEP_ALL);  ALL_REJVSACCEP$gene  <- rownames(res_REJvACCEP_ALL)

write.csv(ALL_REJVSACCEP, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Rejection_Vs_Acceptance_aCD40L/DESEQResults_RejvsAccep_DaysCombined.csv", row.names = TRUE)


## Elastic Net Feature Selection----

sig_genes_REJvACCEP <- ALL_REJVSACCEP$gene[
  !is.na(ALL_REJVSACCEP$pvalue) &
    ALL_REJVSACCEP$pvalue <= 0.05 &
    abs(ALL_REJVSACCEP$log2FoldChange) >= 0.5
]

length(sig_genes_REJvACCEP)
#711
vsd <- vst(dds_IsletTransplant_RejVsAccep, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_RejVsAccep))
design_enet <- model.matrix(~ Group, data = cd)
mat_enet <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  covariates = model.matrix(~ Day, data = cd)[, -1, drop = FALSE],
  design = design_enet
)


# expression matrix for glmnet: samples x genes
x_RejVsAccep <- t(mat_enet[sig_genes_REJvACCEP, , drop = FALSE])

# binary outcome
y_RejVsAccep <- ifelse(cd$Group == "Rejection", 1, 0)

# Find best alpha with LOOCV
set.seed(123)
alpha_grid <- seq(0, 1, by = 0.1)
cv_summary <- data.frame()

for (a in alpha_grid) {
  cvfit <- cv.glmnet(
    x = x_RejVsAccep,
    y = y_RejVsAccep,
    family = "binomial",
    alpha = a,
    foldid = 1:length(y_RejVsAccep),   # LOOCV
    type.measure = "deviance",
    standardize = TRUE
  )
  
  cv_summary <- rbind(
    cv_summary,
    data.frame(
      alpha = a,
      cv_error = min(cvfit$cvm),
      lambda_min = cvfit$lambda.min,
      lambda_1se = cvfit$lambda.1se
    )
  )
}

cv_summary[order(cv_summary$cv_error), ]
cv_summary

best_alpha<-0.3 #0.1 was had lowest with lowest cv error  but 0,3 chosen for sparser gene set 

# Stability selection
set.seed(123)
n_iter <- 1000
n_samples <- nrow(x_RejVsAccep)

selected_list <- vector("list", n_iter)

class1_idx <- which(y_RejVsAccep == unique(y_RejVsAccep)[1])
class2_idx <- which(y_RejVsAccep == unique(y_RejVsAccep)[2])

for (i in 1:n_iter) {
  
  # Subsample ~80% of samples each time
  idx1 <- sample(class1_idx, size = round(0.8 * length(class1_idx)))
  idx2 <- sample(class2_idx, size = round(0.8 * length(class2_idx)))
  idx  <- c(idx1, idx2)
  
  x_sub <- x_RejVsAccep[idx, ]
  y_sub <- y_RejVsAccep[idx]
  
  cvfit <- cv.glmnet(
    x_sub,
    y_sub,
    family = "binomial",
    alpha = best_alpha,
    foldid = 1:length(y_sub),
    type.measure = "deviance",
    standardize = TRUE
  )
  
  coef_mat <- coef(cvfit, s = "lambda.1se")
  genes <- rownames(coef_mat)[coef_mat[,1] != 0]
  genes <- setdiff(genes, "(Intercept)")
  
  selected_list[[i]] <- genes
}

#Select Genes which appear atleast 70% of times
all_genes <- colnames(x_RejVsAccep)

freq <- sapply(all_genes, function(g) {
  mean(sapply(selected_list, function(s) g %in% s))
})

freq_table <- data.frame(
  gene = names(freq),
  selection_frequency = freq
)
freq_table <- freq_table[order(freq_table$selection_frequency, decreasing = TRUE), ]
stable_genes <- subset(freq_table, selection_frequency >= 0.7) #Select genes appearing atleast 70 percent of time
stable_genes
#25

#Plot Heatmp and PCA using Selected Genes 
stable_gene_names <- stable_genes$gene
mat_stable <- mat_enet[stable_gene_names, , drop = FALSE]
mat_scaled <- t(scale(t(mat_stable)))
annotation_col <- data.frame(
  Day = cd$Day,
  Group = cd$Group
)
rownames(annotation_col) <- colnames(mat_scaled)


# Column annotation
annotation_col <- data.frame(
  Day = factor(cd$Day, levels = c("7", "14")),
  Group = factor(
    cd$Group,
    levels = c("Acceptance", "Rejection")
  )
)

rownames(annotation_col) <- colnames(mat_scaled)

stopifnot(
  identical(rownames(annotation_col), colnames(mat_scaled))
)

# Top annotation
ha <- HeatmapAnnotation(
  df = annotation_col,
  col = list(
    Group = c(
      "Acceptance" = "#064273",
      "Rejection" = "#F28500"
    ),
    Day = c(
      "7" = "#F28E6B",
      "14" = "#6FA287"
    )
  ),
  
  border = TRUE,                # <-- adds borders around annotation cells
  
  gp = grid::gpar(
    col = "grey80",              # border color
    lwd = 0.5                   # border width
  ),
  
  annotation_name_gp = grid::gpar(fontsize = 12)
)


Heatmap(
  mat_scaled,
  name = "Z-score",
  
  top_annotation = ha,
  
  # Separate Day 7 and Day 14
  column_split = annotation_col$Day,
  
  # Cluster independently within each day
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  
  column_title = c("Day 7", "Day 14"),
  
  # Cluster genes
  cluster_rows = TRUE,
  
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  show_column_names = FALSE,
  show_row_names = TRUE,
  
  row_names_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  col = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  
  # Borders around individual heatmap cells
  rect_gp = grid::gpar(
    col = "grey80",
    lwd = 0.5
  ),
  
  # Gap between Day 7 and Day 14
  column_gap = grid::unit(4, "mm"),
  
  heatmap_legend_param = list(
    title = "Expression\nZ-score"
  )
)




#PCA Analysis
mat_pca <- t(mat_stable)
pca <- prcomp(mat_pca, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = cd$Group
)

# Define colors & shapes
group_colors <- c("Rejection" = "#F28500", "Acceptance" = "#064273" )
group_shapes <- c("Rejection" = 15, "Acceptance" = 16)

ggplot(pca_df, aes(PC1, PC2, color = Group,, shape = Group)) +
  geom_point(aes(fill = Group), size = 5, stroke = 1.2) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA of Stable Elastic Net Genes",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1],1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2],1), "%)")
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

### Gene Expressions Over Time----


vsd <- vst(dds_IsletTransplant_RejVsAccep, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_RejVsAccep))
design_temporal <- model.matrix(
  ~ Group * Day,
  data = cd
)

mat_temporal <- limma::removeBatchEffect(
  mat,
  batch = cd$Batch,
  design = design_temporal
)

mat_plot <- mat_temporal[
  stable_gene_names,
  ,
  drop = FALSE
]

# Check that metadata and expression samples align
stopifnot(
  identical(colnames(mat_plot), rownames(cd))
)


#  Convert expression matrix to long format
plot_df <- as.data.frame(mat_plot) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  left_join(
    cd %>%
      rownames_to_column("Sample") %>%
      select(Sample, Day, Group),
    by = "Sample"
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = c("7", "14")
    ),
    Group = factor(
      Group,
      levels = c(
        "Acceptance",
        "Rejection"
      )
    )
  )

# Check for missing metadata
stopifnot(
  !any(is.na(plot_df$Day)),
  !any(is.na(plot_df$Group))
)


gene_trajectory_plot <- ggplot(
  plot_df,
  aes(
    x = Day,
    y = Expression,
    color = Group,
    group = Group
  )
) +
  
  # Individual samples
  geom_point(
    aes(
      shape = Group,
      fill = Group
    ),
    position = position_jitter(
      width = 0.06,
      height = 0
    ),
    size = 2,
    alpha = 0.55,
    stroke = 0.5
  ) +
  
  # Group mean connecting Day 7 and Day 14
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.9
  ) +
  
  # Group mean symbols
  stat_summary(
    aes(
      shape = Group,
      fill = Group
    ),
    fun = mean,
    geom = "point",
    size = 3.3,
    stroke = 0.7
  ) +
  
  # Mean ± SEM
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.10,
    linewidth = 0.6
  ) +
  
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    ncol = 5
  ) +
  
  scale_color_manual(
    values = c(
      "Rejection" = "#F28500", 
      "Acceptance" = "#064273" 
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Rejection" = "#F28500", 
      "Acceptance" = "#064273" 
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Rejection" = 15,  # upward triangle
      "Acceptance" = 16  # downward triangle
    )
  ) +
  
  labs(
    title = "Temporal Expression- Anti-CD40L Allogeneic Rejection vs Acceptance Signature",
    x = "Day post-transplant",
    y = "Variance-stabilized (VST) expression",
    color = "Group",
    fill = "Group",
    shape = "Group"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    strip.text = element_text(
      face = "bold.italic",
      colour = "black",
      size = 12
    ),
    
    strip.background = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 0.6
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    panel.spacing = grid::unit(0.8, "lines")
  )

gene_trajectory_plot


## Auto-Allo Genes Signature Cross Preservation----
# subset samples of D7 and D14
sel <- colData(dds_IsletTransplant)$Group %in% c("Acceptance","Rejection") & colData(dds_IsletTransplant)$Day %in% c(7,14)
dds_IsletTransplant_RejVsAccep <- dds_IsletTransplant[, sel]

dds_IsletTransplant_RejVsAccep$Batch       <- factor(dds_IsletTransplant_RejVsAccep$Batch)
dds_IsletTransplant_RejVsAccep$LibraryPrep <- factor(dds_IsletTransplant_RejVsAccep$LibraryPrep)
dds_IsletTransplant_RejVsAccep$Group <- factor(
  dds_IsletTransplant_RejVsAccep$Group,
  levels = c("Acceptance", "Rejection")  # order sets baseline
)
levels(dds_IsletTransplant_RejVsAccep$Group)

dds_IsletTransplant_RejVsAccep$Day <- factor(dds_IsletTransplant_RejVsAccep$Day,
                                             levels = c(7, 14))  # baseline = 7
levels(dds_IsletTransplant_RejVsAccep$Day)

cd <- as.data.frame(colData(dds_IsletTransplant_RejVsAccep))
table(cd$Group,cd$Day)
#            7 14
#Acceptance 14 13
#Rejection   2  3
#Design Formula
# keep <- rowSums(counts(dds_IsletTransplant_RejVsAccep) >= 10) >= 5  # Total Samples in Rejection group used to account for group imbalance
# dds_IsletTransplant_RejVsAccep <- dds_IsletTransplant_RejVsAccep[keep, ]
design(dds_IsletTransplant_RejVsAccep) <- ~ Batch + Day + Group 
dds_IsletTransplant_RejVsAccep <- DESeq(dds_IsletTransplant_RejVsAccep)
design(dds_IsletTransplant_RejVsAccep)
resultsNames(dds_IsletTransplant_RejVsAccep)


vsd <- vst(dds_IsletTransplant_RejVsAccep, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_RejVsAccep))
design_enet <- model.matrix(~ Group, data = cd)
mat_enet <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  covariates = model.matrix(~ Day, data = cd)[, -1, drop = FALSE],
  design = design_enet
)

#Plot Heatmp and PCA using Selected Genes 
Auto_allo_signature <- c("Rsph10b","S1pr5","Erdr1","Pla2g4b","Clec2g","Ctsg","Syce1","Fndc7",
                         "Epgn","Tmem132e","Stkld1","Osbpl6","Sycp3")
setdiff(Auto_allo_signature, rownames(mat_enet))
mat_stable <- mat_enet[Auto_allo_signature, , drop = FALSE]
mat_scaled <- t(scale(t(mat_stable)))
annotation_col <- data.frame(
  Day = cd$Day,
  Group = cd$Group
)
rownames(annotation_col) <- colnames(mat_scaled)


# Column annotation
annotation_col <- data.frame(
  Day = factor(cd$Day, levels = c("7", "14")),
  Group = factor(
    cd$Group,
    levels = c("Acceptance", "Rejection")
  )
)

rownames(annotation_col) <- colnames(mat_scaled)

stopifnot(
  identical(rownames(annotation_col), colnames(mat_scaled))
)

# Top annotation
ha <- HeatmapAnnotation(
  df = annotation_col,
  col = list(
    Group = c(
      "Acceptance" = "#064273",
      "Rejection" = "#F28500"
    ),
    Day = c(
      "7" = "#F28E6B",
      "14" = "#6FA287"
    )
  ),
  
  border = TRUE,                # <-- adds borders around annotation cells
  
  gp = grid::gpar(
    col = "grey80",              # border color
    lwd = 0.5                   # border width
  ),
  
  annotation_name_gp = grid::gpar(fontsize = 12)
)


Heatmap(
  mat_scaled,
  name = "Z-score",
  
  top_annotation = ha,
  
  # Separate Day 7 and Day 14
  column_split = annotation_col$Day,
  
  # Cluster independently within each day
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  
  column_title = c("Day 7", "Day 14"),
  
  # Cluster genes
  cluster_rows = TRUE,
  
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  show_column_names = FALSE,
  show_row_names = TRUE,
  
  row_names_gp = grid::gpar(fontsize = 12),
  column_title_gp = grid::gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  col = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  
  # Borders around individual heatmap cells
  rect_gp = grid::gpar(
    col = "grey80",
    lwd = 0.5
  ),
  
  # Gap between Day 7 and Day 14
  column_gap = grid::unit(4, "mm"),
  
  heatmap_legend_param = list(
    title = "Expression\nZ-score"
  )
)

#PCA Analysis
mat_pca <- t(mat_stable)
pca <- prcomp(mat_pca, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = cd$Group
)

# Define colors & shapes
group_colors <- c("Rejection" = "#F28500", "Acceptance" = "#064273" )
group_shapes <- c("Rejection" = 15, "Acceptance" = 16)

ggplot(pca_df, aes(PC1, PC2, color = Group,, shape = Group)) +
  geom_point(aes(fill = Group), size = 5, stroke = 1.2) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA of Stable Elastic Net Genes",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1],1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2],1), "%)")
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )


### Gene Expressions Over Time

vsd <- vst(dds_IsletTransplant_RejVsAccep, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_RejVsAccep))
design_temporal <- model.matrix(
  ~ Group * Day,
  data = cd
)

mat_temporal <- limma::removeBatchEffect(
  mat,
  batch = cd$Batch,
  design = design_temporal
)

mat_plot <- mat_temporal[
  Auto_allo_signature,
  ,
  drop = FALSE
]

# Check that metadata and expression samples align
stopifnot(
  identical(colnames(mat_plot), rownames(cd))
)


#  Convert expression matrix to long format
plot_df <- as.data.frame(mat_plot) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  left_join(
    cd %>%
      rownames_to_column("Sample") %>%
      select(Sample, Day, Group),
    by = "Sample"
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = c("7", "14")
    ),
    Group = factor(
      Group,
      levels = c(
        "Acceptance",
        "Rejection"
      )
    )
  )

# Check for missing metadata
stopifnot(
  !any(is.na(plot_df$Day)),
  !any(is.na(plot_df$Group))
)


gene_trajectory_plot <- ggplot(
  plot_df,
  aes(
    x = Day,
    y = Expression,
    color = Group,
    group = Group
  )
) +
  
  # Individual samples
  geom_point(
    aes(
      shape = Group,
      fill = Group
    ),
    position = position_jitter(
      width = 0.06,
      height = 0
    ),
    size = 2,
    alpha = 0.55,
    stroke = 0.5
  ) +
  
  # Group mean connecting Day 7 and Day 14
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.9
  ) +
  
  # Group mean symbols
  stat_summary(
    aes(
      shape = Group,
      fill = Group
    ),
    fun = mean,
    geom = "point",
    size = 3.3,
    stroke = 0.7
  ) +
  
  # Mean ± SEM
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.10,
    linewidth = 0.6
  ) +
  
  facet_wrap(
    ~ Gene,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_color_manual(
    values = c(
      "Rejection" = "#F28500", 
      "Acceptance" = "#064273" 
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Rejection" = "#F28500", 
      "Acceptance" = "#064273" 
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Rejection" = 15,  # upward triangle
      "Acceptance" = 16  # downward triangle
    )
  ) +
  
  labs(
    title = "Temporal Expression- Auto-Allo NOD Signature in Anti-CD40L Allo B6",
    x = "Day post-transplant",
    y = "Variance-stabilized (VST) expression",
    color = "Group",
    fill = "Group",
    shape = "Group"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    strip.text = element_text(
      face = "bold.italic",
      colour = "black",
      size = 12
    ),
    
    strip.background = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 0.6
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    panel.spacing = grid::unit(0.8, "lines")
  )

gene_trajectory_plot


## GSEA Analysis----

# C8
CellTypeMSigDB_gene_sets <- msigdbr(species="Mus musculus", category="C8")
mm_c8_sets <- split(CellTypeMSigDB_gene_sets$gene_symbol, CellTypeMSigDB_gene_sets$gs_name)
mm_c8_df <- data.frame(
  gs_name = rep(names(mm_c8_sets), sapply(mm_c8_sets, length)),
  gene_symbol = unlist(mm_c8_sets)
)


# Hallmark
hallmark <- msigdbr(species = "Mus musculus", category  = "H")
mm_h_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
mm_h_df <- data.frame(
  gs_name = rep(names(mm_h_sets), sapply(mm_h_sets, length)),
  gene_symbol = unlist(mm_h_sets)
)


# KEGG
kegg_all <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_LEGACY")
mm_kegg_sets <- split(kegg_all$gene_symbol, kegg_all$gs_name)
mm_kegg_df <- data.frame(
  gs_name = rep(names(mm_kegg_sets), sapply(mm_kegg_sets, length)),
  gene_symbol = unlist(mm_kegg_sets)
)

# inal combined TERM2GENE data frame 
mm_all_df <- rbind(mm_c8_df, mm_h_df, mm_kegg_df)


# Paths to your saved results
ALL_REJVSACCEP_path  <- "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Rejection_Vs_Acceptance_aCD40L/DESEQResults_RejvsAccep_DaysCombined.csv"

# Import
ALL_REJVSACCEP  <- read.csv(ALL_REJVSACCEP_path,  row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_ALL_REJVSACCEP  <- ALL_REJVSACCEP$stat;  names(lfc_vector_ALL_REJVSACCEP)  <- rownames(ALL_REJVSACCEP)


# Drop NAs
lfc_vector_ALL_REJVSACCEP  <- lfc_vector_ALL_REJVSACCEP[!is.na(lfc_vector_ALL_REJVSACCEP)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_ALL_REJVSACCEP  <- sort(lfc_vector_ALL_REJVSACCEP,  decreasing = TRUE)

gsea_results_ALL_REJVSACCEP <- GSEA(
  geneList      = lfc_vector_ALL_REJVSACCEP,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_ALL_REJVSACCEP <- as.data.frame(gsea_results_ALL_REJVSACCEP)


# Full results Combined Timepoins
write.csv(gsea_results_ALL_REJVSACCEP,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Rejection_Vs_Acceptance_aCD40L/GSEAResults_RejvsAccep_DaysCombined.csv",
          row.names = FALSE)


immune_master_REJVSACCEP <- unique(c(
  #Upregukated
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_ALLOGRAFT_REJECTION",
  "KEGG_GRAFT_VERSUS_HOST_DISEASE",
  "KEGG_CELL_ADHESION_MOLECULES_CAMS",
  "KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY",
  "HE_LIM_SUN_FETAL_LUNG_C2_S100A12_HI_CLASSICAL_MONOCYTE",
  "DESCARTES_FETAL_PANCREAS_CCL19_CCL21_POSITIVE_CELLS",
  "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CELL",
  #Downregulated
  "AIZARANI_LIVER_C6_KUPFFER_CELLS_2",
  "HE_LIM_SUN_FETAL_LUNG_C2_APOE_POS_M2_MACROPHAGE_CELL",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "TRAVAGLINI_LUNG_NONCLASSICAL_MONOCYTE_CELL",
  "DESCARTES_FETAL_LIVER_MYELOID_CELLS",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
  "HAY_BONE_MARROW_DENDRITIC_CELL",
  "HALLMARK_PEROXISOME",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "AIZARANI_LIVER_C1_NK_NKT_CELLS_1"

  
))


plot_df <- gsea_results_ALL_REJVSACCEP %>%
  filter(ID %in% immune_master_REJVSACCEP) %>%
  mutate(
    neglog10_padj = -log10(p.adjust),
    neglog10_padj = ifelse(is.infinite(neglog10_padj), NA, neglog10_padj)
  ) %>%
  arrange(NES) %>%
  mutate(
    ID = factor(ID, levels = ID)
  )


ggplot(plot_df, aes(x = NES, y = ID, size = neglog10_padj, color = NES)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](adjusted~italic(p))), range = c(3, 10)) +
  scale_color_gradient2(
    low = "#064273",
    mid = "white",
    high = "#F28500",
    midpoint = 0,
    name = "NES"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 16) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Rejection vs Acceptance"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

## Islet Single Cell Mapping----

IsletScRNA = readRDS("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Islet ScRNASeq/islet_graft_seurat_v8.rds")
### Elastic Net Score ----
#Downregulated and Upregulated wrt Allo vs Syn

downregulated_RejVsAccep_Signature <- c("Ceacam19","Gp5","Tspoap1","Tymp","Pla2g2d","Treml1","Ltc4s","Fabp3","Sarm1","Cd5l","Ifi30")
upregulated_RejVsAccep_Signature <- c("Cpa6","Snhg11","Ebf4","Dok5","Cyp2j9","Flywch2","Fibin","Nrbp2","Hs3st3a1","Ptprf","Lurap1","Eda","Cmah","Prrg1")

downregulated_RejVsAccep_Signature_use <- intersect(downregulated_RejVsAccep_Signature, rownames(IsletScRNA))
upregulated_RejVsAccep_Signature_use <- intersect(upregulated_RejVsAccep_Signature, rownames(IsletScRNA))

combined_RejVsAccep_signature_use<-c(downregulated_RejVsAccep_Signature_use,upregulated_RejVsAccep_Signature_use)
DotPlot(IsletScRNA, features = combined_RejVsAccep_signature_use) +
  scale_color_gradient(low = "grey", high = "red") +
  scale_size(range = c(0, 8), limits = c(0, 100)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )


# Elastic Net Score
# Subset to allogeneic cells only
IsletScRNA_allo <- subset(IsletScRNA, subset = condition == "Allogeneic")


up_use <- intersect(upregulated_RejVsAccep_Signature_use, rownames(IsletScRNA_allo))
down_use <- intersect(downregulated_RejVsAccep_Signature_use, rownames(IsletScRNA_allo))

# Add module scores
IsletScRNA_allo <- AddModuleScore(
  object = IsletScRNA_allo,
  features = list(up_use),
  name = "RejectionScore",
  assay = DefaultAssay(IsletScRNA_allo)
)

IsletScRNA_allo <- AddModuleScore(
  object = IsletScRNA_allo,
  features = list(down_use),
  name = "AcceptanceScore",
  assay = DefaultAssay(IsletScRNA_allo)
)



plot_df <- IsletScRNA_allo@meta.data %>%
  dplyr::select(celltype, RejectionScore1, AcceptanceScore1) %>%
  filter(!is.na(celltype)) %>%
  pivot_longer(
    cols = c(RejectionScore1, AcceptanceScore1),
    names_to = "signature",
    values_to = "score"
  )

stat_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(score ~ signature)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(score, na.rm = TRUE) + 0.08 * diff(range(plot_df$score, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = paste0(signif(p_adj, 3))
  )



ggplot(plot_df, aes(x = celltype, y = score, fill = signature)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    data = stat_df,
    aes(x = celltype, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4.5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Module Score",
    title = "Distribution of Up- and Down-Signature Scores Across Cell Types"
  ) +
  theme(
    axis.title.x = element_text(size = 18),   # x-axis label
    axis.title.y = element_text(size = 18),   # y-axis label
    axis.text.x  = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 18),
    legend.title = element_blank()
  )+scale_fill_manual(
    values = c(
      "RejectionScore1" = "#F28500",
      "AcceptanceScore1" = "#064273"
    ),
    labels = c(
      "RejectionScore1" = "Rejection Score",
      "AcceptanceScore1" = "Acceptance Score"
    )
  )


# 6. Allo Acceptance Vs Syngeneic-----

sel <- colData(dds_IsletTransplant)$Group %in% c("Control Syngeneic","Acceptance")  & colData(dds_IsletTransplant)$Day %in% c(7,14,28,42,56,70)
#IS- Immunosuppressed
dds_IsletTransplant_AccepVsSyn <- dds_IsletTransplant[, sel]

dds_IsletTransplant_AccepVsSyn$Batch       <- factor(dds_IsletTransplant_AccepVsSyn$Batch)
dds_IsletTransplant_AccepVsSyn$LibraryPrep <- factor(dds_IsletTransplant_AccepVsSyn$LibraryPrep)
dds_IsletTransplant_AccepVsSyn$Group <- factor(
  dds_IsletTransplant_AccepVsSyn$Group,
  levels = c("Control Syngeneic", "Acceptance")  # order sets baseline
)
levels(dds_IsletTransplant_AccepVsSyn$Group)

dds_IsletTransplant_AccepVsSyn$Day <- factor(dds_IsletTransplant_AccepVsSyn$Day,
                                                  levels = c(7, 14, 28, 42, 56, 70))  # baseline = 7
levels(dds_IsletTransplant_AccepVsSyn$Day)

cd <- as.data.frame(colData(dds_IsletTransplant_AccepVsSyn))
table(cd$Group,cd$Day)# Check the minimum number of samples in a group

#                   7 14 28 42 56 70
#Control Syngeneic  5  4  5  5  5  5
#Acceptance        14 13 12 12 13 11

keep <- rowSums(counts(dds_IsletTransplant_AccepVsSyn) >= 10) >= 4  # Smallest number in Syngeneic group
dds_IsletTransplant_AccepVsSyn <- dds_IsletTransplant_AccepVsSyn[keep, ]
boxplot(logIEQ ~ Group, data = colData(dds_IsletTransplant_AccepVsSyn)) #Signficantly different between groups


## DESeq-Combined Timepoints----
#Design Formula
design(dds_IsletTransplant_AccepVsSyn) <- ~ Batch + Day + Group 
dds_IsletTransplant_AccepVsSyn <- DESeq(dds_IsletTransplant_AccepVsSyn)
design(dds_IsletTransplant_AccepVsSyn)
resultsNames(dds_IsletTransplant_AccepVsSyn)
summary(dds_IsletTransplant_AccepVsSyn)
# All day combined
res_ACCEPvSYN <- results(dds_IsletTransplant_AccepVsSyn,
                       name = "Group_Acceptance_vs_Control.Syngeneic")

# Thresholds
q_cut  <- 0.10
fc_cut <- 1


## 0) Prep results as data.frames with a 'gene' column
ALL_ACCEPVSSYN  <- as.data.frame(res_ACCEPvSYN);  ALL_ACCEPVSSYN$gene  <- rownames(res_ACCEPvSYN)

write.csv(ALL_ACCEPVSSYN, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/AcceptanceVsSyngeneic/DESEQResults__AccepVsSyn_DaysCombined.csv", row.names = TRUE)

# Function to create color mapping
make_keyvals_fdr_fc_AccepvsSyn <- function(df, q = 0.10, fc = 1,
                                           col_up = "#064273", col_down = "#2E6F40" , col_ns = "gray70") {
  # valid stats (for coloring); everything else becomes NS
  ok   <- !is.na(df$padj) & !is.na(df$log2FoldChange)
  
  up   <- ok & df$padj < q & df$log2FoldChange >=  fc
  down <- ok & df$padj < q & df$log2FoldChange <= -fc
  
  # default color + label
  key   <- rep(col_ns, nrow(df))
  label <- rep("Not significant", nrow(df))
  
  # overwrite where significant
  key[up]   <- col_up
  key[down] <- col_down
  
  label[up]   <- paste0("Allogeneic Acceptance")
  label[down] <- paste0("Syngeneic")
  
  names(key) <- label        # <- legend labels; no NAs
  key
}

keyvals_ALL  <- make_keyvals_fdr_fc_AccepvsSyn(ALL_ACCEPVSSYN)

#selLab_ALL_ACCEPVSSYN  <- pick_labels(ALL_ACCEPVSSYN,  q_cut, fc_cut, 30)

selLab_ALL_ACCEPVSSYN <- c(
  #Up
  "Cd14",
  "Marco",
  "Nos2",
  "Nlrp3",
  "Tnf",
  "Csf3r",
  "Cxcr2",
  "S100a8", 
  "S100a9", 
  "Mmp8",
  "Il17a",
  "Slamf6",
  "Themis",
  #DOwn
  "Ackr2",
  "Il24",
  "Vnn1"
)
## 3) Axis limits
xmax_ALL  <- max(2, ceiling(max(abs(ALL_ACCEPVSSYN$log2FoldChange),  na.rm=TRUE)))

EnhancedVolcano(
  ALL_ACCEPVSSYN,
  lab           = ALL_ACCEPVSSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Allogeneic Acceptance vs Syngeneic",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(ALL_ACCEPVSSYN$padj<0.1 & abs(ALL_ACCEPVSSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-10, 10),
  ylim = c(0,20),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 50,
  colCustom     = keyvals_ALL,
  legendPosition= "right",
  selectLab     = selLab_ALL_ACCEPVSSYN
)

###PCA Analysis-DE Genes----

sig_genes_ACCEPvSYN <- ALL_ACCEPVSSYN$gene[
  !is.na(ALL_ACCEPVSSYN$padj) &
    ALL_ACCEPVSSYN$padj <= 0.05 &
    abs(ALL_ACCEPVSSYN$log2FoldChange) >= 1
]

length(sig_genes_ACCEPvSYN)


vsd <- vst(dds_IsletTransplant_AccepVsSyn, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_AccepVsSyn))
design_heatmap <- model.matrix(~ Group, data = cd)
mat_bc <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  covariates = model.matrix(~ Day, data = cd)[, -1, drop = FALSE],
  design = design_heatmap
)


anno_col <- colData(dds_IsletTransplant_AccepVsSyn)[, c( "Group","Day")] |> 
  as.data.frame()


annotation_colors <- list(
  Group = c(
    "Acceptance" = "#064273" ,   # Ocean
    "Control Syngeneic" = "#2E6F40"     # Moss green
  ),
  Day = c(
    "7" = "#F28E6B",
    "14" = "#6FA287",
    "28" = "#E9C46A",  # soft muted gold (transition from early)
    "42" = "#F4A261",  # warm amber (still mid-stage)
    "56" = "#7B8CDE",  # muted periwinkle (shift to late)
    "70" = "#4C5D8B"   # deep slate blue (late stage)
  )
)

mat_bc_AccepVsSyn <- mat_bc[sig_genes_ACCEPvSYN,]

mat_pca <- t(mat_bc_AccepVsSyn)
pca <- prcomp(mat_pca, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = cd$Group
)

# Define colors & shapes
group_colors <- c("Control Syngeneic" = "#2E6F40", "Acceptance" = "#064273" )
group_shapes <- c("Control Syngeneic" = 24, "Acceptance" = 16)

ggplot(pca_df, aes(PC1, PC2, color = Group,, shape = Group)) +
  geom_point(aes(fill = Group), size = 5, stroke = 1.2) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA of DE Genes: Allo Acceptance vs Syn",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1],1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2],1), "%)")
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )



## GSEA Analysis----
# Paths to your saved results
ALL_ACCEPVSSYN_path  <-  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/AcceptanceVsSyngeneic/DESEQResults__AccepVsSyn_DaysCombined.csv"

# Import
ALL_ACCEPVSSYN  <- read.csv(ALL_ACCEPVSSYN_path,  row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_ALL_ACCEPVSSYN  <- ALL_ACCEPVSSYN$stat;  names(lfc_vector_ALL_ACCEPVSSYN)  <- rownames(ALL_ACCEPVSSYN)


# Drop NAs
lfc_vector_ALL_ACCEPVSSYN  <- lfc_vector_ALL_ACCEPVSSYN[!is.na(lfc_vector_ALL_ACCEPVSSYN)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_ALL_ACCEPVSSYN  <- sort(lfc_vector_ALL_ACCEPVSSYN,  decreasing = TRUE)

gsea_results_ALL_ACCEPVSSYN <- GSEA(
  geneList      = lfc_vector_ALL_ACCEPVSSYN,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_ALL_ACCEPVSSYN <- as.data.frame(gsea_results_ALL_ACCEPVSSYN)


# Full results Combined Timepoins
write.csv(gsea_results_ALL_ACCEPVSSYN,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/AcceptanceVsSyngeneic/GSEAResults_AccepvsSyn_DaysCombined.csv",
          row.names = FALSE)


immune_master_ACCEPVSSYN <- unique(c(
  #Upregukated
  "HAY_BONE_MARROW_NEUTROPHIL",
  "HAY_BONE_MARROW_IMMATURE_NEUTROPHIL",
  "TRAVAGLINI_LUNG_CLASSICAL_MONOCYTE_CELL",
  "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CELL",
  "TRAVAGLINI_LUNG_IGSF21_DENDRITIC_CELL",
  "FAN_OVARY_CL4_T_LYMPHOCYTE_NK_CELL_1",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_GRAFT_VERSUS_HOST_DISEASE",
  #Downregulated
  #"TRAVAGLINI_LUNG_CD4_NAIVE_T_CELL",
  "HALLMARK_KRAS_SIGNALING_DN",
  "HALLMARK_BILE_ACID_METABOLISM",
  "KEGG_GLUTATHIONE_METABOLISM",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_ANDROGEN_RESPONSE"

  
))


plot_df <- gsea_results_ALL_ACCEPVSSYN %>%
  filter(ID %in% immune_master_ACCEPVSSYN) %>%
  mutate(
    neglog10_padj = -log10(p.adjust),
    neglog10_padj = ifelse(is.infinite(neglog10_padj), NA, neglog10_padj)
  ) %>%
  arrange(NES) %>%
  mutate(
    ID = factor(ID, levels = ID)
  )


ggplot(plot_df, aes(x = NES, y = ID, size = neglog10_padj, color = NES)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](adjusted~italic(p))), range = c(3, 10)) +
  scale_color_gradient2(
    low =  "forestgreen",
    mid = "white",
    high =  "#064273",
    midpoint = 0,
    name = "NES"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 16) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Allogeneic Acceptance vs Syngeneic"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

## GSVA Scores Over Time----
sel <- colData(dds_IsletTransplant)$Group %in% c("Control Syngeneic","Acceptance")  & colData(dds_IsletTransplant)$Day %in% c(7,14,28,42,56,70)
#IS- Immunosuppressed
dds_IsletTransplant_AccepVsSyn <- dds_IsletTransplant[, sel]

dds_IsletTransplant_AccepVsSyn$Batch       <- factor(dds_IsletTransplant_AccepVsSyn$Batch)
dds_IsletTransplant_AccepVsSyn$LibraryPrep <- factor(dds_IsletTransplant_AccepVsSyn$LibraryPrep)
dds_IsletTransplant_AccepVsSyn$Group <- factor(
  dds_IsletTransplant_AccepVsSyn$Group,
  levels = c("Control Syngeneic", "Acceptance")  # order sets baseline
)
dds_IsletTransplant_AccepVsSyn$Day <- factor(
  dds_IsletTransplant_AccepVsSyn$Day,
  levels = c(7,14,28,42,56,70)
)

### Inflammatory Pathway Scores----

# VST normalization
vsd <- vst(dds_IsletTransplant_AccepVsSyn, blind = FALSE)

mat <- assay(vsd)

# Remove only sequencing batch effect
cd <- as.data.frame(colData(dds_IsletTransplant_AccepVsSyn))

cd$Day <- factor(cd$Day, levels = c(7,14,28,42,56,70))
cd$Group <- factor(cd$Group, levels = c("Control Syngeneic", "Acceptance"))
mat_gsva <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch
)

cd$MouseID <- paste(cd$Cohort, cd$Animal, sep = "_")

Inflammatory_Macrophage <- c(
  "Cd14",
  "Cd38",
  "Cd300lf",
  "Bst1",
  "Ccl3",
  "Ccl4",
  "Gpr84",
  "Hcar2",
  "Il1b",
  "Il1a",
  "Nos2",
  "Nlrp3",
  "Marco",
  "Ly6c2",
  "Acod1",
  "Ptgs2",
  "Tnf",
  "Trem1",
  "Lyz2",
  "Itgam",
  "Tyrobp",
  "Fcer1g",
  "Fcgr1",
  "Csf1r",
  "Ctss",
  "Clec7a",
  "Cebpb"
)

Inflammatory_Neutrophil  <- c(
  "Csf3r",
  "Cxcr2",
  "S100a8",
  "S100a9",
  "Lcn2",
  "Mmp8",
  "Trem1",
  "Fpr1",
  "Fpr2",
  "Sell",
  "Mcemp1",
  "Prok2",
  "Slpi",
  "Il1b",
  "Ptgs2",
  "Cxcl1",
  "Cxcl2",
  "Cxcl3"
)

Inflammatory_Tcell <- c(
  "Cd3d",
  "Cd3e",
  "Cd3g",
  "Cd4",
  "Themis",
  "Lef1",
  "Slamf6",
  "Tox",
  "Flicr",
  "Art2b",
  "Kcna3",
  "Gpr15",
  "Tnfsf14",
  "Tnfsf18",
  "Il17a",
  "Il17f",
  "Il22",
  "Csf2",
  "Rag1"
)

gene_sets_custom <- list(
  Inflammatory_Macrophage = Inflammatory_Macrophage,
  Inflammatory_Neutrophil = Inflammatory_Neutrophil,
  Inflammatory_Tcell=Inflammatory_Tcell
)

gene_sets_use <- c(gene_sets_custom, gene_sets_atlas)

gsva_par <- gsvaParam(
  mat_gsva,
  gene_sets_use,
  kcdf    = "Gaussian",   # VST is log-like
  minSize = 5,
  maxSize = 500
)
gsva_scores <- gsva(gsva_par)  # pathway x sample
gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$Sample <- rownames(gsva_df)

gsva_df <- cbind(gsva_df, cd[gsva_df$Sample, ])

pathways_use_way <- c(
  "Inflammatory_Macrophage",
  "Inflammatory_Neutrophil",
  "Inflammatory_Tcell"
)

gsva_long <- gsva_df %>%
  pivot_longer(
    cols = all_of(pathways_use_way),
    names_to = "Pathway",
    values_to = "Score"
  )


pvals_df <- gsva_long %>%
  dplyr::group_by(Pathway, Day) %>%
  dplyr::summarise(
    p_value = {
      df <- dplyr::cur_data()
      df$Group <- droplevels(df$Group)
      
      if (nlevels(df$Group) < 2) {
        NA_real_
      } else {
        tryCatch(
          t.test(Score ~ Group, data = df, var.equal = FALSE)$p.value,
          error = function(e) NA_real_
        )
      }
    },
    .groups = "drop"
  )

pvals_df <- pvals_df %>%
  dplyr::mutate(
    signif_label = dplyr::case_when(
      is.na(p_value) ~ "",
      p_value <= 1e-4 ~ "****",
      p_value <= 1e-3 ~ "***",
      p_value <= 1e-2 ~ "**",
      p_value <= 5e-2 ~ "*",
      TRUE ~ ""
    )
  )

gsva_avg <- gsva_long %>%
  dplyr::group_by(Group, Day, Pathway) %>%
  dplyr::summarise(
    mean_score = mean(Score),
    se = sd(Score) / sqrt(dplyr::n()),
    .groups = "drop"
  )

label_pos_df <- gsva_long %>%
  dplyr::group_by(Pathway, Day) %>%
  dplyr::summarise(
    y_pos = max(Score, na.rm = TRUE) + 0.08 * diff(range(Score, na.rm = TRUE)),
    .groups = "drop"
  )

label_pos_df <- gsva_long %>%
  dplyr::group_by(Pathway, Day) %>%
  dplyr::summarise(
    ymin = min(Score, na.rm = TRUE),
    ymax = max(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    span = dplyr::if_else(ymax > ymin, ymax - ymin, 0.2),
    y_pos = ymax - 0.2 * span
  ) %>%
  dplyr::select(Pathway, Day, y_pos)

sig_df <- pvals_df %>%
  dplyr::left_join(label_pos_df, by = c("Pathway", "Day")) %>%
  dplyr::filter(!is.na(p_value), signif_label != "")

ggplot(gsva_avg,
       aes(x = Day,
           y = mean_score,
           color = Group,
           shape = Group,
           group = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_score - se,
        ymax = mean_score + se),
    width = 0.5
  ) +
  geom_text(
    data = sig_df,
    aes(x = Day, y = y_pos, label = signif_label),
    inherit.aes = FALSE,
    size = 7
  ) +
  facet_wrap(
    ~ Pathway,
    labeller = as_labeller(c(
      "Inflammatory_Macrophage" = "Inflammatory Macrophage",
      "Inflammatory_Neutrophil" = "Inflammatory Neutrophil",
      "Inflammatory_Tcell" = "Inflammatory T Cell"
    ))
  ) +
  coord_cartesian(ylim = c(-0.65,0.65)) +
  scale_color_manual(values = c(
    "Control Syngeneic" = "#2E6F40",
    "Acceptance" = "#064273"
  )) +
  scale_shape_manual(values = c(
    "Control Syngeneic" = 17,
    "Acceptance" = 16
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = "Days Post Transplant",
    y = "GSVA score",
    title = "Immune Inflammatory Genesets"
  )+
  theme(
    strip.text = element_text(size = 12)
  )
# Select relevant columns
gsva_export <- gsva_df %>%
  dplyr::select(
    Sample,
    MouseID,
    Group,
    Day,
    all_of(pathways_use_way)
  )
library(writexl)
# Write to CSV
write_xlsx(gsva_export,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Tab7_ Allo_Acceptance_vs_Syn_LongitudinalGSVAScores.xlsx")

# 7.Allo+Anti-CD40L Rejection vs Allo -----

# subset samples
sel <- colData(dds_IsletTransplant)$Group %in% c("Control Allogeneic","Rejection") & colData(dds_IsletTransplant)$Day %in% c(7,14)
#IS- Immunosuppressed
dds_IsletTransplant_IS_RejVsAllo <- dds_IsletTransplant[, sel]

dds_IsletTransplant_IS_RejVsAllo$Batch       <- factor(dds_IsletTransplant_IS_RejVsAllo$Batch)
dds_IsletTransplant_IS_RejVsAllo$LibraryPrep <- factor(dds_IsletTransplant_IS_RejVsAllo$LibraryPrep)
dds_IsletTransplant_IS_RejVsAllo$Group <- factor(
  dds_IsletTransplant_IS_RejVsAllo$Group,
  levels = c("Control Allogeneic", "Rejection")  # order sets baseline
)
levels(dds_IsletTransplant_IS_RejVsAllo$Group)

dds_IsletTransplant_IS_RejVsAllo$Day <- factor(dds_IsletTransplant_IS_RejVsAllo$Day,
                                                  levels = c(7, 14))  # baseline = 7
levels(dds_IsletTransplant_IS_RejVsAllo$Day)

cd <- as.data.frame(colData(dds_IsletTransplant_IS_RejVsAllo))
table(cd$Group,cd$Day)
#                   7 14
#Control Allogeneic 5  4
#Rejection          2  3

#Check if logIEQ varies systemically with Group or Day
boxplot(logIEQ ~ Group, data = colData(dds_IsletTransplant_IS_RejVsAllo)) #Signficantly different between groups

#Design Formula
keep <- rowSums(counts(dds_IsletTransplant_IS_RejVsAllo) >= 10) >= 2   # Smallest number of samples in a group
dds_IsletTransplant_IS_RejVsAllo <- dds_IsletTransplant_IS_RejVsAllo[keep, ]

design(dds_IsletTransplant_IS_RejVsAllo) <- ~ Batch + Day + Group 
dds_IsletTransplant_IS_RejVsAllo <- DESeq(dds_IsletTransplant_IS_RejVsAllo)
design(dds_IsletTransplant_IS_RejVsAllo)
resultsNames(dds_IsletTransplant_IS_RejVsAllo)

# Genes with statistics
res_IS_REJvALLO <- results(dds_IsletTransplant_IS_RejVsAllo,
                          name = "Group_Rejection_vs_Control.Allogeneic")

summary(res_IS_REJvALLO)


# Thresholds
q_cut  <- 0.05
fc_cut <- 1


ALL_IS_REJVSALLO  <- as.data.frame(res_IS_REJvALLO);  ALL_IS_REJVSALLO$gene  <- rownames(res_IS_REJvALLO)

write.csv(ALL_IS_REJVSALLO, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Rejection_Vs_Allogeneic/DESEQResults_RejectionVsAllogeneic.csv", row.names = TRUE)

# Function to create color mapping
make_keyvals_fdr_fc_IsRejvsAllo <- function(df, q = 0.10, fc = 1,
                                           col_up =  "#F28500", col_down = "#FF2400" , col_ns = "gray70") {
  # valid stats (for coloring); everything else becomes NS
  ok   <- !is.na(df$padj) & !is.na(df$log2FoldChange)
  
  up   <- ok & df$padj < q & df$log2FoldChange >=  fc
  down <- ok & df$padj < q & df$log2FoldChange <= -fc
  
  # default color + label
  key   <- rep(col_ns, nrow(df))
  label <- rep("Not significant", nrow(df))
  
  # overwrite where significant
  key[up]   <- col_up
  key[down] <- col_down
  
  label[up]   <- paste0("Anti-CD40L Rejection")
  label[down] <- paste0("Control Allogeneic")
  
  names(key) <- label        # <- legend labels; no NAs
  key
}

keyvals_ALL  <- make_keyvals_fdr_fc_IsRejvsAllo(res_IS_REJvALLO)


selLab_ALL_IS_REJVSALLO  <- pick_labels(ALL_IS_REJVSALLO,  q_cut, fc_cut, 40)


## 3) Axis limits
xmax_ALL  <- max(2, ceiling(max(abs(res_IS_REJvALLO$log2FoldChange),  na.rm=TRUE)))

EnhancedVolcano(
  res_IS_REJvALLO,
  lab           = ALL_IS_REJVSALLO$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.05,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Anti-CD40L Rejection vs Control Rejection",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(res_IS_REJvALLO$padj<0.05 & abs(res_IS_REJvALLO$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_ALL, xmax_ALL),
  ylim = c(0,10),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 35,
  colCustom     = keyvals_ALL,
  legendPosition= "right",
  selectLab     = selLab_ALL_IS_REJVSALLO
)

## #PCA Analysis-DE Genes----

sig_genes_IS_REJVSALLO<- ALL_IS_REJVSALLO$gene[
  !is.na(ALL_IS_REJVSALLO$padj) &
    ALL_IS_REJVSALLO$padj <= 0.05 &
    abs(ALL_IS_REJVSALLO$log2FoldChange) >= 1
]

length(sig_genes_IS_REJVSALLO)
#19

vsd <- vst(dds_IsletTransplant_IS_RejVsAllo, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_IS_RejVsAllo))
design_heatmap <- model.matrix(~ Group, data = cd)
mat_bc <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  covariates = model.matrix(~ Day, data = cd)[, -1, drop = FALSE],
  design = design_heatmap
)


anno_col <- colData(dds_IsletTransplant_IS_RejVsAllo)[, c( "Group","Day")] |> 
  as.data.frame()


annotation_colors <- list(
  Group = c(
    "Rejection" = "#F28500" ,   # Tangerine
    "Control Allogeneic" = "#FF2400"    # Scarlet
  ),
  Day = c(
    "7" = "#F28E6B",
    "14" = "#6FA287"
  )
)

mat_bc_RejVsAllo <- mat_bc[sig_genes_IS_REJVSALLO,]

mat_pca <- t(mat_bc_RejVsAllo)
pca <- prcomp(mat_pca, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = cd$Group
)

# Define colors & shapes
group_colors <- c("Control Allogeneic" ="#FF2400" , "Rejection" = "#F28500" )
group_shapes <- c("Control Allogeneic" = 25, "Rejection" = 15)



ggplot(pca_df, aes(PC1, PC2, color = Group,, shape = Group)) +
  geom_point(aes(fill = Group), size = 5, stroke = 1.2) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA of DE Genes: Allo Rejection vs Control Allo",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1],1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2],1), "%)")
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

## GSEA Analysis----

# C8
CellTypeMSigDB_gene_sets <- msigdbr(species="Mus musculus", category="C8")
mm_c8_sets <- split(CellTypeMSigDB_gene_sets$gene_symbol, CellTypeMSigDB_gene_sets$gs_name)
mm_c8_df <- data.frame(
  gs_name = rep(names(mm_c8_sets), sapply(mm_c8_sets, length)),
  gene_symbol = unlist(mm_c8_sets)
)


# Hallmark
hallmark <- msigdbr(species = "Mus musculus", category  = "H")
mm_h_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
mm_h_df <- data.frame(
  gs_name = rep(names(mm_h_sets), sapply(mm_h_sets, length)),
  gene_symbol = unlist(mm_h_sets)
)


# KEGG
kegg_all <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_LEGACY")
mm_kegg_sets <- split(kegg_all$gene_symbol, kegg_all$gs_name)
mm_kegg_df <- data.frame(
  gs_name = rep(names(mm_kegg_sets), sapply(mm_kegg_sets, length)),
  gene_symbol = unlist(mm_kegg_sets)
)

mm_all_df <- rbind(mm_c8_df, mm_h_df, mm_kegg_df)


# Paths to your saved results
ALL_IS_REJVSALLO_path  <- "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Rejection_Vs_Allogeneic/DESEQResults_RejectionVsAllogeneic.csv"

# Import
ALL_IS_REJVSALLO  <- read.csv(ALL_IS_REJVSALLO_path,  row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_ISREJVSALLO  <- ALL_IS_REJVSALLO$stat;  names(lfc_vector_ISREJVSALLO)  <- rownames(ALL_IS_REJVSALLO)

# Drop NAs
lfc_vector_ISREJVSALLO  <- lfc_vector_ISREJVSALLO[!is.na(lfc_vector_ISREJVSALLO)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_ISREJVSALLO  <- sort(lfc_vector_ISREJVSALLO,  decreasing = TRUE)

gsea_results_ISREJVSALLO <- GSEA(
  geneList      = lfc_vector_ISREJVSALLO,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_ISREJVSALLO <- as.data.frame(gsea_results_ISREJVSALLO)


# Full results Combined Timepoints
write.csv(gsea_results_ISREJVSALLO,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Rejection_Vs_Allogeneic/GSEAResults_ISRejectionVsAllogeneic_DaysCombined.csv",
          row.names = FALSE)




immune_master_ISREJVSALLO <- unique(c(
  #DOWNREGULATED
  "HE_LIM_SUN_FETAL_LUNG_C5_PRO_B_CELL",
  "TRAVAGLINI_LUNG_PROLIFERATING_NK_T_CELL",
  "DESCARTES_FETAL_PANCREAS_LYMPHOID_CELLS",
  "DESCARTES_FETAL_LIVER_LYMPHOID_CELLS",
  "TRAVAGLINI_LUNG_CD8_NAIVE_T_CELL",
  "AIZARANI_LIVER_C1_NK_NKT_CELLS_1",
  "TRAVAGLINI_LUNG_CD4_NAIVE_T_CELL",
  "AIZARANI_LIVER_C12_NK_NKT_CELLS_4",
  "TRAVAGLINI_LUNG_CD8_MEMORY_EFFECTOR_T_CELL",
  "TRAVAGLINI_LUNG_NATURAL_KILLER_CELL",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "TRAVAGLINI_LUNG_NATURAL_KILLER_T_CELL",
  #UPREGULATED
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CELL",
  "HAY_BONE_MARROW_NEUTROPHIL",
  "HAY_BONE_MARROW_IMMATURE_NEUTROPHIL",
  "TRAVAGLINI_LUNG_CLASSICAL_MONOCYTE_CELL",
  "AIZARANI_LIVER_C23_KUPFFER_CELLS_3",
  "DESCARTES_FETAL_PANCREAS_MYELOID_CELLS",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CELL",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HE_LIM_SUN_FETAL_LUNG_C2_S100A12_HI_CLASSICAL_MONOCYTE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "HALLMARK_APOPTOSIS",
  "KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY",
  "TRAVAGLINI_LUNG_IGSF21_DENDRITIC_CELL",
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY"
))

# Helper to standardize clusterProfiler GSEA columns
coerce_gsea_tbl <- function(df, day_label){
  gs  <- if ("Description" %in% names(df)) df$Description else if ("ID" %in% names(df)) df$ID else if ("pathway" %in% names(df)) df$pathway else if ("setName" %in% names(df)) df$setName else rownames(df)
  pad <- if ("p.adjust"   %in% names(df)) df$p.adjust   else if ("padj" %in% names(df)) df$padj else df$pval
  tibble(
    gs_name = as.character(gs),
    NES     = as.numeric(df$NES),
    padj    = as.numeric(pad),
    Day     = day_label
  )
}

plot_df <- gsea_results_ISREJVSALLO %>%
  filter(ID %in% immune_master_ISREJVSALLO) %>%
  mutate(
    neglog10_padj = -log10(p.adjust),
    neglog10_padj = ifelse(is.infinite(neglog10_padj), NA, neglog10_padj)
  ) %>%
  arrange(NES) %>%
  mutate(
    ID = factor(ID, levels = ID)
  )


ggplot(plot_df, aes(x = NES, y = ID, size = neglog10_padj, color = NES)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](adjusted~italic(p))), range = c(3, 10)) +
  scale_color_gradient2(
    low = "#FF2400",
    mid = "white",
    high = "#F28500",
    midpoint = 0,
    name = "NES"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 16) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Anti-CD40L Allogeneic Rejection vs Allogeneic"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )




#8. IN vs Liver Flow Surrogate Analysis----

### Allo vs Syn Transplant----
library(readxl)
raw_df_syn <- read_excel(
  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Table3_LiverVsScaf.xlsx",
  sheet = "Syngeneic"
)
raw_df_allo <- read_excel(
  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Table3_LiverVsScaf.xlsx",
  sheet = "Allogeneic (without aCD40L)"
)

# Keep only the columns present in the Syngeneic sheet
raw_df_allo <- raw_df_allo %>%
  dplyr::select(all_of(names(raw_df_syn)))

# Combine
raw_df <- dplyr::bind_rows(
  raw_df_syn,
  raw_df_allo
)

# Check
dim(raw_df)
table(raw_df$Group)

# Inspect imported column names
names(raw_df)
marker_cols <- names(raw_df)[9:14]

# Rename columns for easier handling
names(raw_df)[1:16] <- c(
  "Samples",
  "Mouse",
  "Batch",
  "Group",
  "Time",
  "ExperimentDate",
  "RecipientSpecies",
  "Tissue",
  "Immune cells",
  "Neutrophils",
  "Macrophages",
  "T cells",
  "CD4 T cells",
  "CD8 T cells",
  "B cells",
  "Dendritic cells"
)


df <- raw_df %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "Syngeneic",
        "Allogeneic"
      )
    ),
    Tissue = factor(
      Tissue,
      levels = c(
        "Scaffold",
        "Liver"
      )
    )
  )

marker_names <- c(
  "Immune cells",
  "Neutrophils",
  "Macrophages",
  "T cells",
  "CD4 T cells",
  "CD8 T cells",
  "B cells",
  "Dendritic cells"
)

paired_long <- df %>%
  pivot_longer(
    cols = all_of(marker_names),
    names_to = "Marker",
    values_to = "Frequency"
  ) %>%
  select(
    Mouse,
    Group,
    Time,
    Batch,
    Tissue,
    Marker,
    Frequency
  ) %>%
  pivot_wider(
    names_from = Tissue,
    values_from = Frequency
  ) %>%
  drop_na(
    Scaffold,
    Liver
  )

paired_long %>%
  dplyr::count(Marker, Group)


marker_axis_labels <- c(
  "Immune cells" = "CD45+ immune cells (% of live cells)",
  "Neutrophils" = "Ly6G+ neutrophils (% of CD45+ cells)",
  "Macrophages" = "F4/80+ macrophages (% of CD45+ cells)",
  "T cells" = "CD3+ T cells (% of CD45+ cells)",
  "CD4 T cells" = "CD4+ T cells (% of CD3+ cells)",
  "CD8 T cells" = "CD8+ T cells (% of CD3+ cells)",
  "B cells" = "CD19+ B cells (% of CD45+ cells)",
  "Dendritic cells" = "CD11c+MHC-II+ dendritic cells (% of CD45+ cells)"
)


make_correlation_plot <- function(marker_name, data = paired_long) 
  {
  marker_df <- data %>%
    dplyr::filter(
      Marker == marker_name,
      !is.na(Liver),
      !is.na(Scaffold)
    )
  
  if (nrow(marker_df) < 3) {
    stop(
      paste0(
        "Not enough paired observations for marker: ",
        marker_name
      )
    )
  }
  
  # Calculate Pearson correlation separately for each group
  cor_stats <- marker_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      cor_test = list(
        if (dplyr::n() >= 3) {
          cor.test(
            Liver,
            Scaffold,
            method = "pearson"
          )
        } else {
          NULL
        }
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      r = purrr::map_dbl(
        cor_test,
        ~ if (is.null(.x)) NA_real_ else unname(.x$estimate)
      ),
      p = purrr::map_dbl(
        cor_test,
        ~ if (is.null(.x)) NA_real_ else .x$p.value
      ),
      cor_label = dplyr::if_else(
        is.na(r),
        paste0(Group, ": insufficient observations"),
        paste0(
          Group,
          ": r = ",
          sprintf("%.2f", r),
          ", P = ",
          format.pval(
            p,
            digits = 2,
            eps = 0.001
          )
        )
      )
    )
  
  # Combine group-specific statistics into one annotation
  stats_label <- paste(
    cor_stats$cor_label,
    collapse = "\n"
  )
  
  axis_label <- marker_axis_labels[[marker_name]]
  
  if (is.null(axis_label)) {
    axis_label <- marker_name
  }
  
  ggplot(
    marker_df,
    aes(
      x = Liver,
      y = Scaffold,
      color = Group,
      shape = Group,
      fill = Group
    )
  ) +
    
    # Separate regression line for each group
    geom_smooth(
      aes(
        group = Group,
        color = Group,
        fill = Group
      ),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 1,
      alpha = 0.20
    ) +
    
    geom_point(
      size = 4.5,
      color = "black",
      stroke = 0.8
    ) +
    
    scale_color_manual(
      values = c(
        "Syngeneic" = "#2E6F40",
        "Allogeneic" = "#FF2400"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Syngeneic" = "#2E6F40",
        "Allogeneic" = "#FF2400"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "Syngeneic" = 24,
        "Allogeneic" = 25
      )
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = stats_label,
      hjust = 1.05,
      vjust = 1.2,
      size = 4.5,
      lineheight = 1.2
    ) +
    
    labs(
      title = marker_name,
      x = paste0("Liver ", axis_label),
      y = paste0("IN ", axis_label),
      color = "Group",
      fill = "Group",
      shape = "Group"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 14,
        face = "bold"
      ),
      
      axis.text = element_text(
        size = 13,
        colour = "black"
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      
      axis.ticks.length = grid::unit(
        0.25,
        "cm"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      )
    )
}
markers_to_plot <- unique(paired_long$Marker)

correlation_plots <- setNames(
  lapply(
    markers_to_plot,
    make_correlation_plot
  ),
  markers_to_plot
)

correlation_plots


### Allo With vs Without-Anti-CD40L Transplant----
library(readxl)

raw_df_allo <- read_excel(
  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Table3_LiverVsScaf.xlsx",
  sheet = "Allogeneic (without aCD40L)"
)

raw_df_alloaCD40L <- read_excel(
  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Table3_LiverVsScaf.xlsx",
  sheet = "Allogeneic (aCD40L)"
)

# Keep only the columns present in the Syngeneic sheet
raw_df_allo <- raw_df_allo %>%
  dplyr::select(all_of(names(raw_df_alloaCD40L)))

# Combine
raw_df <- dplyr::bind_rows(
  raw_df_alloaCD40L,
  raw_df_allo
)

# Check
dim(raw_df)
table(raw_df$Group)

# Inspect imported column names
names(raw_df)
marker_cols <- names(raw_df)[9:18]

# Rename columns for easier handling
names(raw_df)[1:18] <- c(
  "Samples",
  "Mouse",
  "Batch",
  "Group",
  "Time",
  "ExperimentDate",
  "RecipientSpecies",
  "Tissue",
  "Immune cells",
  "Neutrophils",
  "Macrophages",
  "CX3CR1+ Macrophages",
  "T cells",
  "CD4 T cells",
  "CD8 T cells",
  "B cells",
  "Dendritic cells",
  "CX3CR1+ DCs"
)


df <- raw_df %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "Allogeneic_aCD40L",
        "Allogeneic"
      )
    ),
    Tissue = factor(
      Tissue,
      levels = c(
        "Scaffold",
        "Liver"
      )
    )
  )

marker_names <- c(
  "Immune cells",
  "Neutrophils",
  "Macrophages",
  "CX3CR1+ Macrophages",
  "T cells",
  "CD4 T cells",
  "CD8 T cells",
  "B cells",
  "Dendritic cells",
  "CX3CR1+ DCs"
)

paired_long <- df %>%
  pivot_longer(
    cols = all_of(marker_names),
    names_to = "Marker",
    values_to = "Frequency"
  ) %>%
  select(
    Mouse,
    Group,
    Time,
    Batch,
    Tissue,
    Marker,
    Frequency
  ) %>%
  pivot_wider(
    names_from = Tissue,
    values_from = Frequency
  ) %>%
  drop_na(
    Scaffold,
    Liver
  )

paired_long %>%
  dplyr::count(Marker, Group)


marker_axis_labels <- c(
  "Immune cells" = "CD45+ immune cells (% of live cells)",
  "Neutrophils" = "Ly6G+ neutrophils (% of CD45+ cells)",
  "Macrophages" = "F4/80+ macrophages (% of CD45+ cells)",
  "CX3CR1+ Macrophages" = "F4/80+ macrophages (% of CD45+ cells)",
  "T cells" = "CD3+ T cells (% of CD45+ cells)",
  "CD4 T cells" = "CD4+ T cells (% of CD3+ cells)",
  "CD8 T cells" = "CD8+ T cells (% of CD3+ cells)",
  "B cells" = "CD19+ B cells (% of CD45+ cells)",
  "Dendritic cells" = "CD11c+MHC-II+ dendritic cells (% of CD45+ cells)",
  "CX3CR1+ DCs" = "CX3CR1+ DCs (% of CD45+ cells)"
)


make_correlation_plot <- function(marker_name, data = paired_long) 
{
  marker_df <- data %>%
    dplyr::filter(
      Marker == marker_name,
      !is.na(Liver),
      !is.na(Scaffold)
    )
  
  if (nrow(marker_df) < 3) {
    stop(
      paste0(
        "Not enough paired observations for marker: ",
        marker_name
      )
    )
  }
  
  # Calculate Pearson correlation separately for each group
  cor_stats <- marker_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      cor_test = list(
        if (dplyr::n() >= 3) {
          cor.test(
            Liver,
            Scaffold,
            method = "pearson"
          )
        } else {
          NULL
        }
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      r = purrr::map_dbl(
        cor_test,
        ~ if (is.null(.x)) NA_real_ else unname(.x$estimate)
      ),
      p = purrr::map_dbl(
        cor_test,
        ~ if (is.null(.x)) NA_real_ else .x$p.value
      ),
      cor_label = dplyr::if_else(
        is.na(r),
        paste0(Group, ": insufficient observations"),
        paste0(
          Group,
          ": r = ",
          sprintf("%.2f", r),
          ", P = ",
          format.pval(
            p,
            digits = 2,
            eps = 0.001
          )
        )
      )
    )
  
  # Combine group-specific statistics into one annotation
  stats_label <- paste(
    cor_stats$cor_label,
    collapse = "\n"
  )
  
  axis_label <- marker_axis_labels[[marker_name]]
  
  if (is.null(axis_label)) {
    axis_label <- marker_name
  }
  
  ggplot(
    marker_df,
    aes(
      x = Liver,
      y = Scaffold,
      color = Group,
      shape = Group,
      fill = Group
    )
  ) +
    
    # Separate regression line for each group
    geom_smooth(
      aes(
        group = Group,
        color = Group,
        fill = Group
      ),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 1,
      alpha = 0.20
    ) +
    
    geom_point(
      size = 4.5,
      color = "black",
      stroke = 0.8
    ) +
    
    scale_color_manual(
      values = c(
        "Allogeneic_aCD40L" = "#064273",
        "Allogeneic" = "#FF2400"
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Allogeneic_aCD40L" = "#064273",
        "Allogeneic" = "#FF2400"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Allogeneic_aCD40L" = 21,
        "Allogeneic" = 25
      )
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = stats_label,
      hjust = 1.05,
      vjust = 1.2,
      size = 4.5,
      lineheight = 1.2
    ) +
    
    labs(
      title = marker_name,
      x = paste0("Liver ", axis_label),
      y = paste0("IN ", axis_label),
      color = "Group",
      fill = "Group",
      shape = "Group"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 14,
        face = "bold"
      ),
      
      axis.text = element_text(
        size = 13,
        colour = "black"
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      
      axis.ticks.length = grid::unit(
        0.25,
        "cm"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      )
    )
}



markers_to_plot <- unique(paired_long$Marker)

correlation_plots <- setNames(
  lapply(
    markers_to_plot,
    make_correlation_plot
  ),
  markers_to_plot
)

correlation_plots




### Auto-Allo Transplant----
library(readxl)
raw_df <- read_excel(
  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Table3_LiverVsScaf.xlsx",
  sheet = "NOD Allogeneic"
)

# Inspect imported column names
names(raw_df)
marker_cols <- names(raw_df)[9:15]

# Rename columns for easier handling
names(raw_df)[1:15] <- c(
  "SampleID",
  "Mouse",
  "Group",
  "Time",
  "Batch",
  "ExperimentDate",
  "RecipientSpecies",
  "Tissue",
  "Immune cells",
  "Myeloid cells",
  "Neutrophils",
  "Macrophages",
  "T cells",
  "B cells",
  "Dendritic cells"
)


df <- raw_df %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "Early Rejection",
        "Late Rejection"
      )
    ),
    Tissue = factor(
      Tissue,
      levels = c(
        "Scaffold",
        "Liver"
      )
    )
  )

marker_names <- c(
  "Immune cells",
  "Myeloid cells",
  "Neutrophils",
  "Macrophages",
  "T cells",
  "B cells",
  "Dendritic cells"
)

paired_long <- df %>%
  pivot_longer(
    cols = all_of(marker_names),
    names_to = "Marker",
    values_to = "Frequency"
  ) %>%
  select(
    Mouse,
    Group,
    Time,
    Batch,
    Tissue,
    Marker,
    Frequency
  ) %>%
  pivot_wider(
    names_from = Tissue,
    values_from = Frequency
  ) %>%
  drop_na(
    Scaffold,
    Liver
  )

paired_long %>%
  dplyr::count(Marker, Group)


marker_axis_labels <- c(
  "Immune cells" = "CD45+ immune cells (% of live cells)",
  "Myeloid cells" = "CD11b+ myeloid cells (% of CD45+ cells)",
  "Neutrophils" = "Ly6G+ neutrophils (% of CD45+ cells)",
  "Macrophages" = "F4/80+ macrophages (% of CD45+ cells)",
  "T cells" = "CD3+ T cells (% of CD45+ cells)",
  "B cells" = "CD19+ B cells (% of CD45+ cells)",
  "Dendritic cells" = "CD11c+MHC-II+ dendritic cells (% of CD45+ cells)"
)


make_correlation_plot <- function(marker_name, data = paired_long) {
  
  marker_df <- data %>%
    dplyr::filter(Marker == marker_name) %>%
    dplyr::filter(
      !is.na(Liver),
      !is.na(Scaffold)
    )
  
  if (nrow(marker_df) < 3) {
    stop(
      paste0(
        "Not enough paired observations for marker: ",
        marker_name
      )
    )
  }
  
  cor_result <- cor.test(
    marker_df$Liver,
    marker_df$Scaffold,
    method = "pearson"
  )
  
  cor_label <- paste0(
    "Pearson r = ",
    sprintf("%.2f", unname(cor_result$estimate)),
    "\nP = ",
    format.pval(
      cor_result$p.value,
      digits = 2,
      eps = 0.001
    ),
    "\nn = ",
    nrow(marker_df)
  )
  
  axis_label <- marker_axis_labels[[marker_name]]
  
  if (is.null(axis_label)) {
    axis_label <- marker_name
  }
  
  ggplot(
    marker_df,
    aes(
      x = Liver,
      y = Scaffold,
      color = Group,
      shape = Group
    )
  ) +
    
    # One overall regression line
    geom_smooth(
      data = marker_df,
      aes(
        x = Liver,
        y = Scaffold,
        group = 1
      ),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "black",
      fill = "grey80",
      linewidth = 1,
      alpha = 0.3
    ) +
    
    geom_point(
      size = 4.5
    ) +
    
    scale_color_manual(
      values = c(
        "Early Rejection" = "#B23A48",
        "Late Rejection" = "#2A6F97"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "Early Rejection" = 15,
        "Late Rejection" = 16
      )
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = cor_label,
      hjust = 1.1,
      vjust = 1.2,
      size = 5
    ) +
    
    labs(
      title = marker_name,
      x = paste0("Liver ", axis_label),
      y = paste0("IN ", axis_label),
      color = "Group",
      shape = "Group"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 14,
        face = "bold"
      ),
      
      axis.text = element_text(
        size = 13,
        colour = "black"
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      
      axis.ticks.length = grid::unit(
        0.25,
        "cm"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      )
    )
}
markers_to_plot <- unique(paired_long$Marker)

correlation_plots <- setNames(
  lapply(
    markers_to_plot,
    make_correlation_plot
  ),
  markers_to_plot
)

correlation_plots


### Baseline----
library(readxl)
raw_df <- read_excel(
  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Table3_LiverVsScaf.xlsx",
  sheet = "Baseline C57BL6"
)

# Inspect imported column names
names(raw_df)
marker_cols <- names(raw_df)[6:15]

# Rename columns for easier handling
names(raw_df)[1:15] <- c(
  "SampleID",
  "Mouse",
  "Group",
  "Time",
  "Tissue",
  "Immune cells",
  "Neutrophils",
  "Macrophages",
  "CX3CR1+ Macrophages",
  "T cells",
  "CD4 T cells",
  "CD8 T cells",
  "B cells",
  "Dendritic cells",
  "CX3CR1+ Dendritic cells"
)


df <- raw_df %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "Baseline (W/o STZ)",
        "Baseline (With STZ)"
      )
    ),
    Tissue = factor(
      Tissue,
      levels = c(
        "Scaffold",
        "Liver"
      )
    )
  )

marker_names <- c(
  "Immune cells",
  "Neutrophils",
  "Macrophages",
  "CX3CR1+ Macrophages",
  "T cells",
  "CD4 T cells",
  "CD8 T cells",
  "B cells",
  "Dendritic cells",
  "CX3CR1+ Dendritic cells"
)

paired_long <- df %>%
  pivot_longer(
    cols = all_of(marker_names),
    names_to = "Marker",
    values_to = "Frequency"
  ) %>%
  select(
    Mouse,
    Group,
    Time,
    Tissue,
    Marker,
    Frequency
  ) %>%
  pivot_wider(
    names_from = Tissue,
    values_from = Frequency
  ) %>%
  drop_na(
    Scaffold,
    Liver
  )

paired_long %>%
  dplyr::count(Marker, Group)


marker_axis_labels <- c(
  "Immune cells" = "CD45+ immune cells (% of live cells)",
  "Neutrophils" = "Ly6G+ neutrophils (% of CD45+ cells)",
  "Macrophages" = "F4/80+ macrophages (% of CD45+ cells)",
  "CX3CR1+ Macrophages" = "CX3CR1+ macrophages (% of CD45+ cells)",
  "T cells" = "CD3+ T cells (% of CD45+ cells)",
  "CD4 T cells" = "CD4+ T cells (% of CD3+ T cells)",
  "CD8 T cells" = "CD8+ T cells (% of CD3+ T cells)",
  "B cells" = "CD19+ B cells (% of CD45+ cells)",
  "Dendritic cells" = "CD11c+MHC-II+ dendritic cells (% of CD45+ cells)",
  "CX3CR1+ Dendritic cells" = "CX3CR1+ dendritic cells (% of CD45+ cells)"
)

make_correlation_plot <- function(marker_name, data = paired_long) 
{
  marker_df <- data %>%
    dplyr::filter(
      Marker == marker_name,
      !is.na(Liver),
      !is.na(Scaffold)
    )
  
  if (nrow(marker_df) < 3) {
    stop(
      paste0(
        "Not enough paired observations for marker: ",
        marker_name
      )
    )
  }
  
  # Calculate Pearson correlation separately for each group
  cor_stats <- marker_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      cor_test = list(
        if (dplyr::n() >= 3) {
          cor.test(
            Liver,
            Scaffold,
            method = "pearson"
          )
        } else {
          NULL
        }
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      r = purrr::map_dbl(
        cor_test,
        ~ if (is.null(.x)) NA_real_ else unname(.x$estimate)
      ),
      p = purrr::map_dbl(
        cor_test,
        ~ if (is.null(.x)) NA_real_ else .x$p.value
      ),
      cor_label = dplyr::if_else(
        is.na(r),
        paste0(Group, ": insufficient observations"),
        paste0(
          Group,
          ": r = ",
          sprintf("%.2f", r),
          ", P = ",
          format.pval(
            p,
            digits = 2,
            eps = 0.001
          )
        )
      )
    )
  
  # Combine group-specific statistics into one annotation
  stats_label <- paste(
    cor_stats$cor_label,
    collapse = "\n"
  )
  
  axis_label <- marker_axis_labels[[marker_name]]
  
  if (is.null(axis_label)) {
    axis_label <- marker_name
  }
  
  ggplot(
    marker_df,
    aes(
      x = Liver,
      y = Scaffold,
      color = Group,
      shape = Group,
      fill = Group
    )
  ) +
    
    # Separate regression line for each group
    geom_smooth(
      aes(
        group = Group,
        color = Group,
        fill = Group
      ),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 1,
      alpha = 0.20
    ) +
    
    geom_point(
      size = 4.5,
      color = "black",
      stroke = 0.8
    ) +
    
    scale_color_manual(
      values = c(
        "Baseline (W/o STZ)" = "black",
        "Baseline (With STZ)" = "#7A3E00"   # dark brown
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Baseline (W/o STZ)" = "black",
        "Baseline (With STZ)" = "#7A3E00"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "Baseline (W/o STZ)" = 21,  # filled circle
        "Baseline (With STZ)" = 22   # filled square
      )
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = stats_label,
      hjust = 1.05,
      vjust = 1.2,
      size = 4.5,
      lineheight = 1.2
    ) +
    
    labs(
      title = marker_name,
      x = paste0("Liver ", axis_label),
      y = paste0("IN ", axis_label),
      color = "Group",
      fill = "Group",
      shape = "Group"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 14,
        face = "bold"
      ),
      
      axis.text = element_text(
        size = 13,
        colour = "black"
      ),
      
      axis.line = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.8
      ),
      
      axis.ticks.length = grid::unit(
        0.25,
        "cm"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      )
    )
}



markers_to_plot <- unique(paired_long$Marker)

correlation_plots <- setNames(
  lapply(
    markers_to_plot,
    make_correlation_plot
  ),
  markers_to_plot
)

correlation_plots

#9. Long Term aCD40L-Allo Flow Analysis----
getwd()
library(readr)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)

LT_AlloFlow <- read_csv(
  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Islet Rejection Flow/AdvScienceRevision_Flow/Long Term Anti-CD40L Allo Flow/LongTerm_aCD40LAllogeneic_Flow.csv",
  show_col_types = FALSE
) %>%
  select(where(~ !all(is.na(.))))


LT_AlloFlow <- LT_AlloFlow %>%
  filter(!(Time == 70 & MouseID %in% c("R96(4L)", "R97(NL)")))

LT_AlloFlow$MouseID <- factor(LT_AlloFlow$MouseID)

LT_AlloFlow$Time <- factor(
  LT_AlloFlow$Time,
  levels = c(14,28,42,70,106,175)
)
colnames(LT_AlloFlow)
ocean_blue <- "#1F77B4"
markers <- c(
  "CD45+ Immune Cells | Freq. of Live Cells (%)",
  "Ly6G+ Neutrophils | Freq. of CD45+ Immune Cells (%)",
  "F4 80+ Macrophages | Freq. of CD45+ Immune Cells (%)",
  "CCR2+ Macrophages | Freq. of CD45+ Immune Cells (%)",
  "CD80+ Macrophages | Freq. of CD45+ Immune Cells (%)",
  "CD80+ Macrophages | Freq. of F4 80+ Macrophages (%)" ,
  "CD86+ Macrophages | Freq. of CD45+ Immune Cells (%)",
  "CX3CR1+ Macrophages | Freq. of CD45+ Immune Cells (%)",
  "CX3CR1+ Macrophages | Freq. of F4 80+ Macrophages (%)",
  "Ly6C+ Monocytes | Freq. of CD45+ Immune Cells (%)",
  "CD3+ T Cells | Freq. of CD45+ Immune Cells (%)",
  "CD4+ T Cells | Freq. of CD45+ Immune Cells (%)",
  "CD43+ CD4 T Cells | Freq. of CD45+ Immune Cells (%)",
  "CD8+ T Cells | Freq. of CD45+ Immune Cells (%)",
  "CD43+ CD8 T Cells | Freq. of CD45+ Immune Cells (%)",
  "CD19+ B Cells | Freq. of CD45+ Immune Cells (%)",
  "CD43+ B Cells | Freq. of CD45+ Immune Cells (%)",
  "CD11c+MHC-II+ Dendritic Cells | Freq. of CD45+ Immune Cells (%)",
  "CCR2+ Dendritic Cells | Freq. of CD45+ Immune Cells (%)",
  "CD80+ Dendritic Cells | Freq. of CD45+ Immune Cells (%)",
  "CD80+ Dendritic Cells | Freq. of Dendritic Cells (%)" ,
  "CD86+ Dendritic Cells | Freq. of CD45+ Immune Cells (%)",
  "CX3CR1+ Dendritic Cells | Freq. of CD45+ Immune Cells (%)",
  "CX3CR1+ Dendritic Cells | Freq. of Dendritic Cells (%)"
)

plot_longitudinal <- function(marker){
  plot_title <- trimws(sub("\\|.*", "", marker))
  y_label <- trimws(sub(".*\\|", "", marker))
  df <- LT_AlloFlow %>%
    select(MouseID, Time, all_of(marker)) %>%
    rename(Value = all_of(marker)) %>%
    drop_na(Value)
  
  model <- lmer(
    Value ~ Time + (1 | MouseID),
    data = df
  )
  
  p_time <- anova(model)["Time", "Pr(>F)"]
  
  summary_df <- df %>%
    group_by(Time) %>%
    summarise(
      Mean = mean(Value),
      SEM = sd(Value) / sqrt(n()),
      .groups = "drop"
    )
  
  g <- ggplot(
    df,
    aes(
      x = Time,
      y = Value,
      group = MouseID
    )
  ) +
    
    # Individual mouse trajectories
    geom_line(
      colour = "#8EC7E8",
      linewidth = 0.8,
      alpha = 0.8
    ) +
    
    geom_point(
      colour = "#5FA9D3",
      fill = "#B9DDF1",
      shape = 21,
      stroke = 0.8,
      size = 3.2,
      alpha = 0.9
    ) +
    
    # Mean ± SEM
    geom_errorbar(
      data = summary_df,
      aes(
        x = Time,
        ymin = Mean - SEM,
        ymax = Mean + SEM
      ),
      width = 0.15,
      colour = ocean_blue,
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    
    geom_line(
      data = summary_df,
      aes(
        x = Time,
        y = Mean,
        group = 1
      ),
      colour = ocean_blue,
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    
    geom_point(
      data = summary_df,
      aes(
        x = Time,
        y = Mean
      ),
      shape = 21,
      fill = ocean_blue,
      colour = ocean_blue,
      size = 4.5,
      inherit.aes = FALSE
    ) +
    
    annotate(
      "text",
      x = 5.7,
      y = max(df$Value, na.rm = TRUE) * 1.08,
      label = paste0(
        "Mixed model\nP = ",
        signif(p_time, 3)
      ),
      hjust = 1,
      size = 6.5,
      fontface = "bold"
    ) +
    
    labs(
      x = "Days post-transplant",
      y = y_label,
      title = plot_title
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.18))
    ) +
    
    theme_classic(base_size = 18) +
    
    theme(
      axis.title.x = element_text(
        face = "bold",
        size = 22,
        margin = ggplot2::margin(t = 12)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 22,
        margin = ggplot2::margin(r = 12)
      ),
      axis.text.x = element_text(
        colour = "black",
        size = 18
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 18
      ),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      plot.title = element_text(
        face = "bold",
        size = 24,
        hjust = 0.5,
        margin = ggplot2::margin(b = 12)
      )
    )
  
  return(g)
}
dir.create("Longitudinal_Flow_Plots",
           showWarnings = FALSE)

for(i in markers){
  
  p <- plot_longitudinal(i)
  
  ggsave(
    filename=paste0(
      "Longitudinal_Flow_Plots/",
      make.names(i),
      ".pdf"
    ),
    plot=p,
    width=6.5,
    height=5.65,
    dpi=600
  )
}


#10. Baseline Transcriptomics----
# Allo Transplant Metadata Importing
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch3/Metadata_Batch3.csv", sep=",", header=T) # Metadata file
meta_batch4 <- read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch4/Metadata_Batch4.csv", sep=",", header=T) # Metadata file
meta_baseline <- read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Baseline IN/Metadata_Baseline.csv", sep=",", header=T) # Metadata file

meta_batch3 <- as.data.frame(meta_batch3)
meta_batch4 <- as.data.frame(meta_batch4)
meta_baseline <- as.data.frame(meta_baseline)
# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
meta_combined <- rbind(meta_batch3,meta_batch4,meta_baseline)

# Preview the combined metadata
head(meta_combined)
unique(meta_combined$Group)

#All Transplant Counts Data Importing

counts_batch3 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch3/IsTx_gene_expected_count_annot_batch3.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch3 <- na.omit(counts_batch3)

counts_batch4 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Batch4/IsTx_gene_expected_count_annot_batch4.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch4 <- na.omit(counts_batch4)

counts_baseline <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Baseline IN/BaselineIN_counts.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_baseline <- na.omit(counts_baseline)



counts_batch3 <- counts_batch3[!duplicated(counts_batch3[, 1]), ]
genes <- counts_batch3[, 1]
rownames(counts_batch3) <- genes
counts_batch3 <- counts_batch3[, -1]

counts_batch4 <- counts_batch4[!duplicated(counts_batch4[, 1]), ]
genes <- counts_batch4[, 1]
rownames(counts_batch4) <- genes
counts_batch4 <- counts_batch4[, -1]

counts_baseline <- counts_baseline[!duplicated(counts_baseline[, 1]), ]
genes <- counts_baseline[, 1]
rownames(counts_baseline) <- genes
counts_baseline <- counts_baseline[, -1]
#Combine counts data
# First merge counts_batch1 and counts_batch2
combined_counts <- merge(counts_batch3, counts_batch4, by = "row.names", all = TRUE)
# Rename the Row.names column back
rownames(combined_counts) <- combined_counts$Row.names
combined_counts$Row.names <- NULL
# Then merge the result with counts_batch3
combined_counts <- merge(combined_counts, counts_baseline, by = "row.names", all = TRUE)
# Rename the Row.names column back
rownames(combined_counts) <- combined_counts$Row.names
combined_counts$Row.names <- NULL


# Identify the samples to keep (not "Technical Rejection")
samples_to_keep <- meta_combined$Samples[meta_combined$Group != "Technical Rejection"]

# Filter the meta_combined data
meta_combined <- meta_combined[meta_combined$Group != "Technical Rejection", ]

# Filter the combined_counts data to keep only columns corresponding to samples_to_keep
combined_counts <- combined_counts[, colnames(combined_counts) %in% samples_to_keep]


# Preview the combined dataset
head(combined_counts)

#.Preprocessing and Cleaning

# Remove zero and low expressed genes
combined_counts <- combined_counts[, meta_combined$Samples]  # Ensure Sample_IDs match column names in combined_counts

IsletTransplantCounts_Baseline <- flexiDEG.function1(combined_counts, meta_combined, # Run Function 1
                                            convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                            batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0, 0

#Remove undefined and pseudogenes
remove_pattern <- "^Gm[0-9]|^AC[0-9]|^AL[0-9]|^AI[0-9]|^AW[0-9]|^AF[0-9]|^BB[0-9]|^BC[0-9]|^CT[0-9]|^CAAA|^BX[0-9]|^CN[0-9]|^CR[0-9]|^C[0-9]{4,}|^Olfr"
rows_to_remove <- grep(remove_pattern, rownames(IsletTransplantCounts_Baseline)) #Remove Gm genes
IsletTransplantCounts_Baseline <- IsletTransplantCounts_Baseline[-rows_to_remove, ]
# connect to Ensembl mouse database

options(timeout = 120)



pseudo_genes <- gene_info$mgi_symbol[grep("pseudogene", gene_info$gene_biotype)]
# remove them
IsletTransplantCounts_Baseline <- IsletTransplantCounts_Baseline[
  !(rownames(IsletTransplantCounts_Baseline) %in% pseudo_genes), ]


# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_combined)
unique_colSide <- unique(colSide)

# 3.Create DESqEQ Object 

IsletTransplantCounts_Baseline <- as.matrix(IsletTransplantCounts_Baseline)
storage.mode(IsletTransplantCounts_Baseline) <- "integer"

dds_IsletTransplantBaseline <- DESeqDataSetFromMatrix(IsletTransplantCounts_Baseline, meta_combined,
                                              design = ~ 1)   # dummy design for now

sel <- colData(dds_IsletTransplantBaseline)$Group %in% c("Control Allogeneic","Control Syngeneic","Baseline C57BL/6", "Baseline  C57BL/6+STZ" ) & colData(dds_IsletTransplantBaseline)$Day %in% c(14)
dds_IsletTransplant_BaselineD14<- dds_IsletTransplantBaseline[, sel]
dds_IsletTransplant_BaselineD14$Group <- droplevels(
  factor(dds_IsletTransplant_BaselineD14$Group)
)
# Confirm sample numbers
table(dds_IsletTransplant_BaselineD14$Group)


AlloSyn_Signature <- c("Myo18b","Cd59b","Rep15","Shisa9","Rragb","Gdf3",
                       "Lncpint","Ifitm5","Malat1","Myo3b","Hipk4","S1pr5",
                       "Trpm6","Tctn2","Lhfpl4","Ido2","Fsip1","Emx2os",
                       "Hlf","Gfra1","Bmp5","Col6a5","Myh3","Cbs","Vtn","Dnah8","Sybu")

vsd_D14 <- vst(
  dds_IsletTransplant_BaselineD14,
  blind = TRUE
)

vst_mat <- assay(vsd_D14)

# Retain signature genes present in the dataset
signature_present <- intersect(AlloSyn_Signature, rownames(vst_mat))
signature_missing <- setdiff(AlloSyn_Signature, rownames(vst_mat))

message("Signature genes retained: ", length(signature_present))
message(
  "Missing signature genes: ",
  ifelse(
    length(signature_missing) == 0,
    "None",
    paste(signature_missing, collapse = ", ")
  )
)

mat_signature <- vst_mat[signature_present, , drop = FALSE]

# Remove genes with zero variance across samples
gene_variance <- apply(mat_signature, 1, var, na.rm = TRUE)
mat_signature <- mat_signature[
  is.finite(gene_variance) & gene_variance > 0,
  ,
  drop = FALSE
]

# Sample metadata in the same order as the expression matrix
cd <- as.data.frame(
  colData(dds_IsletTransplant_BaselineD14)[
    colnames(mat_signature),
    ,
    drop = FALSE
  ]
)
groups_D14 <- c(
  "Control Allogeneic",
  "Control Syngeneic",
  "Baseline C57BL/6",
  "Baseline  C57BL/6+STZ"
)
cd$Group <- factor(
  cd$Group,
  levels = groups_D14
)

stopifnot(identical(rownames(cd), colnames(mat_signature)))


# 4. Heatmap


# Row-wise Z-score for visualization
mat_scaled <- t(scale(t(mat_signature)))

# Remove any rows that could not be scaled
mat_scaled <- mat_scaled[
  apply(mat_scaled, 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]

annotation_col <- data.frame(
  Group = cd$Group,
  row.names = rownames(cd)
)

stopifnot(
  identical(rownames(annotation_col), colnames(mat_scaled))
)

group_colors <- c(
  "Control Allogeneic" = "#E60000",
  "Control Syngeneic" = "#2E6F40",
  "Baseline C57BL/6" = "#000000",
  "Baseline  C57BL/6+STZ" = "#B8860B"
)



# 5. PCA using signature genes

# Samples in rows and genes in columns
mat_pca <- t(mat_signature)

# Remove any genes with zero variance
pca_gene_sd <- apply(mat_pca, 2, sd, na.rm = TRUE)
mat_pca <- mat_pca[
  ,
  is.finite(pca_gene_sd) & pca_gene_sd > 0,
  drop = FALSE
]

pca <- prcomp(
  mat_pca,
  center = TRUE,
  scale. = TRUE
)

percent_variance <- 100 * (
  pca$sdev^2 / sum(pca$sdev^2)
)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = cd[rownames(pca$x), "Group"]
)

# Requested symbols:
# Baseline C57BL/6 = circle
# Baseline C57BL/6+STZ = square
group_shapes <- c(
  "Control Allogeneic" = 25,       # downward triangle
  "Control Syngeneic" = 24,        # upward triangle
  "Baseline C57BL/6" = 21,         # circle
  "Baseline  C57BL/6+STZ" = 22     # square
)

pca_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    fill = Group,
    shape = Group
  )
) +
  geom_point(
    size = 5,
    stroke = 1.2
  ) +
  stat_ellipse(
    geom = "polygon",
    aes(group = Group),
    alpha = 0.15,
    level = 0.70,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_shape_manual(values = group_shapes, drop = FALSE) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA of Allogeneic–Syngeneic Signature",
    x = paste0(
      "PC1 (",
      round(percent_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(percent_variance[2], 1),
      "%)"
    )
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title = element_text(
      size = 20,
      face = "bold"
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(
      color = "black",
      linewidth = 1
    ),
    panel.grid = element_blank(),
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    )
  )

pca_plot

