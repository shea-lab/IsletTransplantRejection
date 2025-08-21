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
kegg_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:KEGG_LEGACY") # Large df w/ categories
pwl_kegg <- split(kegg_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                  kegg_gene_sets$gs_name) # Pathway names
biocarta_gene_sets <- msigdbr(species="Mus musculus", subcategory="CP:BIOCARTA") # Large df w/ categories
pwl_biocarta <- split(biocarta_gene_sets$gene_symbol, # Genes to split into pathways, by ensembl
                      biocarta_gene_sets$gs_name) # Pathway names
pwl_msigdbr <- c(pwl_hallmark, pwl_reactome, pwl_kegg, pwl_biocarta) # Compile them all
length(pwl_msigdbr)



# 1. Load the Data ----
# Organize Data
setwd("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/")
getwd()

#Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch1/Metadata_Batch1.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch2/Metadata_Batch2.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch3/Metadata_Batch3.csv", sep=",", header=T) # Metadata file
meta_batch4 <- read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch4/Metadata_Batch4.csv", sep=",", header=T) # Metadata file

meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)
meta_batch4 <- as.data.frame(meta_batch4)
# Merge metadata by columns (i.e., add samples from Batch 2 to Batch 1)
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3,meta_batch4)

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

counts_batch4 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Batch4/IsTx_gene_expected_count_annot_batch4.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
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
# Then merge the result with counts_batch4
combined_counts <- merge(combined_counts, counts_batch4, by = "row.names", all = TRUE)
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
# Create new column with high-level analysis groups
meta_combined$Outcomes <- meta_combined$Group
# Change Late Rejection to Tolerance for Initial Analysis 
meta_combined$Group[meta_combined$Group == "Late Rejection"] <- "Tolerance"

# 2.Preprocessing and Cleaning ----

# Select columns in combined_counts that match the remaining sample names in meta_combined
combined_counts <- combined_counts[, meta_combined$Samples]  # Ensure Sample_IDs match column names in combined_counts

IsletTransplantCounts <- flexiDEG.function1(combined_counts, meta_combined, # Run Function 1
                                            convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                            batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0, 0
rows_to_remove <- grep("^Gm[0-9]", rownames(IsletTransplantCounts))
# Remove those rows from case1_f1
IsletTransplantCounts <- IsletTransplantCounts[-rows_to_remove, ]

# Saving case1_f1 dataframe as a CSV file
write.csv(IsletTransplantCounts, file = "/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Combined/Combined_Raw_Counts_ITx_Filtered.csv", row.names = TRUE)

# Saving meta_combined dataframe as a CSV file
write.csv(meta_combined, file = "/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Combined/Combined_Metadata_ITx.csv", row.names = FALSE)

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

saveRDS(dds_IsletTransplant, file = "/Users/jyotirmoyroy/Desktop/Current Projects/Transplant Rejection Sensor Paper/Data/Robjects/dds_IsletTransplant_master.rds")

# 4. Allogeneic Vs Syngeneic Signature----
# subset samples
sel <- colData(dds_IsletTransplant)$Group %in% c("Control Rejected","Control Accepted") & colData(dds_IsletTransplant)$Day %in% c(7,14)
dds_IsletTransplant_AlloVsSyn <- dds_IsletTransplant[, sel]

dds_IsletTransplant_AlloVsSyn$Batch       <- factor(dds_IsletTransplant_AlloVsSyn$Batch)
dds_IsletTransplant_AlloVsSyn$LibraryPrep <- factor(dds_IsletTransplant_AlloVsSyn$LibraryPrep)
dds_IsletTransplant_AlloVsSyn$Group <- factor(
  dds_IsletTransplant_AlloVsSyn$Group,
  levels = c("Control Accepted", "Control Rejected")  # order sets baseline
)
levels(dds_IsletTransplant_AlloVsSyn$Group)

dds_IsletTransplant_AlloVsSyn$Day <- factor(dds_IsletTransplant_AlloVsSyn$Day,
                                            levels = c(7, 14))  # baseline = 7
levels(dds_IsletTransplant_AlloVsSyn$Day)

# Check for Collinearuty
cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
# Basic sanity
lapply(cd[, c("Batch","LibraryPrep","Day","Group")], function(x) table(x, useNA="ifany"))
# Check for NAs
sapply(cd[, c("Batch","LibraryPrep","Day","Group")], function(x) any(is.na(x)))
# Model matrix rank
mm <- model.matrix(~ Batch + LibraryPrep + Day + Group, data = cd)
qr(mm)$rank; ncol(mm)             # if rank < ncol(mm), not full rank


#Design Formula
design(dds_IsletTransplant_AlloVsSyn) <- ~ Batch + Day + Group + Day:Group
dds_IsletTransplant_AlloVsSyn <- DESeq(dds_IsletTransplant_AlloVsSyn)
design(dds_IsletTransplant_AlloVsSyn)
resultsNames(dds_IsletTransplant_AlloVsSyn)

# Day 7
res_ALLOvSYN_D7 <- results(dds_IsletTransplant_AlloVsSyn,
                           name = "Group_Control.Rejected_vs_Control.Accepted")

# Day 14
res_ALLOvSYN_D14 <- results(dds_IsletTransplant_AlloVsSyn,
                            list(c("Group_Control.Rejected_vs_Control.Accepted",
                                   "Day14.GroupControl.Rejected")))

# Interactions
res_interaction <- results(dds_IsletTransplant_AlloVsSyn,
                           name = "Day14.GroupControl.Rejected")
res_interaction_df <- as.data.frame(res_interaction)
res_interaction_df$gene <- rownames(res_interaction_df)

summary(res_ALLOvSYN_D7)
summary(res_ALLOvSYN_D14)

library(EnhancedVolcano)

# Thresholds
q_cut  <- 0.10
fc_cut <- 1

library(EnhancedVolcano)

## 0) Prep results as data.frames with a 'gene' column
d7  <- as.data.frame(res_ALLOvSYN_D7);  d7$gene  <- rownames(res_ALLOvSYN_D7)
d14 <- as.data.frame(res_ALLOvSYN_D14); d14$gene <- rownames(res_ALLOvSYN_D14)

write.csv(d7, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/DESEQResults_Day7_Allo_vs_Syn.csv", row.names = TRUE)
write.csv(d14, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/DESEQResults_Day14_Allo_vs_Syn.csv", row.names = TRUE)
write.csv(res_interaction_df, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/DESEQResults_InteractionTimeGroup_Allo_vs_Syn.csv", row.names = TRUE)

# Function to create color mapping
make_keyvals_fdr_fc_AllovsSyn <- function(df, q = 0.10, fc = 1,
                                          col_up = "#D62728", col_down = "#1F77B4", col_ns = "gray70") {
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
  
  label[up]   <- paste0("Allogeneic-Upregulated")
  label[down] <- paste0("Syngeneic-Upregulated")
  
  names(key) <- label        # <- legend labels; no NAs
  key
}
keyvals_d7  <- make_keyvals_fdr_fc_AllovsSyn(d7)
keyvals_d14 <- make_keyvals_fdr_fc_AllovsSyn(d14)

library(dplyr)
pick_labels <- function(df, q = 0.10, fc = 1, topN = 30) {
  idx <- which(!is.na(df$padj) & !is.na(df$log2FoldChange) &
                 df$padj < q & abs(df$log2FoldChange) >= fc)
  if (length(idx) == 0) return(character(0))
  ord <- order(df$padj[idx], -abs(df$log2FoldChange[idx]), na.last = NA)  # tie-break by |LFC|
  labs <- df$gene[idx][ord]
  labs <- make.unique(labs)  # avoid dup labels
  labs[seq_len(min(topN, length(labs)))]
}

selLab_d7  <- pick_labels(d7,  q_cut, fc_cut, 30)
selLab_d14 <- pick_labels(d14, q_cut, fc_cut, 30)


## 3) Axis limits
xmax_d7  <- max(2, ceiling(max(abs(d7$log2FoldChange),  na.rm=TRUE)))
xmax_d14 <- max(2, ceiling(max(abs(d14$log2FoldChange), na.rm=TRUE)))

## Volcano: Day 7----
EnhancedVolcano(
  d7,
  lab           = d7$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Allogeneic vs Syngeneic — Day 7",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d7$padj<0.10 & abs(d7$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d7, xmax_d7),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_d7,
  legendPosition= "right",
  selectLab     = selLab_d7
)

## Volcano: Day 14----
EnhancedVolcano(
  d14,
  lab           = d14$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Allogeneic vs Syngeneic — Day 14",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d14$padj<0.10 & abs(d14$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d14, xmax_d14),
  ylim = c(0,10),
  boxedLabels   = TRUE,
  pointSize     =3,
  labSize       = 6,
  colAlpha      = 0.8,
  max.overlaps = 15,
  drawConnectors= TRUE,
  colCustom     = keyvals_d14,
  legendPosition= "right",
  selectLab     = selLab_d14
)

## Interactions- Day x Group----
EnhancedVolcano(res_interaction_df,
                lab = rownames(res_interaction_df),
                pCutoff       = 0.10,
                FCcutoff      = 1,
                ylim = c(0,8),
                boxedLabels   = TRUE,
                pointSize     =3,
                labSize       = 6,
                colAlpha      = 0.8,
                max.overlaps = 15,
                drawConnectors= TRUE,
                legendPosition= "right",
                x = "log2FoldChange",
                y = "padj",
                xlab          = expression("log"[2]*"(Fold Change)"),
                ylab          = expression("-log"[10]*"(FDR)"),
                title = "Interaction: (Day 14 vs 7) × (Allo vs Syn)")






## PLSDA Analysis ----

library(DESeq2)
library(limma)
library(mixOmics)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# --- Build the 4-group label from your fitted object ---
CountsAlloVsSyn <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))

# Short group names
shortGroup <- ifelse(CountsAlloVsSyn$Group == "Control Rejected", "Allogeneic", "Syngeneic")
CountsAlloVsSyn$GD4 <- factor(paste0(shortGroup, "_", CountsAlloVsSyn$Day),
                 levels = c("Allogeneic_7","Allogeneic_14","Syngeneic_7","Syngeneic_14"))
# --- Variance-stabilized expression ---
vsd <- vst(dds_IsletTransplant_AlloVsSyn, blind = FALSE)
mat <- assay(vsd)  # genes x samples
# --- Remove batch (for visualization only) ---
mat_corr <- removeBatchEffect(mat, batch = colData(vsd)$Batch)

# --- PLS-DA input: samples x genes ---
X <- t(mat_corr)   # samples in rows
Y <- CountsAlloVsSyn$GD4        # factor with 4 levels

set.seed(123)
plsda_model <- mixOmics::plsda(X, Y, ncomp = 2)

# Scores for LV1 & LV2
scores <- plsda_model$variates$X
plot_df <- data.frame(
  LV1   = scores[, 1],
  LV2   = scores[, 2],
  Group = Y
)


group_colors <- c(
  "Allogeneic_14" = "#8B1A1A",  # darker shade of red (deep crimson)
  "Allogeneic_7"  = "#D62728",  # base red
  "Syngeneic_14"  = "#08306B",  # darker shade of blue (navy/steel blue)
  "Syngeneic_7"   = "#1F77B4"   # base blue
)
group_shapes <- c(
  "Allogeneic_7"  = 16,
  "Allogeneic_14" = 17,
  "Syngeneic_7"   = 15,
  "Syngeneic_14"  = 18
)


expl_var <- round(plsda_model$prop_expl_var$X * 100, 1)  # % explained variance for X
xlab <- paste0("PLS Component 1 (", expl_var[1], "%)")
ylab <- paste0("PLS Component 2 (", expl_var[2], "%)")

# 2D PLS-DA plot (LV1 vs LV2) with filled ellipses
ggplot(plot_df, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  # points
  geom_point(size = 5, alpha = 0.9) +
  # filled confidence ellipses (70% or 68% CI both common)
  stat_ellipse(
    geom  = "polygon", 
    aes(fill = Group), 
    level = 0.70, 
    alpha = 0.3, 
    show.legend = FALSE
  ) +
  # custom colors, shapes
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values  = group_colors) +
  scale_shape_manual(values = group_shapes) +
  # labels
  labs(
    title = "PLS-DA: Allogeneic vs Syngeneic Across Days",
    x = xlab,
    y = ylab
  ) +
  # minimal but with axis lines
  theme_minimal(base_size = 18) +
  theme(
    legend.position   = "top",
    legend.title      = element_blank(),
    axis.title        = element_text(size = 20, face = "bold"),
    axis.text         = element_text(size = 16, color = "black"),
    axis.line         = element_line(color = "black", linewidth = 0.8),
    axis.ticks        = element_line(color = "black"),
    panel.grid        = element_blank(),
    plot.title        = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.margin       = unit(c(10, 40, 10, 10), "pt")  # top, right, bottom, left
  )

# Get VIP scores
vip_scores <- vip(plsda_model)
# Extract only Component 1
vip_axis1 <- vip_scores[, 1]
vip_df_axis1 <- data.frame(
  Gene = rownames(vip_scores),
  VIP_Axis1 = vip_axis1
)
#Sort by descending VIP
vip_df_axis1 <- vip_df_axis1[order(-vip_df_axis1$VIP_Axis1), ]
#Save to CSV
write.csv(vip_df_axis1, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/VIP_scores_AlloVsSyn_PLS1.csv", row.names = FALSE)
# Rank genes by VIP (descending)
top100_genes <- names(sort(vip_scores[, 1], decreasing = TRUE))[1:100]
mat_top100 <- mat[top100_genes, ]
anno <- data.frame(
  Group = Y,        # group labels (factor)
  row.names = colnames(mat_top100)
)

anno_colors <- list(
  Group = group_colors   # same color scheme you used before
)

Y <- factor(Y, levels = c("Allogeneic_7","Syngeneic_7", "Allogeneic_14", "Syngeneic_14"))

# Reorder the columns of the matrix based on group
sample_order <- order(Y)
mat_top100_ordered <- mat_top100[, sample_order]
anno_ordered <- anno[sample_order, , drop = FALSE]
pheatmap(
  mat_top100_ordered,
  scale = "row",                  # row-wise z-score
  annotation_col = anno_ordered,          # add group annotation
  annotation_colors = anno_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50) # publication-style
)

## Immune Markers ----

library(DESeq2)
library(pheatmap)
library(limma)


# Define immune gene list (same as before)
immune_genes <- c("Cd3d","Cd3e","Cd3g","Cd247","Trac","Trbc1","Trbc2","Cd4","Cd8a","Cd8b1",
                  "Cd28","Ctla4","Icos","Pdcd1","Lag3","Havcr2","Tnfrsf9","Tnfrsf18","Tnfrsf4","Cd40lg",
                  "Gzma","Gzmb","Gzmk","Prf1","Ifng","Tnf","Il2",
                  "Foxp3","Il2ra","Ikzf2","Tnfrsf18","Ly6G",
                  "Tbx21","Stat4","Ifng","Gata3","Il4","Il5","Il13","Rorc","Il17a","Il17f","Il21","Il22","Bcl6","Cxcr5","Tox",
                  "Itgax","Zbtb46","Cd80","Cd86","Cd40","H2-Ab1","H2-Aa","H2-Eb1",
                  "Adgre1","Cd68","Cd14","Itgam","Csf1r","Mrc1","Nos2","Arg1",
                  "Cd19","Cd79a","Cd79b","Ms4a1","Cd22",
                  "Ccl2","Ccl5","Cxcl9","Cxcl10","Cxcl11","Cxcr3","Ccr7","Spn")

# immune_genes <- c(
#   "Cd163","Csf2","C4b","C8a","Cfi",
#   "Ifi202b","Il24","Ly6c1","Masp2","Mcpt4",
#   "Plch2","S1pr5","Tac1","Zbtb16"
# )
vsd <- vst(dds_IsletTransplant_AlloVsSyn, blind = FALSE)
mat  <- assay(vsd)
matc <- removeBatchEffect(mat, batch = colData(vsd)$Batch)

# keep only immune genes from earlier step
mat_filt   <- matc[rownames(matc) %in% immune_genes, ]
mat_scaled <- t(scale(t(mat_filt)))

heat_colors <- colorRampPalette(c("navy","white","firebrick3"))(100)

day_vec <- as.character(colData(vsd)$Day)
group_vec <- as.character(colData(vsd)$Group)

ann_colors <- list(
  Group = c("Control Rejected" = "#D62728", "Control Accepted" = "#1F77B4"),
  Day   = c("7" = "gray40", "14" = "gray70")
)

for (day_sel in c("7","14")) {
  cols_idx <- which(day_vec == day_sel)
  
  # reorder columns: first Control Accepted, then Control Rejected
  order_idx <- cols_idx[order(group_vec[cols_idx])]
  mat_day   <- mat_scaled[, order_idx, drop = FALSE]
  
  annotation_col <- data.frame(
    Group = group_vec[order_idx],
    Day   = day_vec[order_idx]
  )
  rownames(annotation_col) <- colnames(mat_day)
  
  print(
    pheatmap(mat_day,
             color = heat_colors,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             cluster_rows = TRUE,
             cluster_cols = FALSE,   # no column clustering
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 8,
             main = paste("Immune Gene Expression (Day", day_sel, ")"))
  )
}


## GSEA Analysis----
# Paths to your saved results
d7_path  <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/DESEQResults_Day7_Allo_vs_Syn.csv"
d14_path <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/DESEQResults_Day14_Allo_vs_Syn.csv"

# Spot-check a known rejection gene (e.g., Ifng) if present:
d7["Ifng", c("log2FoldChange","stat","padj")]
d14["Ifng", c("log2FoldChange","stat","padj")]

# Import
d7  <- read.csv(d7_path,  row.names = 1)
d14 <- read.csv(d14_path, row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_d7  <- d7$stat;  names(lfc_vector_d7)  <- rownames(d7)
lfc_vector_d14 <- d14$stat; names(lfc_vector_d14) <- rownames(d14)

# Drop NAs
lfc_vector_d7  <- lfc_vector_d7[!is.na(lfc_vector_d7)]
lfc_vector_d14 <- lfc_vector_d14[!is.na(lfc_vector_d14)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_d7  <- sort(lfc_vector_d7,  decreasing = TRUE)
lfc_vector_d14 <- sort(lfc_vector_d14, decreasing = TRUE)

library(msigdbr)
library(dplyr)

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

# --- Final combined TERM2GENE data frame ---
mm_all_df <- rbind(mm_c8_df, mm_h_df, mm_kegg_df)

library(clusterProfiler)
library(msigdbr)
# Day 7
gsea_results_d7 <- GSEA(
  geneList      = lfc_vector_d7,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  #keyType       = "SYMBOL",       # <- tell it explicitly
  TERM2GENE     = mm_all_df
)
gsea_results_d7_df <- as.data.frame(gsea_results_d7)

# Day 14
gsea_results_d14 <- GSEA(
  geneList      = lfc_vector_d14,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_d14_df <- as.data.frame(gsea_results_d14)


# Full results Day 7
write.csv(gsea_results_d7_df,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/GSEAResults_AlloVsSyn_Day7.csv",
          row.names = FALSE)

# Full results Day 14
write.csv(gsea_results_d14_df,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/GSEAResults_AlloVsSyn_Day14.csv",
          row.names = FALSE)


library(dplyr)
library(stringr)
library(ggplot2)

# Your curated list EXACTLY as provided (de-dup just in case)
immune_master <- unique(c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V2",
  "KEGG_TRYPTOPHAN_METABOLISM",
  "KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM",
  "KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "KEGG_GRAFT_VERSUS_HOST_DISEASE",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_APOPTOSIS",
  "CUI_DEVELOPING_HEART_C8_MACROPHAGE",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "JONES_OVARY_T_CELL",
  "JONES_OVARY_MACROPHAGE",
  "DESCARTES_FETAL_INTESTINE_LYMPHOID_CELLS",
  "DESCARTES_FETAL_LIVER_LYMPHOID_CELLS",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY"
))

# Helper to standardize clusterProfiler GSEA columns
coerce_gsea_tbl <- function(df, day_label){
  # try common column names: Description/ID/pathway/setName; p.adjust/padj
  gs  <- if ("Description" %in% names(df)) df$Description else if ("ID" %in% names(df)) df$ID else if ("pathway" %in% names(df)) df$pathway else if ("setName" %in% names(df)) df$setName else rownames(df)
  pad <- if ("p.adjust"   %in% names(df)) df$p.adjust   else if ("padj" %in% names(df)) df$padj else df$pval
  tibble(
    gs_name = as.character(gs),
    NES     = as.numeric(df$NES),
    padj    = as.numeric(pad),
    Day     = day_label
  )
}

# Build tidy tables for each day
d7_tbl  <- coerce_gsea_tbl(gsea_results_d7_df,  "Day 7")
d14_tbl <- coerce_gsea_tbl(gsea_results_d14_df, "Day 14")

# Keep ONLY immune_master pathways
d7_sel  <- d7_tbl  %>% filter(gs_name %in% immune_master)
d14_sel <- d14_tbl %>% filter(gs_name %in% immune_master)

# Combine and ensure both days present for every pathway
plot_df <- bind_rows(d7_sel, d14_sel) %>%
  mutate(Day = factor(Day, levels = c("Day 7","Day 14"))) %>%
  # create full grid of (gs_name x Day) and fill missing with neutral values
  complete(gs_name = immune_master, Day,
           fill = list(NES = 0, padj = 1)) %>%
  mutate(logp = -log10(padj))

# Order rows nicely (keep your immune_master order or order by category if you prefer)
plot_df$gs_name <- factor(plot_df$gs_name, levels = rev(immune_master))

# Plot: rows = pathways, columns = Day; size = −log10(padj); color = NES
p <- ggplot(plot_df, aes(x = Day, y = gs_name)) +
  geom_point(aes(size = logp, color = NES)) +
  scale_size_continuous(name = "−log10(padj)", range = c(2, 8)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  labs(x = "", y = "", title = "Enriched Pathways") +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  )
print(p)



## GSVA Analysis----
library(GSVA)
library(msigdbr)
library(limma)
library(BiocParallel)

vsd  <- vst(dds_IsletTransplant_AlloVsSyn, blind = TRUE)
expr <- assay(vsd)
meta <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))

# Make sure factors are clean
meta$Group <- factor(meta$Group, levels = c("Control Accepted", "Control Rejected"))
meta$Batch <- factor(meta$Batch)


c8   <- msigdbr(species="Mus musculus", category="C8")
hall <- msigdbr(species="Mus musculus", category="H")
kegg <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_LEGACY")

sets <- c(
  split(c8$gene_symbol,   c8$gs_name),
  split(hall$gene_symbol, hall$gs_name),
  split(kegg$gene_symbol, kegg$gs_name)
)


pick_day <- function(day) which(meta$Day %in% c(day, paste0("Day",day),"D",day))

expr_d7  <- expr[, pick_day(7)]
expr_d14 <- expr[, pick_day(14)]
meta_d7  <- meta[pick_day(7),]
meta_d14 <- meta[pick_day(14),]


param <- MulticoreParam(workers = max(1, parallel::detectCores() - 1))
gsva_par_d7 <- gsvaParam(
  expr_d7,           # VST matrix
  sets,            # your gene sets
  kcdf        = "Gaussian",        # since it's log-like VST data
  minSize     = 5,
  maxSize     = 500
)
es_d7 <- gsva(gsva_par_d7)

gsva_par_d14 <- gsvaParam(
  expr_d14,           # VST matrix
  sets,            # your gene sets
  kcdf        = "Gaussian",        # since it's log-like VST data
  minSize     = 5,
  maxSize     = 500
)
es_d14 <- gsva(gsva_par_d14)


run_limma <- function(es, meta) {
  design <- model.matrix(~ 0 + Group + Batch, data = meta)
  colnames(design) <- make.names(colnames(design))
  
  fit  <- lmFit(es, design)
  contr <- makeContrasts(Rejected_vs_Accepted = GroupControl.Rejected - GroupControl.Accepted,
                         levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contr))
  topTable(fit2, number=Inf, sort.by="P")
}

GSVA_res_d7  <- run_limma(es_d7,  meta_d7)
GSVA_res_d14 <- run_limma(es_d14, meta_d14)

# # Save results
write.csv(es_d7,  "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/GSVAScores_AlloVsSyn_Day7.csv")
write.csv(es_d14, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/GSVAScores_AlloVsSyn_Day14.csv")
write.csv(GSVA_res_d7,  "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/GSVAStats_AlloVsSyn_Day7.csv")
write.csv(GSVA_res_d14, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/GSVAStats_AlloVsSyn_Day14.csv")


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
