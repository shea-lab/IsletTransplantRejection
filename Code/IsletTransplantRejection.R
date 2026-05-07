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
samples_to_keep <- meta_combined$Sample[meta_combined$Group != "Technical Rejection"]

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
annotation_col <- data.frame(
  Day = cd$Day,
  Group = cd$Group
)
rownames(annotation_col) <- colnames(mat_scaled)

rownames(mat_scaled)
annotation_colors <- list(
  Group = c(
    "Control Allogeneic" = "#FF2400",   # Scarlet
    "Control Syngeneic" = "#2E6F40"     # green
  ),
  Day = c(
    "7" = "#F28E6B",
    "14" = "#6FA287"
  )
)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 14,
  fontsize_row = 12,
  fontsize_col = 10,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  main = "Allogeneic Vs Syngeneic Signature"
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
keep <- rowSums(counts(dds_IsletTransplant_RejVsAccep) >= 10) >= 5  # Samples in Rejection Group
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


annotation_colors <- list(
  Group = c(
    "Rejection" = "#F28500" ,   # Tangerine
    "Acceptance" =  "#064273"   # Ocean
  ),
  Day = c(
    "7" = "#F28E6B",
    "14" = "#6FA287"
  )
)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 14,
  fontsize_row = 12,
  fontsize_col = 10,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  main = "Rejection Vs Acceptance Signature"
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

summary(ALL_ACCEPVSSYN)

## Heatmap and PCA-DE Genes----

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
# row-scale genes
mat_bc_AccepVsSyn_scaled <- t(scale(t(mat_bc_AccepVsSyn)))
# remove genes with zero variance if any
mat_bc_AccepVsSyn_scaled <- mat_bc_AccepVsSyn_scaled[complete.cases(mat_bc_AccepVsSyn_scaled), , drop = FALSE]

samples_AccepVsSyn<-colnames(dds_IsletTransplant_AccepVsSyn)
# annotation
anno_col <- anno_col[samples_AccepVsSyn, , drop = FALSE]

immune_genes_plot <- c(
  "Acod1",
  "Arg1",
  "Arg2",
  "Art2b",
  "Bst1",
  "Ccl3",
  "Ccl4",
  "Cd14",
  "Cd274",
  "Cd300lf",
  "Cd38",
  "Clec4e",
  "Csf2",
  "Csf3",
  "Csf3r",
  "Cxcl1",
  "Cxcl2",
  "Cxcl3",
  "Cxcr2",
  "Flicr",
  "Fpr1",
  "Fpr2",
  "Gbp5",
  "Gpr15",
  "Gpr84",
  "Hcar2",
  "Hp",
  "Ier3",
  "Ifitm1",
  "Il12a",
  "Il17a",
  "Il17f",
  "Il1a",
  "Il1b",
  "Il18rap",
  "Il22",
  "Il23a",
  "Kcna3",
  "Lcn2",
  "Lef1",
  "Ly6c2",
  "Marco",
  "Mcemp1",
  "Mmp8",
  "Nlrp3",
  "Nos2",
  "Prok2",
  "Slamf6",
  "Themis",
  "Tnf",
  "Tox",
  "Trem1"
)
library(ComplexHeatmap)
library(circlize)
at_idx <- which(rownames(mat_bc_AccepVsSyn_scaled) %in% immune_genes_plot)
lab <- rownames(mat_bc_AccepVsSyn_scaled)[at_idx]

ha_col <- HeatmapAnnotation(
  df = anno_col,
  col = annotation_colors
)
ra <- rowAnnotation(
  mark = anno_mark(
    at = at_idx,
    labels = lab,
    side = "right",
    labels_gp = gpar(fontsize = 10, fontface = "italic"),
    link_width = unit(8, "mm"),
    link_height = unit(2, "mm")
  )
)
col_fun <- colorRamp2(
  seq(-2, 2, length.out = 100),
  colorRampPalette(c("navy", "white", "firebrick3"))(100)
)

ht <- Heatmap(
  mat_bc_AccepVsSyn_scaled,
  name = "Z-score",
  top_annotation = ha_col,
  right_annotation = ra,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = paste0(
    "Allogeneic Acceptance vs Syngeneic (n = ",
    nrow(mat_bc_AccepVsSyn_scaled), ")"
  ),
  heatmap_legend_param = list(
    title = "Expression",
    legend_direction = "vertical"
  )
)

pdf("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/AcceptanceVsSyngeneic/Figures/SigDE_AccepVsSyn_Heatmap.pdf", width = 7, height = 9)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

#PCA Analysis
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
  #"SU_HO_CONV_CENT_CHONDROSARCOMA_LEUKOCYTE_C0_M1_MACROPHAGE",
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

# Remove only sequencing batch
cd <- as.data.frame(colData(dds_IsletTransplant_AccepVsSyn))

cd$Day <- factor(cd$Day, levels = c(7,14,28,42,56,70))
cd$Group <- factor(cd$Group, levels = c("Control Syngeneic", "Acceptance"))
design_gsva <- model.matrix(~ Day, data = cd)
mat_gsva <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch#,
  #design = design_gsva
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

### Allo Accep vs Syn Score----
ALL_ACCEPVSSYN<-read.csv("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/AcceptanceVsSyngeneic/DESEQResults__AccepVsSyn_DaysCombined.csv", row.names = 1)


# VST normalization
vsd <- vst(dds_IsletTransplant_AccepVsSyn, blind = FALSE)

mat <- assay(vsd)

# Remove only sequencing batch
cd <- as.data.frame(colData(dds_IsletTransplant_AccepVsSyn))
design_gsva <- model.matrix(~ Day + Group + Day:Group, data = cd)
mat_enet <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  design = design_gsva
)
up_genes_AccepVsSyn <- ALL_ACCEPVSSYN$gene[
  !is.na(ALL_ACCEPVSSYN$padj) &
    ALL_ACCEPVSSYN$padj <= 0.05 &
    ALL_ACCEPVSSYN$log2FoldChange >= 1
]

down_genes_AccepVsSyn <- ALL_ACCEPVSSYN$gene[
  !is.na(ALL_ACCEPVSSYN$padj) &
    ALL_ACCEPVSSYN$padj <= 0.05 &
    ALL_ACCEPVSSYN$log2FoldChange <= -1
]

gsva_par_sig <- gsvaParam(
  mat_enet,
  list(
    Up = up_genes_AccepVsSyn,
    Down = down_genes_AccepVsSyn
  ),
  kcdf = "Gaussian",
  minSize = 5
)

sig_scores <- gsva(gsva_par_sig)

signature_df <- data.frame(
  Sample = colnames(sig_scores),
  GSVA_Up = as.numeric(sig_scores["Up", ]),
  GSVA_Down = as.numeric(sig_scores["Down", ])
)

signature_df$SignatureScore <- signature_df$GSVA_Up - signature_df$GSVA_Down
signature_df <- cbind(signature_df, cd[signature_df$Sample, ])

score_avg <- signature_df %>%
  dplyr::group_by(Group, Day) %>%
  dplyr::summarise(
    mean_score = mean(SignatureScore, na.rm = TRUE),
    se = sd(SignatureScore, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

pvals_score <- signature_df %>%
  dplyr::group_by(Day) %>%
  dplyr::summarise(
    p_value = {
      df <- dplyr::cur_data()
      if (length(unique(df$Group)) < 2) {
        NA_real_
      } else {
        tryCatch(
          t.test(SignatureScore ~ Group, data = df, var.equal = FALSE)$p.value,
          error = function(e) NA_real_
        )
      }
    },
    .groups = "drop"
  ) %>%
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
 
label_pos_score <- signature_df %>%
  dplyr::group_by(Day) %>%
  dplyr::summarise(
    ymin = min(SignatureScore, na.rm = TRUE),
    ymax = max(SignatureScore, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    span = dplyr::if_else(ymax > ymin, ymax - ymin, 0.2),
    y_pos = ymax - 0.05 * span
  ) %>%
  dplyr::select(Day, y_pos)

sig_score_df <- pvals_score %>%
  dplyr::left_join(label_pos_score, by = "Day")

ggplot(score_avg,
       aes(x = Day,
           y = mean_score,
           color = Group,
           shape = Group,
           group = Group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_score - se,
        ymax = mean_score + se),
    width = 0.5
  ) +
  geom_text(
    data = sig_score_df,
    aes(x = Day, y = y_pos, label = signif_label),
    inherit.aes = FALSE,
    size = 7
  ) +
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
    y = "Allo Acceptance vs Syn Score",
    title = "Allogeneic Acceptance vs Syngeneic Score Over Time"
  )

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
keep <- rowSums(counts(dds_IsletTransplant_IS_RejVsAllo) >= 10) >= 2  # Smallest number of samples in a group
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

## Heatmap and PCA-DE Genes----

sig_genes_IS_REJVSALLO<- ALL_IS_REJVSALLO$gene[
  !is.na(ALL_IS_REJVSALLO$padj) &
    ALL_IS_REJVSALLO$padj <= 0.05 &
    abs(ALL_IS_REJVSALLO$log2FoldChange) >= 1
]

length(sig_genes_IS_REJVSALLO)


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
# row-scale genes
mat_bc_RejVsAllo_scaled <- t(scale(t(mat_bc_RejVsAllo)))
# remove genes with zero variance if any
mat_bc_RejVsAllo_scaled <- mat_bc_RejVsAllo_scaled[complete.cases(mat_bc_RejVsAllo_scaled), , drop = FALSE]

samples_RejVsAllo<-colnames(dds_IsletTransplant_IS_RejVsAllo)
# annotation
anno_col <- anno_col[samples_RejVsAllo, , drop = FALSE]


col_fun <- colorRamp2(
  seq(-2, 2, length.out = 100),
  colorRampPalette(c("navy", "white", "firebrick3"))(100)
)
ha_col <- HeatmapAnnotation(
  df = anno_col,
  col = annotation_colors
)
ht <- Heatmap(
  mat_bc_RejVsAllo_scaled,
  name = "Z-score",
  top_annotation = ha_col,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = paste0(
    "Allogeneic Rejection vs Control Allogeneic (n = ",
    nrow(mat_bc_RejVsAllo_scaled), ")"
  ),
  heatmap_legend_param = list(
    title = "Expression",
    legend_direction = "vertical"
  )
)

draw(ht)

#PCA Analysis
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
  #keyType       = "SYMBOL",       # <- tell it explicitly
  TERM2GENE     = mm_all_df
)
gsea_results_ISREJVSALLO <- as.data.frame(gsea_results_ISREJVSALLO)


# Full results Combined Timepoins
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

