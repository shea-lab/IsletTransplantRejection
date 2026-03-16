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
mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "useast"
)
# get gene biotypes
gene_info <- getBM(
  attributes = c("mgi_symbol", "gene_biotype"),
  filters = "mgi_symbol",
  values = rownames(IsletTransplantCounts),
  mart = mart
)
#setdiff(rownames(IsletTransplantCounts), gene_info$mgi_symbol)
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

dds_IsletTransplant_AlloVsSyn$Batch       <- factor(dds_IsletTransplant_AlloVsSyn$Batch)
dds_IsletTransplant_AlloVsSyn$LibraryPrep <- factor(dds_IsletTransplant_AlloVsSyn$LibraryPrep)
dds_IsletTransplant_AlloVsSyn$Group <- factor(
  dds_IsletTransplant_AlloVsSyn$Group,
  levels = c("Control Syngeneic", "Control Allogeneic")  # order sets baseline
)

dds_IsletTransplant_AlloVsSyn$Day <- factor(dds_IsletTransplant_AlloVsSyn$Day,
                                            levels = c(7, 14))  # baseline = 7


# Check for Collinearuty
cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
# Check sample distribution
lapply(cd[, c("Batch","LibraryPrep","logIEQ","Day","Group")], function(x) table(x, useNA="ifany"))
# Model matrix rank
mm <- model.matrix(~ Batch + LibraryPrep + logIEQ+ Day + Group, data = cd)
qr(mm)$rank; ncol(mm)             # if rank < ncol(mm), not full rank and there are collinear columns
#Here bacth and LibraryPrep are collinear

#Check if logIEQ varies systemically with Group or Day
boxplot(logIEQ ~ Group, data = colData(dds_IsletTransplant_AlloVsSyn)) #Signficantly different between groups


#Design Formula
design(dds_IsletTransplant_AlloVsSyn) <- ~ Batch +  Day + Group + Day:Group #IEQ not included since it is higher in allo vs syn group and a function of group
table(cd$Group,cd$Day)# Check the minimum number of samples in a group
keep <- rowSums(counts(dds_IsletTransplant_AlloVsSyn) >= 10) >= 4  # Smallest number of samples in a group
dds_IsletTransplant_AlloVsSyn <- dds_IsletTransplant_AlloVsSyn[keep, ]
dds_IsletTransplant_AlloVsSyn <- DESeq(dds_IsletTransplant_AlloVsSyn)
resultsNames(dds_IsletTransplant_AlloVsSyn)

# Day 7
res_ALLOvSYN_D7 <- results(dds_IsletTransplant_AlloVsSyn,
                           name = "Group_Control.Allogeneic_vs_Control.Syngeneic")

# Day 14
res_ALLOvSYN_D14 <- results(dds_IsletTransplant_AlloVsSyn, 
                            list(c("Group_Control.Allogeneic_vs_Control.Syngeneic", "Day14.GroupControl.Allogeneic")))


summary(res_ALLOvSYN_D7)
summary(res_ALLOvSYN_D14)


# Thresholds
q_cut  <- 0.10
fc_cut <- 1


## Prep results as data.frames with a 'gene' column
d7_ALLOvSYN  <- as.data.frame(res_ALLOvSYN_D7);  d7_ALLOvSYN$gene  <- rownames(res_ALLOvSYN_D7)
d14_ALLOvSYN <- as.data.frame(res_ALLOvSYN_D14); d14_ALLOvSYN$gene <- rownames(res_ALLOvSYN_D14)


write.csv(d7_ALLOvSYN, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/DESEQResults_Day7_Allo_vs_Syn.csv", row.names = TRUE)
write.csv(d14_ALLOvSYN, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/DESEQResults_Day14_Allo_vs_Syn.csv", row.names = TRUE)

# Function to create color mapping
make_keyvals_fdr_fc <- function(df, q = 0.10, fc = 1,
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

keyvals_d7_ALLOvSYN  <- make_keyvals_fdr_fc(d7_ALLOvSYN,fc=fc_cut)
keyvals_d14_ALLOvSYN <- make_keyvals_fdr_fc(d14_ALLOvSYN,fc=fc_cut)

# Select top 30 genes to plot
pick_labels <- function(df, q = 0.10, fc = 1, topN = 30) {
  idx <- which(!is.na(df$padj) & !is.na(df$log2FoldChange) &
                 df$padj < q & abs(df$log2FoldChange) >= fc)
  if (length(idx) == 0) return(character(0))
  ord <- order(df$padj[idx], -abs(df$log2FoldChange[idx]), na.last = NA)  # tie-break by |LFC|
  labs <- df$gene[idx][ord]
  labs <- make.unique(labs)  # avoid dup labels
  labs[seq_len(min(topN, length(labs)))]
}

selLab_d7_ALLOvSYN  <- pick_labels(d7_ALLOvSYN,  q_cut, fc_cut, 30)
selLab_d14_ALLOvSYN <- pick_labels(d14_ALLOvSYN, q_cut, fc_cut, 30)

##  Axis limits
xmax_d7_ALLOvSYN  <- max(2, ceiling(max(abs(d7_ALLOvSYN$log2FoldChange),  na.rm=TRUE)))
xmax_d14_ALLOvSYN <- max(2, ceiling(max(abs(d14_ALLOvSYN$log2FoldChange), na.rm=TRUE)))

## Volcano Plots----
#Day 7
EnhancedVolcano(
  d7_ALLOvSYN,
  lab           = d7_ALLOvSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Allogeneic vs Syngeneic — Day 7",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d7_ALLOvSYN$padj<0.10 & abs(d7_ALLOvSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-5, 5),
  ylim = c(0,4),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_d7_ALLOvSYN,
  legendPosition= "right",
  selectLab     = selLab_d7_ALLOvSYN
)

## Day 14
EnhancedVolcano(
  d14_ALLOvSYN,
  lab           = d14_ALLOvSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Allogeneic vs Syngeneic — Day 14",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d14_ALLOvSYN$padj<0.10 & abs(d14_ALLOvSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d14_ALLOvSYN, xmax_d14_ALLOvSYN),
  ylim = c(0,10),
  boxedLabels   = TRUE,
  pointSize     =3,
  labSize       = 6,
  colAlpha      = 0.8,
  max.overlaps = 15,
  drawConnectors= TRUE,
  colCustom     = keyvals_d14_ALLOvSYN,
  legendPosition= "right",
  selectLab     = selLab_d14_ALLOvSYN
)

## Heatmap-----
library(pheatmap)
# significant genes for each day
sig_d7_genes <- d7_ALLOvSYN$gene[
  !is.na(d7_ALLOvSYN$padj) &
    d7_ALLOvSYN$padj <= q_cut &
    abs(d7_ALLOvSYN$log2FoldChange) >= fc_cut
]

sig_d14_genes <- d14_ALLOvSYN$gene[
  !is.na(d14_ALLOvSYN$padj) &
    d14_ALLOvSYN$padj <= q_cut &
    abs(d14_ALLOvSYN$log2FoldChange) >= fc_cut
]

length(sig_d7_genes)
length(sig_d14_genes)

vsd <- vst(dds_IsletTransplant_AlloVsSyn, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
design_heatmap <- model.matrix(~ Day + Group + Day:Group, data = cd)
mat_bc <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  design = design_heatmap
)

anno_col <- as.data.frame(
  colData(dds_IsletTransplant_AlloVsSyn)[, "Group", drop = FALSE]
)
annotation_colors <- list(
  Group = c(
    "Control Allogeneic" = "#FF2400",   # Scarlet
    "Control Syngeneic" = "#2E6F40"     # Moss green
  )
)
# Day 7 samples
samples_d7 <- rownames(colData(dds_IsletTransplant_AlloVsSyn))[colData(dds_IsletTransplant_AlloVsSyn)$Day == 7]
# matrix for Day 7
mat_d7 <- mat_bc[sig_d7_genes, samples_d7, drop = FALSE]
# row-scale genes
mat_d7_scaled <- t(scale(t(mat_d7)))
# remove genes with zero variance if any
mat_d7_scaled <- mat_d7_scaled[complete.cases(mat_d7_scaled), , drop = FALSE]
# annotation
anno_col_d7 <- anno_col[samples_d7, , drop = FALSE]

pheatmap(
  mat_d7_scaled,
  annotation_col = anno_col_d7,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = paste0("Day 7 DE genes: Allogeneic vs Syngeneic (n = ", nrow(mat_d7_scaled), ")"),
  scale = "none",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)

# Day 14 samples
samples_d14 <- rownames(colData(dds_IsletTransplant_AlloVsSyn))[colData(dds_IsletTransplant_AlloVsSyn)$Day == 14]
# matrix for Day 7
mat_d14 <- mat_bc[sig_d14_genes, samples_d14, drop = FALSE]
# row-scale genes
mat_d14_scaled <- t(scale(t(mat_d14)))
# remove genes with zero variance if any
mat_d14_scaled <- mat_d14_scaled[complete.cases(mat_d14_scaled), , drop = FALSE]
# annotation
anno_col_d14 <- anno_col[samples_d14, , drop = FALSE]

pheatmap(
  mat_d14_scaled,
  annotation_col = anno_col_d14,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = paste0("Day 14 DE genes: Allogeneic vs Syngeneic (n = ", nrow(mat_d14_scaled), ")"),
  scale = "none",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)


## Elastic Net Feature Selection----

# DO analysis for timepoint combined
design(dds_IsletTransplant_AlloVsSyn) <- ~ Batch +  Day + Group  #IEQ not included since it is higher in allo vs syn group and a function of group
keep <- rowSums(counts(dds_IsletTransplant_AlloVsSyn) >= 10) >= 4  # Smallest number of samples in a group
dds_IsletTransplant_AlloVsSyn <- dds_IsletTransplant_AlloVsSyn[keep, ]
dds_IsletTransplant_AlloVsSyn <- DESeq(dds_IsletTransplant_AlloVsSyn)
resultsNames(dds_IsletTransplant_AlloVsSyn)

# All Days combined
res_ALLOvSYN <- results(dds_IsletTransplant_AlloVsSyn,
                        name = "Group_Control.Allogeneic_vs_Control.Syngeneic")

genes_ALLOvSYN  <- as.data.frame(res_ALLOvSYN);  genes_ALLOvSYN$gene  <- rownames(res_ALLOvSYN)

write.csv(genes_ALLOvSYN, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/DESEQResults_DayCombined_Allo_vs_Syn.csv", row.names = TRUE)

sig_genes_ALLOvSYN <- genes_ALLOvSYN$gene[
  !is.na(genes_ALLOvSYN$pvalue) &
    genes_ALLOvSYN$pvalue <= 0.05 &
    abs(genes_ALLOvSYN$log2FoldChange) >= 0.5
]

length(sig_genes_ALLOvSYN)

vsd <- vst(dds_IsletTransplant_AlloVsSyn, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_IsletTransplant_AlloVsSyn))
design_enet <- model.matrix(~ Group, data = cd)
mat_enet <- limma::removeBatchEffect(
  mat,
  batch = vsd$Batch,
  covariates = model.matrix(~ Day, data = cd)[, -1, drop = FALSE],
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

# Stability selection
set.seed(123)
n_iter <- 1000
n_samples <- nrow(x_AlloVsSyn)

selected_list <- vector("list", n_iter)

for (i in 1:n_iter) {
  
  # Subsample ~80% of samples each time
  idx <- sample(1:n_samples, size = round(0.8 * n_samples))
  
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
d7_path  <- "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/DESEQResults_Day7_Allo_vs_Syn.csv"
d14_path <- "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/DESEQResults_Day14_Allo_vs_Syn.csv"
day_combined<-"/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/DESEQResults_DayCombined_Allo_vs_Syn.csv"

# Import
d7_ALLOvSYN  <- read.csv(d7_path,  row.names = 1)
d14_ALLOvSYN <- read.csv(d14_path, row.names = 1)
daycombined_ALLOvSYN<- read.csv(day_combined, row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_d7  <- d7_ALLOvSYN$stat;  names(lfc_vector_d7)  <- rownames(d7_ALLOvSYN)
lfc_vector_d14 <- d14_ALLOvSYN$stat; names(lfc_vector_d14) <- rownames(d14_ALLOvSYN)
lfc_vector_AlloSyn <- daycombined_ALLOvSYN$stat; names(lfc_vector_AlloSyn) <- rownames(daycombined_ALLOvSYN)

# Drop NAs
lfc_vector_d7  <- lfc_vector_d7[!is.na(lfc_vector_d7)]
lfc_vector_d14 <- lfc_vector_d14[!is.na(lfc_vector_d14)]
lfc_vector_AlloSyn<- lfc_vector_AlloSyn[!is.na(lfc_vector_AlloSyn)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_d7  <- sort(lfc_vector_d7,  decreasing = TRUE)
lfc_vector_d14 <- sort(lfc_vector_d14, decreasing = TRUE)
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

# Day 7
gsea_results_d7 <- GSEA(
  geneList      = lfc_vector_d7,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
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


# Full results Day 7
write.csv(gsea_results_d7_df,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/GSEAResults_AlloVsSyn_Day7.csv",
          row.names = FALSE)

# Full results Day 14
write.csv(gsea_results_d14_df,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/GSEAResults_AlloVsSyn_Day14.csv",
          row.names = FALSE)

# Full results Day Combined
write.csv(gsea_results_daycombined_allosyn_df,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/Allogeneic_Vs_Syngeneic/GSEAResults_AlloVsSyn_DayCombined.csv",
          row.names = FALSE)


# Per Day Plot

# Your curated list EXACTLY as provided (de-dup just in case)
immune_master <- unique(c(
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  #"HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  #"HALLMARK_TGF_BETA_SIGNALING",
  "AIZARANI_LIVER_C5_NK_NKT_CELLS_3",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY",
  "AIZARANI_LIVER_C6_KUPFFER_CELLS_2",
  "DESCARTES_FETAL_LIVER_LYMPHOID_CELLS",
  "DESCARTES_FETAL_PANCREAS_LYMPHOID_CELLS",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  "KEGG_PEROXISOME",
  "KEGG_SPHINGOLIPID_METABOLISM",
  "KEGG_TRYPTOPHAN_METABOLISM",
  #Down
  "HALLMARK_TGF_BETA_SIGNALING"
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
  scale_size_continuous(name = "−log10(padj)", range = c(1, 10)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  labs(x = "", y = "", title = "Allogeneic vs Syngeneic") +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  )
print(p)


# Combined Day Plot

library(forcats)
pathways_of_interest <- unique(c(
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
  #Down
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "KEGG_INSULIN_SIGNALING_PATHWAY",
  "AIZARANI_LIVER_C23_KUPFFER_CELLS_3",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
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

IsletScRNA = readRDS("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/Islet ScRNASeq/islet_graft_seurat_v7.rds")
IsletScRNA$celltype<-Idents(IsletScRNA)
IsletScRNA$condition[IsletScRNA$condition == "Allogenic"] <- "Allogeneic"
DimPlot(IsletScRNA,label = TRUE, label.box = T,label.size = 8,repel = T,pt.size = 0.9)+
  NoAxes() +NoLegend()
DimPlot(IsletScRNA,label = TRUE, label.box = T,label.size = 8,repel = T,group.by = "condition",pt.size = 0.9)+
  NoAxes() +NoLegend()
#Downregulated and Upregulated wrt Allo vs Syn
downregulated_Signature <- c("Cd59b","Rep15","Shisa9","Gdf3","Xist","Eif2s3x","Kdm6a","Lncpint","Ifitm5","Malat1","Myo3b")
upregulated_Signature <- c("Hipk4","S1pr5","Trpm6","Tctn2","Lhfpl4","Ido2","Fsip1","Emx2os","Hlf","Gfra1","Bmp5","Col6a5","Myh3","Cbs","Vtn")

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


# Elastic Net Score
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




# 5. Allo+Anti-CD40L-Rejected vs Tolerance-----

# subset samples
sel <- colData(dds_IsletTransplant)$Group %in% c("Tolerance","Rejected") & colData(dds_IsletTransplant)$Day %in% c(7,14)
dds_IsletTransplant_RejVsTol <- dds_IsletTransplant[, sel]

dds_IsletTransplant_RejVsTol$Batch       <- factor(dds_IsletTransplant_RejVsTol$Batch)
dds_IsletTransplant_RejVsTol$LibraryPrep <- factor(dds_IsletTransplant_RejVsTol$LibraryPrep)
dds_IsletTransplant_RejVsTol$Group <- factor(
  dds_IsletTransplant_RejVsTol$Group,
  levels = c("Tolerance", "Rejected")  # order sets baseline
)
levels(dds_IsletTransplant_RejVsTol$Group)

dds_IsletTransplant_RejVsTol$Day <- factor(dds_IsletTransplant_RejVsTol$Day,
                                            levels = c(7, 14))  # baseline = 7
levels(dds_IsletTransplant_RejVsTol$Day)

cd <- as.data.frame(colData(dds_IsletTransplant_RejVsTol))
mm <- model.matrix(~ Batch + LibraryPrep + Day + Group, data = cd)
qr_mm <- qr(mm)

kept    <- colnames(mm)[qr_mm$pivot[seq_len(qr_mm$rank)]]
dropped <- if (qr_mm$rank < ncol(mm)) colnames(mm)[qr_mm$pivot[(qr_mm$rank+1):ncol(mm)]] else character()

kept
dropped

#Design Formula
design(dds_IsletTransplant_RejVsTol) <- ~ Batch + Day + Group + Day:Group
dds_IsletTransplant_RejVsTol <- DESeq(dds_IsletTransplant_RejVsTol)
design(dds_IsletTransplant_RejVsTol)
resultsNames(dds_IsletTransplant_RejVsTol)

# Day 7
res_REJvTOL_D7 <- results(dds_IsletTransplant_RejVsTol,
                           name = "Group_Rejected_vs_Tolerance")

# Day 14
res_REJvTOL_D14 <- results(dds_IsletTransplant_RejVsTol,
                            list(c("Group_Rejected_vs_Tolerance",
                                   "Day14.GroupRejected")))

# Interactions
res_interaction_df_REJvTOL <- results(dds_IsletTransplant_RejVsTol,
                           name = "Day14.GroupRejected")
res_interaction_df_REJvTOL <- as.data.frame(res_interaction_df_REJvTOL)
res_interaction_df_REJvTOL$gene <- rownames(res_interaction_df_REJvTOL)

summary(res_REJvTOL_D7)
summary(res_REJvTOL_D14)

library(EnhancedVolcano)

# Thresholds
q_cut  <- 0.10
fc_cut <- 1

library(EnhancedVolcano)

## 0) Prep results as data.frames with a 'gene' column
d7_REJVSTOL  <- as.data.frame(res_REJvTOL_D7);  d7_REJVSTOL$gene  <- rownames(res_REJvTOL_D7)
d14_REJVSTOL <- as.data.frame(res_REJvTOL_D14); d14_REJVSTOL$gene <- rownames(res_REJvTOL_D14)

write.csv(d7_REJVSTOL, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsRejection_AntiCD40LTreated/DESEQResults_Day7_Rejection_vs_Tolerance.csv", row.names = TRUE)
write.csv(d14_REJVSTOL, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsRejection_AntiCD40LTreated/DESEQResults_Day14_Rejection_vs_Tolerancen.csv", row.names = TRUE)
write.csv(res_interaction_df_REJvTOL, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsRejection_AntiCD40LTreated/DESEQResults_InteractionTimeGroup_Rejection_vs_Tolerance.csv", row.names = TRUE)

# Function to create color mapping
make_keyvals_fdr_fc_TolvsRej <- function(df, q = 0.10, fc = 1,
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
  
  label[up]   <- paste0("Rejection-Upregulated")
  label[down] <- paste0("Tolerance-Upregulated")
  
  names(key) <- label        # <- legend labels; no NAs
  key
}
keyvals_d7  <- make_keyvals_fdr_fc_TolvsRej(d7_REJVSTOL)
keyvals_d14 <- make_keyvals_fdr_fc_TolvsRej(d14_REJVSTOL)

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

selLab_d7_REJVSTOL  <- pick_labels(d7_REJVSTOL,  q_cut, fc_cut, 30)
selLab_d14_REJVSTOL <- pick_labels(d14_REJVSTOL, q_cut, fc_cut, 30)


## 3) Axis limits
xmax_d7  <- max(2, ceiling(max(abs(d7_REJVSTOL$log2FoldChange),  na.rm=TRUE)))
xmax_d14 <- max(2, ceiling(max(abs(d14_REJVSTOL$log2FoldChange), na.rm=TRUE)))

## Volcano: Day 7----
EnhancedVolcano(
  d7_REJVSTOL,
  lab           = d7_REJVSTOL$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Rejection vs Tolerance — Day 7",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d7_REJVSTOL$padj<0.10 & abs(d7_REJVSTOL$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d7, xmax_d7),
  ylim = c(0,6),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_d7,
  legendPosition= "right",
  selectLab     = selLab_d7_REJVSTOL
)

## Volcano: Day 14----
EnhancedVolcano(
  d14_REJVSTOL,
  lab           = d14_REJVSTOL$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Rejection vs Tolerance — Day 14",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d14_REJVSTOL$padj<0.10 & abs(d14_REJVSTOL$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d14, xmax_d14),
  ylim = c(0,9),
  boxedLabels   = TRUE,
  pointSize     =3,
  labSize       = 6,
  colAlpha      = 0.8,
  max.overlaps = 35,
  drawConnectors= TRUE,
  colCustom     = keyvals_d14,
  legendPosition= "right",
  selectLab     = selLab_d14_REJVSTOL
)



## Timepoints Combined----

design(dds_IsletTransplant_RejVsTol) <- ~ Batch + Day + Group
dds_IsletTransplant_RejVsTol <- DESeq(dds_IsletTransplant_RejVsTol)
resultsNames(dds_IsletTransplant_RejVsTol)

# Combined Timrpoints
res_REJvTOL_ALL <- results(dds_IsletTransplant_RejVsTol,
                          name = "Group_Rejected_vs_Tolerance")

summary(res_REJvTOL_ALL)


library(EnhancedVolcano)

# Thresholds
q_cut  <- 0.10
fc_cut <- 1

library(EnhancedVolcano)

## 0) Prep results as data.frames with a 'gene' column
ALL_REJVSTOL  <- as.data.frame(res_REJvTOL_ALL);  ALL_REJVSTOL$gene  <- rownames(res_REJvTOL_ALL)

write.csv(ALL_REJVSTOL, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsRejection_AntiCD40LTreated/DESEQResults_Day7_14_Combined_Rejection_vs_Tolerance.csv", row.names = TRUE)


keyvals_ALL  <- make_keyvals_fdr_fc_TolvsRej(ALL_REJVSTOL)

library(dplyr)

selLab_ALL_REJVSTOL  <- pick_labels(ALL_REJVSTOL,  q_cut, fc_cut, 30)


## 3) Axis limits
xmax_ALL  <- max(2, ceiling(max(abs(ALL_REJVSTOL$log2FoldChange),  na.rm=TRUE)))

EnhancedVolcano(
  ALL_REJVSTOL,
  lab           = ALL_REJVSTOL$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Rejection vs Tolerance — Combined Timepoints",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(ALL_REJVSTOL$padj<0.10 & abs(ALL_REJVSTOL$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_ALL, xmax_ALL),
  ylim = c(0,8),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 35,
  colCustom     = keyvals_ALL,
  legendPosition= "right",
  selectLab     = selLab_ALL_REJVSTOL
)

## PLSDA Analysis-Combined Timepoints ----

library(DESeq2)
library(limma)
library(mixOmics)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# --- Build the  label from your fitted object ---
CountsRejVsTol <- as.data.frame(colData(dds_IsletTransplant_RejVsTol))

# --- Variance-stabilized expression ---
vsd <- vst(dds_IsletTransplant_RejVsTol, blind = FALSE)
mat <- assay(vsd)  # genes x samples
# --- Remove batch (for visualization only) ---
mat_corr <- removeBatchEffect(mat, batch = colData(vsd)$Batch)

# --- PLS-DA input: samples x genes ---
X <- t(mat_corr)   # samples in rows
Y <- CountsRejVsTol$Group        # factor with 4 levels

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
  "Rejected"  = "#D62728",  # base red
  "Tolerance"  = "#08306B"  # darker shade of blue (navy/steel blue)
)
group_shapes <- c(
  "Rejected"  = 16,
  "Tolerance"   = 15
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
    title = "PLS-DA: Rejection vs Tolerance",
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
write.csv(vip_df_axis1, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsRejection_AntiCD40LTreated/VIP_scores_RejVsTol_PLS1.csv", row.names = FALSE)
# Rank genes by VIP (descending)
top100_genes <- names(sort(vip_scores[, 1], decreasing = TRUE))[1:100]
mat_top100 <- mat_corr[top100_genes, ]
anno <- data.frame(
  Group = Y,        # group labels (factor)
  row.names = colnames(mat_top100)
)

anno_colors <- list(
  Group = group_colors   # same color scheme you used before
)

Y <- factor(Y, levels = c("Tolerance","Rejected"))

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
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50) # publication-style
)

### Error Rate in PLSDA----
# Ensure Y is a factor (important)
Y <- factor(Y)

set.seed(123)
perf_plsda <- perf(plsda_model,
                   validation = "Mfold",
                   folds = 5,
                   nrepeat = 10,          # increase if you want more stable estimates
                   dist = "max.dist",     # common choice; try "centroids.dist" too
                   progressBar = TRUE)

## 1) Overall classification error rate (per component, per repeat)
perf_plsda$error.rate$overall        # matrix: repeats x ncomp

## Mean overall error across repeats (per component)
colMeans(perf_plsda$error.rate$overall, na.rm = TRUE)

## 2) Balanced error rate (handles class imbalance)
perf_plsda$error.rate$BER            # matrix: repeats x ncomp
colMeans(perf_plsda$error.rate$BER, na.rm = TRUE)

##) 95% CI for error (quick & dirty via repeats)
er <- perf_plsda$error.rate$overall
mean_er <- colMeans(er, na.rm = TRUE)
se_er   <- apply(er, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
ci_low  <- mean_er - 1.96 * se_er
ci_high <- mean_er + 1.96 * se_er
data.frame(comp = seq_along(mean_er), mean_er, ci_low, ci_high)

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
vsd <- vst(dds_IsletTransplant_RejVsTol, blind = FALSE)
mat  <- assay(vsd)
matc <- removeBatchEffect(mat, batch = colData(vsd)$Batch)

# keep only immune genes from earlier step
mat_filt   <- matc[rownames(matc) %in% immune_genes, ]
mat_scaled <- t(scale(t(mat_filt)))

heat_colors <- colorRampPalette(c("navy","white","firebrick3"))(100)

day_vec <- as.character(colData(vsd)$Day)
group_vec <- as.character(colData(vsd)$Group)

ann_colors <- list(
  Group = c("Rejected" = "#D62728", "Tolerance" = "#1F77B4"),
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
### Combined Timepoints----
# Paths to your saved results
ALL_REJVSTOL_path  <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsRejection_AntiCD40LTreated/DESEQResults_Day7_14_Combined_Rejection_vs_Tolerance.csv"

# Import
ALL_REJVSTOL  <- read.csv(ALL_REJVSTOL_path,  row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_ALL_REJVSTOL  <- ALL_REJVSTOL$stat;  names(lfc_vector_ALL_REJVSTOL)  <- rownames(ALL_REJVSTOL)


# Drop NAs
lfc_vector_ALL_REJVSTOL  <- lfc_vector_ALL_REJVSTOL[!is.na(lfc_vector_ALL_REJVSTOL)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_ALL_REJVSTOL  <- sort(lfc_vector_ALL_REJVSTOL,  decreasing = TRUE)

gsea_results_ALL_REJVSTOL <- GSEA(
  geneList      = lfc_vector_ALL_REJVSTOL,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  #keyType       = "SYMBOL",       # <- tell it explicitly
  TERM2GENE     = mm_all_df
)
gsea_results_ALL_REJVSTOL <- as.data.frame(gsea_results_ALL_REJVSTOL)


# Full results Combined Timepoins
write.csv(gsea_results_ALL_REJVSTOL,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsRejection_AntiCD40LTreated/GSEAResults_RejectionVsTolerance_CombinedTimepoints.csv",
          row.names = FALSE)


library(dplyr)
library(stringr)
library(ggplot2)

# Your curated list EXACTLY as provided (de-dup just in case)
immune_master_ALLREJVSTOL <- unique(c(
  #UPREGULATED
  # Immune proliferation / effector programs
  # Immune/inflammatory
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING","KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_ALLOGRAFT_REJECTION","KEGG_GRAFT_VERSUS_HOST_DISEASE",
  "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CELL",
  "TRAVAGLINI_LUNG_NATURAL_KILLER_CELL","TRAVAGLINI_LUNG_NATURAL_KILLER_T_CELL",
  "HE_LIM_SUN_FETAL_LUNG_C4_TH17_CELL","HE_LIM_SUN_FETAL_LUNG_C4_CD8_T_CELL","HE_LIM_SUN_FETAL_LUNG_C4_CD4_T_CELL",
  # Endothelium/vascular
  "HALLMARK_ANGIOGENESIS",
   # ECM/stromal activation
  "KEGG_ECM_RECEPTOR_INTERACTION",
  # Stress/fibrosis axis
  "HALLMARK_HYPOXIA","HALLMARK_TGF_BETA_SIGNALING","KEGG_TGF_BETA_SIGNALING_PATHWAY",
 
   #DOWNREGUL:ATED
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION","KEGG_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM","KEGG_FATTY_ACID_METABOLISM","KEGG_PPAR_SIGNALING_PATHWAY",
  "HALLMARK_PEROXISOME","KEGG_PEROXISOME","KEGG_GLUTATHIONE_METABOLISM",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  # Regulatory/naïve lymphoid
  "HALLMARK_IL2_STAT5_SIGNALING","TRAVAGLINI_LUNG_CD4_NAIVE_T_CELL","TRAVAGLINI_LUNG_CD8_NAIVE_T_CELL",
  # Tolerogenic myeloid
  "HE_LIM_SUN_FETAL_LUNG_C2_APOE_POS_M2_MACROPHAGE_CELL",
  "TRAVAGLINI_LUNG_PLASMACYTOID_DENDRITIC_CELL",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  # Islet preservation
  "HALLMARK_PANCREAS_BETA_CELLS","VANGURP_PANCREATIC_BETA_CELL"
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
ALL_REJVSTOL_tbl  <- coerce_gsea_tbl(gsea_results_ALL_REJVSTOL,  "Combined Timepoints")

# Keep ONLY immune_master pathways
ALL_REJVSTOL_tbl  <- ALL_REJVSTOL_tbl  %>% filter(gs_name %in% immune_master_ALLREJVSTOL)

library(dplyr)
library(ggplot2)

plot_df <- ALL_REJVSTOL_tbl %>%
  mutate(
    padj = pmax(padj, .Machine$double.eps),
    logp = -log10(padj)
  ) %>%
  arrange((NES)) %>%                     # order by NES (descending)
  mutate(gs_name = factor(gs_name, levels = unique(gs_name)))

# Dotplot: x = NES, y = pathways (ordered by NES)
p_dot <- ggplot(plot_df, aes(x = NES, y = gs_name)) +
  geom_point(aes(size = logp, color = NES)) +
  scale_size_continuous(name = expression(-log[10]("padj")), range = c(2, 8)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "NES") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL,
       title = "GSEA (Rejection vs Tolerance, Combined Timepoints)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )

print(p_dot)

## GSVA Analysis----

library(GSVA)
library(msigdbr)
library(limma)
library(BiocParallel)

# VST
vsd  <- vst(dds_IsletTransplant_RejVsTol, blind = TRUE)
expr <- assay(vsd)
meta <- as.data.frame(colData(dds_IsletTransplant_RejVsTol))

# Factors
meta$Group <- factor(meta$Group, levels = c("Tolerance", "Rejected"))
if ("Batch" %in% names(meta)) meta$Batch <- factor(meta$Batch)

# --- Gene sets (same as yours) ---
c8   <- msigdbr(species = "Mus musculus", category = "C8")
hall <- msigdbr(species = "Mus musculus", category = "H")
kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY")
sets <- c(
  split(c8$gene_symbol,   c8$gs_name),
  split(hall$gene_symbol, hall$gs_name),
  split(kegg$gene_symbol, kegg$gs_name)
)
# Optional: restrict to your curated list
# sets <- sets[names(sets) %in% immune_master_ALLREJVSTOL]

# --- Select Day 7 + Day 14 samples ---
pick_day <- function(day) which(meta$Day %in% c(day, paste0("Day", day), paste0("Day ", day), paste0("D", day)))
idx <- sort(unique(c(pick_day(7), pick_day(14))))

expr_all <- expr[, idx, drop = FALSE]
meta_all <- meta[idx, , drop = FALSE]

# Normalize Day labels to two levels (Day7, Day14)
norm_day <- function(x) {
  x <- as.character(x)
  out <- ifelse(grepl("(^7$|Day ?7|^D7$)", x, ignore.case = TRUE), "Day7",
                ifelse(grepl("(^14$|Day ?14|^D14$)", x, ignore.case = TRUE), "Day14", NA))
  # Fallback for numeric
  out[is.na(out) & x %in% c("7")]  <- "Day7"
  out[is.na(out) & x %in% c("14")] <- "Day14"
  out
}
meta_all$Day <- factor(norm_day(meta_all$Day), levels = c("Day7","Day14"))

param <- MulticoreParam(workers = max(1, parallel::detectCores() - 1))
gsva_par <- gsvaParam(
  expr_all,
  sets,
  kcdf    = "Gaussian",   # VST is log-like
  minSize = 5,
  maxSize = 500
)
es_all <- gsva(gsva_par, BPPARAM = param)  # pathway x sample

# --- LIMMA on GSVA scores: Group main effect, adjusting for Day (+ Batch if present) ---
run_limma_combined <- function(es, meta) {
  has_day   <- "Day"   %in% names(meta) && nlevels(droplevels(meta$Day))   > 1
  has_batch <- "Batch" %in% names(meta) && nlevels(droplevels(meta$Batch)) > 1
  
  # Build design: ~ 0 + Group + (Day) + (Batch)
  terms <- c("Group", if (has_day) "Day", if (has_batch) "Batch")
  form  <- as.formula(paste("~ 0 +", paste(terms, collapse = " + ")))
  design <- model.matrix(form, data = meta)
  colnames(design) <- make.names(colnames(design))
  
  fit  <- lmFit(es, design)
  contr <- makeContrasts(Rejected_vs_Tolerance = GroupRejected - GroupTolerance,
                         levels = colnames(design))
  fit2 <- eBayes(contrasts.fit(fit, contr))
  topTable(fit2, number = Inf, sort.by = "P")
}

GSVA_res_combined <- run_limma_combined(es_all, meta_all)

# --- Save outputs ---
out_dir <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsRejection_AntiCD40LTreated"
write.csv(es_all,             file.path(out_dir, "GSVAScores_RejVsTol_Combined_D7_D14.csv"))
write.csv(GSVA_res_combined,  file.path(out_dir, "GSVAStats_RejVsTol_Combined_D7_D14.csv"))

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# --- Step 1: Select significant pathways
# Rename first column as "ID"
GSVA_res_combined$ID<-rownames(GSVA_res_combined)

sig_paths <- GSVA_res_combined %>%
  filter(P.Value < 0.05, abs(logFC) >= 0.2) %>%
  arrange(P.Value) %>%
  pull(ID) %>%
  unique()

# --- Step 2: Subset GSVA scores (es_all is pathways x samples)
mat <- es_all[sig_paths, , drop = FALSE]

# --- Step 3: Z-score per pathway for visualization
Z <- t(scale(t(mat)))
Z[is.na(Z)] <- 0
min(Z)
min_separators()# --- Step 4: Build sample annotations
meta_sub <- meta_all[colnames(Z), , drop = FALSE]

ann_cols <- list(
  Group = c(Tolerance = "#1f77b4", Rejected = "#d62728"),
  Day   = c(Day7 = "#66c2a5", Day14 = "#fc8d62")
)
top_anno <- HeatmapAnnotation(
  Group = meta_sub$Group,
  Day   = meta_sub$Day,
  col   = ann_cols
)

# --- Step 5: Heatmap
col_fun <- colorRamp2(c(min(Z), 0,max(Z)), c("#2b6cb0", "white", "#b91c1c"))

ht <- Heatmap(
  Z, name = "GSVA Z-score",
  col = col_fun,
  top_annotation = top_anno,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = meta_sub$Group,   # split by Group for clarity
  heatmap_legend_param = list(
    title = "GSVA Z-score"
  )
)


draw(
  ht,
  heatmap_legend_side = "left",
  annotation_legend_side = "left",
  padding = unit(c(5, 20, 5, 70), "mm")  # top, right, bottom, left
)

# 5.Allo Rejected vs Allo+Anti-CD40L Rejected -----

# subset samples
sel <- colData(dds_IsletTransplant)$Group %in% c("Control Rejected","Rejected") & colData(dds_IsletTransplant)$Day %in% c(7,14)
#IS- Immunosuppressed
Celldds_IsletTransplant_IS_RejVsRej <- dds_IsletTransplant[, sel]

Celldds_IsletTransplant_IS_RejVsRej$Batch       <- factor(Celldds_IsletTransplant_IS_RejVsRej$Batch)
Celldds_IsletTransplant_IS_RejVsRej$LibraryPrep <- factor(Celldds_IsletTransplant_IS_RejVsRej$LibraryPrep)
Celldds_IsletTransplant_IS_RejVsRej$Group <- factor(
  Celldds_IsletTransplant_IS_RejVsRej$Group,
  levels = c("Control Rejected", "Rejected")  # order sets baseline
)
levels(Celldds_IsletTransplant_IS_RejVsRej$Group)

Celldds_IsletTransplant_IS_RejVsRej$Day <- factor(Celldds_IsletTransplant_IS_RejVsRej$Day,
                                           levels = c(7, 14))  # baseline = 7
levels(Celldds_IsletTransplant_IS_RejVsRej$Day)

cd <- as.data.frame(colData(Celldds_IsletTransplant_IS_RejVsRej))
mm <- model.matrix(~ Batch + LibraryPrep + Day + Group, data = cd)
qr_mm <- qr(mm)

kept    <- colnames(mm)[qr_mm$pivot[seq_len(qr_mm$rank)]]
dropped <- if (qr_mm$rank < ncol(mm)) colnames(mm)[qr_mm$pivot[(qr_mm$rank+1):ncol(mm)]] else character()

kept
dropped

#Design Formula
design(Celldds_IsletTransplant_IS_RejVsRej) <- ~ Batch + Day + Group 
Celldds_IsletTransplant_IS_RejVsRej <- DESeq(Celldds_IsletTransplant_IS_RejVsRej)
design(Celldds_IsletTransplant_IS_RejVsRej)
resultsNames(Celldds_IsletTransplant_IS_RejVsRej)

# Da7+Day14 Combined
res_IS_REJvREJ <- results(Celldds_IsletTransplant_IS_RejVsRej,
                          name = "Group_Rejected_vs_Control.Rejected")


summary(res_IS_REJvREJ)

## DESeq-Combined Timepoints----
library(EnhancedVolcano)

# Thresholds
q_cut  <- 0.10
fc_cut <- 1

library(EnhancedVolcano)

## 0) Prep results as data.frames with a 'gene' column
ALL_IS_REJVSREJ  <- as.data.frame(res_IS_REJvREJ);  ALL_IS_REJVSREJ$gene  <- rownames(res_IS_REJvREJ)

write.csv(ALL_IS_REJVSREJ, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/AntiCD40LRejection_Rejection/DESEQResults_Day7_14_Combined_AntiCD40L_Rejection_vs_Rejection.csv", row.names = TRUE)

# Function to create color mapping
make_keyvals_fdr_fc_IsRejvsRej <- function(df, q = 0.10, fc = 1,
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
  
  label[up]   <- paste0("Anti-CD40L Rejection-Upregulated")
  label[down] <- paste0("Rejection-Upregulated")
  
  names(key) <- label        # <- legend labels; no NAs
  key
}

keyvals_ALL  <- make_keyvals_fdr_fc_IsRejvsRej(ALL_IS_REJVSREJ)

library(dplyr)

selLab_ALL_IS_REJVSREJ  <- pick_labels(ALL_IS_REJVSREJ,  q_cut, fc_cut, 40)


## 3) Axis limits
xmax_ALL  <- max(2, ceiling(max(abs(ALL_IS_REJVSREJ$log2FoldChange),  na.rm=TRUE)))

EnhancedVolcano(
  ALL_IS_REJVSREJ,
  lab           = ALL_IS_REJVSREJ$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Anti-CD40L Rejection vs Control Rejection — Combined Timepoints",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(ALL_IS_REJVSREJ$padj<0.10 & abs(ALL_IS_REJVSREJ$log2FoldChange)>=1, na.rm=TRUE), ")"),
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
  selectLab     = selLab_ALL_IS_REJVSREJ
)

## PLSDA Analysis-Combined Timepoints ----

library(DESeq2)
library(limma)
library(mixOmics)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# --- Build the  label from your fitted object ---
CountsISRejVsRej <- as.data.frame(colData(Celldds_IsletTransplant_IS_RejVsRej))

# --- Variance-stabilized expression ---
vsd_ISRejVSRej <- vst(Celldds_IsletTransplant_IS_RejVsRej, blind = FALSE)
mat_ISRejVSRej <- assay(vsd_ISRejVSRej)  # genes x samples
# --- Remove batch (for visualization only) ---
mat_corr_ISRejVSRej <- removeBatchEffect(mat_ISRejVSRej, batch = colData(vsd_ISRejVSRej)$Batch)

# --- PLS-DA input: samples x genes ---
X <- t(mat_corr_ISRejVSRej)   # samples in rows
Y <- CountsISRejVsRej$Group        

set.seed(123)
plsda_model_ISRejVsRej <- mixOmics::plsda(X, Y, ncomp = 2)

# Scores for LV1 & LV2
scores <- plsda_model_ISRejVsRej$variates$X
plot_df <- data.frame(
  LV1   = scores[, 1],
  LV2   = scores[, 2],
  Group = Y
)


group_colors <- c(
  "Rejected"  = "#D62728",  # base red
  "Control Rejected"  = "#08306B"  # darker shade of blue (navy/steel blue)
)
group_shapes <- c(
  "Rejected"  = 16,
  "Control Rejected"   = 15
)


expl_var <- round(plsda_model_ISRejVsRej$prop_expl_var$X * 100, 1)  # % explained variance for X
xlab <- paste0("PLS Component 1 (", expl_var[1], "%)")
ylab <- paste0("PLS Component 2 (", expl_var[2], "%)")
unique(plot_df$Group)
# 2D PLS-DA plot (LV1 vs LV2) with filled ellipses
ggplot(plot_df, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  # points
  geom_point(size = 5, alpha = 0.9) +
  # filled confidence ellipses (70% CI both common)
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
    title = "PLS-DA: Anti-CD40L Rejection vs Control Rejection",
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
vip_scores_ISRejVsRej <- vip(plsda_model_ISRejVsRej)
# Extract only Component 1
vip_axis1_ISRejVsRej <- vip_scores_ISRejVsRej[, 1]
vip_df_axis1_ISRejVsRej <- data.frame(
  Gene = rownames(vip_scores_ISRejVsRej),
  VIP_Axis1 = vip_axis1_ISRejVsRej
)
#Sort by descending VIP
vip_df_axis1_ISRejVsRej <- vip_df_axis1_ISRejVsRej[order(-vip_df_axis1_ISRejVsRej$VIP_Axis1), ]
#Save to CSV
write.csv(vip_df_axis1_ISRejVsRej, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/AntiCD40LRejection_Rejection/VIP_scores_ISRejVsRej_PLS1.csv", row.names = FALSE)
# Rank genes by VIP (descending)
top100_genes_ISRejVsRej <- names(sort(vip_scores_ISRejVsRej[, 1], decreasing = TRUE))[1:100]
mat_top100_ISRejVsRej <- mat_corr_ISRejVSRej[top100_genes_ISRejVsRej, ]
anno <- data.frame(
  Group = Y,        # group labels (factor)
  row.names = colnames(mat_top100_ISRejVsRej)
)

anno_colors <- list(
  Group = group_colors   # same color scheme you used before
)

Y <- factor(Y, levels = c("Control Rejected","Rejected"))

# Reorder the columns of the matrix based on group
sample_order <- order(Y)
mat_top100_ordered_ISRejVsRej <- mat_top100_ISRejVsRej[, sample_order]
anno_ordered <- anno[sample_order, , drop = FALSE]
pheatmap(
  mat_top100_ordered_ISRejVsRej,
  scale = "row",                  # row-wise z-score
  annotation_col = anno_ordered,          # add group annotation
  annotation_colors = anno_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50) # publication-style
)

### Error Rate in PLSDA----
# Ensure Y is a factor (important)
Y <- factor(Y)

set.seed(123)
perf_plsda <- perf(plsda_model_ISRejVsRej,
                   validation = "Mfold",
                   folds = 5,
                   nrepeat = 10,          # increase if you want more stable estimates
                   dist = "max.dist",     # common choice; try "centroids.dist" too
                   progressBar = TRUE)

## 1) Overall classification error rate (per component, per repeat)
perf_plsda$error.rate$overall        # matrix: repeats x ncomp

## Mean overall error across repeats (per component)
colMeans(perf_plsda$error.rate$overall, na.rm = TRUE)

## 2) Balanced error rate (handles class imbalance)
perf_plsda$error.rate$BER            # matrix: repeats x ncomp
colMeans(perf_plsda$error.rate$BER, na.rm = TRUE)

##) 95% CI for error (quick & dirty via repeats)
er <- perf_plsda$error.rate$overall
mean_er <- colMeans(er, na.rm = TRUE)
se_er   <- apply(er, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
ci_low  <- mean_er - 1.96 * se_er
ci_high <- mean_er + 1.96 * se_er
data.frame(comp = seq_along(mean_er), mean_er, ci_low, ci_high)

## GSEA Analysis----

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
### Combined Timepoints----
# Paths to your saved results
ISREJVSREJ_path  <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/AntiCD40LRejection_Rejection/DESEQResults_Day7_14_Combined_AntiCD40L_Rejection_vs_Rejection.csv"

# Import
ISREJVSREJ_path  <- read.csv(ISREJVSREJ_path,  row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_ISREJVSREJ  <- ISREJVSREJ_path$stat;  names(lfc_vector_ISREJVSREJ)  <- rownames(ISREJVSREJ_path)

# Drop NAs
lfc_vector_ISREJVSREJ  <- lfc_vector_ISREJVSREJ[!is.na(lfc_vector_ISREJVSREJ)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_ISREJVSREJ  <- sort(lfc_vector_ISREJVSREJ,  decreasing = TRUE)

gsea_results_ISREJVSREJ <- GSEA(
  geneList      = lfc_vector_ISREJVSREJ,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  #keyType       = "SYMBOL",       # <- tell it explicitly
  TERM2GENE     = mm_all_df
)
gsea_results_ISREJVSREJ <- as.data.frame(gsea_results_ISREJVSREJ)


# Full results Combined Timepoins
write.csv(gsea_results_ISREJVSREJ,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/AntiCD40LRejection_Rejection/GSEAResults_AntiCD40LRejectionVsControlRejection_CombinedTimepoints.csv",
          row.names = FALSE)



# Your curated list EXACTLY as provided (de-dup just in case)
immune_master_ISREJVSREJ <- unique(c(
  #UPREGULATED
  # Innate myeloid / neutrophil
  "TRAVAGLINI_LUNG_NEUTROPHIL_CELL",
  "HAY_BONE_MARROW_NEUTROPHIL",
  "TRAVAGLINI_LUNG_CLASSICAL_MONOCYTE_CELL",
  "HE_LIM_SUN_FETAL_LUNG_C2_S100A12_HI_CLASSICAL_MONOCYTE",
  "SU_HO_CONV_CENT_CHONDROSARCOMA_LEUKOCYTE_C0_M1_MACROPHAGE",
  "TRAVAGLINI_LUNG_MACROPHAGE_CELL",
  "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CELL",
  # Inflammatory signaling
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  # NK/NKT
  "AIZARANI_LIVER_C18_NK_NKT_CELLS_5",
  # Injury / stress / thrombosis
  "HALLMARK_HYPOXIA",
  "HALLMARK_APOPTOSIS", "KEGG_APOPTOSIS", "HALLMARK_P53_PATHWAY",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
  "HE_LIM_SUN_FETAL_LUNG_C2_PLATELET_CELL",
  # Remodeling
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "KEGG_REGULATION_OF_AUTOPHAGY",
    #DOWNREGULATED
  # Proliferation / biosynthesis
  "HALLMARK_E2F_TARGETS","HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2","HALLMARK_G2M_CHECKPOINT",
  "KEGG_CELL_CYCLE","KEGG_DNA_REPLICATION","KEGG_RIBOSOME","KEGG_PROTEASOME",
  "KEGG_SPLICEOSOME","KEGG_RNA_DEGRADATION","KEGG_AMINOACYL_TRNA_BIOSYNTHESIS",
  "KEGG_PURINE_METABOLISM","KEGG_PYRIMIDINE_METABOLISM","HALLMARK_MTORC1_SIGNALING",
  # B cell / APC
  "HE_LIM_SUN_FETAL_LUNG_C5_PRO_B_CELL","HAY_BONE_MARROW_PRO_B",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  # Effector T/NK
  "TRAVAGLINI_LUNG_CD8_MEMORY_EFFECTOR_T_CELL","TRAVAGLINI_LUNG_CD4_MEMORY_EFFECTOR_T_CELL",
  "SU_HO_CONV_CENT_CHONDROSARCOMA_LEUKOCYTE_C2_T_CELL",
  # "TRAVAGLINI_LUNG_NATURAL_KILLER_CELL",     # <- include if consistent across contrasts
  # Metabolic support
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION","KEGG_OXIDATIVE_PHOSPHORYLATION","KEKK_CITRATE_CYCLE_TCA_CYCLE",
  # "HALLMARK_FATTY_ACID_METABOLISM",          # <- optional
  # IFNα / damage
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES","HALLMARK_COAGULATION"
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
ISREJVSREJ_tbl  <- coerce_gsea_tbl(gsea_results_ISREJVSREJ,  "Combined Timepoints")

# Keep ONLY immune_master pathways
ISREJVSREJ_tbl  <- ISREJVSREJ_tbl  %>% filter(gs_name %in% immune_master_ISREJVSREJ)

library(dplyr)
library(ggplot2)

plot_df <- ISREJVSREJ_tbl %>%
  mutate(
    padj = pmax(padj, .Machine$double.eps),
    logp = -log10(padj)
  ) %>%
  arrange((NES)) %>%                     # order by NES (descending)
  mutate(gs_name = factor(gs_name, levels = unique(gs_name)))

# Dotplot: x = NES, y = pathways (ordered by NES)
p_dot <- ggplot(plot_df, aes(x = NES, y = gs_name)) +
  geom_point(aes(size = logp, color = NES)) +
  scale_size_continuous(name = expression(-log[10]("padj")), range = c(2, 8)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "NES") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL,
       title = "GSEA (Anti-CD40L-Rejection vs Control Rejection)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )

print(p_dot)

## GSVA Analysis----

library(GSVA)
library(msigdbr)
library(limma)
library(BiocParallel)

# VST
vsd  <- vst(Celldds_IsletTransplant_IS_RejVsRej, blind = TRUE)
expr <- assay(vsd)
meta <- as.data.frame(colData(Celldds_IsletTransplant_IS_RejVsRej))

# Factors
meta$Group <- factor(meta$Group, levels = c("Control Rejected", "Rejected"))
if ("Batch" %in% names(meta)) meta$Batch <- factor(meta$Batch)

# --- Gene sets (same as yours) ---
c8   <- msigdbr(species = "Mus musculus", category = "C8")
hall <- msigdbr(species = "Mus musculus", category = "H")
kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG_LEGACY")
sets <- c(
  split(c8$gene_symbol,   c8$gs_name),
  split(hall$gene_symbol, hall$gs_name),
  split(kegg$gene_symbol, kegg$gs_name)
)
# Optional: restrict to your curated list
# sets <- sets[names(sets) %in% immune_master_ALLREJVSTOL]

# --- Select Day 7 + Day 14 samples ---
pick_day <- function(day) which(meta$Day %in% c(day, paste0("Day", day), paste0("Day ", day), paste0("D", day)))
idx <- sort(unique(c(pick_day(7), pick_day(14))))

expr_all <- expr[, idx, drop = FALSE]
meta_all <- meta[idx, , drop = FALSE]

# Normalize Day labels to two levels (Day7, Day14)
norm_day <- function(x) {
  x <- as.character(x)
  out <- ifelse(grepl("(^7$|Day ?7|^D7$)", x, ignore.case = TRUE), "Day7",
                ifelse(grepl("(^14$|Day ?14|^D14$)", x, ignore.case = TRUE), "Day14", NA))
  # Fallback for numeric
  out[is.na(out) & x %in% c("7")]  <- "Day7"
  out[is.na(out) & x %in% c("14")] <- "Day14"
  out
}
meta_all$Day <- factor(norm_day(meta_all$Day), levels = c("Day7","Day14"))

param <- MulticoreParam(workers = max(1, parallel::detectCores() - 1))
gsva_par <- gsvaParam(
  expr_all,
  sets,
  kcdf    = "Gaussian",   # VST is log-like
  minSize = 5,
  maxSize = 500
)
es_all <- gsva(gsva_par, BPPARAM = param)  # pathway x sample

# --- LIMMA on GSVA scores: Group main effect, adjusting for Day (+ Batch if present) ---
run_limma_combined <- function(es, meta) {
  has_day   <- "Day"   %in% names(meta) && nlevels(droplevels(meta$Day))   > 1
  has_batch <- "Batch" %in% names(meta) && nlevels(droplevels(meta$Batch)) > 1
  
  # Build design: ~ 0 + Group + (Day) + (Batch)
  terms <- c("Group", if (has_day) "Day", if (has_batch) "Batch")
  form  <- as.formula(paste("~ 0 +", paste(terms, collapse = " + ")))
  design <- model.matrix(form, data = meta)
  colnames(design) <- make.names(colnames(design))
  
  fit  <- lmFit(es, design)
  contr <- makeContrasts(Rejected_vs_Tolerance = GroupRejected - GroupControl.Rejected,
                         levels = colnames(design))
  fit2 <- eBayes(contrasts.fit(fit, contr))
  topTable(fit2, number = Inf, sort.by = "P")
}

GSVA_res_combined <- run_limma_combined(es_all, meta_all)

# --- Save outputs ---
out_dir <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/AntiCD40LRejection_Rejection"
write.csv(es_all,             file.path(out_dir, "GSVAScores_ISRejVsRej_Combined_D7_D14.csv"))
write.csv(GSVA_res_combined,  file.path(out_dir, "GSVAStats_ISRejVsRej_Combined_D7_D14.csv"))

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# --- Step 1: Select significant pathways
# Rename first column as "ID"
GSVA_res_combined$ID<-rownames(GSVA_res_combined)

sig_paths <- GSVA_res_combined %>%
  filter(P.Value < 0.05, abs(logFC) >= 0.2) %>%
  arrange(P.Value) %>%
  pull(ID) %>%
  unique()

# --- Step 2: Subset GSVA scores (es_all is pathways x samples)
mat <- es_all[sig_paths, , drop = FALSE]

# --- Step 3: Z-score per pathway for visualization
Z <- t(scale(t(mat)))
Z[is.na(Z)] <- 0
min(Z)
min_separators()# --- Step 4: Build sample annotations
meta_sub <- meta_all[colnames(Z), , drop = FALSE]

ann_cols <- list(
  Group = c("Control Rejected"= "#1f77b4", "Rejected" = "#d62728"),
  Day   = c(Day7 = "#66c2a5", Day14 = "#fc8d62")
)
top_anno <- HeatmapAnnotation(
  Group = meta_sub$Group,
  Day   = meta_sub$Day,
  col   = ann_cols
)

# --- Step 5: Heatmap
col_fun <- colorRamp2(c(min(Z), 0,max(Z)), c("#2b6cb0", "white", "#b91c1c"))

ht <- Heatmap(
  Z, name = "GSVA Z-score",
  col = col_fun,
  top_annotation = top_anno,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = meta_sub$Group,   # split by Group for clarity
  heatmap_legend_param = list(
    title = "GSVA Z-score"
  )
)


draw(
  ht,
  heatmap_legend_side = "left",
  annotation_legend_side = "left",
  padding = unit(c(5, 20, 5, 70), "mm")  # top, right, bottom, left
)

# 6. Allo Tolerance Vs Syngeneic-----

library(biomaRt)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "start_position", "end_position"),
  mart = mart
)

gene_info$gene_length <- gene_info$end_position - gene_info$start_position

sel <- colData(dds_IsletTransplant)$Group %in% c("Control Accepted","Tolerance")  & colData(dds_IsletTransplant)$Day %in% c(7,14,28,42,56,70)
#IS- Immunosuppressed
dds_IsletTransplant_TolVsSyn <- dds_IsletTransplant[, sel]

dds_IsletTransplant_TolVsSyn$Batch       <- factor(dds_IsletTransplant_TolVsSyn$Batch)
dds_IsletTransplant_TolVsSyn$LibraryPrep <- factor(dds_IsletTransplant_TolVsSyn$LibraryPrep)
dds_IsletTransplant_TolVsSyn$Group <- factor(
  dds_IsletTransplant_TolVsSyn$Group,
  levels = c("Control Accepted", "Tolerance")  # order sets baseline
)
levels(dds_IsletTransplant_TolVsSyn$Group)

dds_IsletTransplant_TolVsSyn$Day <- factor(dds_IsletTransplant_TolVsSyn$Day,
                                                  levels = c(7, 14, 28, 42, 56, 70))  # baseline = 7
levels(dds_IsletTransplant_TolVsSyn$Day)

cd <- as.data.frame(colData(dds_IsletTransplant_TolVsSyn))
mm <- model.matrix(~ Batch + LibraryPrep + Day + Group, data = cd)
qr_mm <- qr(mm)

kept    <- colnames(mm)[qr_mm$pivot[seq_len(qr_mm$rank)]]
dropped <- if (qr_mm$rank < ncol(mm)) colnames(mm)[qr_mm$pivot[(qr_mm$rank+1):ncol(mm)]] else character()

kept
dropped


summary(res_TOLvSYN)

## DESeq-Combined Timepoints----
#Design Formula
design(dds_IsletTransplant_TolVsSyn) <- ~ Batch + Day + Group 
dds_IsletTransplant_TolVsSyn <- DESeq(dds_IsletTransplant_TolVsSyn)
design(dds_IsletTransplant_TolVsSyn)
resultsNames(dds_IsletTransplant_TolVsSyn)

# All day combined
res_TOLvSYN <- results(dds_IsletTransplant_TolVsSyn,
                       name = "Group_Tolerance_vs_Control.Accepted")



library(EnhancedVolcano)

# Thresholds
q_cut  <- 0.10
fc_cut <- 1

library(EnhancedVolcano)

## 0) Prep results as data.frames with a 'gene' column
ALL_TOLVSSYN  <- as.data.frame(res_TOLvSYN);  ALL_TOLVSSYN$gene  <- rownames(res_TOLvSYN)

write.csv(ALL_TOLVSSYN, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Days_Combined_Tolerance_vs_Syngeneic.csv", row.names = TRUE)

# Function to create color mapping
make_keyvals_fdr_fc_TolvsSyn <- function(df, q = 0.10, fc = 1,
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
  
  label[up]   <- paste0("Tolerance-Upregulated")
  label[down] <- paste0("Syngeneic-Upregulated")
  
  names(key) <- label        # <- legend labels; no NAs
  key
}

keyvals_ALL  <- make_keyvals_fdr_fc_TolvsSyn(ALL_TOLVSSYN)

library(dplyr)

selLab_ALL_TOLVSSYN  <- pick_labels(ALL_TOLVSSYN,  q_cut, fc_cut, 50)


## 3) Axis limits
xmax_ALL  <- max(2, ceiling(max(abs(ALL_TOLVSSYN$log2FoldChange),  na.rm=TRUE)))

EnhancedVolcano(
  ALL_TOLVSSYN,
  lab           = ALL_TOLVSSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Tolerance vs Syngeneic — Combined Timepoints",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(ALL_TOLVSSYN$padj<0.10 & abs(ALL_TOLVSSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_ALL, xmax_ALL),
  #ylim = c(0,10),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 50,
  colCustom     = keyvals_ALL,
  legendPosition= "right",
  selectLab     = selLab_ALL_TOLVSSYN
)

summary(res_TOLvSYN)

## DESeq-Individual Timepoints----
#Design Formula
design(dds_IsletTransplant_TolVsSyn) <- ~ Batch + Day + Group + Day:Group
dds_IsletTransplant_TolVsSyn <- DESeq(dds_IsletTransplant_TolVsSyn)
design(dds_IsletTransplant_TolVsSyn)
resultsNames(dds_IsletTransplant_TolVsSyn)

# All day combined
res_TOLvSYN_d7 <- results(dds_IsletTransplant_TolVsSyn,
                       name = "Group_Tolerance_vs_Control.Accepted")
# Other days: main group effect + interaction term for that day
res_TOLvSYN_d14 <- results(dds_IsletTransplant_TolVsSyn, contrast = list(c("Group_Tolerance_vs_Control.Accepted","Day14.GroupTolerance")))
res_TOLvSYN_d28 <- results(dds_IsletTransplant_TolVsSyn, contrast = list(c("Group_Tolerance_vs_Control.Accepted","Day28.GroupTolerance")))
res_TOLvSYN_d42 <- results(dds_IsletTransplant_TolVsSyn, contrast = list(c("Group_Tolerance_vs_Control.Accepted","Day42.GroupTolerance")))
res_TOLvSYN_d56 <- results(dds_IsletTransplant_TolVsSyn, contrast = list(c("Group_Tolerance_vs_Control.Accepted","Day56.GroupTolerance")))
res_TOLvSYN_d70 <- results(dds_IsletTransplant_TolVsSyn, contrast = list(c("Group_Tolerance_vs_Control.Accepted","Day70.GroupTolerance")))

## 0) Prep results as data.frames with a 'gene' column
d7_TOLVSYN  <- as.data.frame(res_TOLvSYN_d7);  d7_TOLVSYN$gene  <- rownames(res_TOLvSYN_d7)
d14_TOLVSYN <- as.data.frame(res_TOLvSYN_d14); d14_TOLVSYN$gene <- rownames(res_TOLvSYN_d14)
d28_TOLVSYN <- as.data.frame(res_TOLvSYN_d28); d28_TOLVSYN$gene <- rownames(res_TOLvSYN_d28)
d42_TOLVSYN <- as.data.frame(res_TOLvSYN_d42); d42_TOLVSYN$gene <- rownames(res_TOLvSYN_d42)
d56_TOLVSYN <- as.data.frame(res_TOLvSYN_d56); d56_TOLVSYN$gene <- rownames(res_TOLvSYN_d56)
d70_TOLVSYN <- as.data.frame(res_TOLvSYN_d70); d70_TOLVSYN$gene <- rownames(res_TOLvSYN_d70)

write.csv(d7_TOLVSYN, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day7_ToleranceVsSyngeneic.csv", row.names = TRUE)
write.csv(d14_TOLVSYN, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day14_ToleranceVsSyngeneic.csv", row.names = TRUE)
write.csv(d28_TOLVSYN, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day28_ToleranceVsSyngeneic.csv", row.names = TRUE)
write.csv(d42_TOLVSYN, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day42_ToleranceVsSyngeneic.csv", row.names = TRUE)
write.csv(d56_TOLVSYN, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day56_ToleranceVsSyngeneic.csv", row.names = TRUE)
write.csv(d70_TOLVSYN, "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day70_ToleranceVsSyngeneic.csv", row.names = TRUE)


keyvals_d7  <- make_keyvals_fdr_fc_TolvsSyn(d7_TOLVSYN)
keyvals_d14 <- make_keyvals_fdr_fc_TolvsSyn(d14_TOLVSYN)
keyvals_d28 <- make_keyvals_fdr_fc_TolvsSyn(d28_TOLVSYN)
keyvals_d42 <- make_keyvals_fdr_fc_TolvsSyn(d42_TOLVSYN)
keyvals_d56 <- make_keyvals_fdr_fc_TolvsSyn(d56_TOLVSYN)
keyvals_d70 <- make_keyvals_fdr_fc_TolvsSyn(d70_TOLVSYN)

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

selLab_d7_TOLVSSYN  <- pick_labels(d7_TOLVSYN,  q_cut, fc_cut, 40)
selLab_d14_TOLVSSYN <- pick_labels(d14_TOLVSYN, q_cut, fc_cut, 40)
selLab_d28_TOLVSSYN <- pick_labels(d28_TOLVSYN, q_cut, fc_cut, 40)
selLab_d42_TOLVSSYN <- pick_labels(d42_TOLVSYN, q_cut, fc_cut, 40)
selLab_d56_TOLVSSYN <- pick_labels(d56_TOLVSYN, q_cut, fc_cut, 40)
selLab_d70_TOLVSSYN <- pick_labels(d70_TOLVSYN, q_cut, fc_cut, 40)


## 3) Axis limits
xmax_d7  <- max(2, ceiling(max(abs(d7_TOLVSYN$log2FoldChange),  na.rm=TRUE)))
xmax_d14 <- max(2, ceiling(max(abs(d14_TOLVSYN$log2FoldChange), na.rm=TRUE)))
xmax_d28 <- max(2, ceiling(max(abs(d28_TOLVSYN$log2FoldChange), na.rm=TRUE)))
xmax_d42 <- max(2, ceiling(max(abs(d42_TOLVSYN$log2FoldChange), na.rm=TRUE)))
xmax_d56 <- max(2, ceiling(max(abs(d56_TOLVSYN$log2FoldChange), na.rm=TRUE)))
xmax_d70 <- max(2, ceiling(max(abs(d70_TOLVSYN$log2FoldChange), na.rm=TRUE)))


## Volcano: Day 7----
EnhancedVolcano(
  d7_TOLVSYN,
  lab           = d7_TOLVSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Tolerance vs Syngeneic — Day 7",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d7_TOLVSYN$padj<0.10 & abs(d7_TOLVSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d7, xmax_d7),
  ylim = c(0,14),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 35,
  colCustom     = keyvals_d7,
  legendPosition= "right",
  selectLab     = selLab_d7_TOLVSSYN
)

## Volcano: Day 14----
EnhancedVolcano(
  d14_TOLVSYN,
  lab           = d14_TOLVSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Tolerance vs Syngeneic — Day 14",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d14_TOLVSYN$padj<0.10 & abs(d14_TOLVSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d14, xmax_d14),
  ylim = c(0,13),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 35,
  colCustom     = keyvals_d14,
  legendPosition= "right",
  selectLab     = selLab_d14_TOLVSSYN
)


## Volcano: Day 28----
EnhancedVolcano(
  d28_TOLVSYN,
  lab           = d28_TOLVSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Tolerance vs Syngeneic — Day 28",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d28_TOLVSYN$padj<0.10 & abs(d28_TOLVSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d28, xmax_d28),
  #ylim = c(0,6),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 35,
  colCustom     = keyvals_d28,
  legendPosition= "right",
  selectLab     = selLab_d28_TOLVSSYN
)


## Volcano: Day 42----
EnhancedVolcano(
  d42_TOLVSYN,
  lab           = d42_TOLVSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Tolerance vs Syngeneic — Day 42",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d42_TOLVSYN$padj<0.10 & abs(d42_TOLVSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d42, xmax_d42),
  ylim = c(0,15),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 35,
  colCustom     = keyvals_d42,
  legendPosition= "right",
  selectLab     = selLab_d42_TOLVSSYN
)

## Volcano: Day 56----
EnhancedVolcano(
  d56_TOLVSYN,
  lab           = d56_TOLVSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Tolerance vs Syngeneic — Day 56",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d56_TOLVSYN$padj<0.10 & abs(d56_TOLVSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d56, xmax_d56),
  #ylim = c(0,6),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 35,
  colCustom     = keyvals_d56,
  legendPosition= "right",
  selectLab     = selLab_d56_TOLVSSYN
)

## Volcano: Day 70----
EnhancedVolcano(
  d70_TOLVSYN,
  lab           = d70_TOLVSYN$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Tolerance vs Syngeneic — Day 70",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d70_TOLVSYN$padj<0.10 & abs(d70_TOLVSYN$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d70, xmax_d70),
  #ylim = c(0,6),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 35,
  colCustom     = keyvals_d70,
  legendPosition= "right",
  selectLab     = selLab_d70_TOLVSSYN
)


## GSEA Analysis----
# Paths to your saved results
d7_path_TolVsSyn  <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day7_ToleranceVsSyngeneic.csv"
d14_path_TolVsSyn <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day14_ToleranceVsSyngeneic.csv"
d28_path_TolVsSyn <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day28_ToleranceVsSyngeneic.csv"
d42_path_TolVsSyn <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day42_ToleranceVsSyngeneic.csv"
d56_path_TolVsSyn <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day56_ToleranceVsSyngeneic.csv"
d70_path_TolVsSyn <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/DESEQResults_Day70_ToleranceVsSyngeneic.csv"


# Import
d7_TolVsSyn  <- read.csv(d7_path_TolVsSyn,  row.names = 1)
d14_TolVsSyn <- read.csv(d14_path_TolVsSyn, row.names = 1)
d28_TolVsSyn <- read.csv(d28_path_TolVsSyn, row.names = 1)
d42_TolVsSyn <- read.csv(d42_path_TolVsSyn, row.names = 1)
d56_TolVsSyn <- read.csv(d56_path_TolVsSyn, row.names = 1)
d70_TolVsSyn <- read.csv(d70_path_TolVsSyn, row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_d7_TolVsSyn   <- d7_TolVsSyn$stat;  names(lfc_vector_d7_TolVsSyn )  <- rownames(d7_TolVsSyn)
lfc_vector_d14_TolVsSyn  <- d14_TolVsSyn$stat; names(lfc_vector_d14_TolVsSyn ) <- rownames(d14_TolVsSyn)
lfc_vector_d28_TolVsSyn  <- d28_TolVsSyn$stat; names(lfc_vector_d28_TolVsSyn ) <- rownames(d28_TolVsSyn)
lfc_vector_d42_TolVsSyn  <- d42_TolVsSyn$stat; names(lfc_vector_d42_TolVsSyn ) <- rownames(d42_TolVsSyn)
lfc_vector_d56_TolVsSyn  <- d56_TolVsSyn$stat; names(lfc_vector_d56_TolVsSyn ) <- rownames(d56_TolVsSyn)
lfc_vector_d70_TolVsSyn  <- d70_TolVsSyn$stat; names(lfc_vector_d70_TolVsSyn ) <- rownames(d70_TolVsSyn)

# Drop NAs
lfc_vector_d7_TolVsSyn  <- lfc_vector_d7_TolVsSyn[!is.na(lfc_vector_d7_TolVsSyn)]
lfc_vector_d14_TolVsSyn <- lfc_vector_d14_TolVsSyn[!is.na(lfc_vector_d14_TolVsSyn)]
lfc_vector_d28_TolVsSyn <- lfc_vector_d28_TolVsSyn[!is.na(lfc_vector_d28_TolVsSyn)]
lfc_vector_d42_TolVsSyn <- lfc_vector_d42_TolVsSyn[!is.na(lfc_vector_d42_TolVsSyn)]
lfc_vector_d56_TolVsSyn <- lfc_vector_d56_TolVsSyn[!is.na(lfc_vector_d56_TolVsSyn)]
lfc_vector_d70_TolVsSyn <- lfc_vector_d70_TolVsSyn[!is.na(lfc_vector_d70_TolVsSyn)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_d7_TolVsSyn  <- sort(lfc_vector_d7_TolVsSyn,  decreasing = TRUE)
lfc_vector_d14_TolVsSyn <- sort(lfc_vector_d14_TolVsSyn, decreasing = TRUE)
lfc_vector_d28_TolVsSyn <- sort(lfc_vector_d28_TolVsSyn, decreasing = TRUE)
lfc_vector_d42_TolVsSyn <- sort(lfc_vector_d42_TolVsSyn, decreasing = TRUE)
lfc_vector_d56_TolVsSyn <- sort(lfc_vector_d56_TolVsSyn, decreasing = TRUE)
lfc_vector_d70_TolVsSyn <- sort(lfc_vector_d70_TolVsSyn, decreasing = TRUE)

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
gsea_results_d7_TolVsSyn <- GSEA(
  geneList      = lfc_vector_d7_TolVsSyn,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  #keyType       = "SYMBOL",       # <- tell it explicitly
  TERM2GENE     = mm_all_df
)
gsea_results_d7_TolVsSyn_df <- as.data.frame(gsea_results_d7_TolVsSyn)

# Day 14
gsea_results_d14_TolVsSyn <- GSEA(
  geneList      = lfc_vector_d14_TolVsSyn,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_d14_TolVsSyn_df <- as.data.frame(gsea_results_d14_TolVsSyn)


# Day 28
gsea_results_d28_TolVsSyn <- GSEA(
  geneList      = lfc_vector_d28_TolVsSyn,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_d28_TolVsSyn_df <- as.data.frame(gsea_results_d28_TolVsSyn)

# Day 42
gsea_results_d42_TolVsSyn <- GSEA(
  geneList      = lfc_vector_d42_TolVsSyn,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_d42_TolVsSyn_df <- as.data.frame(gsea_results_d42_TolVsSyn)

# Day 56
gsea_results_d56_TolVsSyn <- GSEA(
  geneList      = lfc_vector_d56_TolVsSyn,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_d56_TolVsSyn_df <- as.data.frame(gsea_results_d56_TolVsSyn)

# Day 70
gsea_results_d70_TolVsSyn <- GSEA(
  geneList      = lfc_vector_d70_TolVsSyn,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_d70_TolVsSyn_df <- as.data.frame(gsea_results_d70_TolVsSyn)


# Full results Day 7
write.csv(gsea_results_d7_TolVsSyn_df,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/Allogeneic_Vs_Syngeneic/GSEAResults_TolVsSyn_Day7.csv",
          row.names = FALSE)

# Full results Day 14
write.csv(gsea_results_d14_TolVsSyn_df,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/GSEAResults_TolVsSyn_Day14.csv",
          row.names = FALSE)

# Full results Day 28
write.csv(gsea_results_d28_TolVsSyn_df,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/GSEAResults_TolVsSyn_Day28.csv",
          row.names = FALSE)

# Full results Day 42
write.csv(gsea_results_d42_TolVsSyn_df,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/GSEAResults_TolVsSyn_Day42.csv",
          row.names = FALSE)


# Full results Day 56
write.csv(gsea_results_d56_TolVsSyn_df,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/GSEAResults_TolVsSyn_Day56.csv",
          row.names = FALSE)

# Full results Day 70
write.csv(gsea_results_d70_TolVsSyn_df,
          "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/GSEAResults_TolVsSyn_Day70.csv",
          row.names = FALSE)


library(dplyr)
library(stringr)
library(ggplot2)

# Your curated list EXACTLY as provided (de-dup just in case)
immune_master <- unique(c(
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "KEGG_GRAFT_VERSUS_HOST_DISEASE",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "CUI_DEVELOPING_HEART_C8_MACROPHAGE",
  "HALLMARK_HYPOXIA",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "TRAVAGLINI_LUNG_CLASSICAL_MONOCYTE_CELL","TRAVAGLINI_LUNG_NONCLASSICAL_MONOCYTE_CELL",
  "TRAVAGLINI_LUNG_NEUTROPHIL_CELL","TRAVAGLINI_LUNG_PROLIFERATING_NK_T_CELL",
  "HE_LIM_SUN_FETAL_LUNG_C4_CD8_T_CELL",
  "HALLMARK_E2F_TARGETS",
  "HE_LIM_SUN_FETAL_LUNG_C5_PRO_B_CELL"
  
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
d7_tbl  <- coerce_gsea_tbl(gsea_results_d7_TolVsSyn_df,  "Day 7")
d14_tbl <- coerce_gsea_tbl(gsea_results_d14_TolVsSyn_df, "Day 14")
d28_tbl <- coerce_gsea_tbl(gsea_results_d28_TolVsSyn_df, "Day 28")
d42_tbl <- coerce_gsea_tbl(gsea_results_d42_TolVsSyn_df, "Day 42")
d56_tbl <- coerce_gsea_tbl(gsea_results_d56_TolVsSyn_df, "Day 56")
d70_tbl <- coerce_gsea_tbl(gsea_results_d70_TolVsSyn_df, "Day 70")

# Keep ONLY immune_master pathways
d7_sel  <- d7_tbl  %>% filter(gs_name %in% immune_master)
d14_sel <- d14_tbl %>% filter(gs_name %in% immune_master)
d28_sel <- d28_tbl %>% filter(gs_name %in% immune_master)
d42_sel <- d42_tbl %>% filter(gs_name %in% immune_master)
d56_sel <- d56_tbl %>% filter(gs_name %in% immune_master)
d70_sel <- d70_tbl %>% filter(gs_name %in% immune_master)

# Combine and ensure both days present for every pathway
plot_df <- bind_rows(d7_sel, d14_sel, d28_sel, d42_sel, d56_sel, d70_sel) %>%
  mutate(Day = factor(Day, levels = c("Day 7","Day 14", "Day 28","Day 42","Day 56","Day 70"))) %>%
  # create full grid of (gs_name x Day) and fill missing with neutral values
  complete(gs_name = immune_master, Day,
           fill = list(NES = 0, padj = 1)) %>%
  mutate(logp = -log10(padj))

# Order rows nicely (keep your immune_master order or order by category if you prefer)
plot_df$gs_name <- factor(plot_df$gs_name, levels = rev(immune_master))

# Plot: rows = pathways, columns = Day; size = −log10(padj); color = NES
p <- ggplot(plot_df, aes(x = Day, y = gs_name)) +
  geom_point(aes(size = logp, color = NES)) +
  scale_size_continuous(name = "−log10(padj)", range = c(0, 10)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  labs(x = "", y = "", title = "Enriched Pathways") +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  )
print(p)

library(dplyr)
library(purrr)
library(tidyr)

# helper: pick NES + FDR for one day and rename columns to include the day
extract_NES_FDR <- function(df, day_label){
  df <- as.data.frame(df)
  if (!"Description" %in% names(df)) df$Description <- df$ID
  # prefer qvalues (FDR); else fall back to BH-adjusted p (p.adjust)
  fdr_col <- if ("qvalues" %in% names(df)) "qvalues" else "p.adjust"
  df %>%
    distinct(ID, Description, .keep_all = TRUE) %>%   # safety: remove any duplicates
    transmute(
      ID, Description,
      !!paste0("NES_", day_label) := NES,
      !!paste0("FDR_", day_label) := .data[[fdr_col]]
    )
}

# build per-day slim tables
d7  <- extract_NES_FDR(gsea_results_d7_TolVsSyn_df,  "D7")
d14 <- extract_NES_FDR(gsea_results_d14_TolVsSyn_df, "D14")
d28 <- extract_NES_FDR(gsea_results_d28_TolVsSyn_df, "D28")
d42 <- extract_NES_FDR(gsea_results_d42_TolVsSyn_df, "D42")
d56 <- extract_NES_FDR(gsea_results_d56_TolVsSyn_df, "D56")
d70 <- extract_NES_FDR(gsea_results_d70_TolVsSyn_df, "D70")
library(purrr)

# full outer join on ID + Description to align rows
gsea_wide <- purrr::reduce(list(d7, d14, d28, d42, d56, d70),
                    full_join, by = c("ID","Description"))

# optional: order rows by best (min) FDR across days, then by max |NES|
fdr_cols <- grep("^FDR_", names(gsea_wide), value = TRUE)
nes_cols <- grep("^NES_", names(gsea_wide), value = TRUE)
gsea_wide <- gsea_wide %>%
  mutate(
    minFDR = do.call(pmin, c(across(all_of(fdr_cols)), list(na.rm = TRUE))),
    maxAbsNES = do.call(pmax, c(across(all_of(nes_cols), ~abs(.x)), list(na.rm = TRUE)))
  ) %>%
  arrange(minFDR, desc(maxAbsNES)) %>%
  select(-minFDR, -maxAbsNES)

# write combined table
out_file <- "/Users/jyotirmoyroy/Desktop/IsletTransplantRejection/ToleranceVsSyngeneic/GSEA_NES_FDR_AllDaysCombined.csv"
write.csv(gsea_wide, out_file, row.names = FALSE)

library(dplyr)
library(tidyr)

# --- choose the timepoints you want to require significance for ---
days5 <- c("D14","D28","D42","D56","D70")
fdr_cols_5 <- paste0("FDR_", days5)

# safety: ensure those columns exist
stopifnot(all(fdr_cols_5 %in% names(gsea_wide)))

# 1) Significant across ALL 5 chosen timepoints (FDR < 0.1 at each)
sig_all5 <- gsea_wide %>%
  filter(if_all(all_of(fdr_cols_5), ~ !is.na(.x) & .x < 0.10)) %>%
  arrange(Description)

# Inspect a compact list
sig_all5_paths <- sig_all5 %>% select(ID, Description)
sig_all5_paths

library(dplyr)

#  define windows (edit these if you want different bins)
early_days <- c("D7","D14")
late_days  <- c("D42","D56","D70")

early_fdr <- paste0("FDR_", early_days)
late_fdr  <- paste0("FDR_", late_days)
early_nes <- paste0("NES_", early_days)
late_nes  <- paste0("NES_", late_days)

# Safety checks
stopifnot(all(c(early_fdr, late_fdr, early_nes, late_nes) %in% names(gsea_wide)))

base_sig <- gsea_wide %>%
  filter(if_all(all_of(c(early_fdr, late_fdr)), ~ !is.na(.x) & .x < 0.10)) %>%
  mutate(
    early_mean_NES = rowMeans(across(all_of(early_nes)), na.rm = TRUE),
    late_mean_NES  = rowMeans(across(all_of(late_nes)),  na.rm = TRUE)
  )


#7. Tolerance Group- Sublcustering----

# Core
library(DESeq2)
library(limma)
library(tidyverse)

# Embeddings & kNN graph clustering
library(RANN)      # fast kNN
library(igraph)    # Louvain clustering
library(uwot)      # UMAP
library(irlba)     # fast PCA

# Optional visualization helpers
library(ggplot2)
library(ggrepel)

hvg_select <- function(mat, n = 2000) {
  # mat: genes x samples (e.g., VST)
  vars <- matrixStats::rowVars(mat)
  ord  <- order(vars, decreasing = TRUE)
  mat[ord[seq_len(min(n, nrow(mat)))], , drop = FALSE]
}

build_knn_graph <- function(X, k = 10) {
  # X: samples x features (e.g., PCs)
  nn  <- RANN::nn2(X, k = k + 1)                     # +1 because 1-NN is itself
  idx <- nn$nn.idx[, -1, drop = FALSE]               # drop self
  # edges: each sample i connects to its k neighbors
  edges <- cbind(rep(seq_len(nrow(idx)), times = ncol(idx)),
                 as.vector(idx))
  g <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
  igraph::simplify(g)
}


louvain_clusters <- function(g) {
  igraph::cluster_louvain(g)$membership |> factor()
}

sel <- colData(dds_IsletTransplant)$Group %in% c("Tolerance")  #& colData(dds_IsletTransplant)$Day %in% c(7,14,28,42,56,70)
#IS- Immunosuppressed
dds_tol <- dds_IsletTransplant[, sel]
vsd  <- vst(dds_tol, blind = FALSE)
mat  <- assay(vsd)
colData(mat)$Batch

# --- Remove batch (for visualization only) ---
mat_bc <- removeBatchEffect(mat, batch = colData(vsd)$Batch)

# mat_bc: genes x samples
mat_hvg <- hvg_select(mat_bc, n = 2000)      # genes x samples (G x N)

# Make rows = samples, cols = genes
X  <- t(mat_hvg)                              # samples x genes (N x G)
Xs <- scale(X, center = TRUE, scale = FALSE)  # center genes across samples

# PCA on samples x genes
pc <- irlba::prcomp_irlba(Xs, n = 30, center = FALSE, scale. = FALSE)
PC <- pc$x                                    # samples x PCs

# Sanity check
stopifnot(nrow(PC) == ncol(mat_bc))           # should be 80 here

# UMAP + clustering on PC space (samples)
set.seed(123)
UMAP <- uwot::umap(PC, n_neighbors = 15, min_dist = 0.2)

g   <- build_knn_graph(PC, k = 10)            # PC has samples in rows
clu <- factor(igraph::cluster_louvain(g)$membership)

# Now lengths match
meta <- as.data.frame(colData(vsd))
stopifnot(nrow(meta) == length(clu))
meta$Cluster <- clu

# assuming UMAP is an N x 2 matrix (samples x coords)
df_plot <- cbind(meta, as.data.frame(UMAP)) |>
  dplyr::rename(UMAP1 = V1, UMAP2 = V2)

# Plot by cluster
library(ggplot2)
ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw() +
  ggtitle("Tolerance samples clustered")

# Plot by timepoint
ggplot(df_plot, aes(x = UMAP1, y = UMAP2, color = Day)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw() +
  ggtitle("Colored by Day")
