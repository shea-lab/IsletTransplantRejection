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
setwd("C:/Users/17343/Desktop/IsletTransplantRejection/Data/")
getwd()

#Metadata Importing
# had issues with read.table so changed to read.csv
meta_batch <- read.csv("C:/Users/17343/Desktop/IsletTransplantRejection/Data/NOD_Data/NODTransplant_Metadata.csv", 
                       header = TRUE, 
                       check.names = FALSE)

meta_batch <- as.data.frame(meta_batch)

# Preview the combined metadata
head(meta_batch)
unique(meta_batch$Group)

#Counts Data Importing
counts_batch <- as.data.frame(read.table(
  "C:/Users/17343/Desktop/IsletTransplantRejection/Data/NOD_Data/NODTransplant_gene_counts_annot.csv"
  , sep=",", header=T,check.names = FALSE))
counts_batch <- na.omit(counts_batch)

#Remove duplicate names
counts_batch <- counts_batch[!duplicated(counts_batch[, 1]), ]
genes <- counts_batch[, 1]
rownames(counts_batch) <- genes
counts_batch <- counts_batch[, -1]

# check to make sure row names are assigned

# Identify the samples to keep (not "Technical Rejection")
# samples_to_keep <- meta_batch$Sample[meta_batch$Group != "Technical Rejection"]
# meta_batch <- meta_batch[meta_batch$Group != "Technical Rejection", ]
# counts_batch <- counts_batch[, colnames(counts_batch) %in% samples_to_keep]
# Not Needed bc no Technical Rejection

# Preview the combined dataset
head(counts_batch)
# Create new column with high-level analysis groups
meta_batch$Outcomes <- meta_batch$Group


# 2.Preprocessing and Cleaning ----
getwd()
setwd("C:/Users/17343/Desktop/IsletTransplantRejection/Code")

# Select columns in combined_counts that match the remaining sample names in meta_combined
counts_batch <- counts_batch[, meta_batch$Samples]  # Ensure Sample_IDs match column names in combined_counts

source("AllFunctions.R") #Added
IsletTransplantCounts <- flexiDEG.function1(counts_batch, meta_batch, # Genes in rows 
                                            convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                            batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0, 0

rows_to_remove <- grep("^Gm[0-9]", rownames(IsletTransplantCounts))

# Remove those rows from case1_f1
IsletTransplantCounts <- IsletTransplantCounts[-rows_to_remove, ]

# Removing any data from beyond day 35 (because Early Rejections do not have data beyond this day)
meta_batch <- meta_batch %>% filter(Day <= 35 & Day != 0)
IsletTransplantCounts <- IsletTransplantCounts[, colnames(IsletTransplantCounts) %in% meta_batch$Samples]

# Saving case1_f1 dataframe as a CSV file
write.csv(IsletTransplantCounts, file = "C:/Users/17343/Desktop/IsletTransplantRejection/Data/NOD/NOD_Raw_Counts_ITx_Filtered.csv", row.names = TRUE)

# Saving meta_combined dataframe as a CSV file
write.csv(meta_batch, file = "C:/Users/17343/Desktop/IsletTransplantRejection/Data/NOD/NOD_Metadata_ITx.csv", row.names = FALSE)

# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_batch)
unique_colSide <- unique(colSide)


# 3.Create Univeral DESqEQ Object ----

IsletTransplantCounts <- as.matrix(IsletTransplantCounts)
storage.mode(IsletTransplantCounts) <- "integer"

dds_IsletTransplant <- DESeqDataSetFromMatrix(IsletTransplantCounts, meta_batch,
                                              design = ~ 1)   # dummy design for now

saveRDS(dds_IsletTransplant, file = "C:/Users/17343/Desktop/IsletTransplantRejection/Data/Robjects/dds_NOD_master.rds")


# 4. Early Rejection Vs Late Rejection Signature----
# subset samples
sel <- colData(dds_IsletTransplant)$Group %in% c("Early Rejection","Late Rejection") # Did not filter for days specific days
dds_IsletTransplant_EarlyVsLate <- dds_IsletTransplant[, sel]

dds_IsletTransplant_EarlyVsLate$Batch       <- factor(dds_IsletTransplant_EarlyVsLate$Batch)
dds_IsletTransplant_EarlyVsLate$LibraryPrep <- factor(dds_IsletTransplant_EarlyVsLate$LibraryPrep)
dds_IsletTransplant_EarlyVsLate$Group <- factor(
  dds_IsletTransplant_EarlyVsLate$Group,
  levels = c("Early Rejection","Late Rejection")  # order sets baseline
)
levels(dds_IsletTransplant_EarlyVsLate$Group)

dds_IsletTransplant_EarlyVsLate$Treatment <- factor(dds_IsletTransplant_EarlyVsLate$Treatment,
                                                  levels = c("Low Dose", "High Dose"))
levels(dds_IsletTransplant_EarlyVsLate$Treatment)

# Add log(IEQ) column
colData(dds_IsletTransplant_EarlyVsLate)$logIEQ <- log(colData(dds_IsletTransplant_EarlyVsLate)$IEQ)

# Check for Collinearuty
cd <- as.data.frame(colData(dds_IsletTransplant_EarlyVsLate))
# Basic sanity
lapply(cd[, c("Batch","LibraryPrep","Day","Group")], function(x) table(x, useNA="ifany"))
# Check for NAs
sapply(cd[, c("Batch","LibraryPrep","Day","Group")], function(x) any(is.na(x)))
# Model matrix rank- Included Group, Treatment, and logIEQ
mm <- model.matrix(~ Group + Treatment + logIEQ, data = cd)
qr(mm)$rank; ncol(mm)             # if rank < ncol(mm), not full rank

# Inital testing if logIEq should be included as covariate
# t test
t.test(logIEQ ~ Group, data = cd) # p = 0.4405 -> indicating it should be covariate

# plot
ggplot(cd, aes(x = Group, y = logIEQ, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2) +  # Shows individual sample points
  theme_minimal() +
  labs(title = "Distribution of IEQ by Rejection Group",
       y = "log10(IEQ)", x = "Group")

# Bucketing the days: 7,14,28 (note: day 35 is being bucketed into 28)
bucketed_cd <- cd %>% mutate(Day = ifelse(Day >= 22, 28, Day))
table(bucketed_cd$Group,bucketed_cd$Day)
#           0   7  14  28
# Early     8   5   6  10
# Late      6   5   5   5
keep <- rowSums(counts(dds_IsletTransplant_EarlyVsLate) >= 10) >= 5  # Smallest number of samples in a group
dds_IsletTransplant_EarlyVsLate <- dds_IsletTransplant_EarlyVsLate[keep, ]


q_cut  <- 0.10
fc_cut <- 1

library(EnhancedVolcano)
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

# Determining if Treatment should be a covariate in the design formula ----
# Design Formula NOT Accounting for Treatment
design(dds_IsletTransplant_EarlyVsLate) <- ~ Day + Group
dds_IsletTransplant_EarlyVsLate <- DESeq(dds_IsletTransplant_EarlyVsLate)
design(dds_IsletTransplant_EarlyVsLate)
resultsNames(dds_IsletTransplant_EarlyVsLate)
res_EarlyVsLate_Without_Treatment <- results(dds_IsletTransplant_EarlyVsLate,
                                             name = "Group_Late.Rejection_vs_Early.Rejection")

result1  <- as.data.frame(res_EarlyVsLate_Without_Treatment);  result1$gene  <- rownames(res_EarlyVsLate_Without_Treatment)
keyvals_result1  <- make_keyvals_fdr_fc(result1)
selLab_result1  <- pick_labels(result1,  q_cut, fc_cut, 30)

## Axis limits
xmax_result1  <- max(2, ceiling(max(abs(result1$log2FoldChange),  na.rm=TRUE)))

# # replace any nas with 1
# result$pvalue[is.na(result$pvalue)] <- 1
# result$padj[is.na(result$padj)] <- 1

dev.new(width = 10, height = 10)
EnhancedVolcano(
  result1,
  lab           = result1$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Early vs Late Progessors",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(result1$padj<0.10 & abs(result1$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_result1, xmax_result1),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_result1,
  legendPosition= "right",
  selectLab     = selLab_result1
)

# Design Formula Accounting for Treatment
design(dds_IsletTransplant_EarlyVsLate) <- ~ Day + Treatment + Group
dds_IsletTransplant_EarlyVsLate <- DESeq(dds_IsletTransplant_EarlyVsLate)
design(dds_IsletTransplant_EarlyVsLate)
resultsNames(dds_IsletTransplant_EarlyVsLate)
res_EarlyVsLate_With_Treatment <- results(dds_IsletTransplant_EarlyVsLate,
                                             name = "Group_Late.Rejection_vs_Early.Rejection")

result2  <- as.data.frame(res_EarlyVsLate_With_Treatment);  result2$gene  <- rownames(res_EarlyVsLate_With_Treatment)
keyvals_result2  <- make_keyvals_fdr_fc(result2)
selLab_result2  <- pick_labels(result2,  q_cut, fc_cut, 30)

## Axis limits
xmax_result2  <- max(2, ceiling(max(abs(result2$log2FoldChange),  na.rm=TRUE)))

dev.new(width = 10, height = 10)
EnhancedVolcano(
  result2,
  lab           = result2$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Early vs Late Progessors",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(result2$padj<0.10 & abs(result2$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_result2, xmax_result2),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_result2,
  legendPosition= "right",
  selectLab     = selLab_result2
)

# Determining if LogIEq should be a covariate in the design forumla ----
# Design Formula NOT Accounting for Treatment and Accounting for LogIEq
design(dds_IsletTransplant_EarlyVsLate) <- ~ Day + logIEQ + Group
dds_IsletTransplant_EarlyVsLate <- DESeq(dds_IsletTransplant_EarlyVsLate)
design(dds_IsletTransplant_EarlyVsLate)
resultsNames(dds_IsletTransplant_EarlyVsLate)
res_EarlyVsLate_Without_logIEq <- results(dds_IsletTransplant_EarlyVsLate,
                                             name = "Group_Late.Rejection_vs_Early.Rejection")

result3  <- as.data.frame(res_EarlyVsLate_Without_logIEq);  result3$gene  <- rownames(res_EarlyVsLate_Without_logIEq)
keyvals_result3  <- make_keyvals_fdr_fc(result3)
selLab_result3  <- pick_labels(result3,  q_cut, fc_cut, 30)

## Axis limits
xmax_result3  <- max(2, ceiling(max(abs(result3$log2FoldChange),  na.rm=TRUE)))

dev.new(width = 10, height = 10)
EnhancedVolcano(
  result3,
  lab           = result3$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Early vs Late Progessors",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(result3$padj<0.10 & abs(result3$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_result3, xmax_result3),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_result3,
  legendPosition= "right",
  selectLab     = selLab_result3
)

# Design Formula Accounting for logIEq and Treatment
design(dds_IsletTransplant_EarlyVsLate) <- ~ Day + Treatment + logIEQ + Group
dds_IsletTransplant_EarlyVsLate <- DESeq(dds_IsletTransplant_EarlyVsLate)
design(dds_IsletTransplant_EarlyVsLate)
resultsNames(dds_IsletTransplant_EarlyVsLate)
res_EarlyVsLate_With_logIEq <- results(dds_IsletTransplant_EarlyVsLate,
                                          name = "Group_Late.Rejection_vs_Early.Rejection")

result4  <- as.data.frame(res_EarlyVsLate_With_logIEq);  result4$gene  <- rownames(res_EarlyVsLate_With_logIEq)
keyvals_result4  <- make_keyvals_fdr_fc(result4)
selLab_result4  <- pick_labels(result4,  q_cut, fc_cut, 30)

## Axis limits
xmax_result4  <- max(2, ceiling(max(abs(result4$log2FoldChange),  na.rm=TRUE)))

dev.new(width = 10, height = 10)
EnhancedVolcano(
  result4,
  lab           = result4$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Early vs Late Progessors",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(result4$padj<0.10 & abs(result4$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_result4, xmax_result4),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_result4,
  legendPosition= "right",
  selectLab     = selLab_result4
)
