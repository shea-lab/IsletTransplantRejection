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
library(ggnewscale)

# MSIGDBR Pathways ----
# Needs msigdbr package: https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
msigdbr_collections() 
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


# 1. Load the Data ----

#Metadata Importing
NOD_meta_batch <- read.csv("/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/NOD_Transplant/NODTransplant_Metadata.csv", 
                       header = TRUE, 
                       check.names = FALSE)

NOD_meta_batch <- as.data.frame(NOD_meta_batch)

# Preview the combined metadata
head(NOD_meta_batch)
unique(NOD_meta_batch$Group)

#Counts Data Importing
NOD_counts_batch <- as.data.frame(read.table(
  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/NOD_Transplant/NODTransplant_gene_counts_annot.csv"
  , sep=",", header=T,check.names = FALSE))
NOD_counts_batch <- na.omit(NOD_counts_batch)

#Remove duplicate names
NOD_counts_batch <- NOD_counts_batch[!duplicated(NOD_counts_batch[, 1]), ]
genes <- NOD_counts_batch[, 1]
rownames(NOD_counts_batch) <- genes
NOD_counts_batch <- NOD_counts_batch[, -1]

# 2. Preprocessing and Cleaning ----

NOD_counts_batch <- NOD_counts_batch[, NOD_meta_batch$Samples]  # Ensure Sample_IDs match column names in combined_counts

NODTransplantCounts <- flexiDEG.function1(NOD_counts_batch, NOD_meta_batch, # Genes in rows 
                                            convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                            batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0, 0

remove_pattern <- "^Gm[0-9]|^AC[0-9]|^AL[0-9]|^AI[0-9]|^AW[0-9]|^AF[0-9]|^BB[0-9]|^BC[0-9]|^CT[0-9]|^CAAA|^BX[0-9]|^CN[0-9]|^CR[0-9]|^C[0-9]{4,}|^Olfr"
rows_to_remove <- grep(remove_pattern, rownames(NODTransplantCounts)) #Remove pseudo genes
NODTransplantCounts <- NODTransplantCounts[-rows_to_remove, ]
mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "www"
)
# get gene biotypes
gene_info <- getBM(
  attributes = c("mgi_symbol", "gene_biotype"),
  filters = "mgi_symbol",
  values = rownames(NODTransplantCounts),
  mart = mart
)
#setdiff(rownames(IsletTransplantCounts), gene_info$mgi_symbol)
pseudo_genes <- gene_info$mgi_symbol[grep("pseudogene", gene_info$gene_biotype)]
# remove them
NODTransplantCounts <- NODTransplantCounts[
  !(rownames(NODTransplantCounts) %in% pseudo_genes), ]

write.csv(NODTransplantCounts, file = "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/NOD_Transplant/NOD_Raw_Counts_ITx_Filtered.csv", row.names = TRUE)


# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(NOD_meta_batch)
unique_colSide <- unique(colSide)


#3. Create Univeral DESqEQ Object ----

NODTransplantCounts <- as.matrix(NODTransplantCounts)
storage.mode(NODTransplantCounts) <- "integer"

dds_NODTransplantCounts<- DESeqDataSetFromMatrix(NODTransplantCounts, NOD_meta_batch,
                                              design = ~ 1)   # dummy design for now


saveRDS(dds_NODTransplantCounts, file = "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Data/NOD_Transplant/dds_NODTransplantCounts.rds")


# 4. Early Rejection Vs Late Rejection Signature----

# subset samples till day 28 since bith groups are present
sel <- colData(dds_NODTransplantCounts)$Day < 30 & 
  colData(dds_NODTransplantCounts)$Day != 0
dds_NODTransplantCounts_EarlyVsLate <- dds_NODTransplantCounts[, sel]


dds_NODTransplantCounts_EarlyVsLate$Batch       <- factor(dds_NODTransplantCounts_EarlyVsLate$Batch)
dds_NODTransplantCounts_EarlyVsLate$LibraryPrep <- factor(dds_NODTransplantCounts_EarlyVsLate$LibraryPrep)
dds_NODTransplantCounts_EarlyVsLate$Group <- factor(
  dds_NODTransplantCounts_EarlyVsLate$Group,
  levels = c("Late Rejection","Early Rejection")  # order sets baseline
)
dds_NODTransplantCounts_EarlyVsLate$Treatment <- factor(dds_NODTransplantCounts_EarlyVsLate$Treatment,
                                                    levels = c("Low Dose", "High Dose"))
#Bin timepoints for analysis
colData(dds_NODTransplantCounts_EarlyVsLate)$Day <-
  ifelse(colData(dds_NODTransplantCounts_EarlyVsLate)$Day >= 22,
         28,
         colData(dds_NODTransplantCounts_EarlyVsLate)$Day)

dds_NODTransplantCounts_EarlyVsLate$Day <- factor(dds_NODTransplantCounts_EarlyVsLate$Day,
                                                        levels = c("7", "14","28"))
levels(dds_NODTransplantCounts_EarlyVsLate$Day)

# Add log(IEQ) column
colData(dds_NODTransplantCounts_EarlyVsLate)$logIEQ <- log(colData(dds_NODTransplantCounts_EarlyVsLate)$IEQ)
cd <- as.data.frame(colData(dds_NODTransplantCounts_EarlyVsLate))
# Basic sanity
lapply(cd[, c("Batch","LibraryPrep","Day","Group")], function(x) table(x, useNA="ifany"))
# Check for NAs
sapply(cd[, c("Batch","LibraryPrep","Day","Group")], function(x) any(is.na(x)))
# plot
boxplot(logIEQ ~ Group, data = colData(dds_NODTransplantCounts_EarlyVsLate)) 
table(cd$Group,cd$Day)
#             7  14  28
# Late        5   5   5
# Early       5   6   8
keep <- rowSums(counts(dds_NODTransplantCounts_EarlyVsLate) >= 10) >= 5  # Smallest number of samples in a group
dds_NODTransplantCounts_EarlyVsLate <- dds_NODTransplantCounts_EarlyVsLate[keep, ]

#Run DESEQ with Treatment and Day as covirate----

design(dds_NODTransplantCounts_EarlyVsLate) <- ~ Treatment+ Day + Group
dds_NODTransplantCounts_EarlyVsLate <- DESeq(dds_NODTransplantCounts_EarlyVsLate)
design(dds_NODTransplantCounts_EarlyVsLate)
resultsNames(dds_NODTransplantCounts_EarlyVsLate)
res_NOD_EarlyVsLate <- results(dds_NODTransplantCounts_EarlyVsLate,
                                             name = "Group_Early.Rejection_vs_Late.Rejection")

result_NOD_EarlyVsLate  <- as.data.frame(res_NOD_EarlyVsLate);  result_NOD_EarlyVsLate$gene  <- rownames(res_NOD_EarlyVsLate)
write.csv(result_NOD_EarlyVsLate, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/NOD Islet Transplant/DESEQResults_EarlyvsLateRejection_NOD_DaysCombined.csv", row.names = TRUE)


# Elastic Net Feature Selection----

sig_genes_EARLYvLATE <- result_NOD_EarlyVsLate$gene[
  !is.na(result_NOD_EarlyVsLate$pvalue) &
    result_NOD_EarlyVsLate$pvalue <= 0.05 &
    abs(result_NOD_EarlyVsLate$log2FoldChange) >= 0.5
]

length(sig_genes_EARLYvLATE) #149 genes

# Regress out attributable to Day and Treatment, Preserve variation assciated with Group
vsd <- vst(dds_NODTransplantCounts_EarlyVsLate, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_NODTransplantCounts_EarlyVsLate))
design_enet <- model.matrix(~ Group, data = cd)
mat_enet <- limma::removeBatchEffect(
  mat,
  covariates = model.matrix(~ Day+Treatment, data = cd)[, -1, drop = FALSE],
  design = design_enet
)


# expression matrix for glmnet: samples x genes
x_EarlyVsLate <- t(mat_enet[sig_genes_EARLYvLATE, , drop = FALSE])

# binary outcome
y_EarlyVsLate <- ifelse(cd$Group == "Early Rejection", 1, 0)

# Find best alpha with LOOCV
set.seed(123)
alpha_grid <- seq(0, 1, by = 0.1)
cv_summary <- data.frame()

for (a in alpha_grid) {
  cvfit <- cv.glmnet(
    x = x_EarlyVsLate,
    y = y_EarlyVsLate,
    family = "binomial",
    alpha = a,
    foldid = 1:length(y_EarlyVsLate),   # LOOCV
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
best_alpha<-0.6 # Second best, but gives more genes for interpretation

# Stability selection
set.seed(123)
n_iter <- 1000
n_samples <- nrow(x_EarlyVsLate)

selected_list <- vector("list", n_iter)

class1_idx <- which(y_EarlyVsLate == unique(y_EarlyVsLate)[1])
class2_idx <- which(y_EarlyVsLate == unique(y_EarlyVsLate)[2])

for (i in 1:n_iter) {
  
  # Subsample ~80% of samples each time
  idx1 <- sample(class1_idx, size = round(0.8 * length(class1_idx)))
  idx2 <- sample(class2_idx, size = round(0.8 * length(class2_idx)))
  idx  <- c(idx1, idx2)
  
  x_sub <- x_EarlyVsLate[idx, ]
  y_sub <- y_EarlyVsLate[idx]
  
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
all_genes <- colnames(x_EarlyVsLate)

freq <- sapply(all_genes, function(g) {
  mean(sapply(selected_list, function(s) g %in% s))
})

freq_table <- data.frame(
  gene = names(freq),
  selection_frequency = freq
)
freq_table <- freq_table[order(freq_table$selection_frequency, decreasing = TRUE), ]
stable_genes <- subset(freq_table, selection_frequency >= 0.7) #Select genes appearing atleast 70 percent of time
stable_genes #13 genes

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
  Day = factor(cd$Day, levels = c("7", "14","28")),
  Group = factor(
    cd$Group,
    levels = c("Early Rejection", "Late Rejection")
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
      "Early Rejection" = "#B23A48",    # Tangerine
      "Late Rejection" =  "#2A6F97"    # Ocean
    ),
    Day = c(
      "7" = "#F28E6B",
      "14" = "#6FA287",
      "28" = "#E9C46A"
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
  
  column_title = c("Day 7", "Day 14", "Day 28"),
  
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
  
  col =circlize::colorRamp2(
    c(min(mat_scaled), 0, max(mat_scaled)),
    c("navy", "white", "firebrick3")
  ),
  
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



#PCA Analysis using IN Signature
mat_pca <- t(mat_stable)
pca <- prcomp(mat_pca, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = cd$Group
)

# Define colors & shapes
group_colors <- c( "Early Rejection" = "#B23A48",    # Tangerine
                   "Late Rejection" =  "#2A6F97" )
group_shapes <- c("Early Rejection" = 15, "Late Rejection" = 16)

ggplot(pca_df, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(aes(fill = Group), size = 5, stroke = 1.2) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA of NOD Early vs Late Rejecyion",
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


vsd <- vst(dds_NODTransplantCounts_EarlyVsLate, blind = FALSE)

mat_temporal <- assay(vsd)

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
      levels = c("7", "14","28")
    ),
    Group = factor(
      Group,
      levels = c(
        "Early Rejection",
        "Late Rejection"
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
      "Early Rejection" = "#B23A48",    
      "Late Rejection" =  "#2A6F97" 
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Early Rejection" = "#B23A48",    
      "Late Rejection" =  "#2A6F97"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Early Rejection" = 15,  # upward triangle
      "Late Rejection" = 16  # downward triangle
    )
  ) +
  
  labs(
    title = "Temporal Expression- Auto-Allo Rejection: Early vs Late Signature",
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



# GSVA Scores Over Time----

## Biological Pathway Scores----
dds_use <- dds_NODTransplantCounts

## 1) Remove sample 15402-JR-8,37,55 (Excluding early rejection samples beyond Day 28 since these are scattered and not enough n to do stats)
dds_use <- dds_use[, !(colnames(dds_use) %in% c("15402-JR-37", "15402-JR-8","15402-JR-55"))]

## 2) Recode Day in metadata
cd <- as.data.frame(colData(dds_use))

cd$Day <- as.numeric(as.character(cd$Day))
unique(cd$Day)
#Binning of Days for analysis
cd$Day[cd$Day %in% c(22, 23, 26)] <- 28
cd$Day[cd$Day == 45] <- 42
cd$Day[cd$Day == 55] <- 56
cd$Day[cd$Day %in% c(61, 63, 71)] <- 70

## Put corrected Day back into dds
colData(dds_use)$Day <- cd$Day
table(cd$Day,cd$Group)
#.    Early Rejection Late Rejection
#0                8              6
#7                5              5
#14               6              5
#28               8              5
#42               0              5
#56               0              4
#70               0              4
## make Day a factor if you want it treated as discrete timepoints
colData(dds_use)$Day <- factor(colData(dds_use)$Day,
                               levels = c(0,7, 14, 28, 42, 56, 70))

dds_use$Group <- factor(
  dds_use$Group,
  levels = c("Late Rejection","Early Rejection")  # order sets baseline
)
levels(dds_use$Group)


# VST normalization
vsd <- vst(dds_use, blind = FALSE)
mat_gsva<- assay(vsd)
cd <- as.data.frame(colData(dds_use))
cd$MouseID <- paste(cd$Cohort, cd$Animal, sep = "_")


genesetsOfinterest<-c(  
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "AIZARANI_LIVER_C1_NK_NKT_CELLS_1",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"
  )

mm_list <- mm_all_df %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop")

# convert to named list
gene_sets_all <- setNames(mm_list$genes, mm_list$gs_name)
# subset
gene_sets_use <- gene_sets_all[genesetsOfinterest]

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

gsva_long <- gsva_df %>%
  pivot_longer(
    cols = all_of(genesetsOfinterest),
    names_to = "Pathway",
    values_to = "Score"
  )

gsva_long <- gsva_long %>%
  dplyr::mutate(
    Pathway = factor(Pathway, levels = genesetsOfinterest)
  )

pvals_df <- gsva_long %>%
  dplyr::group_by(Pathway, Day) %>%
  dplyr::summarise(
    p_value = {
      df <- dplyr::cur_data()
      if(length(unique(df$Group)) < 2) {
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
      p_value <= 1e-1 ~  "+",   
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
    y_pos = ymax - 0.1 * span
  ) %>%
  dplyr::select(Pathway, Day, y_pos)

sig_df <- pvals_df %>%
  dplyr::left_join(label_pos_df, by = c("Pathway", "Day"))



pathway_labels <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNF–NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "IL6–JAK-STAT3",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "IFN-gamma Response",
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammation",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY" = "T Cell Signaling",
  "AIZARANI_LIVER_C1_NK_NKT_CELLS_1" = "NK/NKT Cells"
)

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
    nrow = 3,
    ncol = 3,
    labeller = labeller(Pathway = pathway_labels)
  ) +
  coord_cartesian(ylim = c(-0.6, 0.6)) +
  scale_color_manual(values = c(
    "Early Rejection" =  "#B23A48",
    "Late Rejection" = "#2A6F97" 
  )) +
  scale_shape_manual(values = c(
    "Early Rejection" = 17,
    "Late Rejection" = 16
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = "Days Post Transplant",
    y = "GSVA score",
    title = "NOD Rejection: Immune Inflammatory Genesets"
  )+
  theme(
    strip.text = element_text(size = 16)
  )

library(writexl) 

gsva_NOD_pathways <- gsva_df  
# Write to Excel
write_xlsx(gsva_NOD_pathways, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Table8_NODAllo_Early_vs_Late_Rejection.xlsx")

## EN  Scores----

dds_use <- dds_NODTransplantCounts
## 1) Remove sample 15402-JR-8,37,55 
dds_use <- dds_use[, !(colnames(dds_use) %in% c("15402-JR-37", "15402-JR-8","15402-JR-55"))]
## 2) Recode Day in metadata
cd <- as.data.frame(colData(dds_use))

cd$Day <- as.numeric(as.character(cd$Day))
cd$Day[cd$Day %in% c(22, 23, 26)] <- 28
cd$Day[cd$Day == 45] <- 42
cd$Day[cd$Day == 55] <- 56
cd$Day[cd$Day %in% c(61, 63, 71)] <- 70

## Put corrected Day back into dds
colData(dds_use)$Day <- cd$Day

## make Day a factor if you want it treated as discrete timepoints
colData(dds_use)$Day <- factor(colData(dds_use)$Day,
                               levels = c(0,7, 14, 28, 42, 56, 70))

dds_use$Group <- factor(
  dds_use$Group,
  levels = c("Late Rejection","Early Rejection")  # order sets baseline
)

# VST normalization
vsd <- vst(dds_use, blind = FALSE)
mat_enet <- assay(vsd)

up_genes_EarlyVsLate <- c("Rsph10b","S1pr5","Erdr1","Pla2g4b","Clec2g","Cstg","Syce1","Fndc7")
down_genes_EarlyVsLate <- c("Epgn","Tmem132e","Stkld1","Osbpl6","Sycp3")
gsva_par_sig <- gsvaParam(
  mat_enet, 
  list(
    Up = up_genes_EarlyVsLate,
    Down = down_genes_EarlyVsLate
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
      p_value <= 1e-1 ~  "†", 
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

score_avg <- score_avg %>%
  dplyr::mutate(
    Day = factor(Day, levels = c(0, 7, 14, 28, 42, 56, 70))
  )

sig_score_df <- sig_score_df %>%
  dplyr::mutate(
    Day = factor(Day, levels = c(0, 7, 14, 28, 42, 56, 70))
  )


ggplot(score_avg,
       aes(x = Day,
           y = mean_score,
           color = Group,
           shape = Group,
           group = Group)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 4) +
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
    "Early Rejection" =  "#B23A48",
    "Late Rejection" = "#2A6F97" 
  )) +
  scale_shape_manual(values = c(
    "Early Rejection" = 17,
    "Late Rejection" = 16
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = "Days Post Transplant",
    y = "Early vs Late Rejection GSVA Score",
    title = "Early vs Late Rejection GSVA Score Over Time"
  )


gsva_NOD_Rejscore <- signature_df  
# Write to Excel
write_xlsx(gsva_NOD_Rejscore, "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Supplementary Tables/Supp_Table8_RejScore_NODAllo_Early_vs_Late_Rejection.xlsx")

# GSEA Analysis----
# Paths to your saved results
result_NOD_EarlyVsLate_path  <-  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/NOD Islet Transplant/DESEQResults_EarlyvsLateRejection_NOD_DaysCombined.csv"

# Import
result_NOD_EarlyVsLate  <- read.csv(result_NOD_EarlyVsLate_path,  row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_ALL_EARLYVSLATE  <- result_NOD_EarlyVsLate$stat;  names(lfc_vector_ALL_EARLYVSLATE)  <- rownames(result_NOD_EarlyVsLate)


# Drop NAs
lfc_vector_ALL_EARLYVSLATE  <- lfc_vector_ALL_EARLYVSLATE[!is.na(lfc_vector_ALL_EARLYVSLATE)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_ALL_EARLYVSLATE  <- sort(lfc_vector_ALL_EARLYVSLATE,  decreasing = TRUE)

gsea_results_ALL_EARLYVSLATE<- GSEA(
  geneList      = lfc_vector_ALL_EARLYVSLATE,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_ALL_EARLYVSLATE <- as.data.frame(gsea_results_ALL_EARLYVSLATE)


# Full results Combined Timepoins
write.csv(gsea_results_ALL_EARLYVSLATE,
          "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/NOD Islet Transplant/GSEAResults_EarlyvsLateRejection_NOD_DaysCombined.csv",
          row.names = FALSE)

# GSEA: Auto-allo Vs  Allo  ----

immune_master_EARLYVSLATE <- unique(c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "TRAVAGLINI_LUNG_NEUTROPHIL_CELL",
  "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CEL",
 "HE_LIM_SUN_FETAL_LUNG_C2_S100A12_HI_CLASSICAL_MONOCYTE",
 "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
 "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CELL",
 "HE_LIM_SUN_FETAL_LUNG_C2_APOE_POS_M2_MACROPHAGE_CELL",
 "AIZARANI_LIVER_C1_NK_NKT_CELLS_1"
))

# NOD Result
result_NOD_EarlyVsLate_path  <-  "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Sequencing Results/NOD Islet Transplant/DESEQResults_EarlyvsLateRejection_NOD_DaysCombined.csv"
# Import
result_NOD_EarlyVsLate  <- read.csv(result_NOD_EarlyVsLate_path,  row.names = 1)
# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_ALL_EARLYVSLATE  <- result_NOD_EarlyVsLate$stat;  names(lfc_vector_ALL_EARLYVSLATE)  <- rownames(result_NOD_EarlyVsLate)
# Drop NAs
lfc_vector_ALL_EARLYVSLATE  <- lfc_vector_ALL_EARLYVSLATE[!is.na(lfc_vector_ALL_EARLYVSLATE)]
# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_ALL_EARLYVSLATE  <- sort(lfc_vector_ALL_EARLYVSLATE,  decreasing = TRUE)
gsea_results_ALL_EARLYVSLATE<- GSEA(
  geneList      = lfc_vector_ALL_EARLYVSLATE,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_ALL_EARLYVSLATE <- as.data.frame(gsea_results_ALL_EARLYVSLATE)


# Allogeneic
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
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  #keyType       = "SYMBOL",       # <- tell it explicitly
  TERM2GENE     = mm_all_df
)
gsea_results_ALL_REJVSACCEP <- as.data.frame(gsea_results_ALL_REJVSACCEP)


plot_df_nod <- gsea_results_ALL_EARLYVSLATE %>%
  dplyr::filter(ID %in% immune_master_EARLYVSLATE) %>%
  dplyr::select(ID, Description, NES, p.adjust) %>%
  dplyr::mutate(
    Comparison = "NOD Early vs Late",
    neglog10_padj = -log10(p.adjust + 1e-300)
  )

plot_df_allo <- gsea_results_ALL_REJVSACCEP %>%
  dplyr::filter(ID %in% immune_master_EARLYVSLATE) %>%
  dplyr::select(ID, Description, NES, p.adjust) %>%
  dplyr::mutate(
    Comparison = "Allogeneic Rej vs Accep",
    neglog10_padj = -log10(p.adjust + 1e-300)
  )

plot_df <- bind_rows(plot_df_nod, plot_df_allo)

pathway_order <- plot_df %>%
  filter(Comparison == "Allogeneic Rej vs Accep") %>%
  arrange(NES) %>%   # increasing NES (bottom = most negative, top = most positive)
  pull(ID)

plot_df$ID <- factor(plot_df$ID, levels = pathway_order)

plot_df_nod2 <- plot_df %>%
  filter(Comparison == "NOD Early vs Late") %>%
  mutate(ypos = as.numeric(factor(ID, levels = levels(plot_df$ID))) + 0.01)

plot_df_allo2 <- plot_df %>%
  filter(Comparison == "Allogeneic Rej vs Accep") %>%
  mutate(ypos = as.numeric(factor(ID, levels = levels(plot_df$ID))) - 0.01)

ggplot() +
  geom_point(
    data = plot_df_nod2,
    aes(x = NES, y = ypos, size = neglog10_padj, color = NES),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#2A6F97",
    mid = "white",
    high = "#B23A48",
    midpoint = 0,
    name = "Auto-Allo Rejection"
  ) +
  ggnewscale::new_scale_color() +
  geom_point(
    data = plot_df_allo2,
    aes(x = NES, y = ypos, size = neglog10_padj, color = NES),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#064273",
    mid = "white",
    high = "#F28500",
    midpoint = 0,
    name = "Allo Rejection"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(
    breaks = seq_along(levels(plot_df$ID)),
    labels = levels(plot_df$ID)
  ) +
  scale_size_continuous(
    name = expression(-log[10](adjusted~italic(p))),
    range = c(3, 10)
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


# 5. Allo Genes Signature Cross Preservatiom----
# subset samples of D7 and D14
dds_NODTransplantCounts
dds_NODTransplantCounts<- DESeqDataSetFromMatrix(NODTransplantCounts, NOD_meta_batch,
                                                 design = ~ 1)   # dummy design for now


# subset samples till day 28 since bith groups are present
sel <- colData(dds_NODTransplantCounts)$Day < 30 & 
  colData(dds_NODTransplantCounts)$Day != 0
dds_NODTransplantCounts_EarlyVsLate <- dds_NODTransplantCounts[, sel]


dds_NODTransplantCounts_EarlyVsLate$Batch       <- factor(dds_NODTransplantCounts_EarlyVsLate$Batch)
dds_NODTransplantCounts_EarlyVsLate$LibraryPrep <- factor(dds_NODTransplantCounts_EarlyVsLate$LibraryPrep)
dds_NODTransplantCounts_EarlyVsLate$Group <- factor(
  dds_NODTransplantCounts_EarlyVsLate$Group,
  levels = c("Late Rejection","Early Rejection")  # order sets baseline
)
dds_NODTransplantCounts_EarlyVsLate$Treatment <- factor(dds_NODTransplantCounts_EarlyVsLate$Treatment,
                                                        levels = c("Low Dose", "High Dose"))
#Bin timepoints for analysis
colData(dds_NODTransplantCounts_EarlyVsLate)$Day <-
  ifelse(colData(dds_NODTransplantCounts_EarlyVsLate)$Day >= 22,
         28,
         colData(dds_NODTransplantCounts_EarlyVsLate)$Day)

dds_NODTransplantCounts_EarlyVsLate$Day <- factor(dds_NODTransplantCounts_EarlyVsLate$Day,
                                                  levels = c("7", "14","28"))

design(dds_NODTransplantCounts_EarlyVsLate) <- ~ Treatment+ Day + Group
dds_NODTransplantCounts_EarlyVsLate <- DESeq(dds_NODTransplantCounts_EarlyVsLate)
design(dds_NODTransplantCounts_EarlyVsLate)
resultsNames(dds_NODTransplantCounts_EarlyVsLate)
res_NOD_EarlyVsLate <- results(dds_NODTransplantCounts_EarlyVsLate,
                               name = "Group_Early.Rejection_vs_Late.Rejection")

vsd <- vst(dds_NODTransplantCounts_EarlyVsLate, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_NODTransplantCounts_EarlyVsLate))
design_enet <- model.matrix(~ Group, data = cd)
mat_enet <- limma::removeBatchEffect(
  mat,
  covariates = model.matrix(~ Day+Treatment, data = cd)[, -1, drop = FALSE],
  design = design_enet
)


#Plot Heatmp and PCA using Selected Genes 
AlloaCD40L_signature <- c(
  "Cd5l", "Cpa6", "Ebf4", "Ifi30", "Cmah",
  "Hs3st3a1", "Snhg11", "Tspoap1", "Tymp", "Sarm1",
  "Ceacam19", "Ptprf", "Gp5", "Eda", "Pla2g2d",
  "Flywch2", "Fibin", "Lurap1", "Cyp2j9", "Prrg1",
  "Treml1", "Ltc4s", "Dok5", "Fabp3", "Nrbp2"
)

# AlloSyn_signature <- c(
#   "Hipk4", "Lncpint", "Malat1", "Myh3", "Rep15", "Trpm6", "Cd59b",
#   "Hlf", "Bmp5", "S1pr5", "Cbs", "Tctn2", "Gdf3", "Shisa9", "Ido2",
#   "Vtn", "Ifitm5", "Emx2os", "Fsip1", "Gfra1", "Col6a5", "Lhfpl4",
#   "Myo18b", "Myo3b", "Sybu", "Dnah8", "Rragb"
# )
setdiff(AlloaCD40L_signature, rownames(mat_enet))
mat_stable <- mat_enet[AlloaCD40L_signature, , drop = FALSE]
mat_scaled <- t(scale(t(mat_stable)))
annotation_col <- data.frame(
  Day = cd$Day,
  Group = cd$Group
)
rownames(annotation_col) <- colnames(mat_scaled)


# Column annotation
annotation_col <- data.frame(
  Day = factor(cd$Day, levels = c("7", "14", "28")),
  Group = factor(
    cd$Group,
    levels = c("Early Rejection", "Late Rejection")
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
      "Early Rejection" = "#B23A48",    # Tangerine
      "Late Rejection" =  "#2A6F97"    # Ocean
    ),
    Day = c(
      "7" = "#F28E6B",
      "14" = "#6FA287",
      "28" = "#E9C46A"
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
  
  column_title = c("Day 7", "Day 14", "Day 28"),
  
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
  
  col =circlize::colorRamp2(
    c(min(mat_scaled), 0, max(mat_scaled)),
    c("navy", "white", "firebrick3")
  ),
  
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
group_colors <- c( "Early Rejection" = "#B23A48",    # Tangerine
                   "Late Rejection" =  "#2A6F97" )
group_shapes <- c("Early Rejection" = 15, "Late Rejection" = 16)

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

dds_use <- dds_NODTransplantCounts

## 1) Remove sample 15402-JR-8,37,55 (Excluding early rejection samples beyond Day 28 since these are scattered and not enough n to do stats)
dds_use <- dds_use[, !(colnames(dds_use) %in% c("15402-JR-37", "15402-JR-8","15402-JR-55"))]

## 2) Recode Day in metadata
cd <- as.data.frame(colData(dds_use))

cd$Day <- as.numeric(as.character(cd$Day))
unique(cd$Day)
#Binning of Days for analysis
cd$Day[cd$Day %in% c(22, 23, 26)] <- 28
cd$Day[cd$Day == 45] <- 42
cd$Day[cd$Day == 55] <- 56
cd$Day[cd$Day %in% c(61, 63, 71)] <- 70

## Put corrected Day back into dds
colData(dds_use)$Day <- cd$Day
table(cd$Day,cd$Group)

colData(dds_use)$Day <- factor(colData(dds_use)$Day,
                               levels = c(0,7, 14, 28, 42, 56, 70))

dds_use$Group <- factor(
  dds_use$Group,
  levels = c("Late Rejection","Early Rejection")  # order sets baseline
)
levels(dds_use$Group)

dds_use$Day
# VST normalization
vsd <- vst(dds_use, blind = FALSE)
mat_temporal<- assay(vsd)
mat_plot <- mat_temporal[
#  AlloaCD40L_signature,
  AlloSyn_signature,
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
      levels = c("0","7","14","28","42","56","70")
    ),
    Group = factor(
      Group,
      levels = c(
        "Early Rejection",
        "Late Rejection"
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
      "Early Rejection" = "#B23A48",    # Tangerine
      "Late Rejection" =  "#2A6F97"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Early Rejection" = "#B23A48",    # Tangerine
      "Late Rejection" =  "#2A6F97"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Early Rejection" = 15,  # upward triangle
      "Late Rejection" = 16  # downward triangle
    )
  ) +
  
  labs(
    title = "Temporal Expression- Allo vs Syn Signature in Auto-Allo NOD",
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

# 6. Luminex Analysis--
# getwd()
# library(readr)
# library(dplyr)
# library(tidyverse)
# library(lme4)
# library(lmerTest)
# library(emmeans)
# library(ggplot2)
# 
# Scaf_IsTx_Luminex <- read_csv(
#   "/Users/jyotirmoyroy/Desktop/Islet Transplant Rejection Paper/Luminex Data/IsletTransplantScafLuminex.csv",
#   show_col_types = FALSE
# ) %>%
#   # Remove columns where every value is NA
#   select(where(~ !all(is.na(.)))) %>%
#   # Remove rows where every value is NA
#   filter(!if_all(everything(), is.na))
# 
# Scaf_IsTx_Luminex
# colnames(Scaf_IsTx_Luminex)
# unique(Scaf_IsTx_Luminex$Sample)
# unique(Scaf_IsTx_Luminex$Day)
# unique(Scaf_IsTx_Luminex$Group)
# table(Scaf_IsTx_Luminex$Day,Scaf_IsTx_Luminex$Group)
# Scaf_IsTx_Luminex <- Scaf_IsTx_Luminex %>%
#   rename(
#     Mouse = `Mouse \nNumber`,
#     Total_Protein = `Total protein\n(ng/ul)`,
#     IL2 = `IL-2\r\n pg/ml`,
#     IL6 = `IL-6\r\n pg/ml`,
#     IFNg = `IFNg\r\n pg/ml`,
#     GMCSF = `GM-CSF\r\n pg/ml`,
#     TNFa = `TNFa\r\n pg/ml`
#   )
# 
# Scaf_IsTx_Luminex <- Scaf_IsTx_Luminex %>%
#   mutate(
#     Total_Protein_mg_ml = Total_Protein / 1000, #Convery to mg/ml
#     
#     IL2_norm   = IL2 / Total_Protein_mg_ml, #pg of cytokine/mg protein
#     IL6_norm   = IL6 / Total_Protein_mg_ml,
#     IFNg_norm  = IFNg / Total_Protein_mg_ml,
#     GMCSF_norm = GMCSF / Total_Protein_mg_ml,
#     TNFa_norm  = TNFa / Total_Protein_mg_ml
#   )
# 
# 
# Scaf_IsTx_Luminex <- Scaf_IsTx_Luminex %>%
#   mutate(
#     Day_plot = case_when(
#       Day %in% c(22, 23, 26, 28) ~ 28,
#       Day %in% c(55, 56) ~ 56,
#       Day %in% c(63, 71) ~ 70,
#       TRUE ~ Day
#     )
#   )
# Luminex_long <- Scaf_IsTx_Luminex %>%
#   select(
#     Sample, Mouse, Day, Day_plot, Group,
#     IL2_norm, IL6_norm, IFNg_norm, GMCSF_norm, TNFa_norm
#   ) %>%
#   pivot_longer(
#     cols = ends_with("_norm"),
#     names_to = "Cytokine",
#     values_to = "Normalized_Conc"
#   ) %>%
#   mutate(
#     Cytokine = recode(
#       Cytokine,
#       IL2_norm = "IL-2",
#       IL6_norm = "IL-6",
#       IFNg_norm = "IFN-\u03B3",
#       GMCSF_norm = "GM-CSF",
#       TNFa_norm = "TNF-\u03B1"
#     )
#   )
# 
# # Keep only the four cytokines
# Luminex_plot <- Luminex_long %>%
#   dplyr::filter(Cytokine != "GM-CSF") %>%
#   dplyr::mutate(
#     Day_factor = factor(
#       Day_plot,
#       levels = c(0, 7, 28, 56, 70)
#     ),
#     x_base = as.numeric(Day_factor),
#     
#     # Fixed position for each group
#     x_pos = dplyr::case_when(
#       Group == "Early Rejection" ~ x_base - 0.18,
#       Group == "Late Rejection"  ~ x_base + 0.18
#     )
#   )
# Luminex_summary <- Luminex_plot %>%
#   dplyr::group_by(
#     Cytokine,
#     Day_plot,
#     Group,
#     x_pos
#   ) %>%
#   dplyr::summarise(
#     mean_conc = mean(Normalized_Conc, na.rm = TRUE),
#     n = sum(!is.na(Normalized_Conc)),
#     se = sd(Normalized_Conc, na.rm = TRUE) / sqrt(n),
#     .groups = "drop"
#   )
# 
# ggplot() +
#   
#   # Mean bars
#   geom_col(
#     data = Luminex_summary,
#     aes(
#       x = x_pos,
#       y = mean_conc,
#       fill = Group
#     ),
#     width = 0.32,
#     color = "black",
#     linewidth = 0.7,
#     alpha = 0.55
#   ) +
#   
#   # Mean +/- SEM
#   geom_errorbar(
#     data = Luminex_summary,
#     aes(
#       x = x_pos,
#       ymin = mean_conc - se,
#       ymax = mean_conc + se
#     ),
#     width = 0.10,
#     linewidth = 0.7
#   ) +
#   
#   # Individual samples
#   geom_point(
#     data = Luminex_plot,
#     aes(
#       x = x_pos,
#       y = Normalized_Conc,
#       color = Group
#     ),
#     position = position_jitter(
#       width = 0.04,
#       height = 0
#     ),
#     size = 2.5,
#     alpha = 1
#   ) +
#   
#   facet_wrap(
#     ~ Cytokine,
#     scales = "free_y",
#     ncol = 2
#   ) +
#   
#   scale_x_continuous(
#     breaks = 1:5,
#     labels = c("0", "7", "28", "56", "70")
#   ) +
#   
#   scale_fill_manual(
#     values = c(
#       "Early Rejection" = "#E89A9F",
#       "Late Rejection"  = "#91B9D4"
#     )
#   ) +
#   
#   # Keep individual points darker
#   scale_color_manual(
#     values = c(
#       "Early Rejection" = "#B23A48",
#       "Late Rejection"  = "#2A6F97"
#     )
#   )+
#   
#   labs(
#     x = "Days Post Transplant",
#     y = "Cytokine concentration (pg/mg total protein)",
#     fill = NULL,
#     color = NULL
#   ) +
#   
#   theme_classic(base_size = 16) +
#   
#   theme(
#     strip.text = element_text(
#       size = 16,
#       face = "bold"
#     ),
#     legend.position = "top"
#   )
