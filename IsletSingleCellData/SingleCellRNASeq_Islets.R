library(dplyr)
#install.packages("Seurat")
library(Seurat)
library(SeuratObject)
library(patchwork)
library(BiocManager)
#install.packages('BiocManager')
BiocManager::install('glmGamPoi')
library(glmGamPoi)
#install.packages(c("rlang", "ggplot2", "vctrs"))
library(rlang)
library(ggplot2)
library(vctrs)

# Source Paper: Chen P, Yao F, Lu Y, Peng Y et al. Single-Cell Landscape of Mouse Islet Allograft and Syngeneic Graft. Front Immunol 2022;13:853349. PMID: 35757709
# GEO Source: GSE198865

# Import Data ----
getwd()
allo_data_dir <- "C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/Allogeneic"
allo_counts <- Read10X(data.dir = allo_data_dir)
allo_seurat <- CreateSeuratObject(
  counts = allo_counts,
  project = "Allogenic",
  min.cells = 3,
  min.features = 200)
rm(allo_counts)
syn_data_dir <- "C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/Syngeneic"
syn_counts <- Read10X(data.dir = syn_data_dir)
syn_seurat <- CreateSeuratObject(
  counts = syn_counts,
  project = "Syngeneic",
  min.cells = 3,
  min.features = 200)
rm(syn_counts)

# Quality Control ----
#  Add Mitochondrial Percentage
allo_seurat[["percent.mt"]] <- PercentageFeatureSet(allo_seurat, pattern = "^mt-")
syn_seurat[["percent.mt"]] <- PercentageFeatureSet(syn_seurat, pattern = "^mt-")

# Visualize QC metrics as a violin plot- Allo
VlnPlot(allo_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1_allo <- FeatureScatter(allo_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_allo <- FeatureScatter(allo_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_allo + plot2_allo

# Visualize QC metrics as a violin plot- Syn
VlnPlot(syn_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1_syn <- FeatureScatter(syn_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_syn <- FeatureScatter(syn_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_syn + plot2_syn
s
# Subset the data
allo_seurat <- subset(allo_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)
syn_seurat <- subset(syn_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)

# Create Islet Seurat Object ----
islet_graft_seurat <- merge(allo_seurat, y = syn_seurat, add.cell.ids = c("Allogeneic", "Syngeneic"))
DefaultAssay(islet_graft_seurat)<-"RNA"
# SCT transform to normalize data
islet_graft_seurat<- SCTransform(islet_graft_seurat, vars.to.regress = "percent.mt", verbose = TRUE)
islet_graft_seurat <- RunPCA(islet_graft_seurat, features = VariableFeatures(object = islet_graft_seurat))

print(islet_graft_seurat[["pca"]], dims = 1:5, nfeatures = 5)

# Initial Clustering ----
DimPlot(islet_graft_seurat, reduction = "pca") + NoLegend()
ElbowPlot(islet_graft_seurat)
islet_graft_seurat <- FindNeighbors(islet_graft_seurat, dims = 1:10)
islet_graft_seurat <- FindClusters(islet_graft_seurat, resolution = 0.5)
islet_graft_seurat <- RunUMAP(islet_graft_seurat, dims = 1:10)
islet_graft_seurat$condition <- islet_graft_seurat$orig.ident
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p1<-DimPlot(islet_graft_seurat, reduction = "umap")

p2<-DimPlot(islet_graft_seurat, reduction = "umap",group.by ="condition" )
p1+p2
saveRDS(islet_graft_seurat, "C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/ProcessedSingleCellObjects/islet_graft_seurat_v1.rds")

islet_graft_seurat<-readRDS("C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/ProcessedSingleCellObjects/islet_graft_seurat_v1.rds")
DimPlot(islet_graft_seurat, reduction = "umap",label = T,label.size = 4,label.box = T,pt.size = 1)


# marker_genes <- c(
#   # Immune cells
#   "Cd3","Ptprc",
#   "Cd4", "Cd8a",
#   "Cd68", "Fcgr3a",
#   "Clec9a", "Flt3",
#   "Ncr1", "Klra9",
#   "Cd19", "Ms4a1",
#   
#   # Non-immune cells
#   "Col3a1", "Fbn1",
#   "Pecam1", "Cdh5",
#   "Ins", "Chga"
# )

# Pulling genes linked to types of cells from paper
marker_genes <- list(
  # markers only found in immune cells
  "Immune"=c("Ptprc","Lyz2"),
  
  # figure unclear about what kind of genes these are- checked
  "Unknown" = c("Aldob", "Miox", "Fxyd2", "Spink1", "Gpx3"),
  
  # B Cell- checked
  "B Cell" = c("Cd19", "Cd79a", "Ms4a1", "Ly6d", "Cd79b", "Ebf1"),
  
  # Tconv- checked
  "Tconv" = c("Cd4", "Tnfsf8", "Lat", "Ets1", "Cd28", "Ms4a4b", "Cd3d"),
  
  # CD8+ T- checked
  "CD8+ T" = c("Lat", "Ets1", "Cd28","Cd8a", "Cd8b1", "Ms4a4b","Nkg7", "Cd3d"),
  
  #DC- checked
  "DC" = c("Clec9a", "Xcr1", "Cd24a", "Lsp1", "Flt3","Itgax","H2-Ab1"),
  
  #VEC- checked
  "VEC" = c("Pecam1", "Egfl7", "Plvap", "Emcn", "Ly6c1"),
  
  #Islet Cells- checked
  "Islet Cells" = c("Tmem27", "Chga", "Chgb", "Scg2", "Pcsk2"),
  
  # MO (M-Phi)- checked
  "MO" = c("Cd68", "Adgre1", "Csf1r", "Pla2g7", "Fcgr3", "Cybb"),
  
  # MES- checked
  "MES" = c("Col3a1", "Col1a1", "Col1a2", "Bgn", "Fstl1"),
  
  # NK- checked
  "NK" = c("Gzma", "Cd7", "Klrb1c", "Klrc2", "Klrk1"),
  
  # Treg- checked
  "Treg" = c("Il2ra", "Ctla4", "Cd2", "Tnfrsf4", "Ikzf2")
)

# Dot plot to annotate cell clusters
DotPlot(
  islet_graft_seurat,
  features = unique(unlist(marker_genes))
) +
  RotatedAxis() +
  scale_color_gradient(low = "grey80", high = "red")

# Annotation of Initial Clusters ----
new.cluster.ids <- c(
  "NK/TCell",
  "Macrophage",
  "Macrophage",
  "Macrophage",
  "Macrophage",
  "NK/TCell",
  "MES",
  "NK/TCell",
  "MES",
  "DC",
  "NK/TCell",
  "Islet Cells",
  "MES",
  "MES",
  "MES",
  "VEC",
  "BCell",
  "NK/TCell",
  "MES",
  "DC",
  "Macrophage",
  "DC",
  "Islet Cells"
)
names(new.cluster.ids) <- levels(islet_graft_seurat)
islet_graft_seurat <- RenameIdents(islet_graft_seurat, new.cluster.ids)
DimPlot(islet_graft_seurat, reduction = "umap",label = T,label.size = 4,label.box = T,pt.size = 1)

# v3 has NK/TCell label (Tcell sub-clustering yet to happen)
saveRDS(islet_graft_seurat, "C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/ProcessedSingleCellObjects/islet_graft_seurat_v3.rds")

islet_graft_seurat <- LoadSeuratRds(file = "C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/ProcessedSingleCellObjects/islet_graft_seurat_v3.rds")

# NK and T Cell Sub clustering ----
NK_TCell = subset(x = islet_graft_seurat, idents = "NK/TCell")
# Re-run SCTransform on the subset
NK_TCell <- SCTransform(NK_TCell, verbose = FALSE)
NK_TCell <- RunPCA(NK_TCell, features = VariableFeatures(object = NK_TCell))

ElbowPlot(NK_TCell)
NK_TCell <- FindNeighbors(NK_TCell, dims = 1:15, verbose = F)
NK_TCell <- FindClusters(NK_TCell, res = 0.7)
NK_TCell <- RunUMAP(NK_TCell, dims = 1:15)

NK_TCell$condition <- NK_TCell$orig.ident
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p1<-DimPlot(NK_TCell, reduction = "umap")
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p2<-DimPlot(NK_TCell, reduction = "umap",group.by ="condition" )
p1+p2

NK_TCell$cluster<-Idents(NK_TCell)

DimPlot(NK_TCell, reduction = "umap",label = T,label.size = 4,label.box = T,pt.size = 1)

marker_genes_TNK <- list(

  # B Cell- checked
  "B Cell" = c("Cd19", "Cd79a", "Ms4a1", "Ly6d", "Cd79b"),
  
  # Tconv- checked
  "Tconv" = c("Cd4", "Il7r", "Ikzf2","Cd3d", "Tnfsf8", "Lat", "Ets1", "Cd28", "Ms4a4b"),
  
  # CD8+ T- checked
  "CD8+ T" = c("Cd8a", "Cd8b1", "Trac", "Lat", "Ets1", "Cd28", "Ms4a4b","Nkg7", "Cd3d"),
  
  # NK- checked
  "NK" = c("Ncr1", "Klrb1c", "Klrk1", "Gzma", "Cd7", "Klrc2"),
  
  # Treg- checked
  "Treg" = c("Foxp3","Il2ra", "Ctla4", "Cd2", "Tnfrsf4", "Ikzf2")
)

# Dot plot to annotate cell clusters
DotPlot(
  NK_TCell,
  features = unique(unlist(marker_genes_TNK))
) +
  RotatedAxis() +
  scale_color_gradient(low = "grey80", high = "red")

# Further Cluster Analysis for Annotation- start
features_to_extract <- c("Cd4", "Cd8a", "Foxp3")
p <- DotPlot(
  NK_TCell,
  features = features_to_extract,
  idents = c(4, 7, 9)
)
dotplot_data <- p$data

DefaultAssay(NK_TCell) <- "SCT"
NK_TCell <- PrepSCTFindMarkers(NK_TCell)

cluster4.markers <- FindMarkers(NK_TCell, ident.1 = 4, only.pos = TRUE)
head(cluster4.markers, n = 5)
#                p_val avg_log2FC pct.1 pct.2 p_val_adj
# 2810417H13Rik      0   6.000251 0.976 0.039         0
# Birc5              0   5.358028 0.931 0.042         0
# Cks1b              0   4.707019 0.942 0.064         0
# Stmn1              0   4.854362 0.997 0.131         0
# Mki67              0   5.224568 0.960 0.111         0

cluster7.markers <- FindMarkers(NK_TCell, ident.1 = 7, only.pos = TRUE)
head(cluster7.markers, n = 5)
#         p_val          avg_log2FC  pct.1  pct.2  p_val_adj
# Prkcdbp 2.493809e-20   4.002253    0.106  0.015  3.461906e-16
# Ckap4   1.110139e-16   4.426637    0.086  0.012  1.541095e-12
# Bok     3.351234e-16   4.667646    0.056  0.005  4.652183e-12
# Prlr    5.225256e-16   4.552168    0.040  0.002  7.253700e-12
# Tpm2    8.722613e-16   3.895989    0.126  0.026  1.210873e-11

cluster9.markers <- FindMarkers(NK_TCell, ident.1 = 9, only.pos = TRUE)
head(cluster9.markers, n = 5)
#                p_val avg_log2FC pct.1 pct.2     p_val_adj
# Hivep3 4.407666e-123   3.757219 0.671 0.087 6.118722e-119
# Bdh2    1.685760e-86   5.110816 0.187 0.006  2.340172e-82
# Ptprs   4.864369e-70   3.093951 0.445 0.062  6.752718e-66
# Lag3    1.286647e-69   2.560693 0.652 0.145  1.786124e-65
# Cd200   1.312366e-68   3.019954 0.542 0.099  1.821827e-64

p <- DotPlot(
  NK_TCell,
  features = c("Cd4", "Cd8a", "Cd8b"), # Adjust CD8 names if needed
  idents = c("7", "8", "10")
)

# 2. Extract and print the exact data table
# 'avg.exp' is the non-scaled average expression
# 'pct.exp' is the percentage of cells expressing the gene
print(p$data)

test <- LoadSeuratRds(file = "C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/ProcessedSingleCellObjects/Archive/islet_graft_seurat_v1.rds")
DotPlot(
  test,
  features = unique(unlist(marker_genes))
) +
  RotatedAxis() +
  scale_color_gradient(low = "grey80", high = "red")
DimPlot(test, reduction = "umap",label = T,label.size = 4,label.box = T,pt.size = 1)


# Further Cluster Analysis for Annotation- end

# Annotation of Cells ----
cluster_ids_TNK <- c(
  "Treg",       # 0
  "Tconv",      # 1
  "CD8_T",      # 2
  "CD8_T",      # 3
  "CD8_T",      # 4
  "NK",         # 5
  "Tconv",      # 6
  "CD8_T",      # 7
  "Tconv",      # 8
  "CD8_T",      # 9
  "Tconv",      # 10
  "CD8_T",      # 11
  "Tconv"       # 12
)

# Apply the names
names(cluster_ids_TNK) <- levels(NK_TCell)
NK_TCell <- RenameIdents(NK_TCell, cluster_ids_TNK)
# Save to metadata for future plotting
NK_TCell$cell_type <- Idents(NK_TCell)

DimPlot(NK_TCell, reduction = "umap",label = T,label.size = 4,label.box = T,pt.size = 1)


# Merge NK_TCell with islet_graft_seurat for Final Annotation ----
# Create a copy of the original labels in the islet_graft_seurat object
islet_graft_seurat$final_annotation <- as.character(Idents(islet_graft_seurat))

# Match barcodes between the subset and the original object to determine which cells in islet_graft_seurat belong to sub-clustered NK_TCell object
matching_cells <- match(colnames(NK_TCell), colnames(islet_graft_seurat))

# Overwrite the "NK/TCell" labels with subset annotations
islet_graft_seurat$final_annotation[matching_cells] <- as.character(Idents(NK_TCell))

# Set the new identities as the default for plotting
Idents(islet_graft_seurat) <- "final_annotation"

# Verify the new list of cells
table(Idents(islet_graft_seurat))

# Plot the fully annotated global dataset
DimPlot(islet_graft_seurat, reduction = "umap",label = T,label.size = 4,label.box = T,pt.size = 1) +
ggtitle("Annotated Single Cell Data")
DimPlot(islet_graft_seurat, reduction = "umap",label = T,label.size = 6,label.box = T,pt.size = 1,repel =T) 
  
# Save final version
saveRDS(islet_graft_seurat, "C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/ProcessedSingleCellObjects/islet_graft_seurat_v7.rds")
