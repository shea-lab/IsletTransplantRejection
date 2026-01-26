library(dplyr)
library(Seurat)
library(patchwork)
library(BiocManager)
#install.packages('BiocManager')
#BiocManager::install('glmGamPoi')
library(glmGamPoi)
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

# SCT transform to normalize data
#allo_seurat <- SCTransform(allo_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
#syn_seurat <- SCTransform(syn_seurat, vars.to.regress = "percent.mt", verbose = TRUE)

# Create Islet Seurat Object----
islet_graft_seurat <- merge(allo_seurat, y = syn_seurat, add.cell.ids = c("Allogeneic", "Syngeneic"))
DefaultAssay(islet_graft_seurat)<-"RNA"
islet_graft_seurat<- SCTransform(islet_graft_seurat, vars.to.regress = "percent.mt", verbose = TRUE)
#islet_graft_seurat<- FindVariableFeatures(islet_graft_seurat)
islet_graft_seurat <- RunPCA(islet_graft_seurat, features = VariableFeatures(object = islet_graft_seurat))

print(islet_graft_seurat[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(islet_graft_seurat, reduction = "pca") + NoLegend()
ElbowPlot(islet_graft_seurat)
islet_graft_seurat <- FindNeighbors(islet_graft_seurat, dims = 1:10)
islet_graft_seurat <- FindClusters(islet_graft_seurat, resolution = 0.25)
islet_graft_seurat <- RunUMAP(islet_graft_seurat, dims = 1:10)
islet_graft_seurat$condition <- islet_graft_seurat$orig.ident
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p1<-DimPlot(islet_graft_seurat, reduction = "umap")
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p2<-DimPlot(islet_graft_seurat, reduction = "umap",group.by ="condition" )
saveRDS(islet_graft_seurat, "C:/Users/17343/Desktop/IsletTransplantRejection/IsletSingleCellData/ProcessedSingleCellObjects/islet_graft_seurat_v1.rds")
