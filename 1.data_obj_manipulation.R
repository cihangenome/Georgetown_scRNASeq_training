library(Seurat)
library(tidyverse)
library(matrixStats)

# read in the downsampled seurat object
# data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161918
seurat_obj_ds <- readRDS("GSE161918_ds_obj.rds")

# assays in the data
seurat_obj_ds@assays

head(seurat_obj_ds@meta.data)

table(seurat_obj_ds@meta.data$Class)
table(seurat_obj_ds@meta.data$severity.outcome)

cellid <- sample(colnames(seurat_obj_ds), size = 100, replace = FALSE)
seurat_obj_select <- subset(seurat_obj_ds, cells = cellid)
Mono_Classical <- subset(seurat_obj_ds, subset = celltype == "Mono_Classical")

# clustering on surface protein data
DefaultAssay(seurat_obj_ds) <- "limmaCITE"
seurat_obj_ds <- FindVariableFeatures(seurat_obj_ds, nfeatures = 100)
seurat_obj_ds <- ScaleData(seurat_obj_ds)
seurat_obj_ds <- RunPCA(seurat_obj_ds, reduction.name = "cite.pca", reduction.key = "CITE.PC_")
seurat_obj_ds <- FindNeighbors(seurat_obj_ds, reduction = "cite.pca", dims = 1:10)
seurat_obj_ds <- FindClusters(seurat_obj_ds, resolution = c(0.5))
seurat_obj_ds <- RunUMAP(seurat_obj_ds, 
                         features = seurat_obj_ds[["CITE"]]@var.features,
                         reduction.name = "citeumap", n.neighbors=30L, min.dist = 0.3)

# visualize clusters
DimPlot(seurat_obj_ds, pt.size = 0.2, reduction = "citeumap", label=TRUE)
DimPlot(seurat_obj_ds, pt.size = 0.2, reduction = "citeumap", label=TRUE, group.by = "celltype")+
  NoLegend()


saveRDS(seurat_obj_ds, "GSE161918_ds_obj.rds")


