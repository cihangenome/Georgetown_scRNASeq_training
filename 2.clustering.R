library(Seurat)
library(tidyverse)
library(matrixStats)
library(dsb)

# read in the downsampled seurat object
seurat_obj_ds <- readRDS("GSE161918_ds_obj.rds")

# before dsb
Idents(seurat_obj_ds) <- "orig.ident"
DefaultAssay(seurat_obj_ds) <- "CITE"
FeatureScatter(seurat_obj_ds, feature1="cite_CD4", feature2="cite_CD8", slot = "counts")
# CLR normalization
seurat_obj_ds <- NormalizeData(seurat_obj_ds, normalization.method = "CLR", margin = 2, assay = "CITE")
p_clr <- FeaturePlot(seurat_obj_ds, reduction = "citeumap",
                     features = c( "cite_CD4", "cite_CD16", "cite_CD278", "cite_CD123"), 
                     min.cutoff = "q05", max.cutoff = "q95", ncol = 1, pt.size = 0.05)
p_clr

# dsb
VlnPlot(seurat_obj_ds, features = c( "nCount_CITE"), ncol = 3, pt.size =0, group.by = "orig.ident") + geom_hline(yintercept = c(2500))
quantile(seurat_obj_ds$nCount_CITE, probs = c(0, .1, .2,.25, .5, .75, .95, .99, .999, .9999, 1))

NegativesObj <- subset(seurat_obj_ds, subset = nCount_CITE < 2500 & nCount_RNA < 3000)
isotype.control.name <- c("IgG2aKiso", "IgG2bKiso", 
                          "ratIgG2bKiso", "IgG1Kiso")
background_cite_matrix <- as.matrix(GetAssayData(NegativesObj[["CITE"]], slot = "counts"))
raw_cite_matrix <- as.matrix(GetAssayData(seurat_obj_ds[["CITE"]], slot = "counts"))
dsb_norm_cite <- DSBNormalizeProtein(cell_protein_matrix = raw_cite_matrix,
                                     empty_drop_matrix = background_cite_matrix,
                                     denoise.counts = TRUE, 
                                     use.isotype.control = TRUE,
                                     isotype.control.name.vec = isotype.control.name)

# replace the data slot for assay:CITE
seurat_obj_ds[["CITE"]] <- SetAssayData(seurat_obj_ds[["CITE"]], slot = "data", new.data = dsb_norm_cite)

FeatureScatter(seurat_obj_ds, feature1="cite_CD4", feature2="cite_CD8", slot = "data")

p_dsb <- FeaturePlot(seurat_obj_ds, reduction = "citeumap",
                     features = c( "cite_CD4", "cite_CD16", "cite_CD278", "cite_CD123"), 
                     min.cutoff = "q05", max.cutoff = "q95", ncol = 1, pt.size = 0.05)
CombinePlots(list(p_clr, p_dsb))


# cluster on RNA
DefaultAssay(seurat_obj_ds) <- "RNA"
seurat_obj_ds <- NormalizeData(seurat_obj_ds)
seurat_obj_ds <- FindVariableFeatures(seurat_obj_ds, selection.method = "vst", nfeatures = 2000)
seurat_obj_ds <- ScaleData(seurat_obj_ds, assay = "RNA")
seurat_obj_ds <- RunPCA(seurat_obj_ds, features = VariableFeatures(object = seurat_obj_ds))
ElbowPlot(seurat_obj_ds, reduction = "pca")
seurat_obj_ds <- FindNeighbors(seurat_obj_ds, reduction = "pca", dims = 1:10, nn.eps = 1)
seurat_obj_ds <- FindClusters(seurat_obj_ds, resolution = c(0.5))
seurat_obj_ds <- RunUMAP(seurat_obj_ds, dims = 1:10, n.neighbors=30L, min.dist = 0.3, reduction = "pca", reduction.name = "rnaumap")

DimPlot(seurat_obj_ds, pt.size = 0.2, reduction = "rnaumap", label=TRUE)
DimPlot(seurat_obj_ds, pt.size = 0.2, reduction = "rnaumap", label=TRUE, group.by = "celltype")+
  NoLegend()

# compare with clusters using cell surface protein
p1 <- FeaturePlot(seurat_obj_ds, reduction = "rnaumap",
                  features = c("cite_CD4", "cite_CD8"), 
                  min.cutoff = "q05", max.cutoff = "q95", ncol = 2, pt.size = 0.05)
p2 <- FeaturePlot(seurat_obj_ds, reduction = "citeumap",
                  features = c("cite_CD4", "cite_CD8"), cols = c("lightgrey","darkgreen"),
                  min.cutoff = "q05", max.cutoff = "q95", ncol = 2, pt.size = 0.05)

p1/p2

# Wnn -- combine two modalities
seurat_obj_ds <- FindMultiModalNeighbors(
  seurat_obj_ds, reduction.list = list("cite.pca", "pca"), 
  dims.list = list(1:15, 1:15), modality.weight.name = "RNA.weight"
)

seurat_obj_ds <- RunUMAP(seurat_obj_ds, nn.name = "weighted.nn", reduction.name = "wnnumap", reduction.key = "wnnUMAP_")
seurat_obj_ds <- FindClusters(seurat_obj_ds, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)

DimPlot(seurat_obj_ds, reduction = 'wnnumap', label = TRUE, repel = TRUE) + NoLegend()
DimPlot(seurat_obj_ds, pt.size = 0.2, reduction = "wnnumap", label=TRUE, group.by = "celltype", label.size = 3)+
  NoLegend()

VlnPlot(seurat_obj_ds, features = "RNA.weight", group.by = 'celltype', sort = TRUE, pt.size = 0.1) +
  NoLegend()



