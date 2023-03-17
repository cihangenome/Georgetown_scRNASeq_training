library(Seurat)
library(tidyverse)
library(matrixStats)
library(SingleR)
library(celldex)
library(pheatmap)

# Annotation by findMarkers
seurat_obj_ds <- readRDS("GSE161918_ds_obj.rds")

Idents(seurat_obj_ds) <- "CITE_snn_res.0.5"
cite_markers <- FindAllMarkers(seurat_obj_ds, assay = "CITE", only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
cite_markers_selected <- cite_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = desc(p_val))
View(cite_markers_selected)
p1 <- DimPlot(seurat_obj_ds, pt.size = 0.2, reduction = "citeumap", label=TRUE, group.by = "CITE_snn_res.0.5", label.size = 3)+
  NoLegend()
p2 <- DimPlot(seurat_obj_ds, pt.size = 0.2, reduction = "citeumap", label=TRUE, group.by = "celltype", label.size = 3)+
  NoLegend()
p1+p2

DefaultAssay(seurat_obj_ds) <- "CITE"
DoHeatmap(subset(seurat_obj_ds, downsample = 100), features = cite_markers_selected$gene)


# using some other packages
seurat_obj_ds_sce <- as.SingleCellExperiment(seurat_obj_ds, assay = "RNA")
ref <- HumanPrimaryCellAtlasData()

pred <- SingleR(test=seurat_obj_ds_sce, ref=ref, labels=ref$label.main)

table(pred$labels)
seurat_obj_ds$pred_label <- pred$labels

celltype_overlap <- table(predicted=seurat_obj_ds$pred_label, annotated=seurat_obj_ds$celltype)

pheatmap(log2(celltype_overlap+10), color=colorRampPalette(c("white", "blue"))(101))


