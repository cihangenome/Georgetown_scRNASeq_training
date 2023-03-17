library(Seurat)
library(MAST)
library(tidyverse)
library(matrixStats)

# Differential expression gene analysis at single-cell level
### COVID vs HC
# use Classical monocyte as an example
### using seurat ##############

seurat_obj_ds <- readRDS("GSE161918_ds_obj.rds")

Mono_Classical <- subset(seurat_obj_ds, subset = WCTcoursecelltype == "Mono_Classical")
Idents(Mono_Classical) <- "Class"
markers.COVIDvsHC <- FindMarkers(Mono_Classical, assay = "RNA", ident.1 = "COVID", ident.2 = "HC",
                                 logfc.threshold = 0)

markers.COVIDvsHC.up <- markers.COVIDvsHC %>%
  rownames_to_column("gene") %>%
  filter(avg_log2FC>0, p_val_adj < 0.05)

markers.COVIDvsHC.dn <- markers.COVIDvsHC %>%
  rownames_to_column("gene") %>%
  filter(avg_log2FC<0, p_val_adj < 0.05)

Mono_Classical_avg <- data.frame(log1p(AverageExpression(Mono_Classical, verbose = FALSE)$RNA))
Mono_Classical_avg$gene <- rownames(Mono_Classical_avg)

label1 = paste("COVIDvsHC", sep = ":")
label2 = paste("# of up DE genes", nrow(markers.COVIDvsHC.up), sep = ":")
label3 = paste("# of dn DE genes", nrow(markers.COVIDvsHC.dn), sep = ":")

genes.to.label1 = markers.COVIDvsHC.up$gene
genes.to.label2 = markers.COVIDvsHC.dn$gene

p <- ggplot(Mono_Classical_avg, aes(x = HC, y = COVID)) + 
  geom_point(alpha = 0.2, color="#1b262c") + 
  ggtitle(paste(label1, label2, label3, sep = "\n")) + 
  geom_point(data = Mono_Classical_avg[genes.to.label1, ], color="#fe91ca") +
  geom_point(data = Mono_Classical_avg[genes.to.label2, ], color="#2bb2bb") +
  theme_bw()
p <- LabelPoints(plot = p, points = c(genes.to.label1[1:30],genes.to.label2[1:25]), 
                 repel = TRUE, xnudge = 0, ynudge = 0, size = 2.5)
p

### using MAST ##############
# use MAST linear model and considering other covariates
# Need to Convert objects to SingleCellExperiment objects
Mono_Classical_sce <- as.SingleCellExperiment(Mono_Classical, assay = "RNA")

cdr <-colSums(assay(Mono_Classical_sce)>0)
colData(Mono_Classical_sce)$cngeneson <- scale(cdr)[,1]
Mono_Classical_sca <- SceToSingleCellAssay(Mono_Classical_sce, class = "SingleCellAssay")

## Remove invariant genes and select a subset of genes for illustration
sca <- Mono_Classical_sca[sample(which(freq(Mono_Classical_sca)>0), 4000),]

lmer <- zlm(~ Class + cngeneson + Age + Gender + days_of_symptoms_onset, sca, silent = TRUE)
show(lmer)
summaryCond <- summary(lmer, doLRT='ClassCOVID') 
summarydf <- summaryCond$datatable
foldchange_res <- merge(summarydf[contrast=='ClassCOVID' & component=='H',.(primerid, `Pr(>Chisq)`)], # P values
                  summarydf[contrast=='ClassCOVID' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)],  #logFC coefficients
                  by='primerid')
foldchange_res[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
View(foldchange_res %>% filter(fdr < 0.01, coef > 0))




