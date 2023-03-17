library(fgsea)
library(tidyverse)

# GESA -- geneset enrichment analysis on the DE results

# read in genesets
hallmark_genesets <- readRDS("HALLMARK_genesets.rds")

# read in previous results
mono_classical_de <- read_tsv("Mono_Classical--model@healthy_vs_covid--coef@COVID-Healthy--toptab.tsv")
mono_classical_de
mono_classical_de$rank <- -log10(mono_classical_de$P.Value)
mono_classical_de$rank <- ifelse(mono_classical_de$logFC > 0, mono_classical_de$rank, -mono_classical_de$rank)
rank <- mono_classical_de$rank
names(rank) <- mono_classical_de$gene

Mono_fgseaRes <- fgsea(pathways = hallmark_genesets, 
                       stats = rank,
                       minSize=10,
                       maxSize=610)
topPathwaysUp <- Mono_fgseaRes[ES > 0][head(rev(order(NES)), n=20), pathway]
topPathwaysDown <- Mono_fgseaRes[ES < 0][head(order(NES), n=20), pathway]

topPathwaysUp
View(Mono_fgseaRes)

plotEnrichment(hallmark_genesets[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
               rank) + labs(title="HALLMARK_INTERFERON_ALPHA_RESPONSE")

Mono_fgseaRes_filt <- Mono_fgseaRes %>% 
  filter( abs(NES) > 1 )  %>%
  mutate(n_logp = -log10(padj),
         celltype = "Mono_Classical")

p <- ggplot(Mono_fgseaRes_filt, aes(y = celltype, x = pathway, color = NES, size = n_logp)) + 
  geom_point() +
  scale_color_gradient2(low = "skyblue", mid = "white", high = "red", midpoint = 0) +
  theme_bw() +
  scale_x_discrete(position = "top") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
  theme(axis.title.y = element_blank()) +
  labs(color = 'Normalized Enrichment Score', size = '-log10(adjusted P value)')

p


# plot with other celltypes together
CD4_Mem_de <- read_tsv("CD4_Mem--model@healthy_vs_covid--coef@COVID-Healthy--toptab.tsv")
CD4_Mem_de$rank <- -log10(CD4_Mem_de$P.Value)
CD4_Mem_de$rank <- ifelse(CD4_Mem_de$logFC > 0, CD4_Mem_de$rank, -CD4_Mem_de$rank)
rank <- CD4_Mem_de$rank
names(rank) <- CD4_Mem_de$gene

CD4_fgseaRes <- fgsea(pathways = hallmark_genesets, 
                      stats = rank,
                      minSize=10,
                      maxSize=610)

CD4_fgseaRes_filt <- CD4_fgseaRes %>% 
  filter( abs(NES) > 1 )  %>%
  mutate(n_logp = -log10(padj),
         celltype = "CD4_Mem")

fgseaRes_filt <- rbind(CD4_fgseaRes_filt, Mono_fgseaRes_filt)
p <- ggplot(fgseaRes_filt, aes(y = celltype, x = pathway, color = NES, size = n_logp)) + 
  geom_point() +
  scale_color_gradient2(low = "skyblue", mid = "white", high = "red", midpoint = 0) +
  theme_bw() +
  scale_x_discrete(position = "top") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
  theme(axis.title.y = element_blank()) +
  # theme(legend.position = "bottom") + 
  labs(color = 'Normalized Enrichment Score', size = '-log10(adjusted P value)')

p




