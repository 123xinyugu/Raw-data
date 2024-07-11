##############AUCell#################################
# 
library(tidyverse)
library(dplyr)
library(KEGGREST)
library(clusterProfiler)
library(org.Hs.eg.db)

kegg_hsa = "hsa04010"
# 
kegg_gene = keggGet(kegg_hsa)
gene_id = kegg_gene[[1]]$GENE
gene_id
length(gene_id)
gene = bitr(gene_id, fromType = "ENTREZID",toType = "SYMBOL", OrgDb = "org.Hs.eg.db", drop = T)$SYMBOL
length(gene) # 
gene


# 
library(tidyverse)
library(Seurat)
library(cols4all)
library(glmGamPoi) # 
library(cowplot)
library(clustree)
library(AUCell)
library(clusterProfiler) # 
library(ggsignif)
library(GSVA)

# 
scRNA = readRDS("final_scRNA.rds")
table(scRNA@meta.data$cell_type)

#Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(scRNA@assays$RNA@data))
gene = list("MAPK signaling pathway" = gene)
head(gene)

#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(gene, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)
aucs <- getAUC(cells_AUC)
head(aucs)
aucs = t(aucs) %>%
  as.data.frame() %>%
  write.csv("aucs.csv")

aucs = read.csv("aucs.csv", header = T)
head(aucs)

# 
cell_info = dplyr::select(scRNA@meta.data, cell_type)
cell_info$X = rownames(cell_info)
head(cell_info)

aucs1 = merge(aucs, cell_info, by = "X")
head(aucs1)

colors = c("#ff85a8","#ffd685","#fff385","#85ffd0","#85c2ff","#ff85f3")
ggplot(aucs1, aes(reorder(cell_type, -Foxo.signaling.pathway), Foxo.signaling.pathway, fill = cell_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = "", y = "MAPK signaling pathway")

