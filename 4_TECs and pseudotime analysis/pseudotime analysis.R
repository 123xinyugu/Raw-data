library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
library(monocle) #
library(cols4all)

# 
scRNA =  read_rds("final_tecs.rds")
head(scRNA)

# 
expr_matrix = as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
p_data = scRNA@meta.data
p_data$cell_type = scRNA@active.ident # 
f_data = data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
dim(expr_matrix)
dim(p_data)
dim(f_data)

# 
pd = new('AnnotatedDataFrame', data = p_data)
fd = new('AnnotatedDataFrame', data = f_data)

# 
cds = newCellDataSet(expr_matrix,
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

# 
cds = estimateSizeFactors(cds) # 
cds = estimateDispersions(cds)

# 
cds = detectGenes(cds, min_expr = 0.1) #num_cells_expressed
print(head(fData(cds)))
expressed_genes = row.names(subset(fData(cds),
                                   num_cells_expressed >= 10))

# 
# 
# 
head(scRNA)
# Idents(scRNA) = scRNA@meta.data$Condition # 
deg.cluster = FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
head(deg.cluster)
express_genes = rownames(subset(deg.cluster, p_val_adj < 0.05, abs(avg_log2FC) > 0.25))
head(express_genes)
cds = setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

# 
cds = reduceDimension(cds, max_components = 2, method = "DDRTrees")
# 
cds = orderCells(cds) # 

colors = c("#AEC7E8", "#98DF8A", "#FF9896", "#FFBB78","red")

# 
plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = T)
# 
p1 = plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = T) +
  scale_color_manual(values = colors) +
  labs(title = "TECs")
p1

# 
colors = c("#e22b21","#11e2ff")
p2 = plot_cell_trajectory(cds, color_by = "cell_type", size = 1, show_backbone = T) +
  scale_color_manual(values = colors) +
  labs(title = "TECs")
p2


# 
Time_diff = differentialGeneTest(cds[deg.cluster$gene], cores = 1,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes = Time_diff %>% pull(gene_short_name) %>%as.character()
Time_genes = unique(Time_genes)
genes = c("DUSP6","HSPA1A","HSPA1B","HSPB1","TGFBR2","TNFRSF1A","JUN","JUND","NR4A1","PDGFRB","PGF","RRAS","GADD45B")
p3 = plot_pseudotime_heatmap(cds[genes,], 
                             num_clusters = 2, 
                             show_rownames = T, 
                             return_heatmap = T)
p3

p3 = plot_pseudotime_heatmap(cds[Time_genes,], 
                             num_clusters = 4, 
                             show_rownames = F, 
                             return_heatmap = T)
p3

p3$tree_row
# 
clusters = cutree(p3$tree_row, k=4) # 
clustering = data.frame(clusters)
head(clustering)
clustering[,1] = as.character(clustering[,1])
colnames(clustering) = "Gene_Clusters"
table(clustering)
# write.csv(clustering, "Time_clustering_all.csv")



# 
library(clusterProfiler)
library(org.Hs.eg.db)

# 
bp_go <- function(gene) {
  enrichGO(
    gene,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
}

cluster4_gene = rownames(subset(clustering, Gene_Clusters == "4"))
# 
bp_cluster4 = bp_go(cluster4_gene)
# 
# write.table(bp_cluster4,"bp_cluster4.txt", sep = "\t")
