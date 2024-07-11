# 
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)


# 
deg = read.csv("markers.csv", header = T, row.names = 1)
head(deg)


tecs = subset(deg, cluster %in% 0)
head(tecs)

# 
gene = tecs$gene
length(gene) # 1069


gene_id = bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")$ENTREZID
length(gene_id) # 1023

bp = enrichGO(
  gene_id,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "bp",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2)
# 
write.csv(bp, "bp.csv")

# 
barplot(bp)



bp = read.csv("bp.csv", header = T, row.names = 1)
head(bp)

bp = arrange(bp, -Count)
head(bp)

ggplot(bp, aes(reorder(Description,Count),Count,color = pvalue)) +
  geom_point() +
  theme_bw() +
  coord_flip() +
  scale_color_gradient(low = "#ff00b2", high = "#10d0ee")




















