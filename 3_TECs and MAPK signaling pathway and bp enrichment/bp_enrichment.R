# 
library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(cols4all)
library(glmGamPoi) # 
library(cowplot)
library(clustree)

# 
scRNA = readRDS("final_scRNA.rds")
table(scRNA@meta.data$cell_type)

tecs = subset(scRNA, cell_type == "Tumor-associated endothelial cells")
tecs@meta.data$cell_type = droplevels(tecs@meta.data$cell_type)
table(tecs@meta.data$cell_type)
# saveRDS(tecs, "tecs.rds")


rm(list = ls())
scRNA = readRDS("tecs.rds")
# 
colors = c4a("poly.alphabet2",26)

# 
scRNA = SCTransform(scRNA, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = F)
# 
# saveRDS(scRNA, file = "./sct_scRNA.rds")


# 
scRNA = readRDS("sct_scRNA.rds")
#
scRNA <- RunPCA(scRNA)
colnames(scRNA@meta.data)[1] = c("patient_ID", "nCount_RNA", "nFeature_RNA", "percent.mt") #
# 
scRNA <- RunHarmony(scRNA, group.by.vars="patient_ID", max.iter.harmony=20, lambda=0.5, assay.use = "SCT")
# 
ElbowPlot(scRNA,ndims = 50)
# 
scRNA <- FindNeighbors(scRNA, dims=1:15, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:15, reduction="harmony")
# 
p3 = DimPlot(scRNA, reduction="umap", group.by="patient_ID", pt.size=1, cols = colors)+
  theme(legend.position="right", plot.title=element_blank())
p3

# 
obj = FindClusters(scRNA, resolution = seq(0.1, 1,by=0.1))
p4 = clustree(obj)
p4

# resolution = 0.1
scRNA <- FindClusters(scRNA, resolution=0.1)
p5 = UMAPPlot(scRNA, pt.size=1, label=T)+NoLegend()
p5

# 
scRNA.markers <- FindAllMarkers(scRNA, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# write.csv(scRNA.markers, "tecs_markers.csv")

cell_label = c("TECs 1", "TECs 2")
## 
names(cell_label) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, cell_label)
scRNA[["cell_type"]] = Idents(scRNA)

# 
colors = c("#e22b21","#11e2ff")
p8 = UMAPPlot(scRNA, pt.size=1, label=T, label.size=5, cols = colors)+
  NoLegend()
p8

# 
cell_count1 = scRNA@meta.data %>%
  group_by(cell_type) %>%
  count() %>%
  mutate(Percent = n/sum(.$n))
cell_count1
p12 = ggplot(cell_count1, aes(cell_type,Percent *100, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(Percent * 100,2))), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "no") +
  labs(title = "Cell Ration", x = "", y = "Ratio (%)") +
  scale_fill_manual(values = colors)
p12

# 
# saveRDS(scRNA, "final_tecs.rds")


# 
library(clusterProfiler)
library(org.Hs.eg.db)

deg = read.csv("tecs_markers.csv", header = T, row.names = 1)
head(deg)

# 
tecs1_gene = subset(deg,cluster == 0)$gene
tecs2_gene = subset(deg, cluster == 1)$gene

tecs1_gene_id = bitr(tecs1_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
tecs2_gene_id = bitr(tecs2_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID

bp_tecs1 = enrichGO(tecs1_gene_id,
                    OrgDb = "org.Hs.eg.db",
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2)

bp_tecs2 = enrichGO(tecs2_gene_id,
                    OrgDb = "org.Hs.eg.db",
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2)
barplot(bp_tecs1)
barplot(bp_tecs2)

# write.csv(bp_tecs1,"bp_tecs1.csv")
# write.csv(bp_tecs2,"bp_tecs2.csv")


# 
bp_tecs1 = read.csv("bp_tecs1_select.csv", header = T, row.names = 1)
bp_tecs2 = read.csv("bp_tecs2_select.csv", header = T, row.names = 1)
head(bp_tecs1)
head(bp_tecs2)

p13 = ggplot(bp_tecs1, aes(reorder(Description,Count), Count, fill = pvalue)) +
  geom_col() +
  theme_bw() +
  scale_fill_gradient(low = "#e22b21",high = "#11e2ff") +
  coord_flip()

p14 = ggplot(bp_tecs2, aes(reorder(Description,Count), Count, fill = pvalue)) +
  geom_col() +
  theme_bw() +
  scale_fill_gradient(low = "#e22b21",high = "#11e2ff") +
  coord_flip()
p13;p14
plot_grid(p13,p14, ncol = 2)





# hallmark ssGSEA
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
scRNA = readRDS("final_tecs.rds")
#Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(scRNA@assays$RNA@data))

# 
hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H")
head(hallmark_genes)
hallmark = hallmark_genes %>%
  split(x = .$gene_symbol, f = .$gs_name)
head(hallmark)

#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(hallmark, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)
aucs <- getAUC(cells_AUC)
# t(aucs) %>%
#   as.data.frame() %>%
#   write.csv("aucs.csv")

# 
aucs = read.csv("aucs.csv", header = T, row.names = 1)
aucs$samples = rownames(aucs)
head(aucs)
head(scRNA@meta.data)
cell_info = dplyr::select(scRNA@meta.data, cell_type)
cell_info$samples = rownames(cell_info)
head(cell_info)

data = merge(aucs, cell_info, by = "samples")
head(data)
dim(data)

data_long = pivot_longer(data, cols = colnames(data)[2:51])
head(data_long)

ggplot(data_long, aes(reorder(name,value), value, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("#e22b21","#11e2ff")) +
  stat_compare_means(aes(group = cell_type),label = "p.signif",size=3,method = "t.test")


head(data)
aucs_tec1 = subset(data, cell_type == "TECs 1")
aucs_tec2 = subset(data, cell_type == "TECs 2")
rownames(aucs_tec1) = aucs_tec1$samples
rownames(aucs_tec2) = aucs_tec2$samples

aucs_tec1 = dplyr::select(aucs_tec1, -c(samples,cell_type))
head(aucs_tec1)

aucs_tec2 = dplyr::select(aucs_tec2, -c(samples,cell_type))
head(aucs_tec2)

aucs_tec1_mean = colMeans(aucs_tec1) %>% 
  as.data.frame()
aucs_tec2_mean = colMeans(aucs_tec2) %>% 
  as.data.frame()
names(aucs_tec1_mean) = c("tecs1_mean")
names(aucs_tec2_mean) = c("tecs2_mean")

aucs_tec1_mean$hallmark = rownames(aucs_tec1_mean)
aucs_tec2_mean$hallmark = rownames(aucs_tec2_mean)
head(aucs_tec1_mean)
head(aucs_tec2_mean)

df = merge(aucs_tec1_mean,aucs_tec2_mean)
head(df)
rownames(df) = df$hallmark
df = df %>% 
  dplyr::select(-hallmark)
head(df)

library(pheatmap)
library(RColorBrewer)
pheatmap(df,
         colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         scale = "row",
         cluster_cols = F,
         angle_col = 45)






##############MAPK AUCells############################################
# Foxo
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
genes = bitr(gene_id, fromType = "ENTREZID",toType = "SYMBOL", OrgDb = "org.Hs.eg.db", drop = T)$SYMBOL
length(genes) # 
genes
genes = list("MAPK signaling pathway" = genes)
head(genes)

# AUCells
#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(genes, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)
aucs <- getAUC(cells_AUC)
t(aucs) %>%
  as.data.frame() %>%
  write.csv("mapk_aucs.csv")

# 
aucs = read.csv("mapk_aucs.csv", header = T, row.names = 1)
aucs$samples = rownames(aucs)
head(aucs)
head(scRNA@meta.data)
cell_info = dplyr::select(scRNA@meta.data, cell_type)
cell_info$samples = rownames(cell_info)
head(cell_info)

data = merge(aucs, cell_info, by = "samples")
head(data)

ggplot(data, aes(cell_type, MAPK.signaling.pathway, fill = cell_type)) +
  geom_boxplot() +
  theme_bw() +
  stat_compare_means(aes(group = cell_type),label = "p.format",size=3,method = "t.test")




######################MAPK############################
data = read.table("mapk_gene_CNV_count.txt", header = T, sep = "\t") 
dim(data) # 19  5

table(data$group)
amp_gene = subset(data, group == "Amplification" & percent == 100)$gene
del_gene = subset(data, group == "Deletion" & percent == 100)$gene

colors = c("#e22b21","#11e2ff")
genes = c("CDC42","EFNA1","HSPB1","MAP3K11","RAC1","RAP1A","STMN1","TNFRSF1A","MAP4K4","MAPK3","PDGFD","TAOK2")
p15 = VlnPlot(scRNA, features = genes, group.by = "cell_type", col = colors, pt.size = 0, ncol = 4)
colors = c("#11e2ff","#e22b21")
p16 = DotPlot(scRNA, features = genes, group.by = "cell_type", cols = colors) +
  coord_flip()


plot_grid(p15,p16, ncol = 2, labels = c("A","B"))


head(data)
head(deg)
library(ggvenn)
intersect_gene = intersect(data$gene, deg$gene)
intersect_gene
ggvenn(list("MAPK gene" = data$gene, "TECs markers" = deg$gene))
colors = c("#e22b21","#11e2ff")
p17 = VlnPlot(scRNA, features = intersect_gene, group.by = "cell_type", col = colors, pt.size = 0, ncol = 5)
colors = c("#11e2ff","#e22b21")
p18 = DotPlot(scRNA, features = intersect_gene, group.by = "cell_type", cols = colors) +
  coord_flip()
plot_grid(p17,p18,ncol = 2, labels = c("A","B"))

head(data)
intersect_gene_table = subset(data, gene %in% intersect_gene) %>% 
  arrange(state)
View(intersect_gene_table)


ggplot(intersect_gene_table, aes(reorder(gene, percent), percent, fill = percent)) +
  geom_col()+
  theme_bw() +
  scale_fill_gradient(low = "#11e2ff",high = "#e22b21")+
  coord_flip()+
  facet_grid(.~group)
# write.csv(intersect_gene_table,"intersect_gene_table.csv")
