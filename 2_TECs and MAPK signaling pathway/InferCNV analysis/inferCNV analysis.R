library(Seurat)
library(tidyverse)
library(dplyr)
library(infercnv)

# 
genes = read.table("infercnv_genes.txt", header = T, row.names = 1, sep = "\t")
head(genes)


scRNA = readRDS("final_scRNA.rds")
table(scRNA@meta.data$cell_type)
mydata = subset(scRNA, cell_type %in% c("B cells","Tumor-associated endothelial cells"))

data = as.matrix(mydata@assays$RNA@counts)
write.table(data, "data.txt", col.names=T, row.names=T, quote=F, sep="\t")

annotation = subset(mydata@meta.data, select=c("cell_type"))
write.table(annotation, "annotation.txt", col.names=T, row.names=T, quote=F, sep="\t")
table(annotation$cell_type)

#################################################################################################
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=data,
                                    annotations_file=annotation, delim="\t",
                                    gene_order_file=genes,
                                    ref_group_names=c("B cells")
)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="./try",
                             cluster_by_groups=TRUE, 
                             analysis_mode="subcluster",
                             HMM_type="i3",
                             denoise=TRUE,
                             HMM_report_by="subcluster",
                             HMM=TRUE,
                             output_format="pdf",
                             no_prelim_plot=T
)





##########################################################################
genes_dat = read.delim("./try/HMM_CNV_predictions.HMMi3.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat", header = T, sep = "\t")
head(genes_dat)

cell_number = read.delim("./try/17_HMM_predHMMi3.leiden.hmm_mode-subclusters.cell_groupings", header = T, sep = "\t")
head(cell_number)


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
genes = bitr(gene_id, fromType = "ENTREZID",toType = "SYMBOL", OrgDb = "org.Hs.eg.db", drop = T)$SYMBOL
length(genes) # 
genes

# write.table(genes,"genes of mapk signaling pathway.txt", row.names = F, sep = "\t", quote = F)

dim(genes_dat) # 135019      7

# 
mapk_gene_dat = subset(genes_dat, gene %in% genes)
head(mapk_gene_dat)
dim(mapk_gene_dat) # 2955    7


# 
tumor = mapk_gene_dat$cell_group_name[grep("^Tumor", mapk_gene_dat$cell_group_name)]
filter_data = subset(mapk_gene_dat,cell_group_name %in% tumor)
dim(filter_data) # 2483    8
head(filter_data)

count = filter_data %>% 
  group_by(gene,state) %>% 
  count() %>% 
  group_by(gene) %>% 
  mutate("percent" = n/sum(n) * 100) %>% 
  as.data.frame()
head(count)
count$group = count$gene
count <- count %>%
  mutate(group = case_when(
    state == 1 ~ "Deletion",
    state == 2 ~ "Normal",
    state == 3 ~ "Amplification",
    TRUE ~ count$group
  ))
head(count)
summary(count)

amp = subset(count, state == 3)
p1 = ggplot(amp,aes(reorder(gene,-percent),percent,fill = percent)) +
  geom_col(width = 0.5) +
  theme_bw() +
  scale_fill_gradient(low = "#10d0ee", high = "#ff00b2") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  labs(title = "Amplification",x = "",y = "Percent(%)");p1
  
del = subset(count, state == 1)
p2 = ggplot(del,aes(reorder(gene,-percent),percent,fill = percent)) +
  geom_col(width = 0.5) +
  theme_bw() +
  scale_fill_gradient(low = "#10d0ee", high = "#ff00b2") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  labs(title = "Deletion",x = "",y = "Percent(%)");p2

library(cowplot)
plot_grid(p1,p2,ncol = 1, labels = c("A","B"))

# 
# write.table(count, "mapk_gene_CNV_count.txt", row.names = F, sep = "\t", quote = F)
