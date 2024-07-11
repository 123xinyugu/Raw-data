####1.1 
library(tidyverse)
library(dplyr)
library(TCGAbiolinks)

expr = read.delim("./TCGA-LIHC.htseq_fpkm.tsv/TCGA-LIHC.htseq_fpkm.tsv", quote = "\t")
probe = read.delim("./gencode.v22.annotation.gene.probeMap") %>% 
  select(c(id,gene))
head(probe)
expr[1:5,1:5]

expr = merge(expr, probe, by.x = "Ensembl_ID","id")
expr[1:5,1:5]
expr = expr %>% 
  select(-Ensembl_ID) %>% 
  select(gene, everything())
expr[1:3,1:5]
# 
expr_mean=aggregate(.~gene,mean,data=expr)
# 
# write.table(expr_mean, "expr.txt", row.names = F, quote = F, sep = "\t")

# 
pd = read.delim("./TCGA-LIHC.GDC_phenotype.tsv/TCGA-LIHC.GDC_phenotype.tsv", quote = "\t")
colnames(pd)
View(pd)

# 
pd_info = select(pd, c(submitter_id.samples,pathologic_M,sample_type.samples))
head(pd_info)
table(pd_info$pathologic_M)
# M0  M1  MX 
# 337   6 126
table(pd_info$sample_type.samples)
# Primary Tumor     Recurrent Tumor Solid Tissue Normal 
# 377                   3                  89 

expr = read.table("expr.txt", header = T, row.names = 1, sep = "\t",)
expr[1:5,1:5]
barcode = gsub("\\.","-",colnames(expr))
barcode
length(barcode) # 424
# 
# ；
# 
samplesTP <- TCGAquery_SampleTypes(barcode = barcode,
                                   typesample = "TP")
length(samplesTP) # 371
####
samplesNT <- TCGAquery_SampleTypes(barcode = barcode,
                                   typesample = "NT")
length(samplesNT) # 50

samples = gsub("-",".",c(samplesTP, samplesNT))
samples
expr1 = dplyr::select(expr, samples)
dim(expr1) # 
# 
# write.table(expr1, "filter_expr.txt", row.names = F, quote = F, sep = "\t")
# save(samplesTP, samplesNT, expr1, pd_info, file = "tcga_data.RData")

# 
load(file = "tcga_data.RData")
samplesTP = gsub("-",".",samplesTP)
samplesNT = gsub("-",".",samplesNT)






###############################################
library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
library(readxl)
library(cols4all)
library(tidyr)

deg = read.csv("markers.csv",header = T, row.names = 1)
# 
deg$cell_type[deg$cluster == 0] <- "TECs"
deg$cell_type[deg$cluster == 1] <- "Liver bud hepatic cells"
deg$cell_type[deg$cluster == 2] <- "T cells"
deg$cell_type[deg$cluster == 3] <- "Myofibroblast"
deg$cell_type[deg$cluster == 4] <- "B cells"
deg$cell_type[deg$cluster == 5] <- "TAM"
deg$cell_type[deg$cluster == 6] <- "T cells"
deg$cell_type[deg$cluster == 7] <- "B cells"

list = split(deg$gene, deg$cell_type)
# 
expr = log2(expr1+1)

re <- gsva(as.matrix(expr), list, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)
# head(re)
# re %>%
#   t() %>% 
#   as.data.frame() %>% 
#   write.table("re.txt", sep = "\t", quote = F)


ssgsea = read.table("re.txt", header = T, sep = "\t")
head(ssgsea)
ssgsea$sample = rownames(ssgsea)
ssgsea$group = ifelse(ssgsea$sample %in% samplesTP, "TP", "NT")
head(ssgsea)

dim(ssgsea)
colnames(ssgsea)[1:6]

ssgsea_long = pivot_longer(ssgsea, cols = colnames(ssgsea)[1:6], names_to = "cell_type")
head(ssgsea_long)

ggplot(ssgsea_long, aes(reorder(cell_type, -value), value, fill = group))+
  geom_boxplot()


#####################FoxO信号通路的ssGSEA分析###################################
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
gene = list("MAPK signaling pathway" = gene)

# 
expr = log2(expr1+1)
# 
re <- gsva(as.matrix(expr), gene, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)
head(re)
re %>%
  t() %>%
  as.data.frame() %>%
  write.table("mapk_re.txt", sep = "\t", quote = F)

foxO = read.table("mapk_re.txt", header = T, sep = "\t")
foxO$sample = rownames(foxO)
foxO$group = ifelse(foxO$sample %in% samplesTP, "TP", "NT")
head(foxO)
ggplot(foxO, aes(group, MAPK.signaling.pathway,fill = group))+
  geom_boxplot()

