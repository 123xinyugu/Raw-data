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
scRNA = readRDS("scRNA.rds")
dim(scRNA@meta.data)
head(scRNA)

# 
set1 = read.delim("./GSE125449_Set1_samples.txt/GSE125449_Set1_samples.txt", sep = "\t")
set2 = read.delim("./GSE125449_Set2_samples.txt/GSE125449_Set2_samples.txt", sep = "\t")
head(set1)
head(set2)
set =rbind(set1,set2)
dim(set)
head(set)
hcc = c("S02_P01_LCP21","S07_P02_LCP28","S10_P05_LCP23","S12_P07_LCP30","S15_P09_LCP38",
        "S16_P10_LCP18","S21_P13_LCP37","S351_P10_LCP34","S364_P21_LCP65")
hcc_sample = subset(set, Sample %in% hcc)
dim(hcc_sample)
# 
scRNA@meta.data$Cell.Barcode = substr(rownames(scRNA@meta.data), 6, nchar(rownames(scRNA@meta.data)))
head(scRNA)
scRNA = subset(scRNA, Cell.Barcode %in% hcc_sample$Cell.Barcode)
head(scRNA@meta.data)

#
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

# 
colors = c4a("poly.alphabet2",26)

# 
p1 = VlnPlot(
  scRNA, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  cols = colors,
  pt.size = 0.1, 
  ncol = 3
)
p1

# 
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<5500&percent.mt<15)
p2 = VlnPlot(
  scRNA, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  cols = colors,
  pt.size = 0, 
  ncol = 3
)
p2
dim(scRNA) # 

# 
# scRNA = SCTransform(scRNA, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = F)
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
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:30, reduction="harmony")

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
# scRNA.markers <- FindAllMarkers(scRNA, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# write.csv(scRNA.markers, "markers.csv")

# 
VlnPlot(scRNA, features=c("AKR1C1"), pt.size=0,group.by = "seurat_clusters")+
  NoLegend()+
  theme(axis.title.x=element_blank())

# cluster 0: Tumor-associated endothelial cells: VWF,ENG,CDH5
# cluster 1: Liver bud hepatic cells: CPB2,MASP2,CFHR1,UGT2B4
# cluster 2: T cells:CD3E,CD3G,CD2,CD3D
# cluster 3: Myofibroblast: COL1A1, COL3A1
# cluster 4: B cells:CD79A,SLAMF7,FCRL5,MZB1
# cluster 5: Tumor-associated macrophages: CD163,CD68,CSF1R
# cluster 6: T cells
# cluster 7: B cells

tec = c('VWF','ENG','CDH5')
lbc = c('CPB2','MASP2','CFHR1','UGT2B4')
tcells = c('CD3E','CD3G','CD2','CD3D')
myc = c('COL1A1', 'COL3A1')
bcells = c('CD79A','SLAMF7','FCRL5','MZB1')
tam = c('CD163','CD68','CSF1R')s

anno_markers = c(tec,lbc,tcells,myc,bcells,tam)

# cluster 
cell_label = c("Tumor-associated endothelial cells","Liver bud hepatic cells","T cells","Myofibroblast",
               "B cells","Tumor-associated macrophages","T cells","B cells")
## 
names(cell_label) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, cell_label)
scRNA[["cell_type"]] = Idents(scRNA)


# 
colors = c("#ff85a8","#ffd685","#fff385","#85ffd0","#85c2ff","#ff85f3")
p8 = UMAPPlot(scRNA, pt.size=1, label=T, label.size=5, cols = colors)+
  NoLegend()
p8

# 
p9 = DotPlot(scRNA, features=anno_markers, cols=c("#10d0ee", "#ff00b2"))+coord_flip()+
  theme(
    axis.text.x=element_text(angle=30, hjust=1, size=10, face="bold"), 
    axis.text.y=element_text(face="bold", size=12), 
    axis.title.x=element_blank(), 
    axis.title.y=element_blank()
  )
p9

genes = c("ENG","CFHR1","CD3D","COL3A1","MZB1","CD68")
# 
p10 = VlnPlot(scRNA, features=genes, pt.size=0, ncol = 3, cols = colors)+
  NoLegend()+
  theme(axis.title.x=element_blank())
p10

head(scRNA)

# 
cell_count = scRNA@meta.data %>%
  group_by(cell_type) %>%
  count()
View(cell_count)

p11 = ggplot(cell_count, aes(reorder(cell_type, -n, median), n)) +
  geom_col(aes(fill = cell_type)) +
  geom_text(aes(label = paste0(n)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5)+
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_color_manual() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "no") +
  labs(title = "Cell Count", x = "", y = "Count")
p11

# 
# write.csv(cell_count,"cell_count.csv",row.names = F)
# 
saveRDS(scRNA, file = "./final_scRNA.rds")
