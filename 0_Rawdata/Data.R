# 
library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(cols4all)
library(glmGamPoi) # 
library(cowplot)
library(clustree)

## 
### 
assays <- dir("./HCC/")
dir <- paste0("./HCC/", assays)
# 
samples_name = assays
samples_name

# 
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}

### 
names(scRNAlist) <- samples_name

# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])

# 
saveRDS(scRNA, file = "./scRNA.rds")
