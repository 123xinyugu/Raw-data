library(Seurat)
library(tidyverse)
library(dplyr)

scRNA = readRDS("final_tecs.rds")

head(scRNA)
dim(scRNA) # 14528   726

tecs2 = subset(scRNA, cell_type == "TECs 2")
tecs2@meta.data$cell_type = droplevels(tecs2@meta.data$cell_type)
head(tecs2)

# 
exprMat <- as.matrix(tecs2@assays$RNA@counts)
head(scRNA@meta.data)

cellInfo = tecs2@meta.data[, c("cell_type", "nCount_RNA", "nFeature_RNA")]
head(cellInfo)


# 
saveRDS(exprMat, "./exprMat.rds")
saveRDS(cellInfo, "./cellInfo.rds")







library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
rm(list=ls())
##====##
dir.create("SCENIC")
dir.create("SCENIC/int")


# 
exprMat = readRDS("exprMat.rds")
cellInfo = readRDS("cellInfo.rds")

##
mydbDIR <- "./cisTarget"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather",
           "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")


scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNSCC")
# saveRDS(scenicOptions, "int/scenicOptions.rds")

genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)

exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)



runSCENIC_1_coexNetwork2modules(scenicOptions)

runSCENIC_2_createRegulons(scenicOptions, coexMethod = "top10perTarget")

exprMat_all <- as.matrix(exprMat)
exprMat_all <- log2(exprMat_all+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)

