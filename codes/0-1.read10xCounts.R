library(dplyr)
library(DropletUtils)
library(BiocParallel)
library(SingleCellExperiment)
library(scater)
library(scRNAseq)
library(igraph)
library(scran)
library(SeuratObject)
library(Seurat)



args=commandArgs(T)

## Step 1:  data file ===========
coi <- args[1]

# data.dir <- paste0("/zhanglab/home/caoqi/Proj_01_scRNA/1.cellranger/",coi,"/outs/filtered_feature_bc_matrix/")
# data.dir <- paste0("/datb/Proj_scPLC_hs/",coi,"/outs/filtered_feature_bc_matrix/")
data.dir <- paste0("/zhanglab/home/caoqi/Proj_01_scRNA/1.rename_cellranger/",coi,"/outs/filtered_feature_bc_matrix/")


## Step 2: read 10X file ============
sce <- read10xCounts(data.dir, col.names = T)
# rownames and colnames
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
colnames(sce) <- paste0(coi,"_",1:ncol(sce))
colData(sce)[,'Sample'] <- coi
# save file
saveRDS(sce, paste0('./sce_raw_A124/',coi,'_raw.rds'))
