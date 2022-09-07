library(Seurat)
library(ggplot2)

args = commandArgs(T)

# input and output file 
seurat_file <- 'Seurat_com.counts.rds'
nhvg = 1500
nPC = 10
nres = 1

seurat_file <- args[1]

nhvg = as.numeric(args[2])
nPC = as.numeric(args[3])
nres = as.numeric(args[4])
coi <- gsub('.rds','',seurat_file)


res_file <- paste0(coi,'.seurat_rm.hvg',nhvg,'_PC',nPC,'_res',nres,'.rds') 
res_figure <- paste0(coi,'.seurat_rm.plot_hvg',nhvg,'_PC',nPC,'_res',nres,'.png')


# read rm genes 
GeneMT <- read.table("/genome/human.GRCh38/gene_name_type.MT.txt",stringsAsFactors = F)
GeneHSP <- read.table("/genome/human.GRCh38/gene_name_type.HSP.txt",stringsAsFactors = F)
GeneRP <- read.table("/genome/human.GRCh38/gene_name_type.RP.txt",stringsAsFactors = F)
GeneSC <- read.table("/genome/human.GRCh38/gene_name_type.sc_dissociation_induced_gene.txt",stringsAsFactors = F)

GenesRM <- rbind(GeneMT, GeneHSP, GeneRP,GeneSC)
GenesRM <- GenesRM$V2

## a. read seurat file 
seu_data <- readRDS(seurat_file)

print(paste0(coi,': read input RDS is finished'))


# # seu_data <- readRDS('Seurat_A113.counts.rds')
# 
# seu_data_count <- GetAssayData(seu_data, slot = 'counts')
# seu_data_meta <- seu_data$Sample 
# rm(seu_data)
# 
# seu_data <- CreateSeuratObject(counts = seu_data_count,  
#                                  min.cells = 0, 
#                                  min.features = 100) #  20002 937355
# 
# seu_data[['Sample']] <- seu_data_meta
# # saveRDS(seu_data,'Seurat_A113.counts2.rds')
# print(paste0(coi,': CreateSeuratObject is finished'))



## Step 1: proprocessing ==========

# default: regress mt
# seu_data[["percent.mt"]] <- PercentageFeatureSet(seu_data, pattern = "^MT-")
# print(paste0(coi,': PercentageFeatureSet is finished'))

# VlnPlot(seu_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#         ncol = 1, pt.size = 0, group.by = 'batch')
# p1 <- VlnPlot(seu_data, features = c("nFeature_RNA"),
#               ncol = 1, pt.size = 0, group.by = 'Sample') + ylim(0,10000) + NoLegend()
# p2 <- VlnPlot(seu_data, features = c("nCount_RNA"),
#               ncol = 1, pt.size = 0, group.by = 'Sample') + ylim(0,15000) + NoLegend()
# p3 <- VlnPlot(seu_data, features = c("percent.mt"),
#               ncol = 1, pt.size = 0, group.by = 'Sample') + ylim(0,80) + NoLegend()
# cowplot::plot_grid(p2,p3,ncol=1)
# ggsave('Rplot_neu38_percent.mt.pdf', scale=2)


## filter 
# # according to read counts  
# dim(seu_data)
# length( coi <- rownames( seu_data )[ which(rowSums(seu_data@assays$RNA@counts) > 100) ] ) # 5802
# length( coi <- coi[! coi %in% GenesRM[,2] ]  )  # 5628
# seu_data <- seu_data[coi,]

# gene filter 
seu_data <- seu_data[! rownames(seu_data) %in% GenesRM,]


## Step 2: PCA and cluster =========

seu_data <- NormalizeData(seu_data, normalization.method = 'LogNormalize', scale.factor = 10000)

seu_data <- FindVariableFeatures(seu_data, selection.method = "vst", nfeatures = nhvg)

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# VariableFeatures <- VariableFeatures(object = seu_data)
# VariableFeatures <- VariableFeatures[! VariableFeatures %in% 
#                                        c(GenesRM[,2],'GAPDH',VariableFeaturesSup) ] #,'ZFAND2A','IER5','GBP5'
# sort(VariableFeatures)


seu_data <- ScaleData(seu_data, vars.to.regress = "percent.mt", features = VariableFeatures(seu_data))
# seu_data <- ScaleData(seu_data, vars.to.regress = "percent.mt", features = rownames(seu_data))

print(paste0(coi,': ScateData is finished'))


## run PCA 
seu_data <- RunPCA(seu_data, features = VariableFeatures(seu_data))

print(paste0(coi,': RunPCA is finished'))
# ElbowPlot(seu_data)

# Examine and visualize PCA results a few different ways
# print(seu_data[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(seu_data, dims = 1:2, reduction = "pca")
# DimHeatmap(seu_data, dims = 1, cells = 500, balanced = TRUE)
# 
# # NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:15)


seu_data <- FindNeighbors(seu_data, dims = 1:nPC)

seu_data <- FindClusters(seu_data, resolution = nres)

print(paste0(coi,': FindClusters is finished'))

seu_data <- RunUMAP(seu_data, dims = 1:nPC )

print(paste0(coi,': RunUMAP is finished'))

saveRDS(seu_data, res_file)

png(res_figure, width = 2000,height = 900)
p1 <- DimPlot(seu_data, reduction = "umap",label = T)
p2 <- DimPlot(seu_data, reduction = "umap",group.by = 'Sample')
p1 + p2
# ggsave(p1+p2,paste0('Rplot_hvg2000_PC10_res0.6_',coi,'.png'))
dev.off()
