# Rscript findmarker.R seurat|sce seurat.count.rds 

library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(ggplot2)


# test data
input_file <- 'Seurat_tmp.counts.seurat.hvg2000_PC15_res1.rds' # seurat.rds


args = commandArgs(T)
input_format <- args[1]
input_file <- args[2]

coi <- gsub('.rds','',input_file)
output_png <- paste0(coi,'.marker.png')
output_file <- paste0(coi,'.marker.txt')


if(input_format == 'seurat'){
  ## a. seurat 
  seu_data <- readRDS( input_file )
  seu_data
  
  if( ! 'clusters' %in%  colnames(seu_data@meta.data) ) seu_data$clusters <- seu_data$seurat_clusters
  # seu_data@reductions$UMAP <- Reductions(seu_data, 'umap')
}


if(input_format == 'sce'){
  ## b. scanpy.sce 
  sce_data <- readRDS( input_file )
  sce_data
  assayNames(sce_data) <- 'counts'
  
  # sce to seurat 
  seu_data <- as.Seurat(sce_data, data = NULL)
  seu_data$clusters <- seu_data$louvain
  
}

Idents(seu_data) <- seu_data$clusters
seu_data_gene <- rownames(seu_data)
# seu_data@meta.data[,'Type'] <- gsub('A..._([A-Z_]*)[0-9]?$','\\1',seu_data$Sample)



# seu_data$group <- 'N'
# seu_data@meta.data[seu_data$clusters %in% c(1,2,25,28),]$group <- 'Neu'
# Idents(seu_data) <- seu_data$group

markers <- FindAllMarkers(seu_data, slot='data')

markers_group <- markers %>% group_by(cluster)
write.table(markers,output_file)


markers_top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
markers_top5 <- markers_top5[!duplicated(markers_top5$gene),]
DotPlot(seu_data, features = markers_top5$gene) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(output_png, width = 20, height = 10)
