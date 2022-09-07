# command: Rscript script.R fq.list 

library(scran)
library(BiocParallel)
library(scater)
library(batchelor)
library(igraph)
library(BiocNeighbors)
library(SingleCellExperiment)
library(Seurat)


args = commandArgs(T)

## retained genes 
features.IR <- read.table("/genome/human.GRCh38/gene_name_type.TR_IG.txt", stringsAsFactors = F)
features.prot <- read.table("/genome/human.GRCh38/features.prot.tsv", stringsAsFactors = F)
features <- c(features.prot$V2,features.IR$V2) # 20234



## input sce file list  -------
fq.list <- read.table(args[1],stringsAsFactors = F)

data.dir <- './'



## Step 1: read file and build seurat list ========
fq.list2 <- list()
for( i in 1:nrow(fq.list) ){
  # i=1
  coi <- as.character(fq.list[i,1])
  tmp <- readRDS( paste0(data.dir, coi, "_count.rds") )
  
  tmp <- tmp[rownames(tmp) %in% features,]  # select protein and TR_IG genes 
  colData(tmp) <- NULL
  rowData(tmp) <- NULL
  colData(tmp)[,'Sample'] <- coi

  # as.seurat 
  tmp_seu <- CreateSeuratObject(counts = counts(tmp),  
                                            min.cells = 0, 
                                            min.features = 100) #  20002 937355
  tmp_seu[['Sample']] <- coi
  # assign( coi,  tmp)
  
  # as.list
  tmp_list <- list(tmp_seu)
  names(tmp_list) <- coi
  
  # combined to the whole list
  fq.list2 <- c(fq.list2, tmp_list )
  
  print(paste(coi,'is completed at',Sys.time()))
}




## Step 2: seurat combines list ========

combined_seu <- merge(fq.list2[[1]], fq.list2[-1], merge.data = FALSE)

combined_seu <- NormalizeData(combined_seu, normalization.method = 'LogNormalize', scale.factor = 10000)

combined_seu[["percent.mt"]] <- PercentageFeatureSet(combined_seu, pattern = "^MT-")

print(paste( 'combination is completed at',Sys.time() ) )



## Step 3: save rds ========

saveRDS(combined_seu,'Seurat_com.counts.rds')

